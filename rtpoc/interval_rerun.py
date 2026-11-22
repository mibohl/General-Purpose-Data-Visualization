"""Interval metrics for the stacking arms and for fine-tuned DeepLC 4.5.

The earlier runs stored only the mean width of TabPFN's 90% intervals and no
interval for the fine-tuned model. This reruns, with the same splits, seeds
and settings:
- T1_d45 (TabPFN on peptide features + 8 DeepLC 4.5 heads) and T1_d45_B4
  (plus the cross-fitted B4 prediction), saving the 5% / 50% / 95% quantiles.
  rel_mae must reproduce the stored values (checked in the aggregate).
- FT (deeplc.finetune with config D, then calibrate), predicting the queries
  and the 15% validation part of the context. The query predictions must match
  the ones saved by run_ft.py. FT gets a split-conformal 90% interval: pred
  +- the ceil(0.9 (n_val + 1))-th smallest absolute validation residual. The
  validation part was used for early stopping and is inside the calibration
  fit, so this interval is, if anything, a little narrow.
- FT_spread (descriptive, uses query labels, not a method): the central 90%
  range of FT's query residuals, i.e. the width an oracle interval around FT
  would need.
Output: results/finetune/parts/IV_<name>.csv (resumable),
results/finetune/ivpreds/<name>_s<seed>_k<k>.npz
Usage: python rtpoc/interval_rerun.py <name> --label easy|hard
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ft_d45 as FT  # noqa: E402
import methods as M  # noqa: E402
import run as R  # noqa: E402
from heldout_generic import load_d45  # noqa: E402
from run_ft import b4  # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "finetune", "parts")
PDIR = os.path.join(ROOT, "results", "finetune", "ivpreds")
KS = [100, 250, 500, 1000]
CFG = dict(lr=0.01, epochs=150, patience=20)


def iv_metrics(pred, lo, hi, true, rng):
    w = hi - lo
    return dict(rel_mae=float(np.abs(pred - true).mean() / rng),
                coverage90=float(np.mean((true >= lo) & (true <= hi))),
                width90_median=float(np.median(w) / rng), width90_mean=float(np.mean(w) / rng))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--label", required=True, choices=["easy", "hard"])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--ks", type=int, nargs="+", default=KS)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(1)
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(PDIR, exist_ok=True)
    name = a.name
    path = os.path.join(OUT, f"IV_{name}.csv")
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{name}.parquet"))
    D = load_d45(name)
    F = M.peptide_features(df, extra_mods=["TMTpro", "Cysteinyl"])
    cat = [F.columns.get_loc(c) for c in M.CAT_COLS]
    tr = df["tr"].values.astype(float)
    rng = tr.max() - tr.min()
    rows = pd.read_csv(path).to_dict("records") if os.path.exists(path) else []
    done = {(r["seed"], r["n"]) for r in rows}
    for seed in range(a.seeds):
        q, pool = R.split(len(tr), seed)
        for k in a.ks:
            if (seed, k) in done:
                continue
            t0 = time.time()
            ctx = pool[:k]
            base = dict(dataset=name, group=a.label, n=k, seed=seed)
            lo_c, hi_c = tr[ctx].min(), tr[ctx].max()
            sc = hi_c - lo_c
            yc = (tr[ctx] - lo_c) / sc
            back = lambda z: z * sc + lo_c  # noqa: E731
            pq_b4, oof_b4, _, heads = b4(D, ctx, q, yc, seed)
            Xc0 = np.hstack([F.values[ctx].astype(float), np.asarray(D[ctx][:, heads], dtype=float)])
            Xq0 = np.hstack([F.values[q].astype(float), np.asarray(D[q][:, heads], dtype=float)])
            save = {}
            r = []
            for meth, Xc, Xq in (("T1_d45", Xc0, Xq0),
                                 ("T1_d45_B4", np.hstack([Xc0, oof_b4[:, None]]),
                                  np.hstack([Xq0, pq_b4[:, None]]))):
                med, q05, q95, _, t = M.tabpfn_fit_predict(Xc, yc, Xq, seed, cat)
                med, q05, q95 = back(med), back(q05), back(q95)
                save[meth] = np.vstack([med, q05, q95])
                r.append(dict(base, method=meth, seconds=t, **iv_metrics(med, q05, q95, tr[q], rng)))
            _, va = FT.val_split(np.asarray(ctx), seed)
            t1 = time.time()
            pred, _, t_ft, _, ep = FT.finetune_predict(df, ctx, np.concatenate([q, va]), seed, CFG, 1)
            p_q, p_va = pred[:len(q)], pred[len(q):]
            res = np.sort(np.abs(p_va - tr[va]))
            j = min(len(res), math.ceil(0.9 * (len(res) + 1))) - 1
            h = res[j]
            m = iv_metrics(p_q, p_q - h, p_q + h, tr[q], rng)
            z = np.load(os.path.join(ROOT, "results", "finetune", "preds", f"{name}_s{seed}_k{k}.npz"))
            diff = float(np.max(np.abs(z["pred_ft"] - p_q)) / rng)
            r.append(dict(base, method="FT", seconds=time.time() - t1, n_val=len(va),
                          max_abs_diff_vs_main=diff, ft_epochs=ep - 1, **m))
            e = p_q - tr[q]
            r.append(dict(base, method="FT_spread", rel_mae=m["rel_mae"], coverage90=0.9,
                          width90_median=float((np.quantile(e, 0.95) - np.quantile(e, 0.05)) / rng),
                          width90_mean=float(2 * np.quantile(np.abs(e), 0.9) / rng)))
            np.savez(os.path.join(PDIR, f"{name}_s{seed}_k{k}.npz"), q=q, FT=p_q, FT_h=h, **save)
            rows.extend(r)
            pd.DataFrame(rows).to_csv(path, index=False)
            print(f"  IV seed={seed} k={k} {time.time() - t0:.0f}s diff={diff:.2e} "
                  + " ".join(f"{x['method']}:{x['rel_mae']:.4f}/{x['coverage90']:.2f}/"
                             f"{x['width90_median']:.4f}" for x in r), flush=True)


if __name__ == "__main__":
    main()
