"""DeepLC 4.5 fine-tuning vs calibration vs TabPFN stacking on the held-out setups.

Same splits, seeds and metric as run_heldout.py (R.split, contexts are the first
k of the pool, up to 2,000 queries, relative MAE / full RT range of the setup).
k in {100, 250, 500, 1000}, seeds 0-4. Fine-tuning hyperparameters are fixed
once (ft_tune.py on Part 1 datasets) and passed in.

Methods written here (B4 and T1_d45 are re-used from results/parts/H_<name>.csv;
B4 is recomputed here too, as B4_check, to confirm the splits match):
- FT: deeplc.finetune on the context (85% train / 15% validation, early stopping)
  then deeplc.calibrate on the whole context; i.e. finetune_and_predict with our
  validation split.
- T1_d45_FT (phase "stack"): T1_d45 inputs + the fine-tuned DeepLC prediction.
  The context gets 3-fold cross-fitted fine-tuned predictions (each fold
  fine-tuned on the other 2/3 only, with its own 15% validation split), the
  queries get the full-context fine-tune saved by phase "main". So TabPFN never
  sees an in-sample fine-tuned value. 3 folds rather than 5 to fit the compute.
- T1_d45_B4: T1_d45 inputs + the B4 calibrated prediction (context: 5-fold
  cross-fitted, queries: full-context fit).
Output: results/finetune/parts/{FT,ST}_<name>.csv (resumable)
Usage: python rtpoc/run_ft.py <name> --label easy|hard --phase main|stack --lr .. --epochs .. --patience ..
"""
from __future__ import annotations

import argparse
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
from deeplc.calibration import MultiHeadRidgeCalibration  # noqa: E402
from sklearn.model_selection import KFold  # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "finetune", "parts")
KS = [100, 250, 500, 1000]
K_HEADS = 8


def b4(D, ctx, q, yc, seed):
    t = time.time()
    cal = MultiHeadRidgeCalibration()
    cal.fit(target=yc, source=D[ctx])
    pq = np.asarray(cal.transform(D[q]), dtype=float)
    dt = time.time() - t
    oof = np.empty_like(yc)
    for tr_, va in KFold(5, shuffle=True, random_state=seed).split(yc):
        c = MultiHeadRidgeCalibration()
        c.fit(target=yc[tr_], source=D[ctx[tr_]])
        oof[va] = np.asarray(c.transform(D[ctx[va]]), dtype=float)
    return pq, oof, dt, cal._head_idx[:K_HEADS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--label", required=True, choices=["easy", "hard"])
    ap.add_argument("--mods", nargs="*", default=["TMTpro", "Cysteinyl"])
    ap.add_argument("--lr", type=float, required=True)
    ap.add_argument("--epochs", type=int, required=True)
    ap.add_argument("--patience", type=int, required=True)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--ks", type=int, nargs="+", default=KS)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--phase", choices=["main", "stack"], required=True)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.threads)
    cfg = dict(lr=a.lr, epochs=a.epochs, patience=a.patience)
    name = a.name
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{'FT' if a.phase == 'main' else 'ST'}_{name}.csv")
    pdir = os.path.join(ROOT, "results", "finetune", "preds")
    os.makedirs(pdir, exist_ok=True)
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{name}.parquet"))
    D = load_d45(name)
    F = M.peptide_features(df, extra_mods=a.mods)
    cat = [F.columns.get_loc(c) for c in M.CAT_COLS]
    mods = df["modifications"].astype(str)
    is_mod = mods.apply(lambda s: any(x != "Carbamidomethyl" for x in s.split("|")[1::2])
                        if s else False).values
    tr = df["tr"].values.astype(float)
    rt_range = tr.max() - tr.min()
    rows = pd.read_csv(path).to_dict("records") if os.path.exists(path) else []
    done = {(r["seed"], r["n"]) for r in rows}
    for seed in range(a.seeds):
        q, pool = R.split(len(tr), seed)
        for k in a.ks:
            if (seed, k) in done:
                continue
            t0 = time.time()
            ctx = pool[:k]
            base = dict(check="FT", dataset=name, group=a.label, n=k, seed=seed, noise=0.0,
                        threads=a.threads, **{f"ft_{x}": v for x, v in cfg.items()})
            lo_c, hi_c = tr[ctx].min(), tr[ctx].max()
            sc = hi_c - lo_c
            yc = (tr[ctx] - lo_c) / sc
            back = lambda z: z * sc + lo_c  # noqa: E731
            r = []
            pq_b4, oof_b4, dt_b4, heads = b4(D, ctx, q, yc, seed)
            Hc = np.asarray(D[ctx][:, heads], dtype=float)
            Hq = np.asarray(D[q][:, heads], dtype=float)
            Xc0 = np.hstack([F.values[ctx].astype(float), Hc])
            Xq0 = np.hstack([F.values[q].astype(float), Hq])
            npz = os.path.join(pdir, f"{name}_s{seed}_k{k}.npz")
            t_ft = t_oof = 0.0
            ep = -1
            if a.phase == "main":
                r += R.rows_for(base, "B4_check", back(pq_b4), tr[q], rt_range, is_mod[q], dt_b4)
                pred_ft, raw_ft, t_ft, t_ft_all, ep = FT.finetune_predict(df, ctx, q, seed, cfg,
                                                                          a.threads)
                np.savez(npz, pred_ft=pred_ft, raw_ft=raw_ft, t_ft_all=t_ft_all, q=q)
                for x in R.rows_for(base, "FT", pred_ft, tr[q], rt_range, is_mod[q], t_ft_all):
                    x.update(ft_train_s=t_ft, ft_epochs=ep)
                    r.append(x)
                r += R.rows_for(base, "FT_raw", raw_ft, tr[q], rt_range, is_mod[q], t_ft)
                stacks = [("T1_d45_B4", oof_b4, pq_b4, dt_b4)]
            else:
                z = np.load(npz)
                assert np.array_equal(z["q"], q)
                pred_ft = z["pred_ft"]
                t1 = time.time()
                oof_ft = np.empty(k)
                for tr_, va in KFold(3, shuffle=True, random_state=seed).split(ctx):
                    p, _, _, _, _ = FT.finetune_predict(df, ctx[tr_], ctx[va], seed, cfg, a.threads)
                    oof_ft[va] = p
                t_oof = time.time() - t1
                stacks = [("T1_d45_FT", (oof_ft - lo_c) / sc, (pred_ft - lo_c) / sc,
                           float(z["t_ft_all"]) + t_oof)]
            for meth, ec, eq, extra_t in stacks:
                Xc = np.hstack([Xc0, ec[:, None]])
                Xq = np.hstack([Xq0, eq[:, None]])
                med, q05, q95, _, t = M.tabpfn_fit_predict(Xc, yc, Xq, seed, cat)
                for x in R.rows_for(base, meth, back(med), tr[q], rt_range, is_mod[q],
                                    t + extra_t, back(q05), back(q95)):
                    x.update(tabpfn_s=t, oof_ft_s=t_oof)
                    r.append(x)
            rows.extend(r)
            summ = {x["method"]: round(x["rel_mae"], 4) for x in r if x["subset"] == "all"}
            print(f"  FT seed={seed} k={k} {time.time() - t0:.0f}s ft={t_ft:.0f}s ep={ep} "
                  f"oof={t_oof:.0f}s {summ}", flush=True)
            pd.DataFrame(rows).to_csv(path, index=False)


if __name__ == "__main__":
    main()
