"""Held-out comparison against DeepLC 4.5 (setups it never saw).

Same protocol as C2 (nested contexts n = 25..1000, 5 seeds, up to 2,000
queries, all fits on the context only). Methods:
- B1, B2, B2_deeplc, B3_ridge, B3_xgb, T1, T2: as in the main study (legacy
  generic DeepLC prediction as the generic feature).
- B4: DeepLC 4.5 default model + its own MultiHeadRidgeCalibration (80 best
  heads, spline per head, ridge on top), fit on the context. This is what
  `deeplc.predict_and_calibrate(queries, psm_list_reference=context)` does.
- T1_d45: T1 plus the raw predictions of the 8 DeepLC 4.5 heads that B4 ranks
  highest on the context (ranking uses context labels only).
Also C4-style intervals for T1 / T1_d45 and a constant window for B4.

Output: results/parts/H_<name>.csv
Usage: python rtpoc/run_heldout.py <name> --label easy|hard [--mods GG TMT6plex ...]
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
import methods as M  # noqa: E402
import run as R  # noqa: E402
from heldout_generic import load_d45  # noqa: E402
from deeplc.calibration import MultiHeadRidgeCalibration  # noqa: E402
from sklearn.model_selection import KFold  # noqa: E402

ROOT = os.path.dirname(HERE)
K_HEADS = 8


def b4_fit_predict(D, ctx, q, y_ctx, seed):
    """DeepLC 4.5 calibration on context (y in context-scaled units)."""
    t = time.time()
    cal = MultiHeadRidgeCalibration()
    cal.fit(target=y_ctx, source=D[ctx])
    pq = np.asarray(cal.transform(D[q]), dtype=float)
    dt = time.time() - t
    # cross-fitted context residuals for a constant 90% window (as for B2)
    oof = np.empty_like(y_ctx)
    for tr_, va in KFold(5, shuffle=True, random_state=seed).split(y_ctx):
        c = MultiHeadRidgeCalibration()
        c.fit(target=y_ctx[tr_], source=D[ctx[tr_]])
        oof[va] = np.asarray(c.transform(D[ctx[va]]), dtype=float)
    h = np.quantile(np.abs(y_ctx - oof), 0.90)
    return pq, h, dt, cal._head_idx[:K_HEADS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--label", required=True, choices=["easy", "hard"])
    ap.add_argument("--mods", nargs="*", default=[])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--ns", type=int, nargs="+", default=R.NS)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.threads)
    name = a.name
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{name}.parquet"))
    D = load_d45(name)
    F = M.peptide_features(df, extra_mods=a.mods)
    Foh = M.onehot_terms(F)
    mods = df["modifications"].astype(str)
    is_mod = mods.apply(lambda s: any(x != "Carbamidomethyl" for x in s.split("|")[1::2])
                        if s else False).values
    tr = df["tr"].values.astype(float)
    R.TRUE = tr
    rt_range = tr.max() - tr.min()
    print(f"{name}: N={len(df)} range={rt_range:.1f} modified={is_mod.mean():.2f}", flush=True)
    rows = []
    for seed in range(a.seeds):
        q, pool = R.split(len(tr), seed)
        for n in a.ns:
            ctx = pool[:n]
            base = dict(check="H", dataset=name, group=a.label, n=n, seed=seed, noise=0.0)
            t0 = time.time()
            r, frac = R.run_context(F, Foh, tr, is_mod, ctx, q, seed, rt_range, base)
            lo_c, hi_c = tr[ctx].min(), tr[ctx].max()
            sc = hi_c - lo_c
            yc = (tr[ctx] - lo_c) / sc
            back = lambda z: z * sc + lo_c  # noqa: E731
            pq, h, dt, heads = b4_fit_predict(D, ctx, q, yc, seed)
            r += R.rows_for(base, "B4", back(pq), tr[q], rt_range, is_mod[q], dt,
                            back(pq - h), back(pq + h))
            Hc = np.asarray(D[ctx][:, heads], dtype=float)
            Hq = np.asarray(D[q][:, heads], dtype=float)
            Xc = np.hstack([F.values[ctx].astype(float), Hc])
            Xq = np.hstack([F.values[q].astype(float), Hq])
            cat = [F.columns.get_loc(c) for c in M.CAT_COLS]
            med, q05, q95, _, t = M.tabpfn_fit_predict(Xc, yc, Xq, seed, cat)
            r += R.rows_for(base, "T1_d45", back(med), tr[q], rt_range, is_mod[q], t + dt,
                            back(q05), back(q95))
            for x in r:
                x["b2_frac"] = frac
            rows.extend(r)
            summ = {x["method"]: round(x["rel_mae"], 4) for x in r if x["subset"] == "all"}
            print(f"  H seed={seed} n={n} {time.time() - t0:.0f}s {summ}", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(R.PARTS, f"H_{name}.csv"), index=False)


if __name__ == "__main__":
    main()
