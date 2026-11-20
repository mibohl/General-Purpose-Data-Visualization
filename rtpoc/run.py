"""Run checks C1-C5 for one dataset; writes results/parts/<check>_<dataset>.csv.

Usage: python rtpoc/run.py <dataset> [--checks C1 C2 C3] [--seeds 5]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from datasets import GROUP  # noqa: E402
import methods as M  # noqa: E402

ROOT = os.path.dirname(HERE)
PROC = os.path.join(ROOT, "data", "processed")
PARTS = os.path.join(ROOT, "results", "parts")
N_QUERY = 2000
NS = [25, 50, 100, 250, 500, 1000]
NOISE = [0.05, 0.10, 0.20]


def load(name):
    df = pd.read_parquet(os.path.join(PROC, name + ".parquet"))
    F = M.peptide_features(df)
    Foh = M.onehot_terms(F)
    mods = df["modifications"].astype(str)
    # "modified" = carries any modification other than fixed Carbamidomethyl
    is_mod = mods.apply(lambda s: any(MOD != "Carbamidomethyl" for MOD in s.split("|")[1::2]) if s else False)
    return df, F, Foh, is_mod.values


def split(N, seed):
    perm = np.random.default_rng(seed).permutation(N)
    nq = min(N_QUERY, N - max(NS))
    return perm[:nq], perm[nq:]


def metrics(pred, true, rt_range, lo=None, hi=None):
    err = np.abs(pred - true)
    d = dict(rel_mae=err.mean() / rt_range, rel_dt95=2 * np.quantile(err, 0.95) / rt_range)
    if lo is not None:
        d["coverage90"] = float(np.mean((true >= lo) & (true <= hi)))
        d["width90"] = float(np.mean(hi - lo)) / rt_range
    return d


def rows_for(base, method, pred, true, rt_range, is_mod_q, runtime, lo=None, hi=None):
    out = []
    for subset, mask in (("all", np.ones_like(is_mod_q, bool)), ("modified", is_mod_q),
                         ("unmodified", ~is_mod_q)):
        if mask.sum() < 20:
            continue
        m = metrics(pred[mask], true[mask], rt_range,
                    None if lo is None else lo[mask], None if hi is None else hi[mask])
        out.append({**base, "method": method, "subset": subset, "n_query": int(mask.sum()),
                    **m, "runtime_s": runtime})
    return out


def run_context(F, Foh, tr, is_mod, ctx, q, seed, rt_range, base, use_tabpfn=True,
                methods=("B1", "B2", "B2_deeplc", "B3_ridge", "B3_xgb", "T1", "T2")):
    """Fit all methods on context `ctx` (tr may be corrupted), evaluate on queries `q`."""
    rows = []
    lo_c, hi_c = tr[ctx].min(), tr[ctx].max()
    scale = hi_c - lo_c
    yc = (tr[ctx] - lo_c) / scale                      # context-only rescaling
    back = lambda z: z * scale + lo_c                   # noqa: E731
    true_q = TRUE[q]                                    # always evaluate on clean labels
    gc, gq = F["g"].values[ctx], F["g"].values[q]
    imq = is_mod[q]

    def add(name, pred_scaled, t, lo=None, hi=None):
        rows.extend(rows_for(base, name, back(pred_scaled), true_q, rt_range, imq, t,
                             None if lo is None else back(lo), None if hi is None else back(hi)))

    if "B1" in methods:
        t = time.time(); p = M.Linear().fit(gc, yc).predict(gq); add("B1", p, time.time() - t)
    t = time.time()
    b2 = M.Lowess(seed=seed).fit(gc, yc)
    b2q = b2.predict(gq)
    t_b2 = time.time() - t
    t = time.time()
    oof = M.crossfit_b2(gc, yc, seed)
    rc = yc - oof                                       # cross-fitted context residuals
    t_cf = time.time() - t
    h = np.quantile(np.abs(rc), 0.90)                   # constant 90% window from context
    if "B2" in methods:
        add("B2", b2q, t_b2, b2q - h, b2q + h)
    if "B2_deeplc" in methods:
        try:
            t = time.time(); p = M.DeepLCSpline().fit(gc, yc).predict(gq)
            add("B2_deeplc", p, time.time() - t)
        except Exception as e:  # noqa: BLE001
            print(f"  B2_deeplc failed (n={len(ctx)}): {type(e).__name__}: {e}", flush=True)
    Xc_oh, Xq_oh = Foh.values[ctx].astype(float), Foh.values[q].astype(float)
    if "B3_ridge" in methods:
        t = time.time(); r = M.ridge_model().fit(Xc_oh, rc).predict(Xq_oh)
        add("B3_ridge", b2q + r, t_b2 + t_cf + time.time() - t)
    if "B3_xgb" in methods:
        t = time.time(); r = M.xgb_model(seed).fit(Xc_oh, rc).predict(Xq_oh)
        add("B3_xgb", b2q + r, t_b2 + t_cf + time.time() - t)
    if use_tabpfn:
        Xc, Xq = F.values[ctx].astype(float), F.values[q].astype(float)
        cat = [F.columns.get_loc(c) for c in M.CAT_COLS]
        if "T1" in methods:
            med, q05, q95, _, t = M.tabpfn_fit_predict(Xc, yc, Xq, seed, cat)
            add("T1", med, t, q05, q95)
        if "T2" in methods:
            med, q05, q95, _, t = M.tabpfn_fit_predict(Xc, rc, Xq, seed, cat)
            add("T2", b2q + med, t + t_b2 + t_cf, b2q + q05, b2q + q95)
    return rows, b2.frac_


def c1(name, F, Foh, tr, is_mod, seeds):
    """Headroom: large context (all but the query set); residual R^2 of B3 on B2 residual."""
    rows = []
    rt_range = tr.max() - tr.min()
    for seed in range(seeds):
        q, pool = split(len(tr), seed)
        lo, hi = tr[pool].min(), tr[pool].max()
        yc, yq = (tr[pool] - lo) / (hi - lo), (tr[q] - lo) / (hi - lo)
        gc, gq = F["g"].values[pool], F["g"].values[q]
        t = time.time()
        b2 = M.Lowess(seed=seed).fit(gc, yc)
        b2q = b2.predict(gq)
        oof = M.crossfit_b2(gc, yc, seed)
        rc, rq = yc - oof, yq - b2q
        tb = time.time() - t
        base = dict(check="C1", dataset=name, group=GROUP[name], n=len(pool), seed=seed, noise=0.0)
        Xc, Xq = Foh.values[pool].astype(float), Foh.values[q].astype(float)
        for mname, model in (("B3_ridge", M.ridge_model()), ("B3_xgb", M.xgb_model(seed))):
            t = time.time()
            rhat = model.fit(Xc, rc).predict(Xq)
            dt = time.time() - t
            for subset, mask in (("all", np.ones(len(q), bool)), ("modified", is_mod[q]),
                                 ("unmodified", ~is_mod[q])):
                if mask.sum() < 20:
                    continue
                pred_rt = (b2q + rhat)[mask] * (hi - lo) + lo
                b2_rt = b2q[mask] * (hi - lo) + lo
                rows.append({**base, "method": mname, "subset": subset, "n_query": int(mask.sum()),
                             "resid_r2": r2_score(rq[mask], rhat[mask]),
                             **metrics(pred_rt, tr[q][mask], rt_range),
                             "b2_rel_mae": np.abs(b2_rt - tr[q][mask]).mean() / rt_range,
                             "runtime_s": dt + tb, "b2_frac": b2.frac_})
        print(f"  C1 seed {seed}: B2 frac={b2.frac_} " + " ".join(
            f"{r['method']} R2={r['resid_r2']:.3f}" for r in rows[-6:] if r["subset"] == "all"
            and r["seed"] == seed), flush=True)
    return rows


def main():
    global TRUE
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--checks", nargs="+", default=["C1", "C2", "C3"])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ns", type=int, nargs="+", default=NS)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.threads)
    os.makedirs(PARTS, exist_ok=True)
    name = a.dataset
    df, F, Foh, is_mod = load(name)
    tr = df["tr"].values.astype(float)
    TRUE = tr
    rt_range = tr.max() - tr.min()
    print(f"{name}: N={len(df)} range={rt_range:.1f} modified={is_mod.mean():.2f}", flush=True)

    if "C1" in a.checks:
        rows = c1(name, F, Foh, tr, is_mod, min(a.seeds, 3))
        pd.DataFrame(rows).to_csv(os.path.join(PARTS, f"C1_{name}.csv"), index=False)

    if "C2" in a.checks:
        rows = []
        for seed in range(a.seeds):
            q, pool = split(len(tr), seed)
            for n in a.ns:
                ctx = pool[:n]
                base = dict(check="C2", dataset=name, group=GROUP[name], n=n, seed=seed, noise=0.0)
                t = time.time()
                r, frac = run_context(F, Foh, tr, is_mod, ctx, q, seed, rt_range, base)
                for x in r:
                    x["b2_frac"] = frac
                rows.extend(r)
                summ = {x["method"]: round(x["rel_mae"], 4) for x in r if x["subset"] == "all"}
                print(f"  C2 seed={seed} n={n} {time.time()-t:.0f}s {summ}", flush=True)
            pd.DataFrame(rows).to_csv(os.path.join(PARTS, f"C2_{name}.csv"), index=False)

    if "C3" in a.checks:
        rows = []
        n = 250
        for seed in range(a.seeds):
            q, pool = split(len(tr), seed)
            ctx = pool[:n]
            rng = np.random.default_rng(1000 + seed)
            for noise in NOISE:
                trn = tr.copy()
                k = int(round(noise * n))
                bad = rng.choice(ctx, size=k, replace=False)
                trn[bad] = rng.uniform(tr[ctx].min(), tr[ctx].max(), size=k)
                base = dict(check="C3", dataset=name, group=GROUP[name], n=n, seed=seed, noise=noise)
                r, frac = run_context(F, Foh, trn, is_mod, ctx, q, seed, rt_range, base,
                                      methods=("B2", "B3_ridge", "B3_xgb", "T2"))
                rows.extend(r)
                summ = {x["method"]: round(x["rel_mae"], 4) for x in r if x["subset"] == "all"}
                print(f"  C3 seed={seed} noise={noise} {summ}", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(PARTS, f"C3_{name}.csv"), index=False)


TRUE = None

if __name__ == "__main__":
    main()
