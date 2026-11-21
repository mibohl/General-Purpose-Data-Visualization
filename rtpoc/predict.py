"""Predict RTs for a new LC setup from a few hundred identified peptides.

    python rtpoc/predict.py --context ctx.csv --query query.csv --out pred.csv \
        [--method T1_d45|T1|T2] [--mods TMTpro ...]

Input CSVs use the PEPREC columns of the study: seq, modifications
("pos|Name|..." with Unimod names, 0 = N-term, -1 = C-term; empty if none);
the context file also needs tr (observed RT, any unit). Output: query rows
with rt_pred (median) and a 90% interval rt_lo90 / rt_hi90, in the context's
RT unit.

Methods (see results/SUMMARY.md):
- T1_d45 (default): TabPFN v2 with peptide features, the legacy DeepLC
  prediction and the raw predictions of the 8 DeepLC 4.5 setup heads that
  correlate best with the context.
- T1: TabPFN v2 with peptide features and the legacy DeepLC prediction only.
- T2: as T1 but predicting the residual of a LOWESS calibration.
Nothing is trained; every fit uses the context only.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generic_pred as G  # noqa: E402
import methods as M  # noqa: E402
from heldout_generic import d45_factors  # noqa: E402

K_HEADS = 8


def read(path, need_tr):
    df = pd.read_csv(path)
    df["modifications"] = df["modifications"].fillna("").astype(str)
    if need_tr and "tr" not in df:
        raise SystemExit(f"{path}: context needs a 'tr' column")
    return df


def predict(ctx: pd.DataFrame, qry: pd.DataFrame, method="T1_d45", extra_mods=(), seed=0):
    allp = pd.concat([ctx[["seq", "modifications"]], qry[["seq", "modifications"]]],
                     ignore_index=True)
    logging.getLogger("feat_extractor").setLevel(logging.ERROR)
    models = [torch.load(os.path.join(G.MODEL_DIR, f), weights_only=False,
                         map_location="cpu").eval() for f in G.LEGACY]
    allp["pred_generic"], _ = G.predict_legacy(G.features(allp, G.FeatExtractor(verbose=False)),
                                               models)
    F = M.peptide_features(allp, extra_mods=extra_mods)
    nc = len(ctx)
    tr = ctx["tr"].values.astype(float)
    lo, hi = tr.min(), tr.max()
    y = (tr - lo) / (hi - lo)
    X = F.values.astype(float)
    cat = [F.columns.get_loc(c) for c in M.CAT_COLS]
    if method == "T1_d45":
        from deeplc._factored import FactoredPredictionMatrix
        from deeplc.calibration.multihead import _rank_heads_by_correlation
        mat = FactoredPredictionMatrix(**d45_factors(allp))
        # same ranking DeepLC's MultiHeadRidgeCalibration uses (context labels only)
        heads = _rank_heads_by_correlation(mat[np.arange(nc)], y)[:K_HEADS]
        H = np.asarray(mat[:, heads], dtype=float)
        X = np.hstack([X, H])
    if method == "T2":
        g = F["g"].values
        b2 = M.Lowess(seed=seed).fit(g[:nc], y)
        rc = y - M.crossfit_b2(g[:nc], y, seed)
        med, q05, q95, _, _ = M.tabpfn_fit_predict(X[:nc], rc, X[nc:], seed, cat)
        off = b2.predict(g[nc:])
        med, q05, q95 = med + off, q05 + off, q95 + off
    else:
        med, q05, q95, _, _ = M.tabpfn_fit_predict(X[:nc], y, X[nc:], seed, cat)
    back = lambda v: v * (hi - lo) + lo  # noqa: E731
    out = qry.copy()
    out["rt_pred"], out["rt_lo90"], out["rt_hi90"] = back(med), back(q05), back(q95)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--context", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--method", default="T1_d45", choices=["T1_d45", "T1", "T2"])
    ap.add_argument("--mods", nargs="*", default=[])
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    out = predict(read(a.context, True), read(a.query, False), a.method, a.mods)
    out.to_csv(a.out, index=False)
    print(f"wrote {len(out)} predictions to {a.out}")


if __name__ == "__main__":
    main()
