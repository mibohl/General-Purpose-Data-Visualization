"""Generic predictions for the held-out datasets (data/heldout/<name>.csv).

Per peptidoform, inference only:
- pred_generic: legacy DeepLC ensemble, uncalibrated (same as the main study)
- DeepLC 4.5 default multitask model (6,543 LC-setup heads): the full
  prediction matrix, stored as its low-rank factors (projections per
  peptide + per-head embedding/scale/shift) so any context can be calibrated
  later with DeepLC's own MultiHeadRidgeCalibration.

Output: data/processed/heldout_<name>.parquet and heldout_<name>_d45.npz

Usage: python rtpoc/heldout_generic.py <name> [<name> ...]
       (also accepts a Zenodo dataset name with --zenodo, as a leakage sanity check)
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
from psm_utils import PSM, PSMList
from psm_utils.io.peptide_record import peprec_to_proforma

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generic_pred as G  # noqa: E402
from deeplc import core  # noqa: E402
from deeplc._factored import FactoredPredictionMatrix  # noqa: E402

ROOT = os.path.dirname(HERE)
HELD = os.path.join(ROOT, "data", "heldout")
OUT = os.path.join(ROOT, "data", "processed")


def load(name, zenodo=False):
    path = os.path.join(G.RAW if zenodo else HELD, name + ".csv")
    df = pd.read_csv(path)
    df["modifications"] = df["modifications"].fillna("").astype(str)
    return df.groupby(["seq", "modifications"], as_index=False, sort=False)["tr"].mean()


def d45_factors(df):
    psms = PSMList(psm_list=[
        PSM(peptidoform=peprec_to_proforma(s, m if m else "-"), spectrum_id=str(i))
        for i, (s, m) in enumerate(zip(df.seq, df.modifications))])
    mat = core.predict(psms, return_matrix=True)
    assert isinstance(mat, FactoredPredictionMatrix), type(mat)
    return dict(projections=mat._projections, embedding=mat._embedding, scale=mat._scale,
                shift=mat._shift)


def load_d45(name):
    z = np.load(os.path.join(OUT, f"heldout_{name}_d45.npz"))
    return FactoredPredictionMatrix(z["projections"], z["embedding"], z["scale"], z["shift"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--zenodo", action="store_true")
    a = ap.parse_args()
    logging.getLogger("feat_extractor").setLevel(logging.ERROR)
    torch.set_num_threads(4)
    models = [torch.load(os.path.join(G.MODEL_DIR, f), weights_only=False,
                         map_location="cpu").eval() for f in G.LEGACY]
    fe = G.FeatExtractor(verbose=False)
    for name in a.names:
        t0 = time.time()
        df = load(name, a.zenodo)
        pred, per_model = G.predict_legacy(G.features(df, fe), models)
        df["pred_generic"] = pred
        f = d45_factors(df)
        df.to_parquet(os.path.join(OUT, f"heldout_{name}.parquet"), index=False)
        np.savez_compressed(os.path.join(OUT, f"heldout_{name}_d45.npz"), **f)
        m = load_d45(name)
        rows = np.random.default_rng(0).permutation(len(df))[:5000]
        dense = np.asarray(m[rows], dtype=np.float64)          # (<=5000, 6543)
        t = df.tr.values[rows]
        dz = (dense - dense.mean(0)) / (dense.std(0) + 1e-12)
        best = np.max(np.abs(dz.T @ ((t - t.mean()) / t.std()) / len(t)))
        print(f"{name:32s} n={len(df):7d} r(tr, legacy)={np.corrcoef(df.tr, pred)[0, 1]:+.3f} "
              f"max|r| over DeepLC 4.5 heads={best:.3f} {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
