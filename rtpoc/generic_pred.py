"""Step 1a: uncalibrated generic RT predictions from the legacy DeepLC ensemble.

The legacy ensemble (full_hc_PXD005573_pub_*, 3 models, trained on PXD005573)
ships inside the deeplc 4.5 wheel as onnx2torch-converted GraphModules. Input
features are built with DeepLC 3.1.13's own FeatExtractor (vendored), exactly
as DeepLC 3.x did: inputs = [matrix, matrix_sum, concat(matrix_all, pos_matrix),
matrix_hc]. With several models and calibrate=False, DeepLC 3.x averages the
models' raw outputs; we do the same. Inference only, CPU.

Output: data/processed/<dataset>.parquet with seq, modifications, tr, pred_generic.
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
sys.path.insert(0, os.path.join(HERE, "legacy_deeplc"))
from feat_extractor import FeatExtractor  # noqa: E402

import deeplc  # noqa: E402

ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "zenodo", "train_eval", "datasets")
OUT = os.path.join(ROOT, "data", "processed")
MODEL_DIR = os.path.join(os.path.dirname(deeplc.__file__), "package_data", "models")
LEGACY = [
    "full_hc_PXD005573_pub_1fd8363d9af9dcad3be7553c39396960.pt",
    "full_hc_PXD005573_pub_8c22d89667368f2f02ad996469ba157e.pt",
    "full_hc_PXD005573_pub_cb975cfdd4105f97efa0b3afffe075cc.pt",
]

from datasets import ALL  # noqa: E402


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(RAW, name + ".csv"))
    df["modifications"] = df["modifications"].fillna("").astype(str)
    df = df.groupby(["seq", "modifications"], as_index=False, sort=False)["tr"].mean()
    return df


def features(df: pd.DataFrame, fe: FeatExtractor):
    psms = PSMList(psm_list=[
        PSM(peptidoform=peprec_to_proforma(s, m if m else "-"), spectrum_id=str(i))
        for i, (s, m) in enumerate(zip(df.seq, df.modifications))
    ])
    X = fe.encode_atoms(psms, list(range(len(psms))))
    Xm = np.stack(list(X["matrix"].values())).astype(np.float32)
    Xs = np.stack(list(X["matrix_sum"].values())).astype(np.float32)
    Xg = np.concatenate([np.stack(list(X["matrix_all"].values())),
                         np.stack(list(X["pos_matrix"].values()))], axis=1).astype(np.float32)
    Xh = np.stack(list(X["matrix_hc"].values())).astype(np.float32)
    return Xm, Xs, Xg, Xh


@torch.no_grad()
def predict_legacy(inputs, models, batch=4096):
    n = inputs[0].shape[0]
    preds = []
    for m in models:
        out = []
        for i in range(0, n, batch):
            xs = [torch.from_numpy(x[i:i + batch]) for x in inputs]
            out.append(m(*xs).reshape(-1).numpy())
        preds.append(np.concatenate(out))
    return np.mean(preds, axis=0), np.stack(preds, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="*", default=ALL)
    args = ap.parse_args()
    logging.getLogger("feat_extractor").setLevel(logging.ERROR)
    torch.set_num_threads(4)
    os.makedirs(OUT, exist_ok=True)
    # trusted source: weights shipped in the deeplc wheel from PyPI
    models = [torch.load(os.path.join(MODEL_DIR, f), weights_only=False, map_location="cpu").eval()
              for f in LEGACY]
    fe = FeatExtractor(verbose=False)
    for name in args.datasets:
        t0 = time.time()
        df = load_dataset(name)
        inp = features(df, fe)
        pred, per_model = predict_legacy(inp, models)
        df["pred_generic"] = pred
        for k in range(per_model.shape[1]):
            df[f"pred_m{k}"] = per_model[:, k]
        df.to_parquet(os.path.join(OUT, name + ".parquet"), index=False)
        r = np.corrcoef(df.tr, df.pred_generic)[0, 1]
        print(f"{name:32s} n={len(df):7d} pearson(tr,pred)={r:+.3f} "
              f"model-agree r={np.corrcoef(per_model.T)[0,1]:.3f} {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
