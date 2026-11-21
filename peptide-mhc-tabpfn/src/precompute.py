"""Precompute k/seed-independent predictions for every held-out peptide:
  - MHCflurry 2.x pan-allele affinity (zero target labels; but trained on IEDB).
  - A pan-allele ridge & logistic model we train ourselves on the high-data
    pool (peptide BLOSUM + MHC pseudoseq BLOSUM), which never sees any held-out
    allele (leakage-free pan-allele reference).

Writes cache/precomputed.parquet with one row per held-out (allele, peptide).
"""
import os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from encoders import blosum_encode
import data as D

CACHE = os.path.join(os.path.dirname(__file__), "..", "cache")


def nm_to_log50k(nm):
    return np.clip(1.0 - np.log(np.maximum(nm, 1e-6)) / np.log(50000.0), 0.0, 1.0)


def pseudoseq_blosum(seqs):
    """BLOSUM-encode the 34-mer MHC pseudosequence -> (N, 34*20)."""
    L = max(len(s) for s in seqs)
    return blosum_encode(seqs, length=L)


def main():
    t0 = time.time()
    df = D.load_raw()
    df9 = D.hla_9mers(df)
    held_out, highdata, g = D.define_splits(df9)
    ho = df9[df9["allele"].isin(held_out)].copy().reset_index(drop=True)
    print(f"held-out rows: {len(ho)} across {len(held_out)} alleles", flush=True)

    # ---- MHCflurry pan-allele affinity ----
    print("loading MHCflurry...", flush=True)
    from mhcflurry import Class1AffinityPredictor
    aff = Class1AffinityPredictor.load()
    mhc_pred = np.full(len(ho), np.nan)
    for a in held_out:
        m = (ho["allele"] == a).values
        nm = aff.predict(peptides=ho.loc[m, "peptide"].tolist(), allele=a)
        mhc_pred[m] = nm_to_log50k(nm)
    ho["mhcflurry"] = mhc_pred
    print(f"MHCflurry done ({time.time()-t0:.0f}s)", flush=True)

    # ---- our own pan-allele linear models (trained on high-data pool) ----
    pool = df9[df9["allele"].isin(highdata)].copy()
    Xp_pool = blosum_encode(pool["peptide"].tolist(), length=9)        # 180
    Xm_pool = pseudoseq_blosum(pool["pseudoseq"].tolist())             # 34*20
    Xpool = np.hstack([Xp_pool, Xm_pool])
    ypool = pool["y"].values.astype(np.float32)
    ybpool = pool["binder"].values

    Xp_ho = blosum_encode(ho["peptide"].tolist(), length=9)
    Xm_ho = pseudoseq_blosum(ho["pseudoseq"].tolist())
    Xho = np.hstack([Xp_ho, Xm_ho])

    from sklearn.linear_model import Ridge, LogisticRegression
    print(f"training pan-allele ridge on {len(pool)} rows x {Xpool.shape[1]} feats...", flush=True)
    r = Ridge(alpha=10.0).fit(Xpool, ypool)
    ho["panallele_ridge"] = np.clip(r.predict(Xho), 0, 1)
    print(f"training pan-allele logistic...", flush=True)
    lr = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs").fit(Xpool, ybpool)
    ho["panallele_logistic"] = lr.predict_proba(Xho)[:, 1]
    print(f"pan-allele models done ({time.time()-t0:.0f}s)", flush=True)

    out = ho[["allele", "peptide", "y", "binder",
              "mhcflurry", "panallele_ridge", "panallele_logistic"]].copy()
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "precomputed.csv")
    out.to_csv(path, index=False)
    print(f"wrote {path} ({len(out)} rows) in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
