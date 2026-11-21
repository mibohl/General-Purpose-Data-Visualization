"""Data loading, filtering, splits, and overlap diagnostics for the
Nielsen IEDB MHC class I dataset (TDC: mhc1_iedb-imgt_nielsen).

Target Y is 1-log50k affinity in [0,1]; binder threshold = 500 nM.
"""
import os
import numpy as np
import pandas as pd

DATA_FILE = os.path.join(os.path.dirname(__file__), "..", "data",
                         "mhc1_iedb-imgt_nielsen.tab")

# 500 nM binder cutoff on the 1-log50k scale
BINDER_THRESHOLD = 1.0 - np.log(500.0) / np.log(50000.0)  # ~0.42558

RARE_MIN, RARE_MAX = 50, 1500          # rare-tail definition (measurement count)
HIGHDATA_MIN = 1500                     # pan-allele training pool (> this)
QUAL_MIN_PER_CLASS = 10                 # held-out data-quality filter (per class)
PEPTIDE_LEN = 9                         # first-pass: 9-mers only


def load_raw():
    df = pd.read_csv(DATA_FILE, sep="\t")
    df.columns = ["idx", "peptide", "allele", "pseudoseq", "y"]
    for c in ("peptide", "allele", "pseudoseq"):
        df[c] = df[c].str.strip('"')
    df["length"] = df["peptide"].str.len()
    df["binder"] = (df["y"] >= BINDER_THRESHOLD).astype(int)
    return df


def hla_9mers(df):
    return df[(df["allele"].str.startswith("HLA")) &
              (df["length"] == PEPTIDE_LEN)].reset_index(drop=True)


def allele_stats(df9):
    g = df9.groupby("allele").agg(
        n=("y", "size"),
        n_binder=("binder", "sum"),
    ).reset_index()
    g["n_nonbinder"] = g["n"] - g["n_binder"]
    g["binder_frac"] = g["n_binder"] / g["n"]
    return g.sort_values("n").reset_index(drop=True)


def define_splits(df9):
    """Return (held_out_alleles, highdata_alleles, stats_df).

    Held-out = rare tail (50<=n<=1500) with >=10 of each class (pre-registered,
    result-independent data-quality filter so AUROC/AUPRC are well defined).
    High-data pool = alleles with n>1500 (the pan-allele training pool). No
    held-out allele ever enters that pool.
    """
    g = allele_stats(df9)
    rare = g[(g["n"] >= RARE_MIN) & (g["n"] <= RARE_MAX)]
    qual = rare[(rare["n_binder"] >= QUAL_MIN_PER_CLASS) &
                (rare["n_nonbinder"] >= QUAL_MIN_PER_CLASS)]
    held_out = sorted(qual["allele"].tolist())
    highdata = sorted(g[g["n"] > HIGHDATA_MIN]["allele"].tolist())
    assert not (set(held_out) & set(highdata)), "held-out leaked into pool"
    return held_out, highdata, g


def peptide_overlap_report(df9, held_out, highdata):
    """Diagnostics on identical peptide sequences shared across alleles, and
    specifically between held-out alleles and the high-data pool."""
    pep_alleles = df9.groupby("peptide")["allele"].nunique()
    shared = pep_alleles[pep_alleles > 1]
    pool_peps = set(df9[df9["allele"].isin(highdata)]["peptide"])
    rows = []
    for a in held_out:
        ap = set(df9[df9["allele"] == a]["peptide"])
        inter = ap & pool_peps
        rows.append({"allele": a, "n_peptides": len(ap),
                     "shared_with_pool": len(inter),
                     "frac_shared_with_pool": len(inter) / max(1, len(ap))})
    overlap_df = pd.DataFrame(rows)
    summary = {
        "unique_peptides": int(df9["peptide"].nunique()),
        "peptides_in_multiple_alleles": int((shared).shape[0]),
        "frac_peptides_multiallele": float(shared.shape[0] / df9["peptide"].nunique()),
        "max_alleles_per_peptide": int(pep_alleles.max()),
        "mean_frac_heldout_shared_with_pool": float(overlap_df["frac_shared_with_pool"].mean()),
    }
    return summary, overlap_df


if __name__ == "__main__":
    df = load_raw()
    df9 = hla_9mers(df)
    held, pool, g = define_splits(df9)
    print(f"HLA 9-mers: {len(df9)} rows, {df9.allele.nunique()} alleles")
    print(f"held-out alleles: {len(held)}")
    print(f"high-data pool alleles: {len(pool)} ({df9[df9.allele.isin(pool)].shape[0]} rows)")
    summ, ov = peptide_overlap_report(df9, held, pool)
    print("overlap summary:", summ)
