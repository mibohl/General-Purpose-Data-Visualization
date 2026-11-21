"""
Annotate every TCR with its germline CDR1/CDR2/CDR2.5 (pmhc) loops for BOTH
chains, derived from the V gene via tcrdist3's own germline tables, plus CDR3.
These are exactly the loops the TCRdist metric uses. Writes:
    positives_matched.parquet, controls_matched.parquet   (+ 8 loop columns)
    matched_colmask.npy   (global zero-variance feature mask; label-free)
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import common

CACHE = os.path.join(common.REPO_ROOT, "results", "cache")
KEY = ["cdr3a", "va", "ja", "cdr3b", "vb", "jb"]


def annotate(df):
    from tcrdist.repertoire import TCRrep
    cell = pd.DataFrame({
        "cdr3_a_aa": df.cdr3a.astype(str), "v_a_gene": df.va.astype(str), "j_a_gene": df.ja.astype(str),
        "cdr3_b_aa": df.cdr3b.astype(str), "v_b_gene": df.vb.astype(str), "j_b_gene": df.jb.astype(str),
        "count": 1})
    tr = TCRrep(cell_df=cell, organism="human", chains=["alpha", "beta"],
                compute_distances=False)
    c = tr.clone_df
    c = c.assign(
        _ka=c.cdr3_a_aa + "|" + c.v_a_gene + "|" + c.j_a_gene,
        _kb=c.cdr3_b_aa + "|" + c.v_b_gene + "|" + c.j_b_gene)
    look = {}
    for _, r in c.iterrows():
        look[(r.cdr3_a_aa, r.v_a_gene, r.j_a_gene, r.cdr3_b_aa, r.v_b_gene, r.j_b_gene)] = (
            r.cdr1_a_aa, r.cdr2_a_aa, r.pmhc_a_aa, r.cdr3_a_aa,
            r.cdr1_b_aa, r.cdr2_b_aa, r.pmhc_b_aa, r.cdr3_b_aa)
    out = {k: [] for k in common.LOOP_COLS}
    miss = 0
    for _, r in df.iterrows():
        key = (str(r.cdr3a), str(r.va), str(r.ja), str(r.cdr3b), str(r.vb), str(r.jb))
        v = look.get(key)
        if v is None:
            miss += 1
            v = ("",) * 8
        for col, val in zip(common.LOOP_COLS, v):
            out[col].append(val)
    for col in common.LOOP_COLS:
        df[col] = out[col]
    return df, miss


def main():
    pos = pd.read_parquet(os.path.join(CACHE, "positives.parquet"))
    ctrl = pd.read_parquet(os.path.join(CACHE, "controls.parquet"))
    print(f"annotating {len(pos)} positives + {len(ctrl)} controls ...")
    pos, mp = annotate(pos)
    ctrl, mc = annotate(ctrl)
    print(f"missing annotations: positives={mp}, controls={mc}")
    # drop any row that failed to annotate (empty loops)
    pos = pos[pos.cdr1_b.astype(str).str.len() > 0].reset_index(drop=True)
    ctrl = ctrl[ctrl.cdr1_b.astype(str).str.len() > 0].reset_index(drop=True)

    # global zero-variance mask (label-free -> no leakage)
    allm = common.matched_matrix(pd.concat([pos, ctrl], ignore_index=True))
    keep = allm.std(axis=0) > 1e-9
    np.save(os.path.join(CACHE, "matched_colmask.npy"), keep)
    print(f"matched feature dim: {allm.shape[1]} raw -> {int(keep.sum())} after dropping constant cols")

    pos.to_parquet(os.path.join(CACHE, "positives_matched.parquet"))
    ctrl.to_parquet(os.path.join(CACHE, "controls_matched.parquet"))
    # show loop widths
    for col in common.LOOP_COLS:
        print(f"  {col:8s} width {pos[col].astype(str).str.len().max()}")
    print("wrote matched caches.")


if __name__ == "__main__":
    main()
