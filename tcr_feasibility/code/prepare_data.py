"""
Build cached, cleaned, deduplicated per-epitope positive tables and the
background control pool. Records provenance + dedup counts.

Outputs (in ../results/cache/):
    positives.parquet   cleaned, exact-deduped, clustered positives
    controls.parquet    cleaned, exact-deduped background negatives
    data_provenance.json
"""
import os
import json
import numpy as np
import pandas as pd

import common
import dedup

OUT = os.path.join(common.REPO_ROOT, "results", "cache")
os.makedirs(OUT, exist_ok=True)

MIN_POS = 100          # keep epitopes with >= this many raw labeled (positive) TCRs
IDENT_THRESHOLD = 0.90


def main():
    print("Loading IMMREP22 positives ...")
    pos = common.load_immrep_positives()
    print("Loading IMMREP22 background controls ...")
    ctrl = common.load_immrep_controls()

    prov = {
        "dataset": "IMMREP22 (Working Group 6) TCR specificity benchmark",
        "source_url": "https://github.com/viragbioinfo/IMMREP_2022_TCRSpecificity",
        "positive_source": "VDJdb-derived, tetramer/dextramer-sorted, paired alpha-beta",
        "background_negative_source": "unpublished healthy-donor CD8+CD96+ TCRs (A. Eugster, 11 donors)",
        "ident_threshold": IDENT_THRESHOLD,
        "min_positives_kept": MIN_POS,
        "feature_note": "CDR3-beta only for primary; epitope is constant within a task so no epitope features used within-task.",
        "epitopes": {},
        "n_raw_positive_rows": int(len(pos)),
        "n_control_rows_raw": int(len(ctrl)),
    }

    # --- clean controls ---
    ctrl["valid"] = ctrl["cdr3b"].map(common.valid_cdr3)
    ctrl = ctrl[ctrl["valid"]].drop(columns=["valid"]).reset_index(drop=True)
    ctrl, c_removed = dedup.collapse_exact(ctrl, keys=("cdr3b", "vb"))
    prov["n_control_valid_unique"] = int(len(ctrl))
    prov["n_control_exact_dups_removed"] = int(c_removed)
    ctrl["cluster"] = dedup.greedy_cluster(ctrl["cdr3b"].tolist(),
                                           threshold=IDENT_THRESHOLD)
    ctrl.to_parquet(os.path.join(OUT, "controls.parquet"))

    # --- per-epitope positives ---
    kept_rows = []
    for ep, g in pos.groupby("epitope"):
        n_raw = len(g)
        g = g.copy()
        g["valid"] = g["cdr3b"].map(common.valid_cdr3)
        n_valid = int(g["valid"].sum())
        g = g[g["valid"]].drop(columns=["valid"])
        g, removed = dedup.collapse_exact(g, keys=("cdr3b", "vb"))
        n_unique = len(g)
        clusters = dedup.greedy_cluster(g["cdr3b"].tolist(),
                                        threshold=IDENT_THRESHOLD)
        g["cluster"] = clusters
        n_clusters = int(len(set(clusters)))
        keep = n_raw >= MIN_POS
        prov["epitopes"][ep] = {
            "n_raw_positive": int(n_raw),
            "n_valid": n_valid,
            "n_exact_dups_removed": int(removed),
            "n_unique_positive": int(n_unique),
            "n_clusters_90pct": n_clusters,
            "largest_cluster": int(pd.Series(clusters).value_counts().iloc[0]) if n_unique else 0,
            "kept": bool(keep),
        }
        if keep:
            kept_rows.append(g)

    allpos = pd.concat(kept_rows, ignore_index=True)
    allpos.to_parquet(os.path.join(OUT, "positives.parquet"))
    prov["kept_epitopes"] = sorted(allpos["epitope"].unique().tolist())
    prov["n_kept_epitopes"] = len(prov["kept_epitopes"])

    with open(os.path.join(OUT, "data_provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)

    print(f"\nKept {prov['n_kept_epitopes']} epitopes (>= {MIN_POS} raw positives):")
    for ep in prov["kept_epitopes"]:
        e = prov["epitopes"][ep]
        print(f"  {ep:12s} raw={e['n_raw_positive']:4d} unique={e['n_unique_positive']:4d} "
              f"clusters90={e['n_clusters_90pct']:4d} maxclust={e['largest_cluster']:3d} "
              f"exact_dups_rm={e['n_exact_dups_removed']}")
    print(f"\nControls: {prov['n_control_valid_unique']} unique "
          f"(removed {prov['n_control_exact_dups_removed']} exact dups)")
    print(f"Cache written to {OUT}")


if __name__ == "__main__":
    main()
