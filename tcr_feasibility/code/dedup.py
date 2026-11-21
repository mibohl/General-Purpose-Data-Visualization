"""
Clonal / near-duplicate control for CDR3-beta sequences.

Two operations:
  1. collapse_exact:   drop exact-duplicate (cdr3b, vb) rows (clonal identicals).
  2. greedy_cluster:   single-representative greedy clustering at a high
     normalized-identity threshold, so that an entire cluster of near-identical
     CDR3s can be kept on one side of a support/query split.

Splitting is done by CLUSTER (not by row), guaranteeing no cluster spans the
support and query sets -> no clonal / near-duplicate leakage across the split.
"""
import numpy as np
import pandas as pd
from rapidfuzz import process, fuzz
from rapidfuzz.distance import Levenshtein


def collapse_exact(df, keys=("cdr3b", "vb")):
    """Drop exact duplicate TCRs; return (deduped_df, n_removed)."""
    n0 = len(df)
    out = df.drop_duplicates(subset=list(keys)).reset_index(drop=True)
    return out, n0 - len(out)


def greedy_cluster(seqs, threshold=0.90):
    """
    Greedy single-linkage-ish clustering on normalized CDR3 identity.
    A sequence joins the first existing cluster whose representative is within
    (1 - threshold) normalized Levenshtein distance; otherwise it starts a new
    cluster. Returns an int cluster-id array aligned to `seqs`.

    threshold=0.90 => sequences sharing >=90% identity cluster together.
    Length-bucketing keeps comparisons cheap and correct (identity across very
    different lengths is already low).
    """
    seqs = list(seqs)
    n = len(seqs)
    cluster_id = np.full(n, -1, dtype=int)
    reps = []          # (rep_seq, cluster_index)
    reps_by_len = {}   # len -> list of (rep_seq, cluster_index)
    next_c = 0
    # process longest-first for stable representatives
    order = sorted(range(n), key=lambda i: -len(seqs[i]))
    for i in order:
        s = seqs[i]
        L = len(s)
        best_c = -1
        # compare against representatives of similar length (+/- 2 residues)
        for dl in (0, -1, 1, -2, 2):
            bucket = reps_by_len.get(L + dl)
            if not bucket:
                continue
            cand = [r[0] for r in bucket]
            # max normalized identity match
            m = process.extractOne(
                s, cand, scorer=fuzz.ratio, score_cutoff=threshold * 100.0
            )
            if m is not None:
                best_c = bucket[m[2]][1]
                break
        if best_c == -1:
            best_c = next_c
            next_c += 1
            reps.append((s, best_c))
            reps_by_len.setdefault(L, []).append((s, best_c))
        cluster_id[i] = best_c
    return cluster_id


def cross_identity_mask(query_seqs, ref_seqs, threshold=0.90):
    """
    Boolean mask over query_seqs: True where a query sequence has >=threshold
    identity to ANY ref sequence. Used to (a) drop other-epitope negatives that
    are near-identical to held-out positives (public/cross-reactive TCRs), and
    (b) audit residual support<->query leakage.
    """
    ref = list(ref_seqs)
    if not ref:
        return np.zeros(len(query_seqs), dtype=bool)
    mask = np.zeros(len(query_seqs), dtype=bool)
    # bucket refs by length for speed
    from collections import defaultdict
    rb = defaultdict(list)
    for r in ref:
        rb[len(r)].append(r)
    for i, s in enumerate(query_seqs):
        L = len(s)
        for dl in (0, -1, 1, -2, 2):
            cand = rb.get(L + dl)
            if not cand:
                continue
            m = process.extractOne(s, cand, scorer=fuzz.ratio,
                                    score_cutoff=threshold * 100.0)
            if m is not None:
                mask[i] = True
                break
    return mask
