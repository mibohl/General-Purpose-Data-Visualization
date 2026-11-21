"""Speed / accuracy of TabPFN (T1) for realistic query counts on CPU.

Context n in {250, 1000}; queries 2k to 100k from hela_hf_psms_aligned (161k
peptides). Grid over n_estimators (4, 8) and fit_mode (fit_preprocessors =
what the main runs used; fit_with_cache = context key/value cache reused
across query blocks). Queries are predicted in blocks of 10,000.
Run on an otherwise idle machine. Output: results/timing.csv

Usage: python rtpoc/timing.py [--threads 4] [--queries 2000 20000 100000]
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

ROOT = os.path.dirname(HERE)
BLOCK = 10_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="hela_hf_psms_aligned")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ns", type=int, nargs="+", default=[250, 1000])
    ap.add_argument("--queries", type=int, nargs="+", default=[2000, 20000, 100000])
    ap.add_argument("--n_est", type=int, nargs="+", default=[4, 8])
    ap.add_argument("--modes", nargs="+", default=["fit_preprocessors", "fit_with_cache"])
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.threads)
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", a.dataset + ".parquet"))
    F = M.peptide_features(df)
    tr = df["tr"].values.astype(float)
    rt_range = tr.max() - tr.min()
    perm = np.random.default_rng(0).permutation(len(df))
    qmax = max(a.queries)
    q_all, pool = perm[:qmax], perm[qmax:]
    cat = [F.columns.get_loc(c) for c in M.CAT_COLS]
    out_path = os.path.join(ROOT, "results", "timing.csv")
    rows = pd.read_csv(out_path).to_dict("records") if os.path.exists(out_path) else []
    for n in a.ns:
        ctx = pool[:n]
        lo, hi = tr[ctx].min(), tr[ctx].max()
        yc = (tr[ctx] - lo) / (hi - lo)
        Xc = F.values[ctx].astype(float)
        for n_est in a.n_est:
            for mode in a.modes:
                for nq in sorted(a.queries):
                    q = q_all[:nq]
                    Xq = F.values[q].astype(float)
                    m = M.tabpfn_model(0, n_estimators=n_est)
                    m.set_params(categorical_features_indices=cat, fit_mode=mode)
                    t0 = time.time()
                    m.fit(Xc, yc)
                    t_fit = time.time() - t0
                    med = []
                    for i in range(0, nq, BLOCK):
                        med.append(np.asarray(m.predict(Xq[i:i + BLOCK], output_type="quantiles",
                                                        quantiles=[0.05, 0.5, 0.95])[1]))
                    t_all = time.time() - t0
                    pred = np.concatenate(med) * (hi - lo) + lo
                    r = dict(dataset=a.dataset, n=n, n_query=nq, n_estimators=n_est,
                             fit_mode=mode, threads=a.threads, fit_s=round(t_fit, 2),
                             total_s=round(t_all, 2),
                             ms_per_query=round(1000 * t_all / nq, 3),
                             rel_mae=np.abs(pred - tr[q]).mean() / rt_range)
                    rows.append(r)
                    print(r, flush=True)
                    pd.DataFrame(rows).to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
