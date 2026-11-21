"""Kumo Tabular (NVIDIA open tabular foundation model, `sdm`) on the SAME
held-out splits as run_experiment.py, so TabPFN and Kumo are compared paired
(identical allele / k / seed / query set / BLOSUM features).

Adds methods `kumo_reg` (regression on 1-log50k; point estimate = mean of the
999 predictive quantiles) and `kumo_clf` (classification on the 500 nM binder
label). Writes results_kumo.csv incrementally; resumable.

num_estimators=1 to match TabPFN's single-estimator setting and keep the CPU
sweep tractable (Kumo `small`, ~28M params, CPU).
"""
import os, sys, time, argparse, warnings
import numpy as np
import pandas as pd
import torch

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import data as D
from encoders import blosum_encode
import baselines as B
# reuse the EXACT split logic so query sets match the TabPFN run
from run_experiment import k_list_for, SEEDS, QUERY_CAP

import sdm

ROOT = os.path.join(os.path.dirname(__file__), "..")
RESULTS = os.path.join(ROOT, "results", "results_kumo.csv")
KUMO_SIZE = "small"
KUMO_N_ESTIMATORS = 1
FCOLS = [f"f{i}" for i in range(180)]

COLS = ["allele", "n_allele", "k", "seed", "query_size", "n_support_binder",
        "method", "auroc", "auprc", "spearman", "runtime_s"]


def run_allele(allele, sub, reg, clf, dev):
    peps = sub["peptide"].tolist()
    X = blosum_encode(peps, length=9)
    y = sub["y"].values.astype(np.float32)
    yb = sub["binder"].values.astype(int)
    n = len(sub)

    # Build per-allele TableTensors once (reg target numerical, clf categorical)
    dfr = pd.DataFrame(X, columns=FCOLS); dfr["tgt"] = y
    dfc = pd.DataFrame(X, columns=FCOLS); dfc["tgt"] = yb
    tr = sdm.TableTensor.from_pandas(
        df=dfr, stypes=sdm.infer_stypes(dfr, overrides={"tgt": "numerical"}), device=dev)
    tc = sdm.TableTensor.from_pandas(
        df=dfc, stypes=sdm.infer_stypes(dfc, overrides={"tgt": "categorical"}), device=dev)
    xr = tr.drop_columns("tgt"); yr = tr[:, "tgt"]
    xc = tc.drop_columns("tgt"); yc = tc[:, "tgt"]

    rows = []
    for k in k_list_for(n):
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            perm = rng.permutation(n)
            sup = perm[:k]
            qpool = perm[k:]
            qry = (rng.choice(qpool, QUERY_CAP, replace=False)
                   if len(qpool) > QUERY_CAP else qpool)
            sup_t = torch.as_tensor(sup); qry_t = torch.as_tensor(qry)
            ycq, ybq = y[qry], yb[qry]
            n_sup_bind = int(yb[sup].sum())

            def rec(method, score, rt):
                m = B.metrics(score, ycq, ybq)
                rows.append(dict(allele=allele, n_allele=n, k=int(k), seed=seed,
                                 query_size=len(qry), n_support_binder=n_sup_bind,
                                 method=method, auroc=m["auroc"], auprc=m["auprc"],
                                 spearman=m["spearman"], runtime_s=rt))

            # --- Kumo regressor (point = mean of predictive quantiles) ---
            t = time.time()
            try:
                out = reg(x_context=xr[sup_t], y_context=yr[sup_t],
                          x_query=xr[qry_t], num_estimators=KUMO_N_ESTIMATORS)
                pr = out.to_pandas().mean(axis=1).values
            except Exception as e:
                pr = None
            rec("kumo_reg", pr, time.time() - t)

            # --- Kumo classifier (positive-class prob) ---
            t = time.time()
            try:
                if len(np.unique(yb[sup])) < 2:
                    pc = np.full(len(qry), float(yb[sup][0]))
                else:
                    out = clf(x_context=xc[sup_t], y_context=yc[sup_t],
                              x_query=xc[qry_t], num_estimators=KUMO_N_ESTIMATORS)
                    pdf = out.to_pandas()
                    pc = pdf["1"].values if "1" in pdf.columns else pdf.iloc[:, 0].values
            except Exception as e:
                pc = None
            rec("kumo_clf", pc, time.time() - t)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alleles", type=int, default=0)
    args = ap.parse_args()

    df = D.load_raw(); df9 = D.hla_9mers(df)
    held_out, highdata, g = D.define_splits(df9)
    if args.alleles > 0:
        held_out = held_out[:args.alleles]

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    done = set()
    if os.path.exists(RESULTS):
        try:
            done = set(pd.read_csv(RESULTS)["allele"].unique())
        except Exception:
            done = set()
    todo = [a for a in held_out if a not in done]
    print(f"KUMO: {len(held_out)} alleles, {len(done)} done, {len(todo)} to run", flush=True)

    dev = torch.device("cpu")
    reg = sdm.models.KumoTabular(task="regression", size=KUMO_SIZE, device=dev)
    clf = sdm.models.KumoTabular(task="classification", size=KUMO_SIZE, device=dev)

    header_needed = not os.path.exists(RESULTS)
    for ai, allele in enumerate(todo):
        t0 = time.time()
        sub = df9[df9["allele"] == allele].reset_index(drop=True)
        rows = run_allele(allele, sub, reg, clf, dev)
        pd.DataFrame(rows)[COLS].to_csv(RESULTS, mode="a", index=False,
                                        header=header_needed)
        header_needed = False
        print(f"[{ai+1}/{len(todo)}] {allele} n={len(sub)} -> {len(rows)} rows "
              f"in {time.time()-t0:.0f}s", flush=True)
    print("KUMO DONE", flush=True)


if __name__ == "__main__":
    main()
