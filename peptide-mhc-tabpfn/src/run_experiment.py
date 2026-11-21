"""Main experiment: for each held-out rare HLA allele, each k, and 5 seeds,
score TabPFN (regressor + classifier) and the cheap baselines on the same query
set, plus the k-independent references (MHCflurry, pan-allele linear). Writes
results.csv incrementally and is resumable (skips alleles already finished).

Usage: python src/run_experiment.py [--alleles N] [--query-cap 500]
"""
import os, sys, time, argparse, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import data as D
from encoders import blosum_encode
import baselines as B

ROOT = os.path.join(os.path.dirname(__file__), "..")
CACHE = os.path.join(ROOT, "cache")
RESULTS = os.path.join(ROOT, "results", "results.csv")
KS = [10, 50, 200, 1000]
SEEDS = [0, 1, 2, 3, 4]
QUERY_CAP = 500
# TabPFN ensemble size. n_estimators=1 was verified accuracy-equivalent to the
# default on this task (Spearman 0.719 vs 0.718, AUROC 0.908 vs 0.912) at ~2x
# speed; used here to keep the full CPU sweep tractable.
TABPFN_N_ESTIMATORS = 1

_REG_CKPT = os.path.join(CACHE, "tabpfn-v2-regressor.ckpt")
_CLF_CKPT = os.path.join(CACHE, "tabpfn-v2-classifier-finetuned-zk73skhh.ckpt")

COLS = ["allele", "n_allele", "k", "seed", "query_size", "n_support_binder",
        "method", "auroc", "auprc", "spearman", "runtime_s"]


def k_list_for(n):
    ks = sorted({min(k, n // 2) for k in KS})
    return [k for k in ks if k >= 5 and (n - k) >= 5]


def build_tabpfn():
    from tabpfn import TabPFNRegressor, TabPFNClassifier
    reg = TabPFNRegressor(device="cpu", model_path=_REG_CKPT, n_estimators=TABPFN_N_ESTIMATORS)
    clf = TabPFNClassifier(device="cpu", model_path=_CLF_CKPT, n_estimators=TABPFN_N_ESTIMATORS)
    return reg, clf


def run_allele(allele, sub, pre_sub, reg, clf):
    """sub: allele rows with peptide,y,binder. pre_sub: same rows aligned with
    mhcflurry/panallele columns. Returns list of result dicts."""
    peps = sub["peptide"].tolist()
    X = blosum_encode(peps, length=9)
    y = sub["y"].values.astype(np.float32)
    yb = sub["binder"].values.astype(int)
    n = len(sub)
    mhc = pre_sub["mhcflurry"].values
    paR = pre_sub["panallele_ridge"].values
    paL = pre_sub["panallele_logistic"].values
    rows = []
    for k in k_list_for(n):
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            perm = rng.permutation(n)
            sup = perm[:k]
            qpool = perm[k:]
            if len(qpool) > QUERY_CAP:
                qry = rng.choice(qpool, QUERY_CAP, replace=False)
            else:
                qry = qpool
            ycq, ybq = y[qry], yb[qry]
            n_sup_bind = int(yb[sup].sum())

            def rec(method, score, rt):
                m = B.metrics(score, ycq, ybq)
                rows.append(dict(allele=allele, n_allele=n, k=int(k), seed=seed,
                                 query_size=len(qry), n_support_binder=n_sup_bind,
                                 method=method, auroc=m["auroc"], auprc=m["auprc"],
                                 spearman=m["spearman"], runtime_s=rt))

            # --- TabPFN regressor ---
            t = time.time()
            try:
                reg.fit(X[sup], y[sup]); pr = reg.predict(X[qry])
            except Exception:
                pr = None
            rec("tabpfn_reg", pr, time.time() - t)

            # --- TabPFN classifier ---
            t = time.time()
            try:
                if len(np.unique(yb[sup])) < 2:
                    pc = np.full(len(qry), float(yb[sup][0]))
                else:
                    clf.fit(X[sup], yb[sup]); pc = clf.predict_proba(X[qry])[:, 1]
            except Exception:
                pc = None
            rec("tabpfn_clf", pc, time.time() - t)

            # --- PWM / PSSM ---
            t = time.time()
            pw = B.pwm_logodds_scores([peps[i] for i in sup], yb[sup],
                                      [peps[i] for i in qry], length=9)
            rec("pwm", pw, time.time() - t)

            # --- kNN BLOSUM ---
            t = time.time()
            kn = B.knn_blosum_predict(X[sup], y[sup], X[qry], k_neigh=5)
            rec("knn_blosum", kn, time.time() - t)

            # --- Ridge (same k) ---
            t = time.time()
            rd = B.ridge_predict(X[sup], y[sup], X[qry], alpha=10.0)
            rec("ridge", rd, time.time() - t)

            # --- Logistic (same k) ---
            t = time.time()
            lg = B.logistic_predict(X[sup], yb[sup], X[qry], C=1.0)
            rec("logistic", lg, time.time() - t)

            # --- MHCflurry (k-independent, scored on this query) ---
            rec("mhcflurry", mhc[qry], 0.0)
            # --- MHCflurry recalibrated on the k support peptides ---
            t = time.time()
            try:
                rc = B.recalibrate(mhc[sup], y[sup], mhc[qry], method="isotonic")
            except Exception:
                rc = None
            rec("mhcflurry_recal", rc, time.time() - t)

            # --- our leakage-free pan-allele linear refs ---
            rec("panallele_ridge", paR[qry], 0.0)
            rec("panallele_logistic", paL[qry], 0.0)
    return rows


def main():
    global QUERY_CAP
    ap = argparse.ArgumentParser()
    ap.add_argument("--alleles", type=int, default=0, help="limit #alleles (0=all)")
    ap.add_argument("--query-cap", type=int, default=QUERY_CAP)
    args = ap.parse_args()
    QUERY_CAP = args.query_cap

    df = D.load_raw(); df9 = D.hla_9mers(df)
    held_out, highdata, g = D.define_splits(df9)
    if args.alleles > 0:
        held_out = held_out[:args.alleles]
    pre = pd.read_csv(os.path.join(CACHE, "precomputed.csv"))

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    done = set()
    if os.path.exists(RESULTS):
        try:
            done = set(pd.read_csv(RESULTS)["allele"].unique())
        except Exception:
            done = set()
    todo = [a for a in held_out if a not in done]
    print(f"{len(held_out)} held-out alleles, {len(done)} done, {len(todo)} to run", flush=True)

    reg, clf = build_tabpfn()
    header_needed = not os.path.exists(RESULTS)
    for ai, allele in enumerate(todo):
        t0 = time.time()
        sub = df9[df9["allele"] == allele].reset_index(drop=True)
        pre_sub = pre[pre["allele"] == allele].reset_index(drop=True)
        # align pre_sub to sub order by peptide
        pre_sub = sub[["peptide"]].merge(pre_sub, on="peptide", how="left")
        rows = run_allele(allele, sub, pre_sub, reg, clf)
        pd.DataFrame(rows)[COLS].to_csv(RESULTS, mode="a", index=False,
                                        header=header_needed)
        header_needed = False
        tp = sum(r["runtime_s"] for r in rows if r["method"].startswith("tabpfn"))
        print(f"[{ai+1}/{len(todo)}] {allele} n={len(sub)} "
              f"-> {len(rows)} rows in {time.time()-t0:.0f}s (tabpfn {tp:.0f}s)",
              flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
