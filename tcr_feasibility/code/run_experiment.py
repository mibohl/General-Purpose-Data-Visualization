"""
Leave-one-epitope-out few-shot protocol.

Task = one epitope. Support = k labeled TCRs for that epitope. Query = held-out
TCRs for the same epitope. Two negative-sampling schemes:
    A) background : benchmark-provided healthy-donor control TCRs
    B) otherep    : TCRs that bind the other epitopes

Anti-leakage: splits are made by high-identity CLUSTER, so no clonal / near-
duplicate cluster spans support and query. For scheme B, other-epitope negatives
within 90% identity of any held-out positive are dropped (public/cross-reactive).

Methods per (epitope, scheme, seed, k):
    tabpfn          TabPFN v2 in-context on the k support TCRs
    tcrdist_nn      TCRdist3 nearest binder distance   (baseline that decides it)
    tcrdist_knn     TCRdist3 mean of k nearest binders
    logreg          logistic regression on support features
    lightgbm        LightGBM on support features
    zeroshot_lr     LR trained on all OTHER epitopes (+epitope features), applied
    zeroshot_lgbm   LightGBM trained on all OTHER epitopes, applied

Query is class-balanced (AUPRC chance = 0.5) and fixed across k within a seed.
Metrics: per-epitope AUROC + AUPRC (never pooled). Runtime logged per method.
"""
import os
import sys
import json
import time
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import lightgbm as lgb

import common
import dedup

CACHE = os.path.join(common.REPO_ROOT, "results", "cache")
RESULTS_CSV = os.path.join(common.REPO_ROOT, "results", "results.csv")

KS = [20, 50, 200]
SEEDS = [0, 1, 2, 3, 4]
IDENT = 0.90
QUERY_POS_FRAC = 0.35
QUERY_POS_MIN, QUERY_POS_MAX = 15, 80
KNN_K = 5
FEATURE_MODE = "beta"   # "beta" | "paired" | "matched" (germline loops, both chains)
ALL_METHODS = ["tabpfn", "kumo", "tcrdist", "logreg", "ridge", "lightgbm", "zeroshot"]
ACTIVE_METHODS = set(ALL_METHODS)   # overridden by --methods
TCRDIST_CHAINS = ["beta"]           # ["alpha","beta"] in matched mode
MATCHED_MASK = None                 # global zero-variance column mask (matched mode)

_FIELDS = ["scheme", "epitope", "seed", "k", "method", "auroc", "auprc",
           "n_support", "n_support_pos", "n_query", "n_query_pos",
           "runtime_s", "note"]


# --------------------------------------------------------------------------- #
# Featurization
# --------------------------------------------------------------------------- #
def tcr_features(df, te_v=None, te_j=None, fit_labels=None):
    """
    TCR features.
      beta   (primary):  100 Atchley(CDR3b) + len + V_te + J_te                = 103
      paired (secondary): + 100 Atchley(CDR3a) + len(CDR3a)                    = 204
    The epitope is constant within a task, so NO epitope features are included.
    Target encoders for V/J are fit on the support labels only (no query leak).
    """
    if FEATURE_MODE == "matched":
        # germline CDR1/CDR2/CDR2.5 + CDR3 for BOTH chains (tcrdist3 loops).
        # No V/J target-encoding: the germline loops already carry exactly the
        # V-gene information the TCRdist metric uses.
        X = common.matched_matrix(df)
        if MATCHED_MASK is not None:
            X = X[:, MATCHED_MASK]
        return X.astype(np.float32), None, None

    atch = common.atchley_matrix(df["cdr3b"].tolist())           # (n,100)
    length = np.array([len(s) for s in df["cdr3b"]], dtype=np.float32)[:, None]
    if fit_labels is not None:
        te_v = common.TargetEncoder().fit(df["vb"].tolist(), fit_labels)
        te_j = common.TargetEncoder().fit(df["jb"].tolist(), fit_labels)
    v = te_v.transform(df["vb"].tolist())[:, None]
    j = te_j.transform(df["jb"].tolist())[:, None]
    blocks = [atch, length, v, j]
    if FEATURE_MODE == "paired":
        atch_a = common.atchley_matrix(df["cdr3a"].tolist())     # (n,100)
        len_a = np.array([len(s) for s in df["cdr3a"]], dtype=np.float32)[:, None]
        blocks += [atch_a, len_a]
    X = np.hstack(blocks).astype(np.float32)
    return X, te_v, te_j


def epitope_features(epitopes, maxlen=11):
    """Atchley on the peptide padded to length 11 + peptide length (56 feats)."""
    rows = []
    for ep in epitopes:
        vec = common.atchley_vector(ep, maxlen=maxlen)  # 55
        rows.append(np.concatenate([vec, [len(ep)]]))
    return np.array(rows, dtype=np.float32)


# --------------------------------------------------------------------------- #
# Split helpers
# --------------------------------------------------------------------------- #
def reserve_by_cluster(df, n_target, rng):
    """
    Shuffle clusters, peel whole clusters into a 'reserved' set until it holds
    >= n_target rows. Returns (reserved_df, remaining_df). No cluster spans both.
    """
    clusters = df["cluster"].unique().tolist()
    rng.shuffle(clusters)
    reserved_idx, got = [], 0
    remaining_clusters = []
    for c in clusters:
        rows = df.index[df["cluster"] == c].tolist()
        if got < n_target:
            reserved_idx.extend(rows)
            got += len(rows)
        else:
            remaining_clusters.append(c)
    reserved = df.loc[reserved_idx]
    remaining = df[df["cluster"].isin(remaining_clusters)]
    return reserved, remaining


# --------------------------------------------------------------------------- #
# Negative pools
# --------------------------------------------------------------------------- #
def build_negative_pool(scheme, held_ep, pos_all, controls):
    """Return a negatives DataFrame with a global 'cluster' key for the scheme."""
    if scheme == "background":
        neg = controls.copy()
        neg["gcluster"] = "ctrl:" + neg["cluster"].astype(str)
        return neg
    # scheme == otherep : other epitopes' positives
    neg = pos_all[pos_all["epitope"] != held_ep].copy()
    neg["gcluster"] = neg["epitope"] + ":" + neg["cluster"].astype(str)
    # drop publics near-identical to held-out positives
    held_seqs = pos_all[pos_all["epitope"] == held_ep]["cdr3b"].tolist()
    mask = dedup.cross_identity_mask(neg["cdr3b"].tolist(), held_seqs, threshold=IDENT)
    neg = neg[~mask].reset_index(drop=True)
    neg.attrs["n_public_removed"] = int(mask.sum())
    return neg


def reserve_neg_by_cluster(df, n_target, rng, col="gcluster"):
    clusters = df[col].unique().tolist()
    rng.shuffle(clusters)
    reserved_idx, got, remaining_clusters = [], 0, []
    for c in clusters:
        rows = df.index[df[col] == c].tolist()
        if got < n_target:
            reserved_idx.extend(rows)
            got += len(rows)
        else:
            remaining_clusters.append(c)
    return df.loc[reserved_idx], df[df[col].isin(remaining_clusters)]


# --------------------------------------------------------------------------- #
# TCRdist
# --------------------------------------------------------------------------- #
_KUMO = {"clf": None}


def kumo_get():
    """Lazily build and cache the Kumo Tabular (small) in-context classifier."""
    if _KUMO["clf"] is None:
        import torch  # noqa
        from sdm.models import KumoTabular
        _KUMO["clf"] = KumoTabular(task="classification", size="small",
                                   device="cpu", pretrained=True)
    return _KUMO["clf"]


def kumo_proba(Xs, ys, Xq):
    """P(class==1) for query rows via Kumo Tabular small, in-context on (Xs,ys)."""
    import torch
    from sdm import CategoricalTensor, StringTensor
    clf = kumo_get()
    ycat = CategoricalTensor(code=torch.tensor(ys).view(-1, 1).long(),
                             categories=(StringTensor.from_list(["0", "1"]),))
    clf.fit(torch.from_numpy(np.ascontiguousarray(Xs)), ycat)
    pred = clf.predict(torch.from_numpy(np.ascontiguousarray(Xq)))
    cols = list(pred.numerical_columns)
    idx = cols.index("1")                       # extract by NAME (order not fixed)
    return pred.numerical[:, idx].detach().float().numpy()


def tcrdist_scores(query_df, support_binder_df):
    """
    Return (nn_score, knn_score) arrays for query rows: higher = more likely a
    binder. Scores are negative distances to the support BINDERS.
    """
    from tcrdist.repertoire import TCRrep
    union = pd.concat([query_df, support_binder_df], ignore_index=True)
    paired = "alpha" in TCRDIST_CHAINS
    cols = {"cdr3_b_aa": union["cdr3b"].astype(str),
            "v_b_gene": union["vb"].astype(str),
            "j_b_gene": union["jb"].astype(str), "count": 1}
    if paired:
        cols.update({"cdr3_a_aa": union["cdr3a"].astype(str),
                     "v_a_gene": union["va"].astype(str),
                     "j_a_gene": union["ja"].astype(str)})
    tr = TCRrep(cell_df=pd.DataFrame(cols), organism="human",
                chains=TCRDIST_CHAINS, compute_distances=True)
    cdf = tr.clone_df
    if paired:
        key = (cdf["cdr3_a_aa"] + "|" + cdf["v_a_gene"] + "|" + cdf["j_a_gene"] + "#"
               + cdf["cdr3_b_aa"] + "|" + cdf["v_b_gene"] + "|" + cdf["j_b_gene"]).tolist()
    else:
        key = (cdf["cdr3_b_aa"] + "|" + cdf["v_b_gene"] + "|" + cdf["j_b_gene"]).tolist()
    idx = {k: i for i, k in enumerate(key)}

    def rowkey(d):
        kb = d["cdr3b"].astype(str) + "|" + d["vb"].astype(str) + "|" + d["jb"].astype(str)
        if paired:
            ka = d["cdr3a"].astype(str) + "|" + d["va"].astype(str) + "|" + d["ja"].astype(str)
            return (ka + "#" + kb).map(idx).values
        return kb.map(idx).values

    q = rowkey(query_df)
    s = rowkey(support_binder_df)
    pw = tr.pw_beta.astype(float) + (tr.pw_alpha.astype(float) if paired else 0.0)
    sub = pw[np.ix_(q, s)]      # query x binders (paired = summed distance)
    sub_sorted = np.sort(sub, axis=1)
    nn = -sub_sorted[:, 0]
    kk = min(KNN_K, sub_sorted.shape[1])
    knn = -sub_sorted[:, :kk].mean(axis=1)
    return nn, knn


# --------------------------------------------------------------------------- #
# Zero-shot models (trained once per held-out epitope x scheme)
# --------------------------------------------------------------------------- #
def train_zeroshot(scheme, held_ep, pos_all, controls, rng):
    """Train LR + LGBM on all OTHER epitopes with epitope features included."""
    other_eps = [e for e in pos_all["epitope"].unique() if e != held_ep]
    frames_pos, frames_neg = [], []
    for ep in other_eps:
        p = pos_all[pos_all["epitope"] == ep]
        frames_pos.append(p)
        npos = len(p)
        if scheme == "background":
            neg = controls.sample(n=min(npos, len(controls)),
                                   random_state=rng.randint(1 << 30)).copy()
            neg["epitope"] = ep
        else:
            pool = pos_all[~pos_all["epitope"].isin([ep, held_ep])]
            neg = pool.sample(n=min(npos, len(pool)),
                              random_state=rng.randint(1 << 30)).copy()
            neg["epitope"] = ep  # epitope feature = the target epitope
        frames_neg.append(neg)
    pos_df = pd.concat(frames_pos, ignore_index=True)
    neg_df = pd.concat(frames_neg, ignore_index=True)
    y = np.r_[np.ones(len(pos_df)), np.zeros(len(neg_df))]
    alldf = pd.concat([pos_df, neg_df], ignore_index=True)
    Xtcr, te_v, te_j = tcr_features(alldf, fit_labels=y)
    Xep = epitope_features(alldf["epitope"].tolist())
    X = np.hstack([Xtcr, Xep])
    lr = make_pipeline(StandardScaler(),
                       LogisticRegression(max_iter=2000, C=1.0))
    lr.fit(X, y)
    gb = lgb.LGBMClassifier(n_estimators=200, num_leaves=31, max_depth=-1,
                            learning_rate=0.05, subsample=0.8,
                            colsample_bytree=0.8, verbose=-1)
    gb.fit(X, y)
    return {"lr": lr, "gb": gb, "te_v": te_v, "te_j": te_j}


def zeroshot_predict(model, query_df, held_ep):
    Xtcr, _, _ = tcr_features(query_df, te_v=model["te_v"], te_j=model["te_j"])
    Xep = epitope_features([held_ep] * len(query_df))
    X = np.hstack([Xtcr, Xep])
    return model["lr"].predict_proba(X)[:, 1], model["gb"].predict_proba(X)[:, 1]


# --------------------------------------------------------------------------- #
# One epitope x scheme
# --------------------------------------------------------------------------- #
def run_epitope_scheme(scheme, held_ep, pos_all, controls, tabpfn_clf,
                       done_keys, writer, logf):
    pos_ep = pos_all[pos_all["epitope"] == held_ep].reset_index(drop=True)
    neg_pool = build_negative_pool(scheme, held_ep, pos_all, controls)
    n_public_removed = neg_pool.attrs.get("n_public_removed", 0)

    # zero-shot models depend on scheme + held_ep only (seed-robust): train once
    zs = None
    if "zeroshot" in ACTIVE_METHODS:
        zs_rng = np.random.RandomState(12345)
        zs = train_zeroshot(scheme, held_ep, pos_all, controls, zs_rng)

    npos = len(pos_ep)
    qp = int(np.clip(round(QUERY_POS_FRAC * npos), QUERY_POS_MIN, QUERY_POS_MAX))

    for seed in SEEDS:
        rng = np.random.RandomState(1000 + seed)
        # reserve balanced query (fixed across k for this seed)
        q_pos, pool_pos = reserve_by_cluster(pos_ep, qp, rng)
        q_pos = q_pos.sample(n=min(qp, len(q_pos)),
                             random_state=seed).reset_index(drop=True)
        q_neg_res, pool_neg = reserve_neg_by_cluster(neg_pool, len(q_pos), rng)
        q_neg = q_neg_res.sample(n=min(len(q_pos), len(q_neg_res)),
                                 random_state=seed).reset_index(drop=True)
        # balance query exactly
        m = min(len(q_pos), len(q_neg))
        q_pos, q_neg = q_pos.iloc[:m], q_neg.iloc[:m]
        query = pd.concat([q_pos, q_neg], ignore_index=True)
        yq = np.r_[np.ones(len(q_pos)), np.zeros(len(q_neg))]
        if len(q_pos) < QUERY_POS_MIN or yq.sum() == 0 or yq.sum() == len(yq):
            continue

        # zero-shot: evaluate once per (epitope,scheme,seed) on this query
        if "zeroshot" in ACTIVE_METHODS:
            for mname, probs in zip(("zeroshot_lr", "zeroshot_lgbm"),
                                    zeroshot_predict(zs, query, held_ep)):
                key = (scheme, held_ep, seed, 0, mname)
                if key in done_keys:
                    continue
                _emit(writer, logf, scheme, held_ep, seed, 0, mname,
                      yq, probs, 0, 0, len(yq), int(yq.sum()), 0.0,
                      f"public_rm={n_public_removed}")

        pool_pos = pool_pos.reset_index(drop=True)
        pool_neg = pool_neg.reset_index(drop=True)

        for k in KS:
            half = k // 2
            if len(pool_pos) < half or len(pool_neg) < half:
                # record a skip note once (method-agnostic) for transparency
                key = (scheme, held_ep, seed, k, "SKIP")
                if key not in done_keys:
                    _emit(writer, logf, scheme, held_ep, seed, k, "SKIP",
                          None, None, k, half, len(yq), int(yq.sum()), 0.0,
                          f"insufficient pool pos={len(pool_pos)} neg={len(pool_neg)}")
                continue
            s_pos = pool_pos.sample(n=half, random_state=seed * 97 + k).reset_index(drop=True)
            s_neg = pool_neg.sample(n=half, random_state=seed * 131 + k).reset_index(drop=True)
            support = pd.concat([s_pos, s_neg], ignore_index=True)
            ys = np.r_[np.ones(half), np.zeros(half)]

            # features (target encoders fit on support only)
            feat_needed = ACTIVE_METHODS & {"tabpfn", "kumo", "logreg", "lightgbm"}
            if feat_needed:
                Xs, te_v, te_j = tcr_features(support, fit_labels=ys)
                Xq, _, _ = tcr_features(query, te_v=te_v, te_j=te_j)

            # ---- TabPFN ----
            key = (scheme, held_ep, seed, k, "tabpfn")
            if "tabpfn" in ACTIVE_METHODS and key not in done_keys:
                t = time.time()
                try:
                    tabpfn_clf.fit(Xs, ys)
                    p = tabpfn_clf.predict_proba(Xq)[:, 1]
                    _emit(writer, logf, scheme, held_ep, seed, k, "tabpfn",
                          yq, p, k, half, len(yq), int(yq.sum()),
                          time.time() - t, "")
                except Exception as e:
                    _emit(writer, logf, scheme, held_ep, seed, k, "tabpfn",
                          None, None, k, half, len(yq), int(yq.sum()),
                          time.time() - t, f"ERR:{type(e).__name__}:{e}")

            # ---- Kumo Tabular (small) ----
            key = (scheme, held_ep, seed, k, "kumo")
            if "kumo" in ACTIVE_METHODS and key not in done_keys:
                t = time.time()
                try:
                    p = kumo_proba(Xs, ys, Xq)
                    _emit(writer, logf, scheme, held_ep, seed, k, "kumo",
                          yq, p, k, half, len(yq), int(yq.sum()),
                          time.time() - t, "")
                except Exception as e:
                    _emit(writer, logf, scheme, held_ep, seed, k, "kumo",
                          None, None, k, half, len(yq), int(yq.sum()),
                          time.time() - t, f"ERR:{type(e).__name__}:{e}")

            # ---- Logistic regression ----
            key = (scheme, held_ep, seed, k, "logreg")
            if "logreg" in ACTIVE_METHODS and key not in done_keys:
                t = time.time()
                lr = make_pipeline(StandardScaler(),
                                   LogisticRegression(max_iter=2000, C=1.0))
                lr.fit(Xs, ys)
                p = lr.predict_proba(Xq)[:, 1]
                _emit(writer, logf, scheme, held_ep, seed, k, "logreg",
                      yq, p, k, half, len(yq), int(yq.sum()), time.time() - t, "")

            # ---- Ridge classifier ----
            key = (scheme, held_ep, seed, k, "ridge")
            if "ridge" in ACTIVE_METHODS and key not in done_keys:
                t = time.time()
                rc = make_pipeline(StandardScaler(), RidgeClassifier(alpha=1.0))
                rc.fit(Xs, ys)
                p = rc.decision_function(Xq)      # score for AUROC
                _emit(writer, logf, scheme, held_ep, seed, k, "ridge",
                      yq, p, k, half, len(yq), int(yq.sum()), time.time() - t, "")

            # ---- LightGBM ----
            key = (scheme, held_ep, seed, k, "lightgbm")
            if "lightgbm" in ACTIVE_METHODS and key not in done_keys:
                t = time.time()
                gb = lgb.LGBMClassifier(n_estimators=200, num_leaves=15,
                                        max_depth=-1, learning_rate=0.05,
                                        subsample=0.8, colsample_bytree=0.8,
                                        min_child_samples=5, verbose=-1)
                gb.fit(Xs, ys)
                p = gb.predict_proba(Xq)[:, 1]
                _emit(writer, logf, scheme, held_ep, seed, k, "lightgbm",
                      yq, p, k, half, len(yq), int(yq.sum()), time.time() - t, "")

            # ---- TCRdist NN / kNN (scored vs support binders) ----
            need_nn = (scheme, held_ep, seed, k, "tcrdist_nn") not in done_keys
            need_knn = (scheme, held_ep, seed, k, "tcrdist_knn") not in done_keys
            if "tcrdist" in ACTIVE_METHODS and (need_nn or need_knn):
                t = time.time()
                try:
                    nn, knn = tcrdist_scores(query, s_pos)
                    dt = time.time() - t
                    if need_nn:
                        _emit(writer, logf, scheme, held_ep, seed, k, "tcrdist_nn",
                              yq, nn, k, half, len(yq), int(yq.sum()), dt, "")
                    if need_knn:
                        _emit(writer, logf, scheme, held_ep, seed, k, "tcrdist_knn",
                              yq, knn, k, half, len(yq), int(yq.sum()), dt, "")
                except Exception as e:
                    if need_nn:
                        _emit(writer, logf, scheme, held_ep, seed, k, "tcrdist_nn",
                              None, None, k, half, len(yq), int(yq.sum()),
                              time.time() - t, f"ERR:{type(e).__name__}:{e}")
            logf.flush()


def _emit(writer, logf, scheme, ep, seed, k, method, yq, probs,
          nsup, nsup_pos, nq, nq_pos, rt, note):
    if probs is None or yq is None:
        auroc = auprc = np.nan
    else:
        try:
            auroc = roc_auc_score(yq, probs)
            auprc = average_precision_score(yq, probs)
        except Exception:
            auroc = auprc = np.nan
    row = dict(scheme=scheme, epitope=ep, seed=seed, k=k, method=method,
               auroc=auroc, auprc=auprc, n_support=nsup, n_support_pos=nsup_pos,
               n_query=nq, n_query_pos=nq_pos, runtime_s=round(rt, 3), note=note)
    writer.writerow([row[c] for c in _FIELDS])
    logf.write(f"{scheme:10s} {ep:11s} s{seed} k={k:3d} {method:14s} "
               f"AUROC={auroc if auroc==auroc else float('nan'):.3f} "
               f"AUPRC={auprc if auprc==auprc else float('nan'):.3f} "
               f"({rt:.1f}s) {note}\n")


def load_done_keys():
    done = set()
    if os.path.exists(RESULTS_CSV):
        try:
            prev = pd.read_csv(RESULTS_CSV)
            for _, r in prev.iterrows():
                done.add((r["scheme"], r["epitope"], int(r["seed"]),
                          int(r["k"]), r["method"]))
        except Exception:
            pass
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schemes", nargs="+", default=["background", "otherep"])
    ap.add_argument("--epitopes", nargs="+", default=None)
    ap.add_argument("--fresh", action="store_true", help="ignore existing results.csv")
    ap.add_argument("--feature-mode", choices=["beta", "paired", "matched"], default="beta")
    ap.add_argument("--out", default=None, help="results csv path override")
    ap.add_argument("--seeds", type=int, default=None, help="limit to first N seeds")
    ap.add_argument("--methods", nargs="+", default=None,
                    help=f"subset of {ALL_METHODS}; default all")
    args = ap.parse_args()

    global FEATURE_MODE, RESULTS_CSV, SEEDS, ACTIVE_METHODS, TCRDIST_CHAINS, MATCHED_MASK
    FEATURE_MODE = args.feature_mode
    if args.out:
        RESULTS_CSV = args.out
    if args.seeds:
        SEEDS = SEEDS[:args.seeds]
    if args.methods:
        ACTIVE_METHODS = set(args.methods)

    import csv
    if FEATURE_MODE == "matched":
        pos_all = pd.read_parquet(os.path.join(CACHE, "positives_matched.parquet"))
        controls = pd.read_parquet(os.path.join(CACHE, "controls_matched.parquet"))
        MATCHED_MASK = np.load(os.path.join(CACHE, "matched_colmask.npy"))
        TCRDIST_CHAINS = ["alpha", "beta"]   # match the tabular models' inputs
        print(f"[matched] dim={int(MATCHED_MASK.sum())}, TCRdist chains={TCRDIST_CHAINS}")
    else:
        pos_all = pd.read_parquet(os.path.join(CACHE, "positives.parquet"))
        controls = pd.read_parquet(os.path.join(CACHE, "controls.parquet"))
    if FEATURE_MODE == "paired":
        # keep only TCRs with a usable CDR3-alpha (paired-chain subset)
        pa = pos_all["cdr3a"].map(common.valid_cdr3)
        ca = controls["cdr3a"].map(common.valid_cdr3)
        print(f"[paired] positives with valid CDR3a: {int(pa.sum())}/{len(pos_all)}; "
              f"controls: {int(ca.sum())}/{len(controls)}")
        pos_all = pos_all[pa].reset_index(drop=True)
        controls = controls[ca].reset_index(drop=True)
    epitopes = args.epitopes or sorted(pos_all["epitope"].unique().tolist())

    done = set() if args.fresh else load_done_keys()
    newfile = args.fresh or not os.path.exists(RESULTS_CSV)
    mode = "w" if newfile else "a"

    tabpfn_clf = None
    if "tabpfn" in ACTIVE_METHODS:
        from tabpfn import TabPFNClassifier
        tabpfn_clf = TabPFNClassifier(device="cpu")
    if "kumo" in ACTIVE_METHODS:
        kumo_get()  # warm the model once up front

    f = open(RESULTS_CSV, mode, newline="")
    writer = csv.writer(f)
    if newfile:
        writer.writerow(_FIELDS)
    logf = sys.stdout

    t0 = time.time()
    for scheme in args.schemes:
        for ep in epitopes:
            tep = time.time()
            logf.write(f"\n===== scheme={scheme} epitope={ep} =====\n")
            logf.flush()
            run_epitope_scheme(scheme, ep, pos_all, controls, tabpfn_clf,
                               done, writer, logf)
            f.flush()
            logf.write(f"  [{ep} {scheme} done in {time.time()-tep:.0f}s, "
                       f"total {time.time()-t0:.0f}s]\n")
            logf.flush()
    f.close()
    logf.write(f"\nALL DONE in {time.time()-t0:.0f}s -> {RESULTS_CSV}\n")


if __name__ == "__main__":
    main()
