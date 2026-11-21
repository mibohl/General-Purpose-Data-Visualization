"""
Q2 (headroom) + Q3 (meta-learnability), matched germline+CDR3 features, both
chains. Everything per-epitope, never pooled in the reporting.

Per epitope e, per negative scheme:
  oracle        : best within-epitope model (LightGBM, ridge) under cluster-aware
                  5-fold CV on ALL of e's labels -> ceiling of within-epitope AUROC.
  zeroshot      : train on all OTHER epitopes WITH epitope features, 0 labels from
                  e, predict e's balanced set  (Q2).
  pooled_noepi  : same transfer but with NO epitope conditioning (Q3).
  headroom      : oracle - zeroshot.
Plus per-epitope ridge coefficient vectors -> pairwise cosine similarity (Q3:
if per-task functions are similar, pooled~=oracle and coefficients align).

Outputs: results_headroom.csv, coef_similarity.csv
"""
import os, warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
import lightgbm as lgb
import common, dedup

CACHE = os.path.join(common.REPO_ROOT, "results", "cache")
RES = os.path.join(common.REPO_ROOT, "results")
IDENT = 0.90
SEED = 0
SCHEMES = ["background", "otherep"]

POS = pd.read_parquet(os.path.join(CACHE, "positives_matched.parquet"))
CTRL = pd.read_parquet(os.path.join(CACHE, "controls_matched.parquet"))
MASK = np.load(os.path.join(CACHE, "matched_colmask.npy"))
EPIS = sorted(POS.epitope.unique())


def feats(df):
    return common.matched_matrix(df)[:, MASK].astype(np.float32)


def epi_feats(peptides):
    return np.array([np.concatenate([common.atchley_vector(p, 11), [len(p)]])
                     for p in peptides], dtype=np.float32)


def neg_pool(scheme, held_ep, rng, n):
    """Balanced negatives for one epitope, matched features, cluster ids."""
    if scheme == "background":
        pool = CTRL.copy(); pool["grp"] = "c" + pool["cluster"].astype(str)
    else:
        pool = POS[POS.epitope != held_ep].copy()
        pool["grp"] = pool["epitope"] + ":" + pool["cluster"].astype(str)
        held = POS[POS.epitope == held_ep]["cdr3b"].tolist()
        m = dedup.cross_identity_mask(pool["cdr3b"].tolist(), held, threshold=IDENT)
        pool = pool[~m]
    return pool.sample(n=min(n, len(pool)), random_state=rng).reset_index(drop=True)


def oracle_cv(X, y, groups):
    """Cluster-aware 5-fold CV AUROC for ridge and LightGBM (out-of-fold)."""
    gkf = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    oof = {"ridge": np.full(len(y), np.nan), "lightgbm": np.full(len(y), np.nan)}
    for tr, te in gkf.split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        rc = make_pipeline(StandardScaler(), RidgeClassifier(alpha=1.0)).fit(X[tr], y[tr])
        oof["ridge"][te] = rc.decision_function(X[te])
        gb = lgb.LGBMClassifier(n_estimators=300, num_leaves=15, learning_rate=0.05,
                                subsample=0.8, colsample_bytree=0.8,
                                min_child_samples=5, verbose=-1).fit(X[tr], y[tr])
        oof["lightgbm"][te] = gb.predict_proba(X[te])[:, 1]
    out = {}
    for m, s in oof.items():
        ok = ~np.isnan(s)
        out[m] = roc_auc_score(y[ok], s[ok]) if len(np.unique(y[ok])) == 2 else np.nan
    return out


def build_epi(scheme, ep, rng):
    P = POS[POS.epitope == ep].reset_index(drop=True)
    N = neg_pool(scheme, ep, rng, len(P))
    X = np.vstack([feats(P), feats(N)])
    y = np.r_[np.ones(len(P)), np.zeros(len(N))].astype(int)
    groups = np.r_[P["cluster"].values, np.arange(len(N)) + 100000]  # neg: own groups
    peptides = [ep] * (len(P) + len(N))
    return X, y, groups, np.array(peptides)


def main():
    rows, coef_vecs = [], {}
    for scheme in SCHEMES:
        # pre-build each epitope's matched set once (for transfer training pools)
        built = {ep: build_epi(scheme, ep, SEED + i) for i, ep in enumerate(EPIS)}
        for ep in EPIS:
            Xte, yte, _, _ = built[ep]
            # ---- oracle: within-epitope cluster CV ----
            Xo, yo, go, _ = built[ep]
            orc = oracle_cv(Xo, yo, go)
            oracle = np.nanmax([orc["ridge"], orc["lightgbm"]])
            # ---- transfer: train on all OTHER epitopes ----
            Xtr = np.vstack([built[e][0] for e in EPIS if e != ep])
            ytr = np.concatenate([built[e][1] for e in EPIS if e != ep])
            peps = np.concatenate([built[e][3] for e in EPIS if e != ep])
            # zero-shot WITH epitope features
            Xtr_e = np.hstack([Xtr, epi_feats(peps)])
            Xte_e = np.hstack([Xte, epi_feats([ep] * len(yte))])
            gb = lgb.LGBMClassifier(n_estimators=300, num_leaves=31, learning_rate=0.05,
                                    subsample=0.8, colsample_bytree=0.8, verbose=-1).fit(Xtr_e, ytr)
            zs = roc_auc_score(yte, gb.predict_proba(Xte_e)[:, 1])
            # pooled, NO epitope conditioning
            gb2 = lgb.LGBMClassifier(n_estimators=300, num_leaves=31, learning_rate=0.05,
                                     subsample=0.8, colsample_bytree=0.8, verbose=-1).fit(Xtr, ytr)
            pooled = roc_auc_score(yte, gb2.predict_proba(Xte)[:, 1])
            rows.append(dict(scheme=scheme, epitope=ep, n_pos=int(yte.sum()),
                             oracle_ridge=orc["ridge"], oracle_lgbm=orc["lightgbm"],
                             oracle=oracle, zeroshot=zs, pooled_noepi=pooled,
                             headroom=oracle - zs))
            # ---- per-epitope ridge coefficients (Q3) ----
            if scheme == "otherep":
                rc = make_pipeline(StandardScaler(), RidgeClassifier(alpha=1.0)).fit(Xo, yo)
                coef_vecs[ep] = rc.named_steps["ridgeclassifier"].coef_.ravel()
            print(f"{scheme:10s} {ep:11s} oracle={oracle:.3f} zeroshot={zs:.3f} "
                  f"pooled={pooled:.3f} headroom={oracle-zs:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "results_headroom.csv"), index=False)

    # ---- Q3: coefficient similarity across epitopes ----
    eps = list(coef_vecs)
    M = np.vstack([coef_vecs[e] for e in eps])
    Mn = M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)
    sim = Mn @ Mn.T
    simdf = pd.DataFrame(sim, index=eps, columns=eps)
    simdf.to_csv(os.path.join(RES, "coef_similarity.csv"))
    off = sim[np.triu_indices(len(eps), 1)]
    print(f"\n[Q3] per-epitope ridge coefficient cosine similarity (other-epitope scheme):")
    print(f"     median={np.median(off):.3f}  mean={off.mean():.3f}  "
          f"min={off.min():.3f}  max={off.max():.3f}")
    summ = {"coef_cosine_median": float(np.median(off)),
            "coef_cosine_mean": float(off.mean()),
            "coef_cosine_min": float(off.min()), "coef_cosine_max": float(off.max())}
    json.dump(summ, open(os.path.join(RES, "coef_similarity_summary.json"), "w"), indent=2)
    print("\nwrote results_headroom.csv, coef_similarity.csv")


if __name__ == "__main__":
    main()
