"""Step 4: prediction-model adaptation arms per held-out study.
Usage: python -I 07_arms.py <work_dir> <out_csv> <scripts_dir> <tabpfn_ckpt>
Query cap 1500 per (study, seed) and TabPFN n_estimators=4: CPU budget, fixed before any Step 4 result.
"""
import sys, os, time, warnings
W, OUTCSV, SD, CKPT = sys.argv[1:5]
sys.path.insert(0, SD)
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import KFold
import lightgbm as lgb
from tabpfn import TabPFNRegressor
from feats import peptide_features
sv = pd.read_parquet(f'{W}/study_values.parquet').merge(pd.read_parquet(f'{W}/im2deep_pred.parquet'), on='precursor')
sv = sv.reset_index(drop=True)
F = peptide_features(sv); F['z1'] = sv['k0_pred'].values
cols = list(F.columns)
KS, SEEDS, QCAP = (20, 50, 200, 1000), range(5), 1500
def gbm(k):
    return lgb.LGBMRegressor(n_estimators=200, learning_rate=0.05, num_leaves=15, min_child_samples=max(5, k // 20),
                             subsample=0.8, subsample_freq=1, colsample_bytree=0.8, verbose=-1, n_jobs=4)
done = set()
if os.path.exists(OUTCSV):
    prev = pd.read_csv(OUTCSV); done = set(zip(prev.study, prev.k, prev.seed))
first = not os.path.exists(OUTCSV)
for s, idx in sorted(sv.groupby('study').indices.items(), key=lambda t: len(t[1])):
    n = len(idx)
    if n < 300:
        print('skip', s, n); continue
    y = sv.k0.values[idx]; z1 = sv.k0_pred.values[idx]; X = F.loc[idx, cols].values
    # O1 oracle: cross-fitted LightGBM residual model on all of this study's precursors
    o1 = np.empty(n)
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        o1[te] = z1[te] + gbm(len(tr)).fit(X[tr], y[tr] - z1[tr]).predict(X[te])
    b = np.linalg.lstsq(np.c_[np.ones(n), z1], y, rcond=None)[0]; o_aff = b[0] + b[1] * z1
    for k in KS:
        if n < k + 200: continue
        for seed in SEEDS:
            if (s, k, seed) in done: continue
            t0 = time.time()
            rng = np.random.default_rng(10_000 + 100 * seed + k)
            perm = rng.permutation(n); a = perm[:k]; q = perm[k:k + QCAP]
            preds = {'Z1': z1[q]}
            preds['C1'] = z1[q] + np.median(y[a] - z1[a])
            ba = np.linalg.lstsq(np.c_[np.ones(k), z1[a]], y[a], rcond=None)[0]; preds['C2'] = ba[0] + ba[1] * z1[q]
            preds['C3'] = IsotonicRegression(out_of_bounds='clip').fit(z1[a], y[a]).predict(z1[q])
            rid = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 13))).fit(X[a], y[a] - z1[a])
            preds['R1_ridge'] = z1[q] + rid.predict(X[q])
            preds['R1_lgbm'] = z1[q] + gbm(k).fit(X[a], y[a] - z1[a]).predict(X[q])
            tp = TabPFNRegressor(model_path=CKPT, device='cpu', n_estimators=4, random_state=seed, ignore_pretraining_limits=True)
            preds['T1'] = z1[q] + tp.fit(X[a], y[a] - z1[a]).predict(X[q])
            preds['O1'] = o1[q]; preds['O1_affine_all'] = o_aff[q]
            rows = []
            for arm, p in preds.items():
                e = np.abs(p - y[q])
                rows.append(dict(study=s, k=k, seed=seed, arm=arm, n_query=len(q), n_study=n, mae_median=np.median(e),
                                 q25=np.quantile(e, .25), q75=np.quantile(e, .75), mean_abs=e.mean()))
            pd.DataFrame(rows).to_csv(OUTCSV, mode='a', header=first, index=False); first = False
            print(s, k, seed, {r['arm']: round(r['mae_median'], 4) for r in rows}, f'{time.time() - t0:.0f}s', flush=True)
print('ALL DONE')
