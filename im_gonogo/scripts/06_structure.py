"""Step 5: is the post-affine residual peptide-dependent, and is the dependence study-specific?
Usage: python -I 06_structure.py <work_dir> <out_dir> <scripts_dir>
"""
import sys
W, OUT, SD = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, SD)
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score
import lightgbm as lgb
from feats import peptide_features
sh = pd.read_parquet(f'{W}/shared_with_consensus.parquet')
F = peptide_features(sh)
cols = list(F.columns)
# residual after per-study affine fit on all shared precursors (in-sample, generous to affine)
sh['resid'] = np.nan
for s, idx in sh.groupby('study').indices.items():
    x = sh.k0.values[idx]; c = sh.cons_k0.values[idx]
    b = np.linalg.lstsq(np.c_[np.ones_like(x), x], c, rcond=None)[0]
    sh.loc[sh.index[idx], 'resid'] = c - (b[0] + b[1] * x)
ridge = lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 13)))
gbm = lambda: lgb.LGBMRegressor(n_estimators=300, learning_rate=0.05, num_leaves=15, min_child_samples=30, subsample=0.8,
                                subsample_freq=1, colsample_bytree=0.8, verbose=-1, n_jobs=4)
rows, coefs, loso_pred = [], {}, {}
for subset in ['all', 'no_multi']:
    d = sh if subset == 'all' else sh[~sh.multi]
    for s, g in d.groupby('study'):
        if len(g) < 300: continue
        X = F.loc[g.index, cols].values; y = g.resid.values
        kf = KFold(5, shuffle=True, random_state=0)
        p_r = cross_val_predict(ridge(), X, y, cv=kf); p_g = cross_val_predict(gbm(), X, y, cv=kf)
        # leave-one-study-out: model trained on all other studies' residuals
        o = d[d.study != s]
        m_lo = gbm().fit(F.loc[o.index, cols].values, o.resid.values); p_lo = m_lo.predict(X)
        r_lo = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 13))).fit(F.loc[o.index, cols].values, o.resid.values).predict(X)
        rows.append(dict(study=s, subset=subset, n=len(g), resid_mad=np.median(np.abs(y)),
                         r2_ridge_cv=r2_score(y, p_r), r2_lgbm_cv=r2_score(y, p_g),
                         r2_lgbm_loso_pooled=r2_score(y, p_lo), r2_ridge_loso_pooled=r2_score(y, r_lo),
                         mad_after_lgbm_cv=np.median(np.abs(y - p_g)), mad_after_loso=np.median(np.abs(y - p_lo))))
        if subset == 'all':
            m = ridge().fit(X, y); coefs[s] = pd.Series(m[-1].coef_, index=cols)
            # per-charge mean residual (interpretable pattern)
            for z in (2, 3, 4):
                mz = g.charge.values == z
                rows[-1][f'mean_resid_z{z}'] = y[mz].mean() if mz.sum() > 20 else np.nan
        print(subset, s, {k: round(v, 3) for k, v in rows[-1].items() if isinstance(v, float)}, flush=True)
res = pd.DataFrame(rows); res.to_csv(f'{OUT}/step5_variance_explained.csv', index=False)
C = pd.DataFrame(coefs); C.to_csv(f'{OUT}/step5_ridge_coefs.csv')
cc = C.corr().values[np.triu_indices(C.shape[1], 1)]
print('median pairwise coef correlation across studies', np.median(cc), 'IQR', np.quantile(cc, [.25, .75]))
pd.Series(dict(median_pairwise_coef_corr=np.median(cc), q25=np.quantile(cc, .25), q75=np.quantile(cc, .75))).to_csv(f'{OUT}/step5_coef_similarity.csv')
sh[['study', 'precursor', 'resid']].to_parquet(f'{W}/affine_resid.parquet')
