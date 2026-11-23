"""Pooled ROI proxy (spec fallback): LightGBM on RDKit descriptors + ECFP4 -> within-dataset RT percentile rank.
Held-out datasets get predictions from a model trained on all non-held-out datasets of the same mode.
Pool datasets get 5-fold out-of-fold (grouped by dataset) predictions, used only as a Z1 feature."""
import time, numpy as np, pandas as pd, lightgbm as lgb
from sklearn.model_selection import GroupKFold
OUT = '/home/user/work/out'
comp = pd.read_parquet(f'{OUT}/compounds.parquet')
meta = pd.read_csv(f'{OUT}/datasets_meta.csv', dtype={'id': str}).set_index('id')
D = pd.read_parquet(f'{OUT}/descriptors.parquet')
E = np.load(f'{OUT}/ecfp4.npy')
row = pd.Series(np.arange(len(D)), index=D.index)
comp['y'] = comp.groupby('dataset').rt.rank(pct=True)
P = dict(n_estimators=600, learning_rate=0.05, num_leaves=63, min_child_samples=20, subsample=0.8,
         subsample_freq=1, colsample_bytree=0.5, verbose=-1, n_jobs=2)
def X(sub):
    r = row.loc[sub.smiles].values
    return np.hstack([D.values[r], E[r]]).astype(np.float32)
out, timing = [], []
for mode in ['RP', 'HILIC']:
    ds = meta.index[meta['mode'] == mode]
    held = set(meta.index[meta.heldout_any & (meta['mode'] == mode)])
    pool = comp[comp.dataset.isin(ds) & ~comp.dataset.isin(held)]
    test = comp[comp.dataset.isin(held)]
    t = time.time()
    m = lgb.LGBMRegressor(**P).fit(X(pool), pool.y)
    tt = time.time() - t
    out.append(test[['dataset', 'cid']].assign(roi=m.predict(X(test))))
    oof = np.zeros(len(pool))
    for tr, te in GroupKFold(5).split(pool, groups=pool.dataset):
        oof[te] = lgb.LGBMRegressor(**P).fit(X(pool.iloc[tr]), pool.y.iloc[tr]).predict(X(pool.iloc[te]))
    out.append(pool[['dataset', 'cid']].assign(roi=oof))
    timing.append(dict(mode=mode, n_train=len(pool), n_pool_datasets=pool.dataset.nunique(), fit_seconds=tt))
    print(mode, 'pool rows', len(pool), 'datasets', pool.dataset.nunique(), 'fit', round(tt, 1), 's', flush=True)
pd.concat(out).to_parquet(f'{OUT}/roi_proxy.parquet')
pd.DataFrame(timing).to_csv(f'{OUT}/timing_roi_proxy.csv', index=False)
