"""Step 3 (+Step 6 flags): cross-study drift D, affine residual A (20/50 anchors), oracle O.
Usage: python -I 04_drift.py <work_dir> <out_dir>
"""
import sys
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import KFold
W, OUT = sys.argv[1], sys.argv[2]
sv = pd.read_parquet(f'{W}/study_values.parquet')
sh = sv[sv.n_studies >= 4].copy()

# ---- Step 6 flags: multi-conformer candidates ----
within = sh.groupby('precursor')['max_within_run_spread'].max() > 0.02
def bimodal(v):
    v = np.sort(v)
    if len(v) < 4: return False
    gaps = np.diff(v) / v[:-1]
    for i in range(1, len(v) - 2):        # at least 2 studies on each side
        if gaps[i] > 0.02: return True
    return False
bim = sh.groupby('precursor')['k0'].apply(lambda s: bimodal(s.values))
flag = (within | bim).rename('multi')
sh = sh.merge(flag.reset_index(), on='precursor')
pd.DataFrame({'within_run_split': within, 'across_study_bimodal': bim}).to_csv(f'{OUT}/multiconformer_flags.csv')
print('multi flags: within', within.sum(), 'bimodal', bim.sum(), 'either', flag.sum(), 'of', len(flag), flush=True)

# ---- leave-one-study-out consensus ----
def loo_consensus(df, col):
    out = np.empty(len(df))
    for p, idx in df.groupby('precursor').indices.items():
        v = df[col].values[idx]
        for j, i in enumerate(idx):
            out[i] = np.median(np.delete(v, j))
    return out
sh = sh.reset_index(drop=True)
for col in ['k0', 'ccs']:
    sh[f'cons_{col}'] = loo_consensus(sh, col)
sh.to_parquet(f'{W}/shared_with_consensus.parquet')

def fit_affine(x, y):
    X = np.c_[np.ones_like(x), x]
    return np.linalg.lstsq(X, y, rcond=None)[0]

rows = []      # every (study, space, method, seed) -> stats
devs = {}      # pooled deviations
SEEDS = range(5)
for space in ['k0', 'ccs']:
    for subset in ['all', 'no_multi']:
        d = sh if subset == 'all' else sh[~sh.multi]
        for s, g in d.groupby('study'):
            x = g[space].values; c = g[f'cons_{space}'].values; n = len(g)
            if n < 60: continue
            def rec(method, seed, dev, rel):
                rows.append(dict(study=s, space=space, subset=subset, method=method, seed=seed, n_eval=len(dev),
                                 med=np.median(dev), q25=np.quantile(dev, .25), q75=np.quantile(dev, .75),
                                 med_rel=np.median(rel)))
                devs.setdefault((space, subset, method, seed), []).append(dev)
            rec('raw', -1, np.abs(x - c), np.abs(x - c) / c)
            # oracle affine on all (in-sample), oracle isotonic in-sample, cross-fitted, per charge
            b = fit_affine(x, c); r = np.abs(b[0] + b[1] * x - c); rec('affine_all_insample', -1, r, r / c)
            iso = IsotonicRegression(out_of_bounds='clip').fit(x, c); r = np.abs(iso.predict(x) - c); rec('oracle_isotonic', -1, r, r / c)
            pred = np.empty(n)
            for tr, te in KFold(5, shuffle=True, random_state=0).split(x):
                pred[te] = IsotonicRegression(out_of_bounds='clip').fit(x[tr], c[tr]).predict(x[te])
            r = np.abs(pred - c); rec('oracle_isotonic_cv5', -1, r, r / c)
            pred = np.empty(n); z = g['charge'].values
            for zz in np.unique(z):
                m = z == zz
                pred[m] = IsotonicRegression(out_of_bounds='clip').fit(x[m], c[m]).predict(x[m]) if m.sum() >= 5 else x[m] + np.median(c[m] - x[m])
            r = np.abs(pred - c); rec('oracle_isotonic_per_charge', -1, r, r / c)
            for k in (20, 50):
                for seed in SEEDS:
                    rng = np.random.default_rng(1000 * seed + k)
                    a = rng.choice(n, k, replace=False); q = np.setdiff1d(np.arange(n), a)
                    b = fit_affine(x[a], c[a]); r = np.abs(b[0] + b[1] * x[q] - c[q]); rec(f'affine_{k}', seed, r, r / c[q])
                    sft = np.median(c[a] - x[a]); r = np.abs(x[q] + sft - c[q]); rec(f'shift_{k}', seed, r, r / c[q])
                    r = np.abs(x[q] - c[q]); rec(f'raw_on_query_{k}', seed, r, r / c[q])
        print(space, subset, 'done', flush=True)
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/step3_per_study.csv', index=False)
# pooled: median of all deviations per seed, then median across seeds
pool = []
for (space, subset, method, seed), L in devs.items():
    v = np.concatenate(L)
    pool.append(dict(space=space, subset=subset, method=method, seed=seed, n=len(v), med=np.median(v),
                     q25=np.quantile(v, .25), q75=np.quantile(v, .75)))
pool = pd.DataFrame(pool)
summ = pool.groupby(['space', 'subset', 'method']).agg(med=('med', 'median'), med_seed_min=('med', 'min'), med_seed_max=('med', 'max'),
                                                         q25=('q25', 'median'), q75=('q75', 'median'), n=('n', 'median')).reset_index()
summ.to_csv(f'{OUT}/step3_pooled.csv', index=False)
pd.set_option('display.width', 200)
print(summ.to_string())
print('shared precursors', sh.precursor.nunique(), 'rows', len(sh), 'studies', sh.study.nunique())
