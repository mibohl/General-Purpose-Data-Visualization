"""Step 1b + Step 2: study-level values, shared-precursor matrix, replicate floor F.
Usage: python -I 03_matrix_floor.py <work_dir> <out_dir>
"""
import sys
import numpy as np, pandas as pd
W, OUT = sys.argv[1], sys.argv[2]
psm = pd.read_parquet(f'{W}/psm_filtered.parquet')
R = pd.read_csv(f'{W}/q1_run_index_map.csv')
R['study'] = R['study'].replace({'PXD003206': 'PXD036206'})
# TIMS method groups: runs of one study share a group unless the scan->1/K0 map jumps
grp = {}
for s, g in R.groupby('study'):
    g = g.sort_values(['k0_at_idx500', 'slope'])
    gid, prev = 0, None
    for _, r in g.iterrows():
        if prev is not None and (abs(r.k0_at_idx500 - prev.k0_at_idx500) > 0.01 or abs(r.slope - prev.slope) > 3e-5):
            gid += 1
        grp[r.run] = f'{s}|m{gid}'; prev = r
psm['mgroup'] = psm['Raw file'].map(grp)
psm = psm[psm['mgroup'].notna()]   # runs with < 200 MS/MS rows had no map; they are dropped
psm.groupby(['study', 'mgroup'])['Raw file'].nunique().to_csv(f'{OUT}/method_groups.csv')

# ---- replicate floor F: leave-one-run-out deviation within method group ----
def loo_dev(sub, col):
    sub = sub[['mgroup', 'precursor', 'Raw file', col]]
    n = sub.groupby(['mgroup', 'precursor'])[col].transform('size')
    sub = sub[n >= 3]
    out = []
    for (mg, p), g in sub.groupby(['mgroup', 'precursor'], sort=False):
        v = g[col].values
        for i in range(len(v)):
            out.append(abs(v[i] - np.median(np.delete(v, i))))
    return np.array(out)
rows = []
for s, g in psm.groupby('study'):
    for col in ['1/K0', 'CCS']:
        d = loo_dev(g, col)
        rel = None
        if len(d):
            # relative version uses the precursor's own value as scale
            rel = np.median(d / g[col].median())
        rows.append(dict(study=s, space=col, n_dev=len(d), F=np.median(d) if len(d) else np.nan,
                         F_iqr_lo=np.quantile(d, .25) if len(d) else np.nan, F_iqr_hi=np.quantile(d, .75) if len(d) else np.nan,
                         n_runs=g['Raw file'].nunique()))
    print(s, rows[-2]['F'], rows[-2]['n_dev'], flush=True)
Ft = pd.DataFrame(rows)
Ft.to_csv(f'{OUT}/F_per_study.csv', index=False)

# pooled F: median over all deviations (recompute pooled)
pooled = {}
for col in ['1/K0', 'CCS']:
    d = loo_dev(psm, col); pooled[col] = (np.median(d), len(d), np.quantile(d, .25), np.quantile(d, .75))
pd.DataFrame(pooled, index=['F_pooled', 'n', 'q25', 'q75']).to_csv(f'{OUT}/F_pooled.csv')
print('pooled', pooled)

# ---- study-level values ----
sv = psm.groupby(['study', 'precursor']).agg(k0=('1/K0', 'median'), ccs=('CCS', 'median'), mz=('m/z', 'median'),
        charge=('Charge', 'first'), seq=('Sequence', 'first'), modseq=('modseq', 'first'),
        n_runs=('Raw file', 'nunique'), k0_min_run=('k0_min', 'min'), k0_max_run=('k0_max', 'max'),
        max_within_run_spread=('k0_max', lambda x: np.nan)).reset_index()
# within-run multi-feature spread (relative), max over runs
psm['wr_spread'] = (psm['k0_max'] - psm['k0_min']) / psm['k0_min']
sv = sv.drop(columns='max_within_run_spread').merge(
    psm.groupby(['study', 'precursor'])['wr_spread'].max().rename('max_within_run_spread').reset_index(), on=['study', 'precursor'])
sv['n_studies'] = sv.groupby('precursor')['study'].transform('size')
sv.to_parquet(f'{W}/study_values.parquet')
psm.to_parquet(f'{W}/psm_grouped.parquet')
sh = sv[sv.n_studies >= 4]
print('studies', sv.study.nunique(), 'precursors', sv.precursor.nunique(), 'shared>=4', sh.precursor.nunique())
