"""F variants on the same populations used for D/A/O, for sensitivity of criterion 1.
Usage: python -I 09_floor_variants.py <work_dir> <out_dir>
"""
import sys
import numpy as np, pandas as pd
W, OUT = sys.argv[1], sys.argv[2]
psm = pd.read_parquet(f'{W}/psm_grouped.parquet', columns=['study', 'mgroup', 'Raw file', 'precursor', '1/K0', 'CCS', 'wr_spread'])
sh = pd.read_parquet(f'{W}/shared_with_consensus.parquet', columns=['precursor', 'multi']).drop_duplicates('precursor')
shared = set(sh.precursor); multi = set(sh.precursor[sh.multi])
psm = psm[psm.groupby(['mgroup', 'precursor'])['1/K0'].transform('size') >= 3].reset_index(drop=True)
for col in ['1/K0', 'CCS']:
    dev = np.empty(len(psm))
    for _, idx in psm.groupby(['mgroup', 'precursor']).indices.items():
        v = psm[col].values[idx]
        dev[idx] = [abs(v[i] - np.median(np.delete(v, i))) for i in range(len(v))]
    psm[f'dev_{col}'] = dev
    psm[f'rel_{col}'] = dev / psm[col]
# run-offset removed (diagnostic: pure per-precursor scatter)
med = psm.groupby(['mgroup', 'precursor'])['1/K0'].transform('median')
signed = psm['1/K0'] - med
psm['dev_nooff'] = np.abs(signed - signed.groupby(psm['Raw file']).transform('median'))
psm['is_shared'] = psm.precursor.isin(shared)
psm['is_multi'] = psm.precursor.isin(multi) | (psm.wr_spread > 0.02)
rows = []
def add(name, d):
    rows.append(dict(variant=name, n=len(d), F_k0=d['dev_1/K0'].median(), F_k0_rel=d['rel_1/K0'].median(),
                     F_ccs=d['dev_CCS'].median(), F_ccs_rel=d['rel_CCS'].median(), F_k0_run_offset_removed=d['dev_nooff'].median(),
                     F_k0_median_of_study_F=d.groupby('study')['dev_1/K0'].median().median()))
add('pooled, all precursors (pre-registered)', psm)
add('pooled, shared precursors only', psm[psm.is_shared])
add('pooled, shared, no multi-conformer', psm[psm.is_shared & ~psm.is_multi])
add('pooled, all, no multi-conformer', psm[~psm.is_multi])
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/F_variants.csv', index=False)
pd.set_option('display.width', 250); print(R.round(5).to_string())
