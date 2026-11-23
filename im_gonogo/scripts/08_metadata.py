"""Step 5b: does recorded metadata predict per-study drift? + provenance table.
Usage: python -I 08_metadata.py <work_dir> <out_dir>
"""
import sys, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
W, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(f'{W}/pride_meta.json'))
meta['PXD036206'] = json.load(open(f'{W}/pride_meta_036206.json'))
sv = pd.read_parquet(f'{W}/study_values.parquet')
sh = pd.read_parquet(f'{W}/shared_with_consensus.parquet')
R = pd.read_csv(f'{W}/q1_run_index_map.csv'); R['study'] = R['study'].replace({'PXD003206': 'PXD036206'})
R['k0_idx0'] = R['k0_at_idx500'] - 500 * R['slope']; R['k0_idxmax'] = R['k0_idx0'] + R['slope'] * R['idx_max']
st5 = pd.read_csv(f'{OUT}/step5_variance_explained.csv'); st5 = st5[st5.subset == 'all'].set_index('study')
rows = []
for s, g in sv.groupby('study'):
    m = meta.get(s, {})
    lab = (m.get('labPIs') or [{}])[0].get('name', 'unknown (synthetic 21-PTM library, Kuster-lab naming; no accession in file)') if m else 'unknown (synthetic 21-PTM library; accession not given in file names)'
    aff = (m.get('labPIs') or [{}])[0].get('affiliation', '') if m else ''
    r = R[R.study == s]
    h = sh[sh.study == s]
    row = dict(study=s, title=m.get('title', 'Synthetic peptide library with 21 PTM classes (IM2Deep Zenodo, evidence_*mod_*.txt)')[:90],
               instrument=', '.join(i['name'] for i in m.get('instruments', [])) or 'timsTOF (model not recorded)',
               lab=lab, affiliation=aff.split(',')[0][:60], organism=', '.join(o['name'] for o in m.get('organisms', []))[:60],
               mode='DDA-PASEF', search_engine='MaxQuant (evidence.txt)', n_runs=g.n_runs.max() if False else int(r.run.nunique()),
               n_precursors=len(g), n_shared=len(h),
               tims_scan_slope=r.slope.median(), k0_hi=r.k0_idx0.median(), k0_lo=r.k0_idxmax.median())
    if len(h) >= 60:
        x, c = h.k0.values, h.cons_k0.values
        b = np.linalg.lstsq(np.c_[np.ones_like(x), x], c, rcond=None)[0]
        row.update(D_study=np.median(np.abs(x - c)), offset_at_1=b[0] + b[1] * 1.0 - 1.0, affine_slope=b[1],
                   signed_median_shift=np.median(c - x))
    if s in st5.index:
        row.update(z4_minus_z2=st5.loc[s, 'mean_resid_z4'] - st5.loc[s, 'mean_resid_z2'], r2_lgbm_cv=st5.loc[s, 'r2_lgbm_cv'])
    rows.append(row)
P = pd.DataFrame(rows)
P.to_csv(f'{OUT}/provenance.csv', index=False)
# associations
P['is_pro2'] = P.instrument.str.contains('Pro 2')
P['k0_span'] = P.k0_hi - P.k0_lo
assoc = []
for tgt in ['signed_median_shift', 'D_study', 'affine_slope', 'z4_minus_z2', 'r2_lgbm_cv']:
    d = P.dropna(subset=[tgt])
    for cov in ['tims_scan_slope', 'k0_hi', 'k0_lo', 'k0_span']:
        rho, p = spearmanr(d[cov], d[tgt]); assoc.append(dict(target=tgt, covariate=cov, test='spearman', stat=rho, p=p, n=len(d)))
    a, b = d[d.is_pro2][tgt], d[~d.is_pro2][tgt]
    if len(a) >= 2:
        u, p = mannwhitneyu(a, b); assoc.append(dict(target=tgt, covariate='timsTOF Pro 2 vs Pro', test='mann-whitney', stat=a.median() - b.median(), p=p, n=len(d)))
    gl = d.lab.str.contains('Gerner')
    u, p = mannwhitneyu(d[gl][tgt], d[~gl][tgt]); assoc.append(dict(target=tgt, covariate='Gerner lab vs others', test='mann-whitney', stat=d[gl][tgt].median() - d[~gl][tgt].median(), p=p, n=len(d)))
A = pd.DataFrame(assoc); A.to_csv(f'{OUT}/step5_metadata_assoc.csv', index=False)
pd.set_option('display.width', 250)
print(P[['study', 'instrument', 'lab', 'n_runs', 'n_precursors', 'n_shared', 'tims_scan_slope', 'k0_hi', 'k0_lo', 'signed_median_shift', 'D_study']].round(4).to_string())
print(A.round(3).to_string())
