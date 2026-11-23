"""Render markdown tables for report.md. Usage: python -I 11_md_tables.py <work_out_dir> <deliverable_dir>"""
import sys
import numpy as np, pandas as pd
O, D = sys.argv[1], sys.argv[2]
def md(df, fmt=None):
    fmt = fmt or {}
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|']
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                cells.append('' if np.isnan(v) else (fmt.get(c, '{:.4f}').format(v)))
            else:
                cells.append(str(v))
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)
out = {}
P = pd.read_csv(f'{O}/provenance.csv')
P['lab'] = P['lab'].str.replace('Professor ', '')
P['lab (PI, affiliation)'] = P['lab'] + (', ' + P['affiliation'].fillna('')).where(P['affiliation'].fillna('') != '', '')
P['TIMS 1/K0 range'] = P.k0_lo.round(2).astype(str) + '-' + P.k0_hi.round(2).astype(str)
out['provenance'] = md(P[['study', 'instrument', 'mode', 'search_engine', 'lab (PI, affiliation)', 'organism', 'n_runs', 'n_precursors', 'n_shared', 'TIMS 1/K0 range']])
t = pd.read_csv(f'{D}/step3_per_study_table.csv')
t['A20/F_study'] = t.affine_20 / t.F_study
t['(A20-O)/A20'] = (t.affine_20 - t.oracle_isotonic) / t.affine_20
t['n_shared'] = t.n_shared.astype(int)
out['step3'] = md(t.rename(columns={'raw': 'D (raw)', 'shift_20': 'shift, 20', 'affine_20': 'A (affine, 20)', 'affine_50': 'affine, 50',
                                    'affine_all_insample': 'affine, all (in-sample)', 'oracle_isotonic': 'O (isotonic, all)'}),
                  {'A20/F_study': '{:.2f}', '(A20-O)/A20': '{:.2f}'})
s5 = pd.read_csv(f'{O}/step5_variance_explained.csv')
a = s5[s5.subset == 'all'].set_index('study'); b = s5[s5.subset == 'no_multi'].set_index('study')
T = pd.DataFrame({'n shared': a.n.astype(int), 'residual MAD': a.resid_mad,
                  'R2 ridge (5-fold CV)': a.r2_ridge_cv, 'R2 LightGBM (5-fold CV)': a.r2_lgbm_cv,
                  'R2 pooled other-study model': a.r2_lgbm_loso_pooled,
                  'R2 LightGBM CV, no multi': b.r2_lgbm_cv, 'R2 pooled other-study, no multi': b.r2_lgbm_loso_pooled,
                  'mean resid z2': a.mean_resid_z2, 'mean resid z3': a.mean_resid_z3, 'mean resid z4': a.mean_resid_z4}).reset_index()
med = T.drop(columns='study').median(); med['study'] = '**median**'
T = pd.concat([T, med.to_frame().T], ignore_index=True)
for c in T.columns[1:]: T[c] = T[c].astype(float)
out['step5'] = md(T, {c: '{:.3f}' for c in T.columns if c.startswith('R2')} | {'n shared': '{:.0f}'})
A = pd.read_csv(f'{O}/step5_metadata_assoc.csv')
out['meta'] = md(A, {'stat': '{:.3f}', 'p': '{:.3f}', 'n': '{:.0f}'})
Fs = pd.read_csv(f'{O}/F_per_study.csv'); Fs = Fs[Fs.space == '1/K0'][['study', 'n_runs', 'n_dev', 'F', 'F_iqr_lo', 'F_iqr_hi']]
Fs['n_dev'] = Fs.n_dev.astype(int)
out['F'] = md(Fs.rename(columns={'n_dev': 'n run-level deviations', 'F_iqr_lo': 'q25', 'F_iqr_hi': 'q75'}), {'n run-level deviations': '{:.0f}'})
Fv = pd.read_csv(f'{O}/F_variants.csv'); Fv['n'] = Fv.n.astype(int)
out['Fv'] = md(Fv, {'F_ccs': '{:.3f}', 'n': '{:.0f}'})
for k, v in out.items():
    open(f'{D}/_tbl_{k}.md', 'w').write(v + '\n')
print('ok', list(out))
