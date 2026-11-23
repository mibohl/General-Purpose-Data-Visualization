"""Render Step 4 markdown (pooled table, head-to-head, per-study x k appendix). Usage: python -I 15_step4_md.py <deliverable_dir>"""
import sys
import numpy as np, pandas as pd
D = sys.argv[1]
F = pd.read_csv(f'{D}/decision_numbers.csv').F.iloc[0]
P = pd.read_csv(f'{D}/step4_pooled_median_over_studies.csv').set_index('k')
arms = ['Z1', 'C1', 'C2', 'C3', 'R1_ridge', 'R1_lgbm', 'T1', 'O1', 'G1b_global_then_affine']
L = ['| k | ' + ' | '.join(arms) + ' |', '|' + '---|' * (len(arms) + 1)]
for k, r in P.iterrows():
    L.append(f'| {k} | ' + ' | '.join(f'{r[a]:.4f} ({r[a] / F:.1f} F)' for a in arms) + ' |')
pooled = '\n'.join(L)
H = pd.read_csv(f'{D}/step4_head_to_head_vs_C2.csv')
H = H[H.arm.isin(['C1', 'C3', 'R1_ridge', 'R1_lgbm', 'T1', 'O1', 'G1_global_resid', 'G1b_global_then_affine'])]
L = ['| arm | k | draws | median relative improvement over C2 | share of draws better than C2 | share better by > 10 % |', '|---|---|---|---|---|---|']
for _, r in H.iterrows():
    L.append(f'| {r.arm} | {r.k} | {r.n} | {100 * r.median_rel_improvement_over_C2:+.1f} % | {100 * r.frac_draws_better_than_C2:.0f} % | {100 * r.frac_better_by_10pct:.0f} % |')
h2h = '\n'.join(L)
S = pd.read_csv(f'{D}/step4_per_study_k.csv')
cols = ['Z1', 'C1', 'C2', 'C3', 'R1_ridge', 'R1_lgbm', 'T1', 'O1', 'G1b_global_then_affine']
L = ['| study | k | F_study | ' + ' | '.join(cols) + ' | C2 / F | T1 / F | T1 vs C2 |', '|' + '---|' * (len(cols) + 6)]
for _, r in S.iterrows():
    fs = '' if np.isnan(r.F_study) else f'{r.F_study:.4f}'
    L.append(f'| {r.study} | {r.k} | {fs} | ' + ' | '.join(f'{r[c]:.4f}' for c in cols) +
             f' | {r.C2 / F:.1f} | {r.T1 / F:.1f} | {100 * (r.C2 - r.T1) / r.C2:+.0f} % |')
per = '\n'.join(L)
open(f'{D}/report_src/_tbl_step4_pooled.md', 'w').write(pooled); open(f'{D}/report_src/_tbl_step4_h2h.md', 'w').write(h2h); open(f'{D}/report_src/_tbl_step4_per.md', 'w').write(per)
print(pooled); print(len(S))
