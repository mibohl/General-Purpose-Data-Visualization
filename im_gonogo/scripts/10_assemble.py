"""Assemble results.csv, Step 3/4 tables, decision numbers and the residual-vs-k plot.
Usage: python -I 10_assemble.py <work_out_dir> <deliverable_dir>
"""
import sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
O, D = sys.argv[1], sys.argv[2]
s3 = pd.read_csv(f'{O}/step3_per_study.csv'); p3 = pd.read_csv(f'{O}/step3_pooled.csv')
s4 = pd.read_csv(f'{O}/step4_runs.csv')
Fv = pd.read_csv(f'{O}/F_variants.csv'); Fs = pd.read_csv(f'{O}/F_per_study.csv')
F = Fv.loc[0, 'F_k0']; Fccs = Fv.loc[0, 'F_ccs']
Fst = Fs[Fs.space == '1/K0'].set_index('study')['F']

# ---------- results.csv: every run ----------
a = s3.rename(columns={'method': 'arm', 'med': 'median_abs_dev'}).assign(step='step3_drift_vs_consensus', k=np.nan)
a['k'] = a.arm.str.extract(r'_(\d+)$')[0].astype(float)
b = s4.rename(columns={'mae_median': 'median_abs_dev'}).assign(step='step4_prediction_arms', space='k0', subset='all')
res = pd.concat([a[['step', 'study', 'space', 'subset', 'arm', 'k', 'seed', 'n_eval', 'median_abs_dev', 'q25', 'q75']],
                 b.rename(columns={'n_query': 'n_eval'})[['step', 'study', 'space', 'subset', 'arm', 'k', 'seed', 'n_eval', 'median_abs_dev', 'q25', 'q75']]])
res['F_pooled'] = np.where(res.space == 'ccs', Fccs, F)
res['median_over_F'] = res.median_abs_dev / res.F_pooled
res.to_csv(f'{D}/results.csv', index=False)

# ---------- decision numbers ----------
def g(space, subset, method):
    return p3[(p3.space == space) & (p3.subset == subset) & (p3.method == method)].med.iloc[0]
dec = []
for space, subset, Fval, label in [('k0', 'all', F, '1/K0, all shared precursors (PRIMARY)'),
                                    ('k0', 'no_multi', Fv.loc[2, 'F_k0'], '1/K0, multi-conformer candidates removed'),
                                    ('ccs', 'all', Fccs, 'CCS, all shared precursors'),
                                    ('ccs', 'no_multi', Fv.loc[2, 'F_ccs'], 'CCS, multi-conformer candidates removed'),
                                    ('k0', 'all', Fv.loc[0, 'F_k0_run_offset_removed'], '1/K0, stricter floor (run offsets removed from F)')]:
    Dv, A, A50, Aall = g(space, subset, 'raw'), g(space, subset, 'affine_20'), g(space, subset, 'affine_50'), g(space, subset, 'affine_all_insample')
    for oname in ['oracle_isotonic', 'oracle_isotonic_cv5', 'oracle_isotonic_per_charge']:
        Ov = g(space, subset, oname)
        c1 = A > 2 * Fval; c2 = (A - Ov) / A >= 0.25
        dec.append(dict(setting=label, oracle=oname, F=Fval, D=Dv, A=A, A50=A50, A_all_insample=Aall, O=Ov,
                        A_over_F=A / Fval, A_minus_O_over_A=(A - Ov) / A, crit1_pass=c1, crit2_pass=c2,
                        verdict_from_rules_1_2=('NO-GO' if not (c1 and c2) else 'needs Step 5')))
dec = pd.DataFrame(dec); dec.to_csv(f'{D}/decision_numbers.csv', index=False)

# ---------- Step 3 per-study table ----------
t = s3[(s3.space == 'k0') & (s3.subset == 'all')]
piv = t.groupby(['study', 'method']).med.median().unstack()
piv = piv[['raw', 'shift_20', 'affine_20', 'affine_50', 'affine_all_insample', 'oracle_isotonic']]
piv.insert(0, 'n_shared', t[t.method == 'raw'].set_index('study').n_eval)
piv.insert(1, 'F_study', Fst.reindex(piv.index))
piv.to_csv(f'{D}/step3_per_study_table.csv')

# ---------- Step 4 tables ----------
arms = ['Z1', 'C1', 'C2', 'C3', 'R1_ridge', 'R1_lgbm', 'T1', 'O1', 'O1_affine_all']
st = s4.groupby(['study', 'k', 'arm']).mae_median.median().unstack()[arms]
st['F_study'] = Fst.reindex(st.index.get_level_values(0)).values
st.to_csv(f'{D}/step4_per_study_k.csv')
pooled = s4.groupby(['k', 'arm']).mae_median.median().unstack()[arms]  # median over studies x seeds
per_sk = s4.groupby(['study', 'k', 'arm']).mae_median.median()
by_study_all = per_sk.groupby(['k', 'arm']).median().unstack()[arms]
by_study_all.to_csv(f'{D}/step4_pooled_median_over_studies_all_eligible.csv')
# fixed study set: studies with results at every k (avoids composition change across k)
full = s4.groupby('study').k.nunique(); full = full[full == 4].index
by_study = per_sk[per_sk.index.get_level_values('study').isin(full)].groupby(['k', 'arm']).median().unstack()[arms]
by_study.to_csv(f'{D}/step4_pooled_median_over_studies.csv')
print('fixed study set for plot/pooled table:', len(full), list(full))
# head-to-head vs C2 on identical (study,k,seed) draws
w = s4.pivot_table(index=['study', 'k', 'seed'], columns='arm', values='mae_median')
h2h = []
for arm in ['C1', 'C3', 'R1_ridge', 'R1_lgbm', 'T1', 'O1']:
    rel = (w['C2'] - w[arm]) / w['C2']
    for k, r in rel.groupby(level='k'):
        h2h.append(dict(arm=arm, k=k, n=len(r), median_rel_improvement_over_C2=r.median(),
                        frac_draws_better_than_C2=(r > 0).mean(), frac_better_by_10pct=(r > 0.10).mean()))
pd.DataFrame(h2h).to_csv(f'{D}/step4_head_to_head_vs_C2.csv', index=False)

# ---------- plot: residual vs k, one line per arm ----------
pal = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
mk = ['o', 's', 'D', '^', 'v', 'P', 'X', '*']
plot_arms = ['Z1', 'C1', 'C2', 'C3', 'R1_ridge', 'R1_lgbm', 'T1', 'O1']
labels = {'Z1': 'Z1 IM2Deep as shipped', 'C1': 'C1 shift', 'C2': 'C2 affine', 'C3': 'C3 isotonic',
          'R1_ridge': 'R1 ridge residual', 'R1_lgbm': 'R1 LightGBM residual', 'T1': 'T1 TabPFN v2 residual', 'O1': 'O1 oracle (ceiling)'}
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e', 'ytick.color': '#52514e'})
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=150, facecolor='#fcfcfb', gridspec_kw={'width_ratios': [1, 1.25]})
for ax in axes:
    ax.set_facecolor('#fcfcfb'); ax.grid(axis='y', color='#e4e3df', lw=0.8); ax.set_axisbelow(True)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
ks = sorted(by_study.index)
for i, arm in enumerate(plot_arms):
    for ax in axes:
        ax.plot(ks, by_study[arm].values, color=pal[i], lw=2, marker=mk[i], ms=7, mec='#fcfcfb', mew=1.2, label=labels[arm])
for ax, (lo, hi) in zip(axes, [(0, 0.028), (0.0095, 0.0135)]):
    ax.axhspan(0.007, 0.010, color='#d6d5d0', alpha=0.5, lw=0)
    ax.axhline(F, color='#52514e', lw=1, ls='--'); ax.axhline(2 * F, color='#52514e', lw=1, ls=':')
    ax.set_xscale('log'); ax.minorticks_off(); ax.set_xticks(ks); ax.set_xticklabels([str(k) for k in ks]); ax.set_ylim(lo, hi)
    ax.set_xlabel('k reference precursors from the held-out study')
axes[0].set_ylabel('median |predicted - observed| 1/K0 (Vs/cm²)\nmedian over held-out studies')
axes[0].text(ks[0], F - 0.0009, f'F = {F:.4f} (replicate floor)', va='center', fontsize=8.5, color='#52514e')
axes[0].text(ks[0], 2 * F - 0.0009, f'2F = {2*F:.4f}', va='center', fontsize=8.5, color='#52514e')
axes[0].text(ks[0], 0.0085, 'published fine-tuned range 0.007-0.010', va='center', fontsize=8.5, color='#52514e')
fig.suptitle(f'Prediction residual vs number of reference precursors ({len(full)} held-out studies with results at every k; 5 anchor seeds each)', x=0.01, ha='left', fontsize=10.5, color='#0b0b0b')
axes[0].set_title('All arms', loc='left', color='#0b0b0b'); axes[1].set_title('Zoom on adapted arms', loc='left', color='#0b0b0b')
axes[1].legend(frameon=False, fontsize=8.5, loc='upper left', bbox_to_anchor=(1.0, 1.0))
fig.tight_layout(); fig.savefig(f'{D}/residual_vs_k.png', facecolor='#fcfcfb')
pd.set_option('display.width', 250)
print(dec.round(4).to_string()); print(by_study.round(4).to_string()); print(pd.DataFrame(h2h).round(3).to_string())
