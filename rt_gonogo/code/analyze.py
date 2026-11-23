"""Aggregate runs -> decision quantities, tables, plots."""
import glob, json, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
OUT = '/home/user/work/out'
DEL = '/home/user/General-Purpose-Data-Visualization/rt_gonogo'
os.makedirs(f'{DEL}/figures', exist_ok=True)
GROUPS = ['RP-2step', 'RP-proxy', 'PUB-2step', 'PUB-proxy', 'HILIC']
ICL = ['T1', 'T2', 'R1']
rng = np.random.default_rng(12345)

res = pd.concat([pd.read_parquet(f) for f in glob.glob(f'{OUT}/runs/*.parquet')
                 if not f.endswith(('_pairs.parquet', '_timing.parquet'))], ignore_index=True)
pairs = pd.concat([pd.read_parquet(f) for f in glob.glob(f'{OUT}/runs/*_pairs.parquet')], ignore_index=True)
timing = pd.concat([pd.read_parquet(f) for f in glob.glob(f'{OUT}/runs/*_timing.parquet')], ignore_index=True)
meta = pd.read_csv(f'{OUT}/datasets_meta.csv', dtype={'id': str}).set_index('id')
prov = pd.read_csv(f'{DEL}/provenance.csv', dtype={'report_id': str}).set_index('report_id')
comp = pd.read_parquet(f'{OUT}/compounds.parquet')
r2 = pd.read_parquet(f'{OUT}/roi_2step.parquet').rename(columns={'roi': 'roi_2step'})
rp = pd.read_parquet(f'{OUT}/roi_proxy.parquet').rename(columns={'roi': 'roi_proxy'})
comp = comp.merge(r2, on=['dataset', 'cid'], how='left').merge(rp, on=['dataset', 'cid'], how='left')
res['name'] = meta.loc[res.dataset, 'name'].values
res.to_csv(f'{DEL}/results.csv', index=False, float_format='%.5g')

# ---- purged eligibility: >= 20 compounds outside pool (n_query + k >= 20 -> use max n_query + k) ----
pur = res[(res.eval_set == 'purged') & (res.arm == 'I1')]
outside = (pur.n_query + pur.k).groupby([pur.group, pur.dataset]).max()
elig = {ev: set() for ev in ['purged', 'scaffold']}
for ev in elig:
    s = res[(res.eval_set == ev) & (res.arm == 'I1')]
    o = (s.n_query + s.k).groupby([s.group, s.dataset]).max()  # = compounds outside pool (anchors may be among them)
    elig[ev] = set(o[o >= 20].index)

def valid(df, ev):
    if ev == 'all':
        return df
    df = df[df.n_query >= 5]
    return df[[(g, d) in elig[ev] for g, d in zip(df.group, df.dataset)]]

# ---- per-run relative quantities ----
W = res.pivot_table(index=['group', 'dataset', 'roi_source', 'scheme', 'k', 'seed', 'eval_set', 'n_query'],
                    columns='arm', values='mae_min').reset_index()
W['H'] = (W.I1 - W.O1) / W.I1
W['H_cv'] = (W.I1 - W.O1_CV) / W.I1
for a in ['I2', 'I3', 'I4', 'Z1', 'T1', 'T2', 'R1']:
    if a in W: W[f'G_{a}'] = (W.I1 - W[a]) / W.I1
# per dataset: mean over seeds
GC = [c for c in W.columns if c.startswith('G_')]
ARMS = [c for c in ['I1', 'I2', 'I3', 'I4', 'O1', 'O1_CV', 'O1_L2', 'Z1', 'T1', 'T2', 'R1'] if c in W]
PD = []
for ev in ['all', 'purged', 'scaffold']:
    w = valid(W[W.eval_set == ev], ev)
    PD.append(w.groupby(['group', 'dataset', 'eval_set', 'scheme', 'k'])[['H', 'H_cv'] + GC + ARMS + ['n_query']].mean().reset_index())
PD = pd.concat(PD)
PD['grad_len'] = meta.loc[PD.dataset, 'grad_len'].values
PD.to_csv(f'{OUT}/per_dataset.csv', index=False)

def boot_median(x, B=10000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return (np.nan, np.nan)
    b = np.median(rng.choice(x, (B, len(x))), axis=1)
    return tuple(np.percentile(b, [2.5, 97.5]))

def verdict(H, G, lo):
    if not np.isfinite(H): return 'n/a'
    if H < 0.15: return 'NO-GO'
    if H < 0.25: return 'INDETERMINATE (rule gap)'
    if np.isfinite(G) and G >= 0.15 and lo > 0: return 'GO'
    return 'CONDITIONAL'

rows = []
for (g, ev, sc, k), d in PD.groupby(['group', 'eval_set', 'scheme', 'k']):
    r = dict(group=g, eval_set=ev, scheme=sc, k=k, n_datasets=d.dataset.nunique(),
             median_H=d.H.median(), median_H_cv=d.H_cv.median(),
             median_MAE_I1_min=d.I1.median(), median_MAE_O1_min=d.O1.median())
    r['H_CI_lo'], r['H_CI_hi'] = boot_median(d.H)
    for a in ['I2', 'I3', 'I4', 'Z1', 'T1', 'T2', 'R1']:
        if f'G_{a}' in d: r[f'median_G_{a}'] = d[f'G_{a}'].median()
    best = max(ICL, key=lambda a: r.get(f'median_G_{a}', -np.inf) if np.isfinite(r.get(f'median_G_{a}', np.nan)) else -np.inf)
    r['best_ICL_arm'] = best; r['G'] = r.get(f'median_G_{best}', np.nan)
    r['G_CI_lo'], r['G_CI_hi'] = boot_median(d[f'G_{best}']) if f'G_{best}' in d else (np.nan, np.nan)
    r['verdict'] = verdict(r['median_H'], r['G'], r['G_CI_lo'])
    rows.append(r)
DEC = pd.DataFrame(rows)
DEC.to_csv(f'{DEL}/decision_table.csv', index=False, float_format='%.4f')

ORDER = {'NO-GO': 0, 'INDETERMINATE (rule gap)': 1, 'CONDITIONAL': 2, 'GO': 3, 'n/a': -1}
def headline(group, ev='purged', scheme='spread'):
    d = DEC[(DEC.group == group) & (DEC.eval_set == ev) & (DEC.scheme == scheme) & DEC.k.isin([15, 30])].set_index('k')
    vs = [d.loc[k, 'verdict'] for k in (15, 30) if k in d.index]
    vs = [v for v in vs if v != 'n/a']
    return d, (min(vs, key=lambda v: ORDER[v]) if vs else 'n/a')

summary = {}
for g in GROUPS:
    for ev in ['purged', 'all', 'scaffold']:
        for sc in ['spread', 'random']:
            d, v = headline(g, ev, sc)
            summary[f'{g}|{ev}|{sc}'] = dict(verdict=v, **{f'k{k}': d.loc[k, ['n_datasets', 'median_H', 'H_CI_lo', 'H_CI_hi', 'median_H_cv', 'G', 'best_ICL_arm', 'G_CI_lo', 'G_CI_hi', 'verdict']].to_dict() for k in d.index})
json.dump(summary, open(f'{OUT}/headline.json', 'w'), indent=1, default=float)

# ---- Step 4: error decomposition, minutes, per held-out dataset ----
dec_rows = []
for (g, d), x in PD[(PD.scheme == 'spread') & PD.k.isin([15, 30])].groupby(['group', 'dataset']):
    r = dict(group=g, dataset=d, name=meta.loc[d, 'name'], column=meta.loc[d, 'column_name'],
             n=int(meta.loc[d, 'n_clean']), grad_len_min=meta.loc[d, 'grad_len'])
    for ev in ['all', 'purged']:
        for k in (15, 30):
            y = x[(x.eval_set == ev) & (x.k == k)]
            if len(y):
                y = y.iloc[0]
                r[f'{ev}_k{k}_total_I1'] = y.I1; r[f'{ev}_k{k}_order_O1'] = y.O1; r[f'{ev}_k{k}_warp'] = y.I1 - y.O1
                r[f'{ev}_k{k}_warp_share'] = y.H; r[f'{ev}_k{k}_nq'] = y.n_query
    dec_rows.append(r)
DECOMP = pd.DataFrame(dec_rows)
DECOMP.to_csv(f'{DEL}/error_decomposition.csv', index=False, float_format='%.4f')

# ---- Step 5: documentation split + Z1 ----
s5 = []
for g in GROUPS:
    for ev in ['all', 'purged']:
        for k in (15, 30):
            x = PD[(PD.group == g) & (PD.eval_set == ev) & (PD.scheme == 'spread') & (PD.k == k)].copy()
            if not len(x): continue
            x['well'] = [bool(prov.loc[d, 't0_present']) and prov.loc[d, 'n_metadata_fields_populated'] >= 9 for d in x.dataset]
            for w, y in x.groupby('well'):
                s5.append(dict(group=g, eval_set=ev, k=k, well_documented=w, n_datasets=len(y), median_H=y.H.median(),
                               median_Z1_over_I1=(y.Z1 / y.I1).median(), median_MAE_I1=y.I1.median(), median_MAE_Z1=y.Z1.median()))
            s5.append(dict(group=g, eval_set=ev, k=k, well_documented='all', n_datasets=len(x), median_H=x.H.median(),
                           median_Z1_over_I1=(x.Z1 / x.I1).median(), median_MAE_I1=x.I1.median(), median_MAE_Z1=x.Z1.median()))
S5 = pd.DataFrame(s5); S5.to_csv(f'{DEL}/step5_documentation_split.csv', index=False, float_format='%.4f')

# ---- Step 6: order flips (dataset level, all compounds) + arm recovery at k=30 spread ----
fl = []
for g in GROUPS:
    src = 'roi_2step' if '2step' in g else 'roi_proxy'
    for d in PD[PD.group == g].dataset.unique():
        s = comp[comp.dataset == d].dropna(subset=[src])
        r, o = s.rt.values, s[src].values
        iu = np.triu_indices(len(r), 1)
        drt, dro = r[iu[0]] - r[iu[1]], o[iu[0]] - o[iu[1]]
        m = (np.abs(drt) >= 10 / 60) & (dro != 0)
        fl.append(dict(group=g, dataset=d, name=meta.loc[d, 'name'], n_pairs=int(m.sum()),
                       flip_fraction=float((np.sign(drt[m]) != np.sign(dro[m])).mean()),
                       spearman=spearmanr(r, o).statistic))
FL = pd.DataFrame(fl)
pp = pairs[pairs.scheme == 'spread'].groupby(['group', 'dataset', 'arm'])[['conc_flipped', 'conc_nonflipped', 'conc_all']].mean().reset_index()
rec = pp.pivot_table(index=['group', 'dataset'], columns='arm', values='conc_flipped')
non = pp.pivot_table(index=['group', 'dataset'], columns='arm', values='conc_nonflipped')
alla = pp.pivot_table(index=['group', 'dataset'], columns='arm', values='conc_all')
FL = FL.set_index(['group', 'dataset'])
for a in ['I1', 'I2', 'I4', 'R1', 'T1', 'T2']:
    if a in rec:
        FL[f'flipped_recovered_{a}'] = rec[a]; FL[f'nonflipped_kept_{a}'] = non[a]; FL[f'pair_acc_{a}'] = alla[a]
FL.reset_index().to_csv(f'{DEL}/step6_order_flips.csv', index=False, float_format='%.4f')

# ---- timing ----
TT = timing.groupby(['group', 'arm']).seconds.sum().unstack(0)
nfits = res[(res.eval_set == 'all')].groupby(['group', 'arm']).size().unstack(0)
TIM = (TT / nfits.reindex_like(TT)).round(4)
TIM.to_csv(f'{DEL}/timing_per_fit_seconds.csv')
TT.round(1).to_csv(f'{DEL}/timing_total_seconds.csv')

# ---- plots ----
def plot_vs_k(ev, fn):
    fig, axes = plt.subplots(1, len(GROUPS), figsize=(4.2 * len(GROUPS), 3.8), sharey=False)
    arms = ['I1', 'I2', 'I3', 'I4', 'O1', 'Z1', 'T1', 'T2', 'R1']
    colors = dict(I1='#1f2937', I2='#9ca3af', I3='#d1d5db', I4='#6b7280', O1='#16a34a', Z1='#a855f7',
                  T1='#2563eb', T2='#0ea5e9', R1='#f97316')
    styles = dict(O1='--', Z1=':')
    for ax, g in zip(axes, GROUPS):
        x = PD[(PD.group == g) & (PD.eval_set == ev) & (PD.scheme == 'spread')]
        if not len(x):
            ax.set_title(f'{g}\n(no eligible datasets)', fontsize=9); ax.axis('off'); continue
        for a in arms:
            if a not in x: continue
            m = x.groupby('k')[a].median()
            ax.plot(m.index, m.values, marker='o', ms=3, lw=2.2 if a in ('I1', 'O1') else 1.4,
                    ls=styles.get(a, '-'), color=colors[a], label=a)
        ax.set_title(f'{g}  (n={x.dataset.nunique()} datasets)', fontsize=9)
        ax.set_xlabel('k anchors'); ax.set_xticks([10, 15, 20, 30, 50]); ax.grid(alpha=.3)
    axes[0].set_ylabel('median over datasets of MAE (min)')
    axes[-1].legend(fontsize=7, loc='upper right', ncol=2)
    fig.suptitle(f'Error in minutes vs k, spread anchors, eval set = {ev}', fontsize=10)
    fig.tight_layout(); fig.savefig(fn, dpi=150); plt.close(fig)
plot_vs_k('all', f'{DEL}/figures/error_vs_k_all.png')
plot_vs_k('purged', f'{DEL}/figures/error_vs_k_purged.png')

def plot_headroom(fn):
    gs = [g for g in ['RP-2step', 'RP-proxy', 'PUB-2step', 'PUB-proxy', 'HILIC'] if (PD.group == g).any()]
    fig, axes = plt.subplots(len(gs), 1, figsize=(12, 2.6 * len(gs)))
    for ax, g in zip(np.atleast_1d(axes), gs):
        x = PD[(PD.group == g) & (PD.scheme == 'spread') & (PD.eval_set == 'all')]
        p15 = x[x.k == 15].set_index('dataset').H; p30 = x[x.k == 30].set_index('dataset').H
        order = p30.reindex(p15.index).fillna(p15).sort_values().index
        xi = np.arange(len(order))
        ax.bar(xi - 0.2, p15.reindex(order), 0.4, label='k=15', color='#94a3b8')
        ax.bar(xi + 0.2, p30.reindex(order), 0.4, label='k=30', color='#1e40af')
        xp = PD[(PD.group == g) & (PD.scheme == 'spread') & (PD.eval_set == 'purged') & (PD.k == 30)].set_index('dataset').H
        ax.scatter(xi[np.isin(order, xp.index)] + 0.2, xp.reindex(order).dropna(), marker='D', s=14, color='#f97316', zorder=3, label='k=30 purged')
        for h, c in [(0.15, '#dc2626'), (0.25, '#16a34a')]:
            ax.axhline(h, color=c, lw=1, ls='--')
        ax.set_xticks(xi); ax.set_xticklabels([str(meta.loc[d, 'name'])[:18] for d in order], rotation=60, ha='right', fontsize=6.5)
        ax.set_ylabel('H'); ax.set_title(f'{g}: headroom H = (I1 - O1) / I1 per dataset (spread anchors; dashed: 0.15 NO-GO, 0.25 GO bars)', fontsize=9)
        ax.set_ylim(min(-0.1, np.nanmin(np.r_[p15.values, p30.values, 0]) - 0.05), 1)
        ax.legend(fontsize=7, loc='upper left')
    fig.tight_layout(); fig.savefig(fn, dpi=150); plt.close(fig)
plot_headroom(f'{DEL}/figures/headroom_per_dataset.png')

print(DEC[(DEC.k.isin([15, 30]))][['group', 'eval_set', 'scheme', 'k', 'n_datasets', 'median_H', 'H_CI_lo', 'H_CI_hi', 'median_H_cv',
                                   'median_G_T1', 'median_G_T2', 'median_G_R1', 'G', 'G_CI_lo', 'G_CI_hi', 'verdict']].round(3).to_string())
print({k: v['verdict'] for k, v in summary.items() if '|spread' in k})
