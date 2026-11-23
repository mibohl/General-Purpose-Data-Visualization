import json, numpy as np, pandas as pd
OUT = '/home/user/work/out'
DEL = '/home/user/General-Purpose-Data-Visualization/rt_gonogo'
comp = pd.read_parquet(f'{OUT}/compounds.parquet')
meta = pd.read_csv(f'{OUT}/datasets_meta.csv', dtype={'id': str}).set_index('id')
surv = pd.read_csv(f'{OUT}/survey_all.csv', dtype={'id': str}).set_index('id')
z1 = {k: set(v) for k, v in json.load(open(f'{OUT}/z1_pools.json')).items()}
train2s = set(json.load(open('/home/user/work/2step_train_sets.json')))
r2 = pd.read_parquet(f'{OUT}/roi_2step.parquet')
ok2 = r2.groupby('dataset').roi.apply(lambda s: s.notna().all())

def roles(d):
    r = []
    for g, c in [('heldout:RP-2step', 'g_rp2step'), ('heldout:RP-proxy', 'g_rpproxy'), ('heldout:PUB', 'g_pub'), ('heldout:HILIC', 'g_hilic')]:
        if meta.loc[d, c]: r.append(g)
    if d in z1.get('RP|2step', ()): r.append('pool:Z1-RP-2step')
    if d in z1.get('RP|proxy', ()): r.append('pool:proxy+Z1-RP')
    if d in z1.get('HILIC|proxy', ()): r.append('pool:proxy+Z1-HILIC')
    if d in train2s: r.append('2step-checkpoint-training')
    return ';'.join(r)

FIELDS = [c for c in surv.columns if c.startswith('f_')]
prov = pd.DataFrame({
    'report_id': meta.index, 'name': meta.name, 'mode': meta['mode'], 'column': meta.column_name,
    'n_compounds_raw': surv.loc[meta.index, 'n_compounds'], 'n_compounds_clean': meta.n_clean,
    'gradient_length_min': meta.grad_len, 'gradient_length_source': meta.grad_len_source.fillna(''),
    't0_present': surv.loc[meta.index, 't0_present'], 't0_used_for_void': meta.t0_eff.round(4), 't0_source': meta.t0_source,
    'n_metadata_fields_populated': surv.loc[meta.index, 'n_meta_fields'], 'n_metadata_fields_tracked': 10,
    'pH_recorded': meta.pH, 'pH_imputed_for_2step': meta.pH_imputed, 'pH_2step_input': meta.pH_2step_input,
    'hsm_tanaka_available': meta.hsm_tanaka_available, 'in_2step_checkpoint_training': meta.in_2step_train,
    'checkpoint_roi_available': [bool(ok2.get(d, False)) for d in meta.index],
    'roles': [roles(d) for d in meta.index]})
for f in FIELDS:
    prov[f] = surv.loc[meta.index, f].values
prov = prov[prov.roles != '']
prov.to_csv(f'{DEL}/provenance.csv', index=False)
print('provenance rows', len(prov))

# metadata completeness table for held-out datasets
ho = prov[prov.roles.str.contains('heldout')].copy()
ho['well_documented'] = ho.t0_present & (ho.n_metadata_fields_populated >= 9)
cols = ['report_id', 'name', 'mode', 'roles', 't0_present', 'n_metadata_fields_populated', 'well_documented'] + FIELDS
ho[cols].to_csv(f'{DEL}/metadata_completeness.csv', index=False)

# compound overlap (InChIKey first block) across all used datasets
used = list(prov.report_id)
sets = {d: set(comp.loc[comp.dataset == d, 'ik14']) for d in used}
n = len(used)
ov = np.zeros((n, n), dtype=np.int32)
for i, a in enumerate(used):
    for j, b in enumerate(used):
        ov[i, j] = len(sets[a] & sets[b])
pd.DataFrame(ov, index=used, columns=used).to_csv(f'{OUT}/overlap_counts.csv')
size = np.diag(ov).astype(float)
frac_elsewhere = []
allk = {}
for d in used:
    for k in sets[d]: allk.setdefault(k, set()).add(d)
for d in used:
    frac_elsewhere.append(np.mean([len(allk[k]) > 1 for k in sets[d]]) if sets[d] else np.nan)
summary = pd.DataFrame({'report_id': used, 'n_unique_ik14': size, 'frac_compounds_in_any_other_dataset': frac_elsewhere,
                        'max_overlap_with_single_other_dataset': [(np.delete(ov[i], i) / max(size[i], 1)).max() for i in range(n)]})
summary.to_csv(f'{DEL}/overlap_summary.csv', index=False)
print('median fraction of a dataset\'s compounds found elsewhere:', round(np.nanmedian(frac_elsewhere), 3))
hs = summary.set_index('report_id').loc[ho.report_id]
print('held-out median fraction elsewhere:', round(hs.frac_compounds_in_any_other_dataset.median(), 3))
