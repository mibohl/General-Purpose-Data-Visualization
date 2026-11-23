"""Clean RepoRT datasets, compute molecular + metadata features, select held-out sets."""
import os, json, warnings
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*'); warnings.filterwarnings('ignore')

R = '/home/user/data/RepoRT/processed_data'
OUT = '/home/user/work/out'
train2s = set(json.load(open('/home/user/work/2step_train_sets.json')))
PUB = {'0002': 'FEM long', '0003': 'IPB Halle', '0018': 'UniToyama Atlantis',
       '0019': 'Eawag XBridgeC18', '0054': 'LIFE old', '0055': 'LIFE new'}
survey = pd.read_csv(f'{OUT}/survey_all.csv', dtype={'id': str}).set_index('id')
hsm = pd.read_csv('/home/user/data/RepoRT/resources/hsm_database/hsm_database.tsv', sep='\t')
tan = pd.read_csv('/home/user/data/RepoRT/resources/tanaka_database/tanaka_database.tsv', sep='\t')
HSM_F = ['H', 'S*', 'A', 'B', 'C (pH 2.8)', 'C (pH 7.0)']
TAN_F = ['kPB', 'αCH2', 'αT/O', 'αC/P', 'αB/P', 'αB/P.1']
hsm_m = hsm.groupby('name_normalized')[HSM_F].mean()
tan_m = tan.groupby('name_normalized')[TAN_F].mean()

DESC = {'MolLogP': Descriptors.MolLogP, 'TPSA': Descriptors.TPSA, 'MolWt': Descriptors.MolWt,
        'NumHDonors': Descriptors.NumHDonors, 'NumHAcceptors': Descriptors.NumHAcceptors,
        'NumRotatableBonds': Descriptors.NumRotatableBonds, 'NumAromaticRings': rdMolDescriptors.CalcNumAromaticRings,
        'FractionCSP3': rdMolDescriptors.CalcFractionCSP3, 'HeavyAtomCount': Descriptors.HeavyAtomCount,
        'FormalCharge': lambda m: Chem.GetFormalCharge(m), 'NumHeteroatoms': Descriptors.NumHeteroatoms,
        'MolMR': Descriptors.MolMR}

def num(v):
    try:
        v = float(v); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan

def gradient_features(d):
    try:
        g = pd.read_csv(f'{R}/{d}/{d}_gradient.tsv', sep='\t')
    except Exception:
        return {}
    if len(g) < 2:
        return {}
    g = g.apply(pd.to_numeric, errors='coerce').dropna(subset=['t [min]', 'A [%]'])
    if len(g) < 2:
        return {}
    t = g['t [min]'].values.astype(float)
    org = 100 - g['A [%]'].values.astype(float)
    L = t.max()
    f = {'grad_len': L, 'grad_org_start': org[0], 'grad_org_max': org.max()}
    # initial isocratic hold
    ch = np.where(np.abs(org - org[0]) > 1e-6)[0]
    f['grad_hold'] = t[ch[0] - 1] if len(ch) else L
    span = org.max() - org[0]
    for q in (0.5, 0.9):
        target = org[0] + q * span
        idx = np.where(org >= target - 1e-9)[0]
        if len(idx) and span > 0:
            i = idx[0]
            if i == 0: f[f'grad_t{int(q*100)}'] = t[0]
            else:
                f[f'grad_t{int(q*100)}'] = t[i-1] + (target - org[i-1]) / max(org[i] - org[i-1], 1e-9) * (t[i] - t[i-1])
        else:
            f[f'grad_t{int(q*100)}'] = np.nan
    f['grad_slope'] = span / max(f.get('grad_t90', L) - f['grad_hold'], 1e-3) if span > 0 else 0.0
    fl = g.get('flow rate [ml/min]')
    f['grad_flow_var'] = float(fl.astype(float).std()) if fl is not None else np.nan
    return f

def metadata_features(d):
    md = pd.read_csv(f'{R}/{d}/{d}_metadata.tsv', sep='\t').iloc[0]
    f = {k.replace('column.', 'col_'): num(md.get(k)) for k in
         ['column.length', 'column.id', 'column.particle.size', 'column.temperature', 'column.flowrate']}
    t0 = num(md.get('column.t0'))
    f['t0_reported'] = t0 if (t0 == t0 and t0 > 0) else np.nan
    est = np.nan
    if all(np.isfinite([f['col_length'], f['col_id'], f['col_flowrate']])) and f['col_flowrate'] > 0:
        est = 0.68 * np.pi * (f['col_id'] / 2) ** 2 * f['col_length'] / 1000 / f['col_flowrate']
    f['t0_eff'] = f['t0_reported'] if np.isfinite(f['t0_reported']) else est
    f['t0_source'] = 'metadata' if np.isfinite(f['t0_reported']) else ('estimated' if np.isfinite(est) else 'none')
    phs = [num(md.get(f'eluent.{p}.pH')) for p in 'AB']
    phs = [p for p in phs if p == p and p > 0]
    f['pH'] = phs[0] if phs else np.nan
    for s in ['acn', 'meoh', 'iproh']:
        f[f'B_{s}'] = num(md.get(f'eluent.B.{s}'))
    for a in ['formic', 'acetic', 'nh4ac', 'nh4form', 'trifluoroacetic', 'nh4bicarb', 'nh4oh' ]:
        f[f'A_{a}'] = float(any(num(md.get(f'eluent.{p}.{a}')) > 0 for p in 'AB'))
    cname = md.get('column.name')
    f['column_name'] = cname if isinstance(cname, str) else None
    for F, db, tag in [(HSM_F, hsm_m, 'hsm'), (TAN_F, tan_m, 'tan')]:
        for k in F:
            f[f'{tag}_{k}'] = db.loc[cname, k] if (isinstance(cname, str) and cname in db.index) else np.nan
    f.update(gradient_features(d))
    return f

rows, mols = [], {}
ids = [d for d in survey.index if survey.loc[d, 'mode'] in ('RP', 'HILIC') and d != '0186']
meta = {}
for d in ids:
    f = metadata_features(d); meta[d] = f
    rt = pd.read_csv(f'{R}/{d}/{d}_rtdata_canonical_success.tsv', sep='\t')
    fn = f'{R}/{d}/{d}_rtdata_isomeric_success.tsv'
    if os.path.exists(fn):
        iso = pd.read_csv(fn, sep='\t').set_index('id')['smiles.std']
        rt['smiles.std'] = rt['id'].map(iso).fillna(rt['smiles.std'])
    rt = rt.dropna(subset=['rt', 'smiles.std'])
    n_raw = len(rt)
    rt = rt.drop_duplicates(['smiles.std', 'rt'])
    rt = rt.drop_duplicates('smiles.std', keep=False)
    n_dedup = len(rt)
    void = 2 * f['t0_eff'] if np.isfinite(f['t0_eff']) else -np.inf
    rt = rt[rt.rt >= void]
    for r in rt.itertuples():
        rows.append(dict(dataset=d, cid=r.id, smiles=r._5, inchikey=r._7, ik14=str(r._7)[:14], rt=r.rt))
    meta[d].update(n_raw=n_raw, n_dedup=n_dedup, n_clean=len(rt), void_threshold=void)
df = pd.DataFrame(rows)
print('rows', len(df), 'datasets', df.dataset.nunique())

# molecular features per unique smiles
smi = df.smiles.unique()
feat, fps, scaf = {}, {}, {}
for s in smi:
    m = Chem.MolFromSmiles(s)
    if m is None:
        continue
    feat[s] = [float(fn(m)) for fn in DESC.values()]
    fps[s] = np.frombuffer(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024).ToBitString().encode(), 'u1') - 48
    try:
        scaf[s] = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        scaf[s] = ''
df = df[df.smiles.isin(feat)].copy()
df['scaffold'] = df.smiles.map(scaf)
D = pd.DataFrame.from_dict(feat, orient='index', columns=list(DESC))
D.index.name = 'smiles'
D.to_parquet(f'{OUT}/descriptors.parquet')
np.save(f'{OUT}/ecfp4.npy', np.stack([fps[s] for s in D.index]))
df.to_parquet(f'{OUT}/compounds.parquet')
M = pd.DataFrame.from_dict(meta, orient='index'); M.index.name = 'id'
M['n_clean'] = df.groupby('dataset').size().reindex(M.index).fillna(0).astype(int)
M['mode'] = survey.loc[M.index, 'mode']
M['name'] = survey.loc[M.index, 'name']
M['in_2step_train'] = survey.loc[M.index, 'in_2step_train']
M['n_meta_fields'] = survey.loc[M.index, 'n_meta_fields']
M['hsm_tanaka_available'] = survey.loc[M.index, 'hsm'] & survey.loc[M.index, 'tanaka']

# ---- held-out selection (rules frozen in prereg_interpretation.md, incl. amendment 1) ----
ADD = ['A_formic', 'A_acetic', 'A_nh4ac', 'A_nh4form', 'A_trifluoroacetic', 'A_nh4bicarb', 'A_nh4oh']
M['additive_pattern'] = M[ADD].astype(int).astype(str).agg(''.join, axis=1)
ph_table = M[(M['mode'] == 'RP') & M.pH.notna()].groupby('additive_pattern').pH.median()
M['pH_imputed'] = M.pH.isna() & M.additive_pattern.isin(ph_table.index)
M['pH_2step_input'] = M.pH.fillna(M.additive_pattern.map(ph_table))
ph_table.to_csv(f'{OUT}/ph_imputation_table.csv')

def select(pool, cap):
    out = []
    key = pool.column_name.fillna(pd.Series('NA:' + pool.index, index=pool.index))
    for cname, g in pool.groupby(key):
        out += list(g.sort_values(['n_clean'], ascending=False).index[:cap])
    return sorted(out)
def pick(mask, need):
    pool = M[mask]
    for cap in (2, 3, 4, 5, 10**6):
        sel = select(pool, cap)
        if len(sel) >= need:
            return sel, cap
    return sel, None
unseen_rp = (M['mode'] == 'RP') & (~M.in_2step_train) & M.grad_len.notna()
proxy_mask = unseen_rp & (M.n_clean >= 80)
s2_base = unseen_rp & M.hsm_tanaka_available & M.pH_2step_input.notna()
s2_sel, s2_cap = pick(s2_base & (M.n_clean >= 80), 25)
s2_floor = 80
if s2_cap is None:
    s2_sel, s2_cap = pick(s2_base & (M.n_clean >= 50), 25); s2_floor = 50
hi_mask = (M['mode'] == 'HILIC') & M.grad_len.notna() & (M.n_clean >= 80)
px_sel, px_cap = pick(proxy_mask, 25)
hi_sel, hi_cap = pick(hi_mask, 10)
M['g_rp2step'] = M.index.isin(s2_sel)
M['g_rpproxy'] = M.index.isin(px_sel)
M['g_pub'] = M.index.isin(list(PUB))
M['g_hilic'] = M.index.isin(hi_sel)
M['heldout_any'] = M[['g_rp2step', 'g_rpproxy', 'g_pub', 'g_hilic']].any(axis=1)
M['pub_name'] = M.index.map(PUB)
PUB_GRAD = {'0002': 60, '0003': 20, '0018': 40, '0019': 30, '0054': 6, '0055': 7}
M['grad_len_source'] = np.where(M.grad_len.notna(), 'RepoRT gradient table', '')
for d, L in PUB_GRAD.items():
    if not np.isfinite(M.loc[d, 'grad_len']):
        M.loc[d, 'grad_len'] = L; M.loc[d, 'grad_len_source'] = '2-step README'
print('RP-2step', len(s2_sel), 'cap', s2_cap, 'floor', s2_floor, 'eligible', (s2_base & (M.n_clean >= s2_floor)).sum())
print('RP-proxy', len(px_sel), 'cap', px_cap, 'eligible', proxy_mask.sum())
print('HILIC', len(hi_sel), 'cap', hi_cap, 'eligible', hi_mask.sum())
M.to_csv(f'{OUT}/datasets_meta.csv')
print(M.loc[M.g_rp2step, ['name', 'column_name', 'n_clean', 'grad_len', 'pH', 'pH_2step_input', 't0_source', 'n_meta_fields', 'g_rpproxy']].to_string())
