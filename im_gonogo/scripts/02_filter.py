"""Step 1a: strict PSM filtering and precursor-identity normalization.
Usage: python -I 02_filter.py <work_dir>
Writes psm_filtered.parquet (one row per precursor x run, most intense MS/MS row) and filter_log.csv.
"""
import sys, re
import numpy as np, pandas as pd
W = sys.argv[1]
cols = ['study', 'file', 'Raw file', 'Type', 'Reverse', 'Potential contaminant', 'PEP', 'Sequence',
        'Modified sequence', 'Charge', 'm/z', 'Intensity', '1/K0', 'Calibrated 1/K0', 'CCS',
        'Ion mobility index', 'Retention time']
df = pd.read_parquet(f'{W}/evidence_all.parquet', columns=cols)
log = []
def step(name, mask):
    global df
    n0 = len(df); df = df[mask]
    log.append(dict(step=name, rows_before=n0, rows_removed=n0 - len(df), rows_after=len(df)))
    print(name, n0, '->', len(df), flush=True)

# PXD019086 (Meier et al. 2021) has no 1/K0 column (CCS only, older MaxQuant); the brief reserves it
# as a single-study reference, so it is excluded from the cross-study analysis.
step('drop PXD019086 (no 1/K0 column; single-study reference)', df['study'] != 'PXD019086')
# the truncated accession in the file name resolves to PXD036206 on PRIDE
df['study'] = df['study'].replace({'PXD003206': 'PXD036206'})
# raw files that appear in two projects are kept only in the earlier submission
shared = set(df.loc[df['study'] == 'PXD035987', 'Raw file']) & set(df.loc[df['study'] == 'PXD043511', 'Raw file'])
step(f'drop {len(shared)} raw files of PXD043511 also deposited in PXD035987',
     ~((df['study'] == 'PXD043511') & df['Raw file'].isin(shared)))
step('MS/MS-identified only (drop match-between-runs)', df['Type'].eq('TIMS-MULTI-MSMS'))
step('drop reverse hits', df['Reverse'].fillna('') != '+')
step('drop contaminants', df['Potential contaminant'].fillna('') != '+')
step('PEP <= 0.01', df['PEP'] <= 0.01)
step('1/K0 present and > 0', df['1/K0'] > 0)
# no evidence file has a 'Labeling state' or reporter-intensity column: no SILAC/TMT searches
step('charge 1..6', df['Charge'].between(1, 6))

# ---- modification normalization ----
ms = df['Modified sequence'].astype(str).str.strip('_')
CANON = {'Oxidation (M)': 'ox', 'ox': 'ox', 'Acetyl (Protein N-term)': 'ac', 'ac': 'ac'}
tok_re = re.compile(r'\(((?:[^()]|\([^()]*\))*)\)')
def norm(s):
    toks = tok_re.findall(s)
    bad = [t for t in toks if t not in CANON]
    if bad:
        return None, bad[0]
    out = tok_re.sub(lambda m: '[' + CANON[m.group(1)] + ']', s)
    return out, None
uniq = ms.unique()
mp = {u: norm(u) for u in uniq}
df['modseq'] = ms.map(lambda u: mp[u][0])
df['badmod'] = ms.map(lambda u: mp[u][1])
bad_counts = df['badmod'].value_counts()
step('only unmodified / Oxidation(M) / Acetyl(Protein N-term)', df['modseq'].notna())
# sanity: stripped modseq must equal Sequence
chk = df['modseq'].str.replace(r'\[[a-z]+\]', '', regex=True)
step('modified sequence consistent with Sequence', chk == df['Sequence'])
# oxidation only on M, acetyl only at N-term
step('ox only after M', ~df['modseq'].str.contains(r'[^M]\[ox\]'))
step('ac only at N-terminus', ~df['modseq'].str.contains(r'.\[ac\]') | df['modseq'].str.match(r'^\[ac\]') & (df['modseq'].str.count(r'\[ac\]') == 1))
# cysteine: fixed alkylation is invisible in MaxQuant modified sequences -> drop
step('drop cysteine-containing peptides (fixed alkylation not verifiable)', ~df['Sequence'].str.contains('C'))
df['precursor'] = df['modseq'] + '/' + df['Charge'].astype(int).astype(str)

# one value per precursor x run: most intense MS/MS-identified row; record conformer spread
df['Intensity'] = df['Intensity'].fillna(0)
g = df.groupby(['study', 'Raw file', 'precursor'])
spread = g['1/K0'].agg(k0_min='min', k0_max='max', n_rows='size')
df = df.sort_values('Intensity', ascending=False).drop_duplicates(['study', 'Raw file', 'precursor'])
df = df.merge(spread.reset_index(), on=['study', 'Raw file', 'precursor'])
log.append(dict(step='collapse to one row per precursor x run (most intense)', rows_after=len(df)))
df = df[['study', 'file', 'Raw file', 'precursor', 'modseq', 'Sequence', 'Charge', 'm/z', 'Intensity', '1/K0',
         'Calibrated 1/K0', 'CCS', 'Ion mobility index', 'Retention time', 'k0_min', 'k0_max', 'n_rows']]
df.to_parquet(f'{W}/psm_filtered.parquet')
pd.DataFrame(log).to_csv(f'{W}/filter_log.csv', index=False)
bad_counts.to_csv(f'{W}/rejected_mod_tokens.csv')
print(bad_counts.head(20))
