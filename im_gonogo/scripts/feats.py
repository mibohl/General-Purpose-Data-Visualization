import numpy as np, pandas as pd
AA = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8,
          P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
def peptide_features(df):
    """df needs seq, modseq, charge, mz. Returns feature frame (no target)."""
    seq = df['seq'].values
    L = np.array([len(s) for s in seq], float)
    comp = np.array([[s.count(a) / len(s) for a in AA] for s in seq])
    gravy = np.array([np.mean([KD[c] for c in s]) for s in seq])
    F = pd.DataFrame(comp, columns=[f'aa_{a}' for a in AA], index=df.index)
    F['length'] = L; F['mz'] = df['mz'].values; F['gravy'] = gravy
    F['mass'] = df['mz'].values * df['charge'].values
    for z in (1, 2, 3, 4):
        F[f'z{z}'] = (df['charge'].values == z).astype(float)
    F['z5plus'] = (df['charge'].values >= 5).astype(float)
    F['n_ox'] = df['modseq'].str.count(r'\[ox\]').values
    F['nterm_ac'] = df['modseq'].str.startswith('[ac]').astype(float).values
    F['n_basic'] = np.array([s.count('K') + s.count('R') + s.count('H') for s in seq], float)
    F['c_term_KR'] = np.array([s[-1] in 'KR' for s in seq], float)
    return F
