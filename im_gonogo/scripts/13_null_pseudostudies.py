"""Noise-only null for A: split one study's runs (single TIMS method) into 6 pseudo-studies and rerun D/A/O.
Usage: python -I 13_null_pseudostudies.py <work_dir> <out_dir>
"""
import sys
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
W, OUT = sys.argv[1], sys.argv[2]
psm = pd.read_parquet(f'{W}/psm_grouped.parquet', columns=['study', 'mgroup', 'Raw file', 'precursor', '1/K0'])
rows = []
for mg in psm.groupby('mgroup')['Raw file'].nunique().loc[lambda s: s >= 24].index:
    g0 = psm[psm.mgroup == mg]
    runs = np.array(sorted(g0['Raw file'].unique()))
    for split_seed in range(3):
        rng = np.random.default_rng(split_seed)
        lab = dict(zip(rng.permutation(runs), np.arange(len(runs)) % 6))
        g = g0.assign(ps=g0['Raw file'].map(lab))
        sv = g.groupby(['ps', 'precursor'])['1/K0'].median().reset_index()
        sv['n'] = sv.groupby('precursor').ps.transform('size'); sv = sv[sv.n >= 4].reset_index(drop=True)
        cons = np.empty(len(sv))
        for p, idx in sv.groupby('precursor').indices.items():
            v = sv['1/K0'].values[idx]
            for j, i in enumerate(idx): cons[i] = np.median(np.delete(v, j))
        sv['c'] = cons
        devD, devA, devO = [], [], []
        for ps, h in sv.groupby('ps'):
            x, c = h['1/K0'].values, h.c.values; n = len(x)
            devD.append(np.abs(x - c)); devO.append(np.abs(IsotonicRegression(out_of_bounds='clip').fit(x, c).predict(x) - c))
            for seed in range(5):
                r = np.random.default_rng(1000 * seed + 20); a = r.choice(n, 20, replace=False); q = np.setdiff1d(np.arange(n), a)
                b = np.linalg.lstsq(np.c_[np.ones(20), x[a]], c[a], rcond=None)[0]; devA.append(np.abs(b[0] + b[1] * x[q] - c[q]))
        # replicate floor in this group (same definition as Step 2)
        rows.append(dict(mgroup=mg, n_runs=len(runs), split_seed=split_seed, n_shared=sv.precursor.nunique(),
                         D_null=np.median(np.concatenate(devD)), A_null=np.median(np.concatenate(devA)), O_null=np.median(np.concatenate(devO))))
        print(rows[-1], flush=True)
R = pd.DataFrame(rows)
Fs = pd.read_csv(f'{OUT}/F_per_study.csv'); Fs = Fs[Fs.space == '1/K0'].set_index('study').F
R['F_study'] = R.mgroup.str.split('|').str[0].map(Fs)
R['A_null_over_F_study'] = R.A_null / R.F_study
R.to_csv(f'{OUT}/null_pseudostudies.csv', index=False)
print(R.groupby('mgroup')[['n_runs', 'n_shared', 'F_study', 'D_null', 'A_null', 'O_null', 'A_null_over_F_study']].median().round(4).to_string())
