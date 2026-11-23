"""Post-hoc control for Step 4: is the T1/R1 gain over C2 study-specific, or a global IM2Deep error a one-time
retrain would remove? G1 = C2 (same anchors, same queries as 07_arms.py) + LightGBM residual model trained once on
all OTHER studies (target: observed - per-study affine-calibrated Z1). No held-out-study data enters G1 beyond the k anchors.
Usage: python -I 14_global_control.py <work_dir> <out_csv> <scripts_dir>
"""
import sys, warnings
W, OUTCSV, SD = sys.argv[1:4]
sys.path.insert(0, SD); warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import lightgbm as lgb
from feats import peptide_features
sv = pd.read_parquet(f'{W}/study_values.parquet').merge(pd.read_parquet(f'{W}/im2deep_pred.parquet'), on='precursor').reset_index(drop=True)
F = peptide_features(sv); F['z1'] = sv['k0_pred'].values; X_all = F.values
y_all = sv.k0.values; z_all = sv.k0_pred.values
# per-study in-sample affine calibration of Z1 for training targets
cal = np.empty(len(sv))
for s, idx in sv.groupby('study').indices.items():
    b = np.linalg.lstsq(np.c_[np.ones(len(idx)), z_all[idx]], y_all[idx], rcond=None)[0]; cal[idx] = b[0] + b[1] * z_all[idx]
resid_all = y_all - cal
rows = []
for s, idx in sorted(sv.groupby('study').indices.items(), key=lambda t: len(t[1])):
    n = len(idx)
    if n < 300: continue
    other = np.setdiff1d(np.arange(len(sv)), idx)
    rng0 = np.random.default_rng(0); tr = rng0.choice(other, min(150_000, len(other)), replace=False)
    g = lgb.LGBMRegressor(n_estimators=400, learning_rate=0.05, num_leaves=31, min_child_samples=50, subsample=0.8, subsample_freq=1,
                          colsample_bytree=0.8, verbose=-1, n_jobs=4).fit(X_all[tr], resid_all[tr])
    y = y_all[idx]; z1 = z_all[idx]; gx = g.predict(X_all[idx])
    for k in (20, 50, 200, 1000):
        if n < k + 200: continue
        for seed in range(5):
            rng = np.random.default_rng(10_000 + 100 * seed + k)       # identical draws to 07_arms.py
            perm = rng.permutation(n); a = perm[:k]; q = perm[k:k + 1500]
            ba = np.linalg.lstsq(np.c_[np.ones(k), z1[a]], y[a], rcond=None)[0]
            c2 = ba[0] + ba[1] * z1[q]
            # G1: global residual added after C2; G1b: C2 refit on anchors after global correction
            g1 = c2 + gx[q]
            zb = z1 + gx; bb = np.linalg.lstsq(np.c_[np.ones(k), zb[a]], y[a], rcond=None)[0]; g1b = bb[0] + bb[1] * zb[q]
            for arm, p in (('C2_check', c2), ('G1_global_resid', g1), ('G1b_global_then_affine', g1b)):
                e = np.abs(p - y[q]); rows.append(dict(study=s, k=k, seed=seed, arm=arm, mae_median=np.median(e)))
    print(s, flush=True)
pd.DataFrame(rows).to_csv(OUTCSV, index=False)
