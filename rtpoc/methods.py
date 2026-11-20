"""Features and methods. Every fit uses context peptides only.

Conventions
- Target y is RT rescaled to [0, 1] with the CONTEXT min/max (scaler fit on
  context only). All methods predict in that space and are mapped back to
  minutes/seconds before metrics are computed.
- g = uncalibrated generic DeepLC prediction (legacy PXD005573 ensemble).
"""
from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.nonparametric.smoothers_lowess import lowess

AA = "ACDEFGHIKLMNPQRSTVWY"
MOD_ALIASES = {"Acetylation": "Acetyl"}
MOD_VOCAB = ["Carbamidomethyl", "Oxidation", "Acetyl", "Phospho", "Methyl", "Dimethyl",
             "Trimethyl", "Formyl", "Propionyl", "Succinyl", "Malonyl", "Crotonyl",
             "Deamidated", "Nitro"]
LOWESS_FRACS = (0.1, 0.2, 0.3, 0.5, 0.8)


# --------------------------------------------------------------------- features
def peptide_features(df: pd.DataFrame) -> pd.DataFrame:
    """Per-peptide features (fixed vocabulary, no data-dependent fitting)."""
    seq = df["seq"].astype(str)
    F = pd.DataFrame(index=df.index)
    F["g"] = df["pred_generic"].astype(float)
    F["length"] = seq.str.len()
    for a in AA:
        F[f"n_{a}"] = seq.str.count(a)
    F["n_KRH"] = F["n_K"] + F["n_R"] + F["n_H"]
    F["n_DE"] = F["n_D"] + F["n_E"]
    F["nterm"] = seq.str[0].map({a: i for i, a in enumerate(AA)}).fillna(-1).astype(int)
    F["cterm"] = seq.str[-1].map({a: i for i, a in enumerate(AA)}).fillna(-1).astype(int)
    counts = {m: np.zeros(len(df), dtype=int) for m in MOD_VOCAB}
    for i, s in enumerate(df["modifications"].astype(str)):
        if s:
            for name in s.split("|")[1::2]:
                name = MOD_ALIASES.get(name, name)
                if name in counts:
                    counts[name][i] += 1
    for m in MOD_VOCAB:
        F[f"mod_{m}"] = counts[m]
    return F


CAT_COLS = ["nterm", "cterm"]


def onehot_terms(F: pd.DataFrame) -> pd.DataFrame:
    """Version for ridge / XGBoost: one-hot N/C-terminal residue (fixed 20-letter vocab)."""
    G = F.drop(columns=CAT_COLS).copy()
    for c in CAT_COLS:
        for i, a in enumerate(AA):
            G[f"{c}_{a}"] = (F[c].values == i).astype(int)
    return G


# --------------------------------------------------------------------- B1 / B2
class Linear:
    def fit(self, g, y):
        self.m = LinearRegression().fit(g[:, None], y)
        return self

    def predict(self, g):
        return self.m.predict(g[:, None])


class Lowess:
    """LOWESS of y ~ g; frac chosen by 5-fold CV on the context (context only).

    Outside the context range of g, extrapolate linearly with the global (B1)
    slope from the curve's end points.
    """

    def __init__(self, frac=None, seed=0):
        self.frac, self.seed = frac, seed

    @staticmethod
    def _fit_curve(g, y, frac):
        delta = 0.005 * (g.max() - g.min()) if len(g) > 2000 else 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = lowess(y, g, frac=frac, it=3, delta=delta, return_sorted=True)
        xs, ys = fit[:, 0], fit[:, 1]
        xs, idx = np.unique(xs, return_index=True)
        ys = ys[idx]
        ok = np.isfinite(ys)
        return xs[ok], ys[ok]

    @staticmethod
    def _eval(xs, ys, slope, g):
        out = np.interp(g, xs, ys)
        lo, hi = g < xs[0], g > xs[-1]
        out[lo] = ys[0] + slope * (g[lo] - xs[0])
        out[hi] = ys[-1] + slope * (g[hi] - xs[-1])
        return out

    def fit(self, g, y):
        self.slope = LinearRegression().fit(g[:, None], y).coef_[0]
        frac = self.frac
        if frac is None:
            kf = KFold(5, shuffle=True, random_state=self.seed)
            best = (np.inf, None)
            for f in LOWESS_FRACS:
                if f * len(g) * 0.8 < 8:  # too few points per local fit
                    continue
                err = 0.0
                for tr, va in kf.split(g):
                    xs, ys = self._fit_curve(g[tr], y[tr], f)
                    sl = LinearRegression().fit(g[tr, None], y[tr]).coef_[0]
                    err += np.abs(self._eval(xs, ys, sl, g[va]) - y[va]).sum()
                if err < best[0]:
                    best = (err, f)
            frac = best[1] if best[1] is not None else LOWESS_FRACS[-1]
        self.frac_ = frac
        self.xs, self.ys = self._fit_curve(g, y, frac)
        return self

    def predict(self, g):
        return self._eval(self.xs, self.ys, self.slope, g)


class DeepLCSpline:
    """DeepLC's own built-in single-series calibration (SplineTransformerCalibration)."""

    def fit(self, g, y):
        from deeplc.calibration.simple import SplineTransformerCalibration
        self.m = SplineTransformerCalibration()
        self.m.fit(target=y.astype(np.float64), source=g.astype(np.float64))
        return self

    def predict(self, g):
        return np.asarray(self.m.transform(g.astype(np.float64)), dtype=float)


def crossfit_b2(g, y, seed, k=5):
    """Out-of-fold B2 predictions on the context (residuals not shrunk by in-sample fit)."""
    oof = np.empty_like(y)
    for tr, va in KFold(k, shuffle=True, random_state=seed).split(g):
        oof[va] = Lowess(seed=seed).fit(g[tr], y[tr]).predict(g[va])
    return oof


# --------------------------------------------------------------------- B3
def ridge_model():
    return make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 13)))


def xgb_model(seed):
    from xgboost import XGBRegressor
    return XGBRegressor(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8,
                        colsample_bytree=0.8, min_child_weight=3, reg_lambda=1.0,
                        n_jobs=1, random_state=seed, verbosity=0)


# --------------------------------------------------------------------- TabPFN
_TABPFN_PATH = None


def tabpfn_model(seed, n_estimators=4):
    import os
    from tabpfn import TabPFNRegressor
    from tabpfn.constants import ModelVersion
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "data", "models", "tabpfn-v2-regressor.ckpt")
    return TabPFNRegressor.create_default_for_version(
        ModelVersion.V2, model_path=path, device="cpu", random_state=seed,
        n_estimators=n_estimators, ignore_pretraining_limits=True,
        categorical_features_indices=None)


def tabpfn_fit_predict(Xc, yc, Xq, seed, cat_idx):
    """Returns (median, q05, q95, mean, seconds)."""
    m = tabpfn_model(seed)
    m.set_params(categorical_features_indices=cat_idx)
    t0 = time.time()
    m.fit(Xc, yc)
    out = m.predict(Xq, output_type="full", quantiles=[0.05, 0.5, 0.95])
    dt = time.time() - t0
    q05, q50, q95 = out["quantiles"]
    return np.asarray(q50), np.asarray(q05), np.asarray(q95), np.asarray(out["mean"]), dt
