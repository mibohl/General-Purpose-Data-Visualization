"""Warp / in-context arms. Every arm: fit on anchors (roi, rt, X) and predict for all compounds."""
import sys, numpy as np, pandas as pd, warnings
from scipy.optimize import lsq_linear
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, '/home/user/data/2-step')
import logging
from mapping import LADModel  # the released 2-step mapping code, unmodified
warnings.filterwarnings('ignore')
logging.getLogger().setLevel(logging.ERROR)


def i1_quadratic(roi_a, rt_a, roi_all):
    """Incumbent: 2-step LADModel as configured in predict.py."""
    m = LADModel(pd.DataFrame({'roi': roi_a, 'rt': rt_a}), ols_after=True,
                 ols_discard_if_negative=True, ols_drop_mode='2*median')
    return m.get_mapping(np.asarray(roi_all, dtype=float))


def i2_isotonic(roi_a, rt_a, roi_all):
    return IsotonicRegression(out_of_bounds='clip').fit(roi_a, rt_a).predict(roi_all)


def i3_linear(roi_a, rt_a, roi_all):
    return LinearRegression().fit(roi_a[:, None], rt_a).predict(np.asarray(roi_all)[:, None])


def i4_monotone_spline(roi_a, rt_a, roi_all):
    """Monotone piecewise-linear regression spline, knots at anchor ROI quantiles,
    non-negative segment slopes, linear extrapolation with the boundary slopes."""
    k = len(roi_a)
    nseg = int(np.clip(k // 8, 2, 5))
    knots = np.unique(np.quantile(roi_a, np.linspace(0, 1, nseg + 1)))
    if len(knots) < 3:
        return i3_linear(roi_a, rt_a, roi_all)
    def basis(x):
        x = np.asarray(x, float)
        cols = [np.ones_like(x)]
        inner = knots[1:-1]
        # first segment extends to -inf, last to +inf
        cols.append(np.minimum(x - knots[1], 0.0) if len(inner) else x - knots[0])  # slope of first segment
        for j in range(1, len(knots) - 2):
            cols.append(np.clip(x - knots[j], 0, knots[j + 1] - knots[j]))
        cols.append(np.maximum(x - knots[-2], 0.0))
        return np.stack(cols, 1)
    B = basis(roi_a)
    lb = np.r_[-np.inf, np.zeros(B.shape[1] - 1)]
    sol = lsq_linear(B, rt_a, bounds=(lb, np.full(B.shape[1], np.inf)))
    return basis(roi_all) @ sol.x


def isotonic_l1(x, y):
    """MAE-optimal monotone (non-decreasing) fit of y on x; tied x share one value. Returns fitted y."""
    order = np.argsort(x, kind='mergesort')
    xs, ys = np.asarray(x)[order], np.asarray(y)[order]
    # group ties
    starts = np.r_[0, np.where(np.diff(xs) != 0)[0] + 1]
    groups = np.split(ys, starts[1:])
    blocks = []  # list of [values array, median, count_groups]
    for g in groups:
        blocks.append([g, np.median(g), 1])
        while len(blocks) > 1 and blocks[-2][1] > blocks[-1][1]:
            v = np.concatenate([blocks[-2][0], blocks[-1][0]])
            n = blocks[-2][2] + blocks[-1][2]
            blocks.pop(); blocks[-1] = [v, np.median(v), n]
    fitted_groups = []
    for b in blocks:
        fitted_groups += [b[1]] * b[2]
    fit_sorted = np.repeat(fitted_groups, [len(g) for g in groups])
    out = np.empty_like(fit_sorted, dtype=float)
    out[order] = fit_sorted
    return out


def isotonic_l1_crossfit(x, y, folds=10, seed=0):
    """Sensitivity variant of O1: each compound predicted by an L1 isotonic map fit on the other folds
    (step function, linear interpolation between fitted points, flat beyond the ends)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    rng = np.random.default_rng(seed); f = rng.permutation(len(x)) % folds
    out = np.empty(len(x))
    for j in range(folds):
        tr, te = f != j, f == j
        fit = isotonic_l1(x[tr], y[tr])
        o = np.argsort(x[tr], kind='mergesort')
        xs, fs = x[tr][o], fit[o]
        xu, idx = np.unique(xs, return_index=True)
        out[te] = np.interp(x[te], xu, fs[idx])
    return out


def isotonic_l2(x, y):
    return IsotonicRegression().fit(x, y).predict(x)


def r1_ridge_residual(roi_a, rt_a, Xd_a, i1_a, roi_all, Xd_all, i1_all):
    F = lambda r, X: np.column_stack([r, r ** 2, X])
    sc = StandardScaler().fit(F(roi_a, Xd_a))
    Fa = np.nan_to_num(sc.transform(F(roi_a, Xd_a)))
    Fall = np.nan_to_num(sc.transform(F(roi_all, Xd_all)))
    m = RidgeCV(alphas=np.logspace(-2, 3, 11)).fit(Fa, rt_a - i1_a)
    return i1_all + m.predict(Fall)


_TABPFN = None
def tabpfn_model(seed):
    global _TABPFN
    if _TABPFN is None:
        from tabpfn import TabPFNRegressor
        from tabpfn.constants import ModelVersion
        _TABPFN = TabPFNRegressor.create_default_for_version(
            ModelVersion.V2, model_path='/home/user/data/tabpfn/tabpfn-v2-regressor.ckpt',
            device='cpu', n_estimators=8)
    _TABPFN.random_state = seed
    return _TABPFN


def tabpfn_arm(Fa, rt_a, Fall, seed):
    m = tabpfn_model(seed)
    m.fit(Fa, rt_a)
    return m.predict(Fall)
