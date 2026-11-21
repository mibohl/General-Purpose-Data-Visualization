"""Cheap, label-light baselines evaluated on the same k support peptides:
PSSM/PWM log-odds, kNN on BLOSUM similarity, ridge, logistic, and
recalibration maps for MHCflurry.
"""
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.isotonic import IsotonicRegression

_AA = "ARNDCQEGHILKMFPSTWYV"
_AA_IDX = {a: i for i, a in enumerate(_AA)}

# Background amino-acid frequencies (UniProt-ish), for PSSM pseudocounts/log-odds.
_BG = np.array([0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026,
                0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013,
                0.032, 0.073])
_BG = _BG / _BG.sum()


def _onehot_counts(peptides, length=9):
    """Per-position AA counts (length x 20)."""
    counts = np.zeros((length, 20), dtype=np.float64)
    for pep in peptides:
        for j, aa in enumerate(pep[:length]):
            idx = _AA_IDX.get(aa)
            if idx is not None:
                counts[j, idx] += 1.0
    return counts


def pwm_logodds_scores(support_peptides, support_binder, query_peptides,
                       length=9, beta=50.0):
    """Classic PSSM/PWM built from the BINDERS among the support set, scored by
    summed per-position log-odds vs background. This is the strong cheap
    baseline.  Pseudocount weight beta mixes background into sparse columns.

    Returns a score per query peptide (higher = more binder-like), or None if
    fewer than 2 binders in the support (PWM undefined).
    """
    binders = [p for p, b in zip(support_peptides, support_binder) if b == 1]
    if len(binders) < 2:
        return None
    counts = _onehot_counts(binders, length)              # L x 20
    n = len(binders)
    # pseudocount-smoothed position frequencies (Henikoff-style simple mix)
    freq = (counts + beta * _BG[None, :]) / (n + beta)
    logodds = np.log(freq / _BG[None, :])                 # L x 20
    scores = np.zeros(len(query_peptides))
    for i, pep in enumerate(query_peptides):
        s = 0.0
        for j, aa in enumerate(pep[:length]):
            idx = _AA_IDX.get(aa)
            if idx is not None:
                s += logodds[j, idx]
        scores[i] = s
    return scores


def knn_blosum_predict(X_sup, y_sup, X_qry, k_neigh=5):
    """kNN regression on BLOSUM-encoded peptides using Euclidean distance in
    BLOSUM space (equivalent ranking to BLOSUM similarity). Distance-weighted
    mean of neighbour 1-log50k values."""
    k_neigh = min(k_neigh, len(X_sup))
    # pairwise squared euclidean
    # ||q-s||^2 = |q|^2 + |s|^2 - 2 q.s
    sq_s = np.sum(X_sup**2, axis=1)
    preds = np.empty(len(X_qry))
    for i in range(len(X_qry)):
        d = sq_s - 2.0 * X_sup.dot(X_qry[i]) + np.sum(X_qry[i]**2)
        nn = np.argpartition(d, k_neigh - 1)[:k_neigh]
        w = 1.0 / (np.sqrt(np.maximum(d[nn], 0)) + 1e-6)
        preds[i] = np.sum(w * y_sup[nn]) / np.sum(w)
    return preds


def ridge_predict(X_sup, y_sup, X_qry, alpha=10.0):
    m = Ridge(alpha=alpha)
    m.fit(X_sup, y_sup)
    return m.predict(X_qry)


def logistic_predict(X_sup, yb_sup, X_qry, C=1.0):
    """Logistic regression on binder label. If support is single-class, return
    a constant equal to that class prob."""
    if len(np.unique(yb_sup)) < 2:
        return np.full(len(X_qry), float(yb_sup[0]))
    m = LogisticRegression(C=C, max_iter=2000)
    m.fit(X_sup, yb_sup)
    return m.predict_proba(X_qry)[:, 1]


def recalibrate(base_sup, y_sup, base_qry, method="isotonic"):
    """Fit a monotone/linear map from base predictions -> observed 1-log50k on
    the k support peptides, apply to query. Uses zero *extra* info beyond the k
    labels already available to every other k-method."""
    if method == "isotonic":
        ir = IsotonicRegression(out_of_bounds="clip")
        ir.fit(base_sup, y_sup)
        return ir.predict(base_qry)
    else:  # linear
        a = np.polyfit(base_sup, y_sup, 1)
        return np.polyval(a, base_qry)


def metrics(score, y_cont, y_bin):
    """AUROC & AUPRC vs binder label; Spearman vs 1-log50k. NaN where undefined."""
    out = {"auroc": np.nan, "auprc": np.nan, "spearman": np.nan}
    if score is None or len(score) == 0:
        return out
    if len(np.unique(y_bin)) == 2:
        out["auroc"] = roc_auc_score(y_bin, score)
        out["auprc"] = average_precision_score(y_bin, score)
    if np.std(score) > 0 and np.std(y_cont) > 0:
        out["spearman"] = spearmanr(score, y_cont).statistic
    return out
