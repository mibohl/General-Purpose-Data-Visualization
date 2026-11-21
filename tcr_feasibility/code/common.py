"""
Common utilities: amino-acid featurization (Atchley factors), IMMREP22 data
loading, and sequence featurization.

The epitope is constant within a task, so the within-task feature vector encodes
the TCR only. Epitope features are added only for the zero-shot transfer model.
"""
import os
import re
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Atchley factors (Atchley et al. PNAS 2005): 5 numeric descriptors per residue.
# ---------------------------------------------------------------------------
ATCHLEY = {
    "A": (-0.591, -1.302, -0.733,  1.570, -0.146),
    "C": (-1.343,  0.465, -0.862, -1.020, -0.255),
    "D": ( 1.050,  0.302, -3.656, -0.259, -3.242),
    "E": ( 1.357, -1.453,  1.477,  0.113, -0.837),
    "F": (-1.006, -0.590,  1.891, -0.397,  0.412),
    "G": (-0.384,  1.652,  1.330,  1.045,  2.064),
    "H": ( 0.336, -0.417, -1.673, -1.474, -0.078),
    "I": (-1.239, -0.547,  2.131,  0.393,  0.816),
    "K": ( 1.831, -0.561,  0.533, -0.277,  1.648),
    "L": (-1.019, -0.987, -1.505,  1.266, -0.912),
    "M": (-0.663, -1.524,  2.219, -1.005,  1.212),
    "N": ( 0.945,  0.828,  1.299, -0.169,  0.933),
    "P": ( 0.189,  2.081, -1.628,  0.421, -1.392),
    "Q": ( 0.931, -0.179, -3.005, -0.503, -1.853),
    "R": ( 1.538, -0.055,  1.502,  0.440,  2.897),
    "S": (-0.228,  1.399, -4.760,  0.670, -2.647),
    "T": (-0.032,  0.326,  2.213,  0.908,  1.313),
    "V": (-1.337, -0.279, -0.544,  1.242, -1.262),
    "W": (-0.595,  0.009,  0.672, -2.128, -0.184),
    "Y": ( 0.260,  0.830,  3.097, -0.838,  1.512),
}
AAS = set(ATCHLEY.keys())
MAXLEN = 20  # CDR3 aligned/padded length -> 20 * 5 = 100 Atchley features

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMMREP_DIR = os.path.join(REPO_ROOT, "data", "immrep22")


def valid_cdr3(s, lo=5, hi=MAXLEN):
    """A CDR3 is usable if it is all-standard-AA and within length bounds."""
    if not isinstance(s, str):
        return False
    s = s.strip().upper()
    if len(s) < lo or len(s) > hi:
        return False
    return all(c in AAS for c in s)


def atchley_vector(seq, maxlen=MAXLEN):
    """
    Encode a CDR3 as a fixed maxlen*5 vector. Sequence is anchored at BOTH
    termini (the conserved C... / ...F/W ends), zero-padding the centre. This
    mirrors the IMGT-style alignment used in the TCR literature and keeps the
    conserved flanks in register across lengths.
    """
    seq = seq.strip().upper()
    mat = np.zeros((maxlen, 5), dtype=np.float32)
    L = len(seq)
    if L >= maxlen:
        # center-truncate: keep first half and last half
        first = maxlen // 2
        last = maxlen - first
        kept = seq[:first] + seq[L - last:]
        for i, c in enumerate(kept):
            mat[i] = ATCHLEY[c]
    else:
        first = (L + 1) // 2           # ceil
        last = L - first               # floor
        for i in range(first):
            mat[i] = ATCHLEY[seq[i]]
        for j in range(last):
            mat[maxlen - last + j] = ATCHLEY[seq[first + j]]
    return mat.reshape(-1)


def atchley_matrix(seqs, maxlen=MAXLEN):
    return np.vstack([atchley_vector(s, maxlen) for s in seqs])


class TargetEncoder:
    """Smoothed target encoder fit on SUPPORT labels only (no query leakage)."""

    def __init__(self, smoothing=10.0):
        self.smoothing = smoothing
        self.map_ = {}
        self.prior_ = 0.5

    def fit(self, cats, y):
        y = np.asarray(y, dtype=float)
        self.prior_ = float(y.mean())
        df = pd.DataFrame({"c": list(cats), "y": y})
        g = df.groupby("c")["y"]
        n = g.count()
        m = g.mean()
        self.map_ = ((m * n + self.prior_ * self.smoothing) /
                     (n + self.smoothing)).to_dict()
        return self

    def transform(self, cats):
        return np.array([self.map_.get(c, self.prior_) for c in cats],
                        dtype=np.float32)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
EPITOPE_FILES_IGNORE = {"README"}


def load_immrep_positives():
    """
    Pool all Label==1 TCRs across training_data + true_set for every epitope.
    Returns a DataFrame with columns:
        epitope, cdr3b, vb, jb, cdr3a, va, ja, source
    The within-task protocol re-splits each epitope into support/query itself,
    so the benchmark's 80:20 split is not used for task construction.
    """
    rows = []
    train_dir = os.path.join(IMMREP_DIR, "training_data")
    true_dir = os.path.join(IMMREP_DIR, "true_set")
    epis = sorted(f[:-4] for f in os.listdir(train_dir)
                  if f.endswith(".txt") and f[:-4] not in EPITOPE_FILES_IGNORE)
    for ep in epis:
        for src, d in (("train", train_dir), ("true", true_dir)):
            fp = os.path.join(d, ep + ".txt")
            if not os.path.exists(fp):
                continue
            df = pd.read_csv(fp, sep="\t", dtype=str)
            if "Label" not in df.columns:
                continue
            pos = df[df["Label"] == "1"]
            for _, r in pos.iterrows():
                rows.append({
                    "epitope": ep,
                    "cdr3b": str(r.get("TRB_CDR3", "")).strip().upper(),
                    "vb": str(r.get("TRBV", "")).strip(),
                    "jb": str(r.get("TRBJ", "")).strip(),
                    "cdr3a": str(r.get("TRA_CDR3", "")).strip().upper(),
                    "va": str(r.get("TRAV", "")).strip(),
                    "ja": str(r.get("TRAJ", "")).strip(),
                    "source": src,
                })
    return pd.DataFrame(rows)


def load_immrep_controls():
    """Background negatives: unpublished healthy-donor CD8 TCRs (A. Eugster)."""
    rows = []
    cdir = os.path.join(IMMREP_DIR, "negativeData_controls")
    for f in sorted(os.listdir(cdir)):
        if not f.endswith(".tsv"):
            continue
        df = pd.read_csv(os.path.join(cdir, f), sep="\t", dtype=str)
        for _, r in df.iterrows():
            rows.append({
                "epitope": "__control__",
                "cdr3b": str(r.get("TRB_CDR3", "")).strip().upper(),
                "vb": str(r.get("TRBV", "")).strip(),
                "jb": str(r.get("TRBJ", "")).strip(),
                "cdr3a": str(r.get("TRA_CDR3", "")).strip().upper(),
                "va": str(r.get("TRAV", "")).strip(),
                "ja": str(r.get("TRAJ", "")).strip(),
                "source": os.path.splitext(f)[0],
            })
    return pd.DataFrame(rows)
