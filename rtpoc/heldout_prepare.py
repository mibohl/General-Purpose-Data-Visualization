#!/usr/bin/env python3
"""Convert held-out PRIDE identification tables to DeepLC-style RT CSVs.

Output: data/heldout/<PXD>.csv with columns seq, modifications, tr, in the same
format as the Zenodo CSVs in data/zenodo/train_eval/datasets/:
  - modifications: PEPREC "pos|Name|pos|Name" with Unimod names, position 0 =
    N-terminus, 1-based residue positions, -1 = C-terminus, "" if unmodified.
    Fixed Carbamidomethyl (and fixed TMTpro on K) is written explicitly.
  - tr: retention time in minutes, the median over all identifications of that
    peptidoform within ONE run (the run with the most unique peptidoforms at
    1% FDR), so no cross-run alignment is involved.

Input files are untrusted downloads in data/heldout/raw/<PXD>/; they are only
parsed as tables. Run with `python3 -I rtpoc/heldout_prepare.py [PXD ...]`
(add --check to verify every mapped peptidoform's mass with psm_utils against
the mass reported by the search engine).

One reader per source format: read_maxquant, read_fragpipe, read_diann,
read_peaks, read_pd. Each returns a DataFrame with columns
run, seq, mods (tuple of (pos, name) or None if unmappable), tr (min),
key (the engine's own peptidoform string, used to count identifications per
run) and mass (engine-reported neutral peptide mass, for --check).
"""
import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HELD = ROOT / "data" / "heldout"
RAW = HELD / "raw"
PROTON = 1.007276467
AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

DATASETS = {
    # easy
    "PXD076643": dict(fmt="maxquant", file="PXD076643/220822_ML_Cd_evidence.txt", fixed_cam=True),
    "PXD079751": dict(fmt="diann", file="PXD079751/report.parquet"),
    "PXD062983": dict(fmt="diann", file="PXD062983/report.parquet"),
    # hard
    "PXD050641": dict(fmt="maxquant", file="PXD050641/txt/evidence.txt", fixed_cam=True),
    "PXD074962": dict(fmt="fragpipe", file="PXD074962/P4096_Set1_PP_psm.tsv"),
    "PXD061731": dict(fmt="peaks", file="PXD061731/db.psms.csv"),
    "PXD084949": dict(fmt="pd", file="PXD084949/070224_hMDM_pathMHCI_WT_DDA_PSMs.txt"),
    "PXD085246": dict(fmt="maxquant", file="PXD085246/evidence.txt", fixed_cam=False),
}

# ---------------------------------------------------------------- mod tables
# MaxQuant modification names (long form, MQ >= 2) and abbreviations (MQ 1.x).
MQ_MODS = {
    "Acetyl (Protein N-term)": "Acetyl", "Acetyl (N-term)": "Acetyl", "Acetyl (K)": "Acetyl", "ac": "Acetyl",
    "Oxidation (M)": "Oxidation", "ox": "Oxidation",
    "Phospho (STY)": "Phospho", "ph": "Phospho",
    "GlyGly (K)": "GG", "gl": "GG",
    "Deamidation (NQ)": "Deamidated", "de": "Deamidated",
    "Crotony(K)": "Crotonyl", "Crotonyl (K)": "Crotonyl", "cr": "Crotonyl",
    "Carbamidomethyl (C)": "Carbamidomethyl", "cam": "Carbamidomethyl",
    "Gln->pyro-Glu": "Gln->pyro-Glu", "Glu->pyro-Glu": "Glu->pyro-Glu",
    "Methyl (KR)": "Methyl", "Dimethyl (KR)": "Dimethyl", "Trimethyl (K)": "Trimethyl",
}
# DIA-NN "(UniMod:<id>)" -> Unimod name.
UNIMOD_ID = {
    1: "Acetyl", 4: "Carbamidomethyl", 5: "Carbamyl", 7: "Deamidated", 21: "Phospho",
    26: "Pyro-carbamidomethyl", 27: "Glu->pyro-Glu", 28: "Gln->pyro-Glu", 34: "Methyl",
    35: "Oxidation", 36: "Dimethyl", 37: "Trimethyl", 121: "GG", 312: "Cysteinyl",
    385: "Ammonia-loss", 737: "TMT6plex", 1363: "Crotonyl", 2016: "TMTpro",
}
# Monoisotopic mass deltas -> Unimod name (FragPipe, PEAKS).
MASS_MODS = [
    (57.021464, "Carbamidomethyl"), (15.994915, "Oxidation"), (79.966331, "Phospho"),
    (304.207146, "TMTpro"), (229.162932, "TMT6plex"), (42.010565, "Acetyl"),
    (0.984016, "Deamidated"), (14.01565, "Methyl"), (28.0313, "Dimethyl"),
    (42.04695, "Trimethyl"), (114.042927, "GG"), (68.026215, "Crotonyl"),
    (119.004099, "Cysteinyl"), (43.005814, "Carbamyl"), (39.994915, "Pyro-carbamidomethyl"),
]
# PD modification names that are already Unimod names.
PD_MODS = {"Oxidation", "Carbamidomethyl", "Phospho", "Acetyl", "Deamidated", "GG", "TMTpro",
           "TMT6plex", "Methyl", "Dimethyl", "Trimethyl", "Cysteinyl", "Gln->pyro-Glu",
           "Glu->pyro-Glu", "Carbamyl"}


def mass_to_name(delta, residue, tol):
    if residue == "Q" and abs(delta + 17.026549) < tol:
        return "Gln->pyro-Glu"
    if residue == "E" and abs(delta + 18.010565) < tol:
        return "Glu->pyro-Glu"
    if residue == "C" and abs(delta + 17.026549) < tol:
        return "Ammonia-loss"
    hits = [n for m, n in MASS_MODS if abs(m - delta) < tol]
    return hits[0] if len(hits) == 1 else None


def add_fixed_cam(seq, mods):
    if mods is None:
        return None
    have = {p for p, _ in mods}
    return tuple(mods) + tuple((i + 1, "Carbamidomethyl") for i, a in enumerate(seq)
                               if a == "C" and i + 1 not in have)


def split_mod_string(s, open_="(", close=")"):
    """Split 'AB(x)C(y(z))' into residues and bracketed groups; yields
    (n_residues_before_group, group_text) and returns stripped sequence."""
    seq, groups, i = [], [], 0
    while i < len(s):
        c = s[i]
        if c == open_:
            depth, j = 1, i + 1
            while depth:
                depth += {open_: 1, close: -1}.get(s[j], 0)
                j += 1
            groups.append((len(seq), s[i + 1:j - 1]))
            i = j
        else:
            seq.append(c)
            i += 1
    return "".join(seq), groups


# ------------------------------------------------------------------ readers
def read_maxquant(path, fixed_cam):
    """MaxQuant evidence.txt. Keeps MS/MS-identified rows only (no MBR
    transfers, no second peptides), drops reverse hits; evidence.txt is
    already filtered at 1% PSM FDR. tr = 'Retention time' (min, per run)."""
    cols = ["Sequence", "Modified sequence", "Raw file", "Type", "Retention time", "Mass", "Reverse"]
    e = pd.read_csv(path, sep="\t", usecols=lambda c: c in cols, low_memory=False)
    e = e[e["Type"].str.endswith("MSMS") & (e["Reverse"].fillna("") != "+")]

    def parse(ms):
        body = ms.strip("_")
        cterm = None
        if "_(" in body:  # C-terminal modification after trailing underscore
            body, cterm = body.split("_(", 1)
            cterm = cterm.rstrip(")")
        seq, groups = split_mod_string(body)
        mods = []
        for n, g in groups:
            name = MQ_MODS.get(g)
            if name is None:
                return seq, None
            mods.append((0 if n == 0 else n, name))
        if cterm is not None:
            if cterm not in MQ_MODS:
                return seq, None
            mods.append((-1, MQ_MODS[cterm]))
        return seq, tuple(mods)

    parsed = {ms: parse(ms) for ms in e["Modified sequence"].unique()}
    out = pd.DataFrame({
        "run": e["Raw file"].values,
        "key": e["Modified sequence"].values,
        "seq": [parsed[m][0] for m in e["Modified sequence"]],
        "mods": [parsed[m][1] for m in e["Modified sequence"]],
        "tr": e["Retention time"].values.astype(float),
        "mass": e["Mass"].values.astype(float),
    })
    if fixed_cam:
        out["mods"] = [add_fixed_cam(s, m) for s, m in zip(out.seq, out.mods)]
    return out


def read_fragpipe(path):
    """FragPipe psm.tsv (already 1% PSM/peptide/protein FDR). Mods come from
    'Assigned Modifications' (includes fixed ones). tr = Retention / 60."""
    cols = ["Spectrum", "Peptide", "Modified Peptide", "Retention", "Assigned Modifications",
            "Calculated Peptide Mass", "Is Decoy"]
    p = pd.read_csv(path, sep="\t", usecols=lambda c: c in cols, low_memory=False)
    if "Is Decoy" in p:
        p = p[~p["Is Decoy"].astype(str).str.lower().eq("true")]

    def parse(seq, am):
        mods = []
        for tok in str(am).split(","):
            tok = tok.strip()
            if not tok or tok == "nan":
                continue
            m = re.match(r"^(N-term|C-term|n-term|c-term|(\d+)([A-Z]))\((-?[\d.]+)\)$", tok)
            if not m:
                return None
            delta = float(m.group(4))
            if m.group(2):
                pos, res = int(m.group(2)), m.group(3)
            else:
                pos, res = (0, seq[0]) if m.group(1).lower() == "n-term" else (-1, seq[-1])
            name = mass_to_name(delta, res, 0.002)
            if name is None:
                return None
            mods.append((pos, name))
        return tuple(mods)

    keys = p["Peptide"] + "|" + p["Assigned Modifications"].fillna("")
    parsed = {k: parse(*k.split("|", 1)) for k in keys.unique()}
    return pd.DataFrame({
        "run": p["Spectrum"].str.rsplit(".", n=3).str[0].values,
        "key": keys.values,
        "seq": p["Peptide"].values,
        "mods": [parsed[k] for k in keys],
        "tr": p["Retention"].values.astype(float) / 60.0,
        "mass": p["Calculated Peptide Mass"].values.astype(float),
    })


def read_diann(path):
    """DIA-NN report (.parquet or .tsv). Precursors at run-level Q.Value <= 0.01
    and Lib.Q.Value <= 0.01 (DIA-NN's recommended filter with MBR), no decoys.
    tr = RT (min, apex in that run)."""
    cols = ["Run", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge", "Precursor.Mz",
            "RT", "Q.Value", "Lib.Q.Value", "Decoy"]
    if str(path).endswith(".parquet"):
        import pyarrow.parquet as pq
        have = pq.ParquetFile(path).schema.names
        r = pd.read_parquet(path, columns=[c for c in cols if c in have])
    else:
        r = pd.read_csv(path, sep="\t", usecols=lambda c: c in cols)
    keep = r["Q.Value"] <= 0.01
    if "Lib.Q.Value" in r:
        keep &= r["Lib.Q.Value"] <= 0.01
    if "Decoy" in r:
        keep &= r["Decoy"] == 0
    r = r[keep]

    def parse(ms):
        seq, groups = split_mod_string(ms)
        mods = []
        for n, g in groups:
            m = re.match(r"^UniMod:(\d+)$", g)
            if not m or int(m.group(1)) not in UNIMOD_ID:
                return seq, None
            mods.append((n, UNIMOD_ID[int(m.group(1))]))
        return seq, tuple(mods)

    parsed = {ms: parse(ms) for ms in r["Modified.Sequence"].unique()}
    return pd.DataFrame({
        "run": r["Run"].values,
        "key": r["Modified.Sequence"].values,
        "seq": r["Stripped.Sequence"].values,
        "mods": [parsed[m][1] for m in r["Modified.Sequence"]],
        "tr": r["RT"].values.astype(float),
        "mass": (r["Precursor.Mz"].values - PROTON) * r["Precursor.Charge"].values,
    })


def read_peaks(path):
    """PEAKS db.psms.csv (exported at 1% peptide FDR). Mods are mass deltas
    with 2 decimals in 'Peptide', e.g. AAK(+42.01)QPP; a leading group is
    N-terminal. tr = RT (min)."""
    d = pd.read_csv(path)

    def parse(pep):
        body = re.sub(r"^[A-Z-]\.|\.[A-Z-]$", "", pep)  # flanking residues, if any
        seq, groups = split_mod_string(body)
        mods = []
        for n, g in groups:
            try:
                delta = float(g)
            except ValueError:
                return seq, None
            res = seq[n - 1] if n > 0 else seq[0]
            name = mass_to_name(delta, res, 0.006)
            if name is None:
                return seq, None
            mods.append((n, name))
        return seq, tuple(mods)

    parsed = {p: parse(p) for p in d["Peptide"].unique()}
    return pd.DataFrame({
        "run": d["Source File"].values,
        "key": d["Peptide"].values,
        "seq": [parsed[p][0] for p in d["Peptide"]],
        "mods": [parsed[p][1] for p in d["Peptide"]],
        "tr": d["RT"].values.astype(float),
        "mass": d["Mass"].values.astype(float),
    })


def read_pd(path):
    """Proteome Discoverer PSMs.txt. PSMs at Percolator q-value <= 0.01.
    Mods like 'M8(Oxidation); N-Term(Acetyl)'. tr = 'RT in min'."""
    p = pd.read_csv(path, sep="\t", low_memory=False)
    if "Percolator q-Value" in p:
        p = p[p["Percolator q-Value"] <= 0.01]
    else:
        p = p[p["Confidence"] == "High"]
    seqs = p["Sequence"].str.upper()

    def parse(seq, ms):
        mods = []
        for tok in str(ms).split(";"):
            tok = tok.strip()
            if not tok or tok == "nan":
                continue
            m = re.match(r"^(N-Term|C-Term)(?:\(Prot\))?\((.+)\)$", tok)
            if m:
                pos, name = (0 if m.group(1) == "N-Term" else -1), m.group(2)
            else:
                m = re.match(r"^([A-Z])(\d+)\((.+)\)$", tok)
                if not m:
                    return None
                pos, name = int(m.group(2)), m.group(3)
            if name not in PD_MODS:
                return None
            mods.append((pos, name))
        return tuple(mods)

    keys = seqs + "|" + p["Modifications"].fillna("")
    parsed = {k: parse(*k.split("|", 1)) for k in keys.unique()}
    return pd.DataFrame({
        "run": p["Spectrum File"].values,
        "key": keys.values,
        "seq": seqs.values,
        "mods": [parsed[k] for k in keys],
        "tr": p["RT in min"].values.astype(float),
        "mass": p["Theo MHplus in Da"].values.astype(float) - PROTON,
    })


READERS = {"maxquant": read_maxquant, "fragpipe": read_fragpipe, "diann": read_diann,
           "peaks": read_peaks, "pd": read_pd}


# ------------------------------------------------------------------ common
def peprec(mods):
    def order(pm):
        p = pm[0]
        return (1e9 if p == -1 else p, pm[1])
    return "|".join(f"{p}|{n}" for p, n in sorted(mods, key=order))


def to_proforma(seq, mods):
    res = [[a, ""] for a in seq]
    nterm = cterm = ""
    for p, n in mods:
        if p == 0:
            nterm += f"[{n}]"
        elif p == -1:
            cterm += f"[{n}]"
        else:
            res[p - 1][1] += f"[{n}]"
    return (nterm + "-" if nterm else "") + "".join(a + m for a, m in res) + ("-" + cterm if cterm else "")


def prepare(acc, check=False):
    cfg = DATASETS[acc]
    kw = {k: v for k, v in cfg.items() if k not in ("fmt", "file")}
    df = READERS[cfg["fmt"]](RAW / cfg["file"], **kw)
    ids_per_run = df.groupby("run")["key"].nunique().sort_values(ascending=False)
    run = ids_per_run.index[0]
    d = df[df["run"] == run].copy()
    n_ids = d["key"].nunique()
    bad_seq = ~d["seq"].str.match(AA)
    unmapped = d["mods"].isna()
    drop_keys_mod = d.loc[unmapped & ~bad_seq, "key"].unique()
    drop_keys_seq = d.loc[bad_seq, "key"].unique()
    dropped_mod_types = Counter(d.loc[unmapped & ~bad_seq, "key"].drop_duplicates())
    d = d[~bad_seq & ~unmapped]
    d["modifications"] = [peprec(m) for m in d["mods"]]
    out = (d.groupby(["seq", "modifications"], as_index=False)["tr"].median()
           .sort_values(["seq", "modifications"]))
    out.to_csv(HELD / f"{acc}.csv", index=False, float_format="%.4f")

    # stats
    mod_lists = [m.split("|")[1::2] if m else [] for m in out["modifications"]]
    pf_with = Counter()
    for ml in mod_lists:
        pf_with.update(set(ml))
    non_cam = np.mean([any(x != "Carbamidomethyl" for x in ml) for ml in mod_lists])
    non_label = np.mean([any(x not in ("Carbamidomethyl", "TMTpro", "TMT6plex") for x in ml)
                         for ml in mod_lists])
    print(f"== {acc} ({cfg['fmt']}, {cfg['file']})")
    print(f"   runs: {len(ids_per_run)}; chosen run: {run} ({n_ids} peptidoforms at 1% FDR; "
          f"next best {ids_per_run.iloc[1] if len(ids_per_run) > 1 else '-'})")
    print(f"   rows in run: {(df['run'] == run).sum()}; dropped peptidoforms: "
          f"{len(drop_keys_mod)} unmappable mods, {len(drop_keys_seq)} non-standard residues")
    if len(drop_keys_mod):
        print(f"   examples of unmappable: {list(dropped_mod_types)[:5]}")
    print(f"   written: {len(out)} peptidoforms -> data/heldout/{acc}.csv")
    print(f"   tr range {out.tr.min():.2f} to {out.tr.max():.2f} min; "
          f"length median {out.seq.str.len().median():.0f} (range {out.seq.str.len().min()}-{out.seq.str.len().max()})")
    print(f"   modified (non-CAM): {non_cam:.1%}; modified (non-CAM, non-TMT label): {non_label:.1%}")
    print(f"   mod types (peptidoforms carrying it): {dict(pf_with.most_common())}")

    if check:
        from psm_utils import Peptidoform
        sample = d.drop_duplicates(["seq", "modifications"])
        sample = sample.sample(min(len(sample), 3000), random_state=0)
        err = np.array([Peptidoform(to_proforma(s, m)).theoretical_mass - ms
                        for s, m, ms in zip(sample.seq, sample.mods, sample.mass)])
        ok = np.abs(err) < 0.02
        print(f"   mass check (psm_utils vs engine, n={len(err)}): {ok.mean():.2%} within 0.02 Da; "
              f"median |err| {np.median(np.abs(err)):.4f} Da")
        if not ok.all():
            bad = sample[~ok].head(3)
            print("   mismatches e.g.:", list(zip(bad.key, np.round(err[~ok][:3], 4))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("datasets", nargs="*", default=list(DATASETS))
    ap.add_argument("--check", action="store_true", help="verify peptidoform masses with psm_utils")
    a = ap.parse_args()
    for acc in a.datasets:
        prepare(acc, a.check)


if __name__ == "__main__":
    sys.exit(main())
