"""Metadata coverage tables for the study (results/metadata_study/metadata_coverage.md + csv).

1. DeepLC 4.5 training collection: 958 PRIDE projects (5,553 of 6,543 heads),
   fields extracted from the PRIDE protocol text by rtpoc/meta_parse.py.
   The other 990 heads (MassIVE MSV accessions, named datasets) were not
   looked up and count as missing.
2. The 28 datasets with observed RTs in this repo (20 Zenodo + 8 held-out),
   curated by hand from PRIDE and publications (data/metadata/curated_lc_metadata.csv).
"""
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "data", "metadata")
OUT = os.path.join(ROOT, "results", "metadata_study")
FIELDS = ["chemistry", "column_length_cm", "inner_diameter_um", "particle_um", "gradient_min",
          "start_pct_b", "flow_nl_min", "temperature_c", "additive"]
LABEL = {"chemistry": "column chemistry", "column_length_cm": "column length",
         "inner_diameter_um": "inner diameter", "particle_um": "particle size",
         "gradient_min": "gradient length", "start_pct_b": "gradient composition (%B)",
         "flow_nl_min": "flow rate", "temperature_c": "column temperature",
         "additive": "additive / ion pairing"}
# hand validation of the parser on random projects (3 samples, see report)
PRECISION = {"chemistry": "~100%", "column_length_cm": "~95%", "inner_diameter_um": "~88%",
             "particle_um": "~88%", "gradient_min": "~95%", "start_pct_b": "~65% (not used)",
             "flow_nl_min": "~85%", "temperature_c": "~85%", "additive": "~85%"}


def main():
    os.makedirs(OUT, exist_ok=True)
    md = pd.read_csv(os.path.join(MD, "pride_lc_metadata.csv"))
    tn = [t for t in open(os.path.join(MD, "deeplc45_task_names.txt")).read().split("\n") if t]
    px = [(re.findall(r"PXD\d{6}", t) or [None])[0] for t in tn]
    heads = pd.DataFrame(dict(pxd=px)).merge(md, on="pxd", how="left")
    single = md.usable & ~md.multi_gradient.fillna(False).astype(bool)
    out = ["# LC metadata coverage\n",
           "## 1. DeepLC 4.5 training collection (PRIDE protocol text, automatic extraction)\n",
           f"{len(md)} PRIDE projects behind {heads.pxd.notna().sum():,} of {len(tn):,} DeepLC 4.5 setup "
           f"heads. The other {heads.pxd.isna().sum()} heads (MassIVE and named datasets) were not "
           "looked up and count as missing below. Publications were not mined for these projects, "
           f"although {md.has_publication.mean():.0%} have one, so these numbers are a lower bound on "
           "what is recorded anywhere.\n"]
    t = pd.DataFrame({"field": [LABEL[f] for f in FIELDS],
                      "projects with value": [f"{md[f].notna().mean():.0%}" for f in FIELDS],
                      "heads with value": [f"{heads[f].notna().mean():.0%}" for f in FIELDS],
                      "extraction precision (hand check)": [PRECISION[f] for f in FIELDS]})
    out.append(t.to_markdown(index=False))
    nf = md.n_fields.value_counts().sort_index()
    out.append("\nFields found per project (of 9): " + ", ".join(f"{k}: {v}" for k, v in nf.items()))
    hu = heads.usable.astype("boolean").fillna(False).astype(bool)
    hs = hu & ~heads.multi_gradient.astype("boolean").fillna(False).astype(bool)
    out.append(f"\n- No LC parameter at all: {int((md.n_fields == 0).sum())} projects "
               f"({(md.n_fields == 0).mean():.0%}).\n"
               f"- Usable (chemistry + gradient length + at least 2 of length, ID, particle, flow; rule "
               f"fixed before counting): {int(md.usable.sum())} projects ({md.usable.mean():.0%}), "
               f"{int(hu.sum()):,} heads ({hu.mean():.0%}).\n"
               f"- Usable and one unambiguous gradient length (the training set of Competitor B): "
               f"{int(single.sum())} projects ({single.mean():.0%}), {int(hs.sum()):,} heads.\n")
    md.assign(usable_single_gradient=single).to_csv(os.path.join(OUT, "metadata_coverage_training.csv"), index=False)
    # curated datasets
    c = pd.read_csv(os.path.join(MD, "curated_lc_metadata.csv"))
    g = pd.to_numeric(c.gradient_min, errors="coerce")
    filled = pd.DataFrame({
        "chemistry": c.chemistry.notna(), "length": c.column_length_cm.notna(),
        "ID": c.inner_diameter_um.notna(), "particle": c.particle_um.notna(),
        "gradient": c.gradient_min.notna(),
        "composition": (c.start_pct_b.notna().astype(float) + c.end_pct_b.notna().astype(float)) / 2,
        "flow": c.flow_nl_min.notna(), "temperature": c.temperature_c.notna(),
        "additive": c.additive.notna()}).astype(float)
    c["fields_filled"] = filled.sum(axis=1)
    c["companion_only"] = c.provenance.fillna("").str.contains("companion", case=False) & c.pxd.eq("PXD074962")
    role = {p: "held-out" for p in ["PXD076643", "PXD079751", "PXD062983", "PXD050641", "PXD074962",
                                    "PXD061731", "PXD084949", "PXD085246"]}
    c["role"] = c.pxd.map(role).fillna("Zenodo (Part 1 / DeepLC paper)")
    c["usable_for_B"] = c.chemistry.notna() & (g.notna() | c.total_run_min.notna()) & (
        filled[["length", "ID", "particle", "flow"]].sum(axis=1) >= 2)
    cols = ["role", "dataset", "pxd", "fields_filled", "usable_for_B", "provenance"]
    c[cols + ["notes"]].to_csv(os.path.join(OUT, "metadata_coverage_curated.csv"), index=False)
    out.append("## 2. Datasets with observed RTs in this repo (curated by hand from PRIDE and publications)\n")
    show = c[cols].copy()
    show["provenance"] = show.provenance.fillna("").str.slice(0, 90)
    show["fields_filled"] = show.fields_filled.map(lambda v: f"{v:g} / 9")
    out.append(show.sort_values(["role", "fields_filled"]).to_markdown(index=False))
    h = c[c.role == "held-out"]
    out.append(f"\nHeld-out setups usable for Competitor B: {int(h.usable_for_B.sum())} of 8. "
               "PXD074962's online-LC values come only from a companion paper of the same core facility "
               "(no source for this dataset itself); PXD050641 has no column information; PXD085246 has no "
               "LC parameter at all, so Competitor B is untestable there.")
    open(os.path.join(OUT, "metadata_coverage.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
