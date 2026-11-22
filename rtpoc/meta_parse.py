"""Extract recorded LC method parameters from PRIDE project records.

Input: data/metadata/pride/<PXD>.json (rtpoc/meta_fetch.py). Text searched:
sampleProcessingProtocol + dataProcessingProtocol (PRIDE truncates long
protocols, so some fields are simply absent). Only sentences that talk about
the LC separation are used, so StageTip desalting (C18) and digestion
temperatures are not picked up.

Fields (the ones the study asked for): chemistry, column_length_cm,
inner_diameter_um, particle_um, gradient_min, start_pct_b / end_pct_b,
flow_nl_min, temperature_c, additive, plus lc_system, evosep_spd and
enrichment / label type. A project with more than one distinct gradient
length is flagged `multi_gradient` (its heads cannot be assigned one value).

"Usable" (fixed before looking at the coverage numbers): chemistry and
gradient length known, and at least two of column length, inner diameter,
particle size, flow rate.

Output: data/metadata/pride_lc_metadata.csv
Usage: python -I rtpoc/meta_parse.py
"""
import glob
import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UM = r"(?:µm|μm|um|micron|microns|µ m|μ m)"
NUM = r"(\d+(?:[.,]\d+)?)"
LC_WORDS = re.compile(r"(column|gradient|nL/min|nl/min|µl/min|μl/min|ul/min|nanoLC|nano-LC|nanoflow|"
                      r"separat|eluted|analytical|EASY-nLC|Easy-nLC|nanoElute|Vanquish|UltiMate|"
                      r"Ultimate|Evosep|Acquity|M-Class|buffer A|solvent A|mobile phase|HPLC|UHPLC|"
                      r"oven|Aurora|PepMap|ReproSil|Reprosil|EASY-Spray|Easy-Spray|µPAC|IonOpticks)", re.I)
DESALT = re.compile(r"(stage ?tip|stagetip|desalt|sep-pak|solid.phase extraction|spe |zip ?tip|"
                    r"spin column|microspin|oasis|digest|incubat|lysis|lysed|reduc|alkylat|"
                    r"fractionat)", re.I)
# Evosep standard methods: samples per day -> gradient minutes (Evosep One / Eno)
EVOSEP = {500: 2.2, 300: 2.3, 200: 5.6, 100: 11.5, 60: 21.0, 30: 44.0, 15: 88.0,
          80: 14.0, 40: 31.0, 20: 58.0}


def fnum(s):
    return float(s.replace(",", "."))


def sentences(text):
    s = re.split(r"(?<=[.;])\s+(?=[A-Z(])", text)
    return [x for x in s if LC_WORDS.search(x)]


def extract(text):
    d = {}
    lc = sentences(text)
    lcj = " ".join(lc)
    sep = [s for s in lc if not DESALT.search(s) or re.search(r"analytical|gradient|nL/min", s, re.I)]
    sepj = " ".join(sep)
    # sentences about the analytical column, not the trap / loading step
    ana = [x for x in sep if not re.search(r"trap|pre-?column|precolumn|load(?:ed|ing)? (?:at|with)|loading", x, re.I)
           or re.search(r"analytical|separation column|resolving|separated on|separated using", x, re.I)]
    anaj = " ".join(ana)
    # chemistry
    if re.search(r"\bHILIC\b|hydrophilic interaction", sepj, re.I):
        d["chemistry"] = "HILIC"
    elif re.search(r"\bSCX\b|strong cation", sepj, re.I) and not re.search(r"C18", sepj):
        d["chemistry"] = "SCX"
    elif re.search(r"\bC18\b|C-18|reverse[d]?[- ]phase|\bRP\b|PepMap|ReproSil|Reprosil|Aurora|"
                   r"Acclaim|BEH|Aeris|Luna|Jupiter|µPAC|uPAC|EASY-Spray|Easy-Spray|IonOpticks|"
                   r"Odyssey|Kinetex|Zorbax|CSH", sepj):
        d["chemistry"] = "RP"
    elif re.search(r"\bC8\b|\bC4\b", sepj):
        d["chemistry"] = "RP"
    elif re.search(r"Evosep", text, re.I):
        d["chemistry"] = "RP"  # Evosep runs only its own C18 columns
    # column length (largest plausible value in column sentences: analytical > trap)
    lens = []
    for s in sep:
        if not re.search(r"column|length|long|packed|emitter|tip", s, re.I):
            continue
        for v, u in re.findall(NUM + r"\s*[- ]?(cm|mm)\b", s):
            x = fnum(v) * (0.1 if u == "mm" else 1)
            if 5 <= x <= 100 and not (u == "mm" and fnum(v) < 30):
                lens.append(x)
    if lens:
        d["column_length_cm"] = max(lens)
    # inner diameter: 20-300 um near x / i.d. / inner; or 1-4.6 mm i.d.
    ids = []
    for m in re.finditer(NUM + r"\s*" + UM + r"\s*(?:i\.?\s?d\.?|ID\b|inner|internal|x|×|X\b|\*)", anaj):
        x = fnum(m.group(1))
        if 20 <= x <= 300:
            ids.append(x)
    for m in re.finditer(r"(?:inner|internal) diameter(?: of)?\s*(?:\(ID\))?\s*(?:of\s*)?" + NUM + r"\s*" + UM, anaj, re.I):
        x = fnum(m.group(1))
        if 20 <= x <= 300:
            ids.append(x)
    for m in re.finditer(r"(?:x|×|by)\s*" + NUM + r"\s*" + UM, anaj):
        x = fnum(m.group(1))
        if 20 <= x <= 300:
            ids.append(x)
    for m in re.finditer(r"(?:x|×)\s*" + NUM + r"\s*mm\b|" + NUM + r"\s*mm\s*(?:i\.?d|ID\b|x|×)", anaj):
        v = m.group(1) or m.group(2)
        x = fnum(v)
        if 0.3 <= x <= 4.6:
            ids.append(x * 1000)
    if ids:
        d["inner_diameter_um"] = ids[0]
    # particle size: 1-5 um, not followed by i.d./x
    parts = []
    for m in re.finditer(NUM + r"\s*[- ]?" + UM + r"(?!\s*(?:i\.?\s?d|ID\b|inner|internal|x\b|×|X\b|\*))", anaj):
        x = fnum(m.group(1))
        if 1.0 <= x <= 5.0:
            parts.append(x)
    if parts:
        d["particle_um"] = parts[0]
    # gradient length
    grads = []
    for s in lc:
        if not (re.search(r"gradient", s, re.I)
                or re.search(r"separat|elut", s, re.I) and re.search(r"over|within|during", s, re.I)
                and not DESALT.search(s)):
            continue
        for m in re.finditer(NUM + r"\s*[- ]?(?:min|mins|minutes|minute)\b", s):
            x = fnum(m.group(1))
            if 5 <= x <= 600:
                grads.append(x)
        for m in re.finditer(NUM + r"\s*[- ]?(?:h|hr|hrs|hour|hours)\b(?:\s*[- ]?(?:gradient|long|linear|LC|run|method))", s):
            x = fnum(m.group(1)) * 60
            if 5 <= x <= 600:
                grads.append(x)
    spd = re.findall(r"(\d{2,3})\s*[- ]?(?:SPD|samples[- ]per[- ]day)", text, re.I)
    if spd and re.search(r"Evosep", text, re.I):
        d["evosep_spd"] = int(spd[0])
        if int(spd[0]) in EVOSEP:
            grads.insert(0, EVOSEP[int(spd[0])])
    if grads:
        # take the longest stated duration in gradient sentences: segment times and
        # equilibration steps are shorter than the full gradient
        sig = sorted({round(g) for g in grads})
        d["gradient_min"] = max(grads) if "evosep_spd" not in d else grads[0]
        big = [g for g in sig if g >= 20]
        d["multi_gradient"] = len(big) > 1 and (max(big) / min(big) > 1.6)
    # %B composition
    gsent = " ".join(x for x in lc if re.search(r"gradient", x, re.I)
                     and not re.search(r"fraction|high[- ]pH|pH ?10", x, re.I))
    m = re.search(NUM + r"\s*%?\s*(?:B|buffer B|solvent B|ACN|acetonitrile)?\s*(?:to|-|–|up to)\s*" + NUM + r"\s*%", gsent)
    if m:
        a, b = fnum(m.group(1)), fnum(m.group(2))
        if a < b <= 100:
            d["start_pct_b"], d["end_pct_b"] = a, b
    if re.search(r"acetonitrile|ACN|MeCN", sepj):
        d["solvent_b"] = "ACN"
    elif re.search(r"methanol", sepj, re.I):
        d["solvent_b"] = "MeOH"
    # flow
    elu = " ".join(x for x in lc if re.search(r"elut|gradient|separat", x, re.I)
                   and not re.search(r"load|fraction|high[- ]pH|pH ?10", x, re.I))
    fl = re.findall(NUM + r"\s*(nl|nL|NL|µl|µL|μl|μL|ul|uL|ml|mL)\s*(?:/|per)\s*min", elu)
    if not fl:
        fl = re.findall(NUM + r"\s*(nl|nL|NL|µl|µL|μl|μL|ul|uL|ml|mL)\s*/\s*min", anaj)
    if not fl:
        fl = re.findall(NUM + r"\s*(nl|nL|NL|µl|µL|μl|μL|ul|uL|ml|mL)\s*/\s*min", lcj)
    if not fl:
        fl = re.findall(NUM + r"\s*(nl|nL|µl|µL|μl|μL|ul|uL)\s*(?:per|/)\s*min", lcj)
    if fl:
        v, u = fl[0]
        x = fnum(v) * (1 if u.lower() == "nl" else 1000 if u.lower() in ("µl", "μl", "ul") else 1e6)
        if 10 <= x <= 2e6:
            d["flow_nl_min"] = x
    # temperature: only sentences about the column / oven
    for s in sep:
        if re.search(r"column|oven|heater|heated", s, re.I) and not DESALT.search(s):
            t = re.findall(r"(\d{2})\s*(?:°|º|˚|Â°|degrees?)\s*C", s)
            t = [int(x) for x in t if 20 <= int(x) <= 80]
            if t:
                d["temperature_c"] = t[0]
                break
    # additive / ion pairing (mobile phase)
    mp = " ".join(s for s in lc if re.search(r"buffer A|solvent A|mobile phase|buffer B|solvent B|"
                                             r"phase A|phase B|\bA:|\bB:", s, re.I)) or sepj
    # FA wins over TFA when both appear: TFA is mostly the loading / reconstitution acid
    if re.search(r"formic acid|\bFA\b|HCO2H", mp, re.I):
        d["additive"] = "FA"
    elif re.search(r"trifluoroacetic|\bTFA\b", mp):
        d["additive"] = "TFA"
    elif re.search(r"acetic acid", mp, re.I):
        d["additive"] = "AcOH"
    elif re.search(r"ammonium formate|ammonium hydroxide|pH 10", mp, re.I):
        d["additive"] = "high-pH buffer"
    # LC system
    for pat, name in ((r"Evosep", "Evosep"), (r"nanoElute", "nanoElute"), (r"Vanquish Neo", "Vanquish Neo"),
                      (r"EASY-?nLC|Easy-?nLC|EASY nLC|Easy nLC", "EASY-nLC"),
                      (r"UltiMate|Ultimate 3000|Dionex", "UltiMate 3000"),
                      (r"nanoAcquity|NanoAcquity|nanoACQUITY|M-Class|M Class", "Waters nano"),
                      (r"Vanquish", "Vanquish"), (r"Eksigent|ekspert", "Eksigent")):
        if re.search(pat, text):
            d["lc_system"] = name
            break
    return d


def enrichment(text, ptms, quant):
    ptms = [x.get("name", "") if isinstance(x, dict) else str(x) for x in ptms]
    t = text + " " + " ".join(ptms) + " " + " ".join(quant)
    tags = []
    for pat, name in ((r"HLA|MHC|immunopeptid", "HLA/MHC peptides"),
                      (r"TiO2|IMAC|Fe-NTA|Fe\(III\)|phosphopeptide enrich|phospho-enrich|phosphoenrich", "phospho-enriched"),
                      (r"\bTMT|tandem mass tag", "TMT"), (r"iTRAQ", "iTRAQ"),
                      (r"diGly|di-Gly|K-ε-GG|ubiquitin remnant", "diGly"),
                      (r"glycopeptide|N-glyco|glycoproteom", "glyco"),
                      (r"crosslink|cross-link|XL-MS", "crosslinked"),
                      (r"acetyl-lysine|anti-acetyl|acetylome", "acetyl-enriched"),
                      (r"peptidom", "peptidome")):
        if re.search(pat, t, re.I):
            tags.append(name)
    return ";".join(tags) if tags else "none detected"


def main(sub="pride", out="pride_lc_metadata.csv"):
    rows = []
    for f in sorted(glob.glob(os.path.join(ROOT, "data", "metadata", sub, "PXD*.json"))):
        j = json.load(open(f))
        sp = j.get("sampleProcessingProtocol") or ""
        dp = j.get("dataProcessingProtocol") or ""
        text = sp + " " + dp
        d = extract(text)
        d["enrichment"] = enrichment(text, j.get("identifiedPTMStrings") or [],
                                     [x.get("name", "") if isinstance(x, dict) else str(x)
                                      for x in (j.get("quantificationMethods") or [])])
        inst = j.get("instruments") or []
        d["instrument"] = ";".join(x.get("name", "") if isinstance(x, dict) else str(x) for x in inst)
        d["pxd"] = j["accession"]
        d["protocol_chars"] = len(sp)
        d["protocol_truncated"] = len(sp) >= 4990
        d["has_publication"] = bool(j.get("references"))
        rows.append(d)
    df = pd.DataFrame(rows)
    core = ["chemistry", "column_length_cm", "inner_diameter_um", "particle_um", "gradient_min",
            "start_pct_b", "flow_nl_min", "temperature_c", "additive"]
    for c in core + ["multi_gradient"]:
        if c not in df:
            df[c] = None
    df["n_fields"] = df[core].notna().sum(axis=1)
    hw = df[["column_length_cm", "inner_diameter_um", "particle_um", "flow_nl_min"]].notna().sum(axis=1)
    df["usable"] = df.chemistry.notna() & df.gradient_min.notna() & (hw >= 2)
    cols = ["pxd"] + core + ["end_pct_b", "solvent_b", "multi_gradient", "evosep_spd", "lc_system",
                             "enrichment", "instrument", "protocol_chars", "protocol_truncated",
                             "has_publication", "n_fields", "usable"]
    df = df[[c for c in cols if c in df]]
    df.to_csv(os.path.join(ROOT, "data", "metadata", out), index=False)
    print(df[core].notna().mean().round(3).to_string())
    print("usable", df.usable.mean().round(3), "n", len(df))


if __name__ == "__main__":
    main()
    main("pride_heldout", "pride_lc_metadata_heldout.csv")
