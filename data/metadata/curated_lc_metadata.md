# Curated LC metadata for the 28 RT datasets

File: `curated_lc_metadata.csv` (one row per dataset). Curated 2026-11-22 from
PRIDE Archive v3 records, open-access full texts (Europe PMC, NCBI efetch, the
PMC open-access S3 bucket), publisher supplementary files on figshare, and the
Zenodo 5171792 record. Values were copied only when a source states them.
Empty cell = not reported; each row's `notes` ends with a "Not reported:" list.
`provenance` names the source (record, DOI + section/table) for the row.

Conventions: `inner_diameter_um` and `flow_nl_min` are in µm and nl/min, so the
analytical-scale first-dimension columns read 3000 / 2100 µm and 300000 nl/min.
`gradient_min` is the active gradient; `arabidopsis` keeps the stated range
"110-120". Values tagged `[companion]` come from a closely related paper by
the same lab and are not stated for that dataset itself (see below).

## How the 20 DeepLC datasets map to their sources

| Dataset(s) | Source | Accession |
|---|---|---|
| hela_hf | Kelstrup 2018 JPR (Q Exactive HF-X), 46 fractions, 15 min runs | PXD006932 |
| hela/plasma_lumos_1h/2h | Li 2019 JASMS | PXD013477 |
| pancreas | Wang 2019 Mol Syst Biol (29 tissues) | PXD010154 |
| yeast_60min/120min | Jarnuczak 2016 JPR | PXD003472 |
| arabidopsis | Mucha 2019 Plant Cell | PXD008812 |
| proteometools_library | Zolg 2017 Nat Methods | PXD004732 |
| prosit_ptm_2020 | Zolg 2018 MCP (ProteomeTools PTM) | PXD009449 |
| dia_fixed_mods | Rosenberger 2014 Sci Data (SWATH pan-human library) | PXD000953 |
| PXD005573_mcp | Bruderer 2017 MCP (DIA HF) | PXD005573 |
| unmod_fixed_mods | Nagaraj 2012 MCP (yeast, LysC, 4 h) | Tranche, no PXD |
| mod_fixed_mods | Sharma 2014 Cell Rep (HeLa phospho) | PXD000612 |
| ATLANTIS/LUNA_SILICA, LUNA_HILIC, Xbridge | Spicer & Krokhin 2018 J Chromatogr A | none |
| SCX | Gussakovsky 2017 Anal Chem | none |

Mapping evidence: the DeepLC paper's Data availability list and references;
the Zenodo 5171792 README table and `dataset_configurations.zip` (PXD plus
raw-file patterns for the nine "psms_aligned" sets); and the DeepRT+ GitHub
`data/README.md`, which names the sources of mod/unmod/dia/SCX/HILIC files.

Corrections to earlier assumptions:
- `unmod_fixed_mods` is not the unmodified half of the phospho study. DeepRT+
  pairs `unmod.txt` with Nagaraj 2012 yeast and `mod.txt` with Sharma 2014 HeLa phospho.
- `Xbridge_fixed_mods` is XBridge **Amide** HILIC (pH 4.5), not high-pH RP.
- `dataset_configurations.zip` has the hela_hf and pancreas configs swapped;
  the paper, the README table and the RT ranges all give hela_hf = PXD006932.

## What was found

- Complete (9 of 9 fields): hela_lumos_1h/2h, plasma_lumos_1h/2h,
  proteometools_library, prosit_ptm_2020, PXD076643 (all PRIDE values confirmed by the paper).
- Near complete (8 to 8.5): pancreas, both yeast sets, arabidopsis,
  dia_fixed_mods, unmod_fixed_mods, PXD005573_mcp, PXD062983, PXD061731.
  The usual gap is column temperature.
- 2D-LC first dimension (silica, HILIC, amide, SCX): column, size, eluents,
  pH, flow, gradient slope and fraction count come from the Krokhin lab's
  Yeung et al. 2020 Anal Chem SI Table S1 (doi 10.1021/acs.analchem.9b05407).
  The cited 2017/2018 papers are paywalled. The XBridge Amide row (38
  fractions, 44,489 peptides) matches Krokhin 2017 Anal Chem SI exactly, so
  the 2020 table reuses these runs. Gradient length, end condition, LC system
  and temperature are not reported. These RTs are 1 min fraction times.
- PXD079751: column and solvents come from PMID 40846850, which PRIDE cites "as
  described before". Flow and temperature are not in either source.
- PXD074962: no paper found. Online LC values are `[companion]`, taken from the
  same lab and core's Nat Cardiovasc Res 2026 paper (doi 10.1038/s44161-026-00818-2),
  whose offline high-pH text is identical to the PRIDE record.

## What was not found

- hela_hf: Kelstrup 2018 full text unreachable (ACS and the KU repository
  return 403; the SI has no LC methods). Only the 15 min gradient (raw file
  names) is recorded.
- mod_fixed_mods: flow, temperature, and whether the 120 or 240 min gradient
  applies (the Cell Reports full text was unreachable).
- PXD050641: no publication found; column not named in PRIDE.
- PXD084949: no publication found; the lab's eLife 2023 paper has the same
  text and gives no column ID, flow or oven setpoint.
- PXD085246: no publication found; PRIDE gives only "EASY-nLC 1200".
- Hosts that blocked access: pubs.acs.org, cell.com, sciencedirect.com,
  researchprofiles.ku.dk, research.manchester.ac.uk, krokhinlab.ca, and PMC
  HTML pages (some). Europe PMC `fullTextXML` returned HTTP 500 throughout.
