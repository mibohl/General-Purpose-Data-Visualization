## Decision table

All rows apply the pre-registered rule unchanged. Rows other than the first are sensitivity checks. None of them is used to pick the verdict.

{{DECISION_TABLE}}

How to read it:

- **Rule 2 fails everywhere.** The gap between a 20-anchor affine map and the best per-study monotone map is 8 to 18 % of A. Most of that gap is anchor-sampling noise rather than warp: an affine map fit on *all* shared precursors (in-sample) already gets to 0.0062, and the isotonic oracle only adds a further 5 %. There is no nonlinear 1/K0 warp worth learning.
- **Calibrated against a noise-only null, rule 2 is unambiguous:** splitting single studies into pseudo-studies gives (A − O)/A = 0.15 with zero real drift, the same as the real 0.14 (Step 3).
- **Rule 1 is a coin flip and should not be leaned on.** A sits right at 2F under the pre-registered floor. It falls clearly below 2F in CCS (A/F = 1.56) and without multi-conformer candidates (1.79). It is clearly above 2F (3.07) if F is replaced by a stricter floor that first removes per-run offsets. The NO-GO does not depend on rule 1.

## Question 1: is the drift already recorded?

**Answer: no. The context is latent.** The reported values are already instrument-calibrated, and drift persists in them.

What I verified and how:

1. **MaxQuant `evidence.txt`, `1/K0` column: it is the instrument calibration applied to the scan number, and nothing else.** Every evidence row carries the raw integer `Ion mobility index` (TIMS scan) next to `1/K0`. Within each of the 745 runs with ≥ 200 MS/MS rows, `1/K0` is an exact smooth function of the scan index: a cubic fit in scan index leaves a 99th-percentile |residual| below 1e-6 Vs/cm² in every run (`q1_run_index_map.csv`). So `1/K0` is a deterministic per-run scan→1/K0 map. In MaxQuant that map comes from Bruker's TDF access library (MaxQuant reads `.tdf` through its Bruker plugin API; Prianichnikov et al., MCP 2020). Bruker's `tims_scannum_to_oneoverk0` applies the calibration stored in `analysis.tdf`; AlphaTims' documentation describes the same call for its `mobility_values`. The per-run maps differ between studies, as recorded calibrations and ramp settings should (1/K0 at scan 500 ranges from 0.96 to 1.29). Within one TIMS method of a study, the map shifts slightly between recalibrations (sd of 1/K0 at scan 500 up to 0.0026).
2. **MaxQuant's own cross-run alignment is *not* in `1/K0`.** It is written to a separate `Calibrated 1/K0` column, which differs from `1/K0` by a per-study shift of −0.019 to +0.030 Vs/cm² in studies where it is active, and is identical where it is not. I used `1/K0`, never `Calibrated 1/K0`. This matters: the aligned column would hide within-study run drift and shrink F.
3. **MaxQuant `CCS` is a fixed transform of `1/K0`.** `CCS / im2ccs(1/K0, m/z, z)` = 0.99999 with zero spread in every study, where `im2ccs` is IM2Deep's Mason-Schamp conversion (N₂, 305 K). So no per-run temperature or pressure enters MaxQuant's CCS. Relative drift is identical in CCS and 1/K0; only the per-study affine fits differ (they are fit in different units).
4. **Not verified:** DIA-NN `report.tsv` and FragPipe `psm.tsv`. Every table available here is MaxQuant. One relevant finding for FragPipe-like pipelines: the open-source `timsrust` reader (used by Sage and others) currently converts scans to 1/K0 with an *uncalibrated* linear map between `OneOverK0AcqRangeLower/Upper` (`UncalibratedScan2ImConverter` in `crates/timsrust-tdf/src/calibration.rs`). Outputs from such readers would not carry the stored TIMS calibration at all.
5. **`.tdf` files were not inspected.** Point 1 already shows directly that the reported value *is* the stored calibration applied to the scan index, so downloading raw `.d` folders would not add evidence.

Evidence that drift persists in the calibrated values:

- Cross-study: D = 0.0114 = 3.3 F. Per-study signed median shifts against the consensus range from −0.020 to +0.060 Vs/cm² (`provenance.csv`).
- **Within a single lab:** the six Gerner-lab projects (Vienna, all timsTOF Pro) shift from −0.016 to +0.035. Drift is not a lab constant.
- **Within a single study, run to run:** in PXD043511 the median per-run offset is 0.030 Vs/cm², while the recorded scan→1/K0 map moved by only sd 0.0008. PXD036206 and PXD035987 show per-run offsets of 0.010 to 0.012. The recorded calibration does not track whatever moves these runs.
- **Recorded metadata predicts little.** Across studies, the per-study shift correlates only weakly with recorded TIMS ramp limits (Spearman |ρ| ≈ 0.4, uncorrected p ≈ 0.04) and instrument generation (Pro 2 vs Pro, p = 0.02 with n = 3 Pro 2 studies). None of the ~30 association tests survives multiple-testing correction (table in Step 5).

## Data and provenance

Source: the MaxQuant `evidence.txt` tables IM2Deep was trained on (Zenodo record 14886716, `MQ_output.zip`, 2.07 GB, 69 files; deleted after extraction). Two files were byte-identical duplicates and were skipped. 36 raw files deposited in both PXD035987 and PXD043511 were kept only in PXD035987. The file named `PXD03206_evidence.txt` carries a truncated accession; I resolved it to **PXD036206** through the PRIDE API (title and raw-file names match the SARS-CoV-2 nsp3 IP study). PXD019086 (Meier et al. 2021) was parsed but excluded from the cross-study analysis: its older MaxQuant export has no `1/K0` column, and the brief reserves it as a single-study reference. All projects are DDA-PASEF: the MaxQuant evidence `Type` column is `TIMS-MULTI-MSMS`/`TIMS-MULTI-MATCH`. Instrument model, lab and organism come from the PRIDE archive API. No file has SILAC `Labeling state` or TMT reporter columns.

{{PROVENANCE}}

`n_shared` = precursors of that study that are among the 34,507 shared precursors. Studies with fewer than 60 shared precursors (PXD035675, PXD036482, PXD036746, PXD038824) are excluded from Step 3. Studies with fewer than 300 precursors are excluded from Step 4.

## Step 1: shared-precursor matrix

Identity = stripped sequence + canonical variable-modification string + charge. Rules (pre-registered in `PREREGISTRATION.md`):

- MS/MS-identified rows only (match-between-runs rows dropped), PEP ≤ 0.01, no reverse or contaminant hits.
- Only unmodified residues, Oxidation (M) and Acetyl (Protein N-term) are allowed. Old MaxQuant notation `(ox)`/`(ac)` (PXD038840, PXD040481) and new notation `(Oxidation (M))`/`(Acetyl (Protein N-term))` are mapped to one token. Each mapped sequence was checked to strip back to the `Sequence` column, with oxidation only on M and acetyl only at the N-terminus (0 violations).
- **All cysteine-containing peptides are dropped.** MaxQuant omits fixed modifications from `Modified sequence`, and the parameter files are not in the archive, so carbamidomethyl vs. other or no alkylation cannot be verified. A sloppy join here would manufacture drift.
- One value per precursor and run: the most intense MS/MS-identified feature. Per study: the median over runs.

Candidate matches rejected (unique precursor identities among 467,490 candidates):

| reason | identities rejected | of which would have been shared (≥ 4 studies) |
|---|---|---|
| modification outside {Ox(M), Acetyl N-term} (phospho, GlyGly, deamidation, 21 synthetic PTM types, ...) | 45,664 | 0 |
| cysteine-containing (alkylation not verifiable) | 83,956 | 8,529 |
| PSM-level filters (MBR, PEP, decoy, contaminant) | rows, see `filter_log.csv`: 6.7 M of 12.3 M rows kept | n/a |

Result: **34,507 shared precursors** (≥ 4 studies; median 5, max 24 studies each), comfortably above the ~2,000 threshold. Charge: z=1 9, z=2 23,758, z=3 9,660, z=4 1,080. m/z: median 709, 5th to 95th percentile 480 to 1,074, range 315 to 1,635.

## Step 2: replicate floor F

Replicate groups are runs of one project acquired with the same TIMS method, i.e. the same scan→1/K0 map within 0.01 Vs/cm² at scan 500 and 3e-5 in slope. Projects mixing methods (PXD035987, PXD043166, SYN21PTM) are split accordingly. For each precursor seen in ≥ 3 runs of a group: |x_run − median(other runs)|.

{{F_TABLE}}

Pooled over 4.73 M run-level deviations: **F = 0.00344** Vs/cm² (IQR 0.0018 to 0.0097). Relative: 0.39 %, identical in CCS (2.04 Å²). Most per-study values (0.0015 to 0.004) are in line with the ~0.002 within-study replicate deviation the brief quotes from Carafe2. Two caveats bear on how far rule 1 can be trusted:

- **Quantization.** MaxQuant's `1/K0` is a function of the *integer* apex scan, so it moves in steps of 0.0003 to 0.0012 Vs/cm² (the TIMS ramp slope). F is only 3 to 10 scan steps.
- **Run offsets inside studies.** Several studies (PXD043511, PXD036206, PXD035987, PXD037622) have per-run offsets of 0.007 to 0.030 that inflate their F. Removing a per-run median offset gives a pure-scatter floor of 0.0022. Under that floor A/F = 3.1. Rule 2 still says NO-GO.

{{FV_TABLE}}

## Step 3: drift and corrections

Held out one study at a time. Consensus = median over the other studies' values for the same precursor. Anchors are drawn from the held-out study's shared precursors; deviations are measured on the remaining shared precursors. 5 seeds.

Pooled (median over all study × precursor deviations, median over seeds):

{{STEP3_POOLED}}

Per study (1/K0, all shared precursors, median over 5 seeds):

{{STEP3_TABLE}}

Readings:

- **The affine map does most of the work.** It takes the pooled deviation from 0.0114 to 0.0068. A one-parameter shift is about as good (0.0066 at 20 anchors, slightly *better* than affine because the slope estimate from 20 noisy anchors costs more than it buys).
- **A flexible monotone warp adds almost nothing,** even with all shared precursors and evaluated in-sample (O = 0.0059 vs 0.0062 for affine on all). Per-charge isotonic maps give 0.0058.
- **Studies with very small shared sets** (PXD038840, PXD039646, PXD042114, n < 100) show large A − O gaps. That is in-sample isotonic overfitting on 60 to 85 points. Cross-fitted isotonic removes it (for example PXD038840: in-sample 0.0011, cross-fitted 0.0055, affine 0.0034).
- **Per study, a few exceed the 0.25 gap at 20 anchors** (PXD035986 0.27, PXD039469 0.27, PXD043382 0.29, SYN21PTM 0.29, plus the three tiny sets). In each, an affine map fit on all shared precursors closes most of it: PXD035986 0.0054 vs O 0.0050 (7 %), PXD043382 0.0047 vs 0.0043 (9 %). The rest come from sets of 229 to 474 precursors, where in-sample isotonic overfits (PXD039469 20 %, SYN21PTM 26 %). It is a 20-anchor sampling effect, not a warp.
- **High-precision studies keep a gap to their own floor.** PXD036191, PXD040521 and PXD043226 have F_study ≈ 0.002 and A ≈ 0.007, i.e. 3.5 to 4 F_study, and their oracle is no better. Part of this is consensus noise: the consensus is a median of other studies that each carry their own drift and noise, so A is an upper bound on what is left in the held-out study. The remainder is not a monotone function of 1/K0 (Step 5 looks for peptide dependence).
- **Noise-only null** (`null_pseudostudies.csv`, `scripts/13_null_pseudostudies.py`). For the 10 largest single-method studies (28 to 116 runs), I split the runs at random into 6 pseudo-studies (3 splits) and reran the identical D/A/O pipeline. Pseudo-studies share lab, instrument and method, so any "drift" between them is replicate noise plus run offsets. Median over the 9 groups with > 1,000 shared precursors: D_null = 0.0034, A_null = 0.0030 (0.7 × that study's F), O_null = 0.0025, and **(A − O)/A = 0.15 under pure noise** (range 0.10 to 0.22). Two consequences:
  1. The real (A − O)/A of 0.14 is *exactly what noise alone produces*. Rule 2's gap is anchor noise plus isotonic fitting noise, not warp structure.
  2. The real A (0.0068, 2.0 F) is about 2× the noise-only A. So after the affine map there *is* genuine study × precursor residual, roughly 0.004 to 0.005 Vs/cm² beyond replicate noise. It is not a monotone function of 1/K0 (oracle) and, per Step 5, barely predictable from peptide features. Candidate causes are study-specific and per-precursor, for example ESI-source-dependent conformer populations, co-eluting interference shifting the apex scan, and feature-detection differences. None of them comes with a support-set-learnable rule.
- **CCS vs 1/K0.** Relative drift is identical by construction (point 3 of Question 1). In CCS units A/F = 1.56 and (A − O)/A = 0.15. The brief's "CCS fixes it" escape hatch does not apply, because the conversion MaxQuant uses has no per-run inputs. It does not matter for the verdict either way.

## Step 4: prediction-model arms

{{STEP4}}

## Step 5: is the residual peptide-dependent?

Residual = consensus − (per-study affine map fit on all shared precursors, in-sample). It was regressed on charge (one-hot), m/z, mass, length, GRAVY hydrophobicity, 20 amino-acid fractions, number of basic residues, C-terminal K/R, oxidation count and N-terminal acetylation. Columns: 5-fold CV R² within the study (Ridge, LightGBM); R² of a LightGBM model trained on all *other* studies' residuals ("pooled other-study model"); the same without multi-conformer candidates; and mean residual by charge.

{{STEP5_TABLE}}

Pattern:

- **Some structure, small, and study-specific.** Within-study CV R² has a median of 0.06 (0.13 without multi-conformer candidates; max 0.29). It is mostly charge-dependent: in several studies z=4 precursors sit 0.003 to 0.007 Vs/cm² off after the affine map (positive in PXD036191 and PXD040521, negative in PXD035986 and PXD035987).
- **The sign and shape differ between studies.** A model trained on the other studies explains nothing (median R² −0.01), and ridge coefficient vectors correlate across studies at a median r = 0.02 (IQR −0.19 to 0.25). This is the qualitative signature an in-context model would exploit.
- **But it is worth very little.** Even the within-study cross-fitted model leaves the median residual MAD unchanged at 0.0058 (all precursors). Without multi-conformer candidates it goes from 0.0056 to 0.0048. That is far less than the 25 % rule 2 requires, and less than the anchor-noise gap between 20-anchor and all-precursor affine maps.
- **Metadata association:** see below. The charge pattern (z4 − z2 residual) correlates with the recorded upper 1/K0 limit of the TIMS ramp (ρ = −0.52, uncorrected p = 0.03, n = 17). That hints that part of the peptide dependence is partially observable from recorded ramp settings. With ~30 tests it is not significant after correction.

{{META_TABLE}}

## Step 6: multi-conformer check

A shared precursor is flagged if, in any single run, two MS/MS-identified features differ in 1/K0 by > 2 %, or if its study values split into two groups of ≥ 2 studies separated by a > 2 % gap. 8,737 of 34,507 shared precursors are flagged (6,968 within-run splits, 2,969 across-study bimodal; overlap counted once). They are enriched in z = 3 (36 % of flagged vs 25 % of unflagged precursors). This is a generous flag: it will also catch co-eluting isobaric features. Headline numbers without them: F = 0.0034, D = 0.0107, A = 0.0061, O = 0.0053; A/F = 1.79, (A − O)/A = 0.13. Still NO-GO on both rules. Removing them helps the within-study structure (Step 5 R² doubles to 0.13) but not the decision.

## Caveats that could change the reading

1. **One search engine, one acquisition mode.** All 25 studies are MaxQuant DDA-PASEF, and 21 are timsTOF Pro. diaPASEF, DIA-NN and FragPipe outputs, and newer instruments (timsTOF HT/Ultra, SCP), are not represented. A pipeline whose reader ignores the stored TIMS calibration (see `timsrust`, Question 1 point 4) would show larger, and more likely nonlinear, drift.
2. **IM2Deep was trained on these same projects.** Z1 is in-distribution and therefore optimistic. That could make the learned arms look *less* useful than on genuinely new studies. Rule 2 and Step 5 do not depend on IM2Deep at all, and they give the verdict.
3. **The consensus is noisy.** A, O and the Step 5 residuals are all measured against a median of 3 to 23 other studies that carry their own drift. That noise affects arms equally, but it caps the R² achievable in Step 5.
4. **Quantization to the scan grid** (0.0003 to 0.0012 Vs/cm² per step) is a property of MaxQuant's reported values. It puts a floor under F that has nothing to do with chemistry or drift.
5. **Run-level drift inside studies** (up to 0.03 Vs/cm²) is real and larger than the per-study residual after affine. A method that adapts *per run* rather than per study is a different question from the one asked here. The cheap answer to that question is again a per-run shift: removing an (oracle) per-run median offset takes F from 0.0034 to 0.0022.

## Deviations from the pre-registration (all made before seeing the affected result)

- Replicate groups: runs within a project were split by TIMS method (scan→1/K0 map) before computing F, because runs acquired with different TIMS methods are not replicates. The rule was decided before any F was computed.
- Step 4 query cap lowered from 3,000 to 1,500 precursors per (study, seed), with TabPFN `n_estimators = 4`, for CPU budget. Decided after timing TabPFN on synthetic data, before any Step 4 result.
- PXD019086 excluded from the cross-study analysis (no 1/K0 column), as the brief anticipated.
- No `Labeling state` or reporter columns exist in any file, so the SILAC/TMT filters were no-ops.

Post-hoc analyses (not pre-registered, added after seeing results, and not used for the verdict): the noise-only pseudo-study null (Step 3), the global-retrain controls G1/G1b (Step 4), the stricter run-offset-free floor, cross-fitted and per-charge oracles, and the metadata association tests.

## Files

| file | content |
|---|---|
| `report.md` | this report |
| `PREREGISTRATION.md` | analysis choices written before opening any evidence table |
| `results.csv` | every run: Step 3 (study × space × subset × method × seed) and Step 4 (study × k × seed × arm), median and IQR of absolute deviation, also as a multiple of F |
| `decision_numbers.csv` | F, D, A, A50, O and both rule checks for every variant |
| `residual_vs_k.png` | Step 4 residual vs k, one line per arm |
| `step3_per_study_table.csv`, `step4_per_study_k.csv`, `step4_head_to_head_vs_C2.csv`, `step4_pooled_median_over_studies*.csv` | Step 3 and 4 tables |
| `step5_variance_explained.csv`, `step5_ridge_coefs.csv`, `step5_metadata_assoc.csv` | Step 5 |
| `provenance.csv`, `F_per_study.csv`, `F_variants.csv`, `q1_run_index_map.csv`, `filter_log.csv`, `multiconformer_flags.csv` | supporting tables |
| `scripts/` | run order: 01, 02, 02b, 03, 04, 05, 06, 07, 08, 09, 13, 14, 10, 15, 11, 12. All paths are arguments; `12_build_report.py` assembles this file from `report_src/` |
| `null_pseudostudies.csv` | noise-only null: one study split into 6 pseudo-studies |
| `report_src/` | report text fragments and rendered tables used by `12_build_report.py` |
