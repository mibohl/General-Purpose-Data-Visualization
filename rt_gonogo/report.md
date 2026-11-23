# Go / No-Go: small-molecule retention time on the absolute time axis

**Not GO. The pre-registered rule returns CONDITIONAL for reversed phase and INDETERMINATE (rule gap) for HILIC.** On reversed-phase systems the released 2-step checkpoint never saw (25 held out; 24 and 23 have enough compounds for the purged evaluation at k = 15 and 30; spread anchors), headroom is H = 0.33 (95% CI 0.28 to 0.40) at k = 15 and H = 0.27 (0.24 to 0.31) at k = 30. Both clear the 0.25 bar. The best in-context arm (T2, TabPFN v2 stacked on the quadratic) gains only G = 0.05 (CI 0.01 to 0.13) at k = 15 and G = 0.10 (CI -0.03 to 0.23) at k = 30, below the 0.15 bar, so both decision points are CONDITIONAL. For HILIC (10 datasets, proxy ROI), H = 0.21 (0.16 to 0.33) at k = 15 and 0.22 (0.16 to 0.29) at k = 30, which falls in the gap between the NO-GO (< 0.15) and GO (>= 0.25) thresholds that the rule does not cover. HILIC G is 0.00 at k = 15 and 0.06 (T1, CI 0.02 to 0.20) at k = 30. RP numbers use the **released 2-step checkpoint** (`2-step0525.pt`, repo commit `41a4233`). The pH fed to it was imputed for all 25 RP systems (see Deviations). HILIC numbers rest on the **pooled ROI proxy**, because the checkpoint is reversed-phase only. In absolute terms, on the RP headline set at k = 30, the incumbent's median error is 2.01 min (10.3% of gradient). The oracle monotone map reaches 1.24 min (6.8%), and the best in-context arm reaches 1.73 min.

**I think the rule overstates the prize, and say so here while reporting its verdict unchanged.** O1 is an isotonic map fitted in-sample on every compound, so it also fits noise in the ROI ordering. A 10-fold cross-fitted version of the same oracle gives RP H = 0.20 at k = 15 and 0.14 at k = 30. Under that ceiling, k = 30 would be NO-GO and k = 15 would fall in the rule gap. More importantly, the experiment's premise does not hold on unseen systems. **Most of the remaining error is in the ordering, not the time axis.** The order error (O1, about 1.2 to 1.4 min) exceeds the warp error (I1 - O1, about 0.5 to 0.7 min) on the median dataset. In-context arms help mainly by fixing order flips that no warp can touch, not by learning a better time axis.

## Decision table (headline cells)

Purged evaluation is the headline. H is the median over datasets of the per-dataset mean over 5 seeds of (MAE_I1 - MAE_O1) / MAE_I1. G is the same for the best of T1/T2/R1. CIs are percentile bootstraps over datasets (10,000 resamples). All CIs are per arm and do not adjust for picking the best of three arms. Full grid: `decision_table.csv`.

| Group | ROI | Eval | Anchors | k | n datasets | H (95% CI) | H, cross-fitted O1 (sensitivity) | G best arm (95% CI) | Rule verdict |
|---|---|---|---|---|---|---|---|---|---|
| **RP-2step (headline)** | checkpoint | purged | spread | 15 | 24 | **0.33** (0.28, 0.40) | 0.20 | **0.05** T2 (0.01, 0.13) | **CONDITIONAL** |
| **RP-2step (headline)** | checkpoint | purged | spread | 30 | 23 | **0.27** (0.24, 0.31) | 0.14 | **0.10** T2 (-0.03, 0.23) | **CONDITIONAL** |
| RP-2step | checkpoint | purged | random | 15 | 24 | 0.24 (0.20, 0.27) | 0.09 | -0.02 R1 (-0.12, 0.04) | INDETERMINATE |
| RP-2step | checkpoint | purged | random | 30 | 23 | 0.19 (0.15, 0.23) | 0.06 | 0.09 T2 (-0.03, 0.21) | INDETERMINATE |
| RP-2step | checkpoint | all | spread | 15 | 25 | 0.36 (0.29, 0.45) | 0.18 | 0.07 T1 (0.00, 0.15) | CONDITIONAL |
| RP-2step | checkpoint | all | spread | 30 | 25 | 0.34 (0.26, 0.41) | 0.15 | 0.09 T2 (0.06, 0.24) | CONDITIONAL |
| RP-2step | checkpoint | scaffold-purged | spread | 30 | 21 | 0.33 (0.28, 0.42) | 0.20 | 0.11 T2 (0.04, 0.21) | CONDITIONAL |
| RP-proxy (robustness) | proxy | purged | spread | 15 | 25 | 0.27 (0.14, 0.32) | 0.14 | 0.02 T1 (-0.04, 0.09) | CONDITIONAL |
| RP-proxy (robustness) | proxy | purged | spread | 30 | 25 | 0.26 (0.19, 0.31) | 0.13 | 0.09 T1 (0.01, 0.13) | CONDITIONAL |
| **HILIC** | proxy | purged | spread | 15 | 10 | **0.21** (0.16, 0.33) | 0.11 | **0.00** T2 (-0.15, 0.08) | **INDETERMINATE (rule gap)** |
| **HILIC** | proxy | purged | spread | 30 | 10 | **0.22** (0.16, 0.29) | 0.14 | **0.06** T1 (0.02, 0.20) | **INDETERMINATE (rule gap)** |
| HILIC | proxy | all | spread | 15 | 11 | 0.29 (0.23, 0.36) | 0.17 | 0.01 T2 (-0.27, 0.11) | CONDITIONAL |
| HILIC | proxy | all | spread | 30 | 11 | 0.26 (0.19, 0.31) | 0.15 | 0.07 T1 (-0.05, 0.20) | CONDITIONAL |
| PUB-2step (flagged: in checkpoint training) | checkpoint | all | spread | 15 | 6 | 0.48 (0.36, 0.66) | 0.27 | 0.04 T1 (-0.12, 0.17) | CONDITIONAL |
| PUB-2step (flagged: in checkpoint training) | checkpoint | all | spread | 30 | 6 | 0.40 (0.28, 0.65) | 0.18 | 0.08 T2 (-0.01, 0.29) | CONDITIONAL |
| PUB-proxy | proxy | purged | spread | 15 | 6 | 0.14 (0.06, 0.50) | 0.01 | 0.25 T2 (0.05, 0.41) | NO-GO |
| PUB-proxy | proxy | purged | spread | 30 | 6 | 0.14 (0.07, 0.43) | -0.02 | 0.15 T1 (-0.02, 0.52) | NO-GO |

Decision-cell rule (frozen in `prereg_interpretation.md` before any error was computed): primary cells are purged, spread anchors, k = 15 and k = 30. If they disagree, the more conservative verdict is the headline. Random anchors are a sensitivity check. With random anchors, RP H drops into the rule gap (0.24 and 0.19), which points the same way: not GO. PUB-2step has no purged evaluation, because the checkpoint trained on those datasets and every compound is in its training pool.

## What was run

- **Data.** RepoRT commit `420b550` (412 RP and HILIC datasets, SMRT excluded). Per dataset: canonical RT table with isomeric SMILES substituted where available; duplicate (SMILES, RT) rows dropped; SMILES with conflicting RTs dropped (2-step `remove_doublets`); compounds eluting before 2 * t0 dropped (2-step `void_factor`). Provenance for every dataset used: `provenance.csv`.
- **Held-out sets (deterministic rules, `prereg_interpretation.md`).**
  - **RP-2step, headline, 25 datasets.** Unseen by the checkpoint, checkpoint-runnable, >= 50 compounds, at most 3 per column name, 14 distinct column names.
  - **RP-proxy, 33 datasets.** Unseen RP, >= 80 compounds, at most 2 per column name. Mostly columns without HSM/Tanaka, where the checkpoint cannot run.
  - **PUB.** The six published targets.
  - **HILIC, 11 datasets.** >= 80 compounds, at most 3 per column name, 5 distinct column names.
  - All held-out sets are excluded from every trained component (proxy ROI, Z1).
- **Leakage.** Compound overlap uses the InChIKey first block (`overlap_summary.csv`). As warned, the median RepoRT dataset has 100% of its compounds measured in some other dataset (held-out sets: 96%). The purged set removes query compounds present anywhere in the training pool. For RP-2step that pool is the checkpoint's 171 training datasets plus the Z1 pool. Purged sets are small: the median at k = 30 is 19 query compounds per RP-2step dataset (RP-proxy 29, HILIC 45). Datasets with fewer than 20 compounds outside the pool are excluded from purged medians. Runs with fewer than 5 purged queries are dropped.
- **Arms.**
  - **I1:** the unmodified 2-step `LADModel` (LAD quadratic, OLS refit, same flags as `predict.py`).
  - **I2:** isotonic. **I3:** linear. **I4:** monotone piecewise-linear spline.
  - **O1:** L1 isotonic on all compounds, the MAE-optimal monotone map. Cross-fitted and L2 variants are kept for reference only.
  - **Z1:** LightGBM on [ROI, 12 RDKit descriptors, 36 metadata and gradient features], trained on the pool, no anchors.
  - **T1 / T2:** TabPFN v2 (official v2 regressor checkpoint, 8 estimators). T2 adds the I1 prediction as a feature.
  - **R1:** RidgeCV residual on top of I1.
  - Every arm except O1 and Z1 sees only the k anchors. Nothing is tuned on query compounds or held-out datasets.
- **Grid.** k in {10, 15, 20, 30, 50} (k = 50 only where n >= 80), 5 seeds, random and spread anchors, 81 dataset-by-ROI units. Every run: `results.csv` (130,680 rows).

## Deviations and caveats (read before citing any number)

1. **The published targets are not held out for the released checkpoint.** Its stored config lists 171 RepoRT training datasets, including 0002, 0003, 0018, 0019, 0054 and 0055. PUB-2step numbers are therefore in-sample for ordering and are flagged everywhere. The **published 2-step benchmark was not reproduced** as a held-out experiment, because that requires retraining the D-MPNN without each target, which is out of scope (no GPU training). PUB-proxy is the clean, but weaker, comparison.
2. **The checkpoint requires recorded pH and HSM + Tanaka column parameters.** If any of its 13 system features is missing, it returns NaN. None of the 73 eligible unseen RP datasets had all three. I added a rule before computing any error (amendment 1, committed separately): pH is imputed, as a 2-step input only, as the median recorded pH among RepoRT RP datasets with the same additive pattern (`ph_imputation_table.csv`; e.g. formic acid only gives pH 3.0, based on 226 datasets). **All 25 RP-2step datasets use imputed pH.** Unseen columns without HSM/Tanaka cannot run the checkpoint; those datasets are covered only by RP-proxy. To reach 25 datasets, the size floor dropped from 80 to 50 compounds, as the amendment allows.
3. **HILIC rests entirely on the proxy ROI** (LightGBM on RDKit descriptors + ECFP4, pooled over 60 non-held-out HILIC datasets, condition-unaware). So HILIC H measures headroom over a weaker ordering than RP-2step.
4. Gradient length for UniToyama Atlantis (0018) comes from the 2-step README (40 min), because RepoRT has no gradient table for it. Elsewhere it is the last time point of the RepoRT gradient table.
5. The published % figures use the authors' own split protocol and anchor sets. The comparison below is indicative only.

## Step 4: error decomposition (minutes)

Order error = O1. Warp error = I1 - O1, the part a better warp could attack. Total = I1. Medians over seeds, spread anchors. Remaining groups (PUB-2step, PUB-proxy, RP-proxy) are in the appendix and in `error_decomposition.csv`.

**RP-2step** (spread anchors, minutes; order = O1, warp = I1 - O1, total = I1; purged columns blank when the dataset has < 20 compounds outside the training pool)

| id | dataset | column | n | grad (min) | k=15 total | k=15 order | k=15 warp | k=30 total | k=30 order | k=30 warp | k=30 warp share | purged k=30 total | purged k=30 order | purged k=30 warp | purged n_q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0263 | SPECIMEN_Neg | Waters ACQUITY UPLC BEH C18 | 725 | 30 | 1.91 | 1.46 | 0.45 | 1.72 | 1.44 | 0.28 | 0.16 | 2.01 | 1.74 | 0.27 | 197 |
| 0010 | FEM_lipids | Merck Supelco Ascentis Express C18 | 54 | 30 | 2.31 | 1.83 | 0.48 | 2.42 | 2.02 | 0.39 | 0.16 | 3.26 | 4.78 | -1.51 | 6 |
| 0414 | ROASMI_67 | Agilent ZORBAX Eclipse Plus C18 | 225 | 28 | 1.63 | 1.36 | 0.27 | 1.68 | 1.35 | 0.32 | 0.19 | 2.74 | 2.24 | 0.50 | 19 |
| 0262 | SPECIMEN_Pos | Waters ACQUITY UPLC BEH C18 | 883 | 30 | 1.48 | 1.16 | 0.32 | 1.50 | 1.15 | 0.35 | 0.23 | 1.94 | 1.60 | 0.34 | 253 |
| 0047 | Acquity HSST3-RP-60mn | Waters ACQUITY UPLC HSS T3 | 159 | 60 | 3.39 | 2.58 | 0.80 | 3.26 | 2.52 | 0.74 | 0.23 | 5.62 | 4.66 | 0.96 | 8 |
| 0415 | ROASMI_83 | Agilent ZORBAX Eclipse Plus C18 | 650 | 22 | 1.33 | 0.97 | 0.35 | 1.24 | 0.95 | 0.29 | 0.23 | 1.66 | 1.24 | 0.42 | 290 |
| 0400 | ROASMI_79 | Waters ACQUITY UPLC HSS C18 | 305 | 5.4 | 0.40 | 0.26 | 0.14 | 0.35 | 0.27 | 0.08 | 0.23 | 0.74 | 0.60 | 0.15 | 19 |
| 0428 | ROASMI_62 | Phenomenex Kinetex C18 | 180 | 25 | 0.94 | 0.61 | 0.33 | 0.84 | 0.62 | 0.22 | 0.26 | 0.89 | 0.66 | 0.23 | 113 |
| 0433 | ROASMI_66 | Merck Supelco Ascentis Express C18 | 119 | 25 | 2.47 | 1.75 | 0.72 | 2.48 | 1.78 | 0.70 | 0.28 | 3.71 | 2.75 | 0.96 | 12 |
| 0401 | ROASMI_80 | Waters ACQUITY UPLC BEH C18 | 332 | 5.4 | 0.49 | 0.34 | 0.14 | 0.48 | 0.34 | 0.14 | 0.28 | 0.64 | 0.48 | 0.16 | 21 |
| 0404 | ROASMI_49 | Thermo Scientific Accucore C18 | 97 | 27.5 | 3.24 | 2.05 | 1.19 | 2.60 | 1.79 | 0.80 | 0.30 | 2.78 | 1.91 | 0.87 | 5 |
| 0394 | ACQUITY_UPLC_HSSC18_QTOF_20min | Waters ACQUITY UPLC HSS T3 | 183 | 17 | 2.63 | 1.67 | 0.96 | 2.38 | 1.62 | 0.76 | 0.32 | 2.81 | 1.93 | 0.88 | 108 |
| 0419 | ROASMI_59 | Waters SunFire C18 | 155 | 8 | 0.83 | 0.55 | 0.27 | 0.82 | 0.54 | 0.28 | 0.34 | 0.80 | 0.54 | 0.25 | 121 |
| 0409 | ROASMI_46 | Waters Atlantis T3 | 91 | 29 | 1.36 | 0.68 | 0.68 | 1.09 | 0.71 | 0.38 | 0.35 | 1.69 | 1.16 | 0.53 | 8 |
| 0418 | ROASMI_41 | Agilent ZORBAX Eclipse XDB-C18 | 78 | 38 | 2.20 | 1.36 | 0.84 | 2.09 | 1.32 | 0.77 | 0.37 | 0.71 | 0.51 | 0.20 | 6 |
| 0042 | CS24 | Phenomenex Kinetex EVO C18 | 51 | 12 | 0.97 | 0.58 | 0.39 | 0.83 | 0.51 | 0.32 | 0.38 | 0.91 | 0.58 | 0.34 | 12 |
| 0412 | ROASMI_61 | Thermo Scientific Hypersil GOLD aQ | 173 | 12 | 1.52 | 0.83 | 0.70 | 1.28 | 0.76 | 0.52 | 0.40 | 1.40 | 1.00 | 0.40 | 15 |
| 0410 | ROASMI_30 | Waters Atlantis T3 | 66 | 29 | 1.71 | 1.06 | 0.65 | 1.69 | 0.98 | 0.71 | 0.41 |  |  |  |  |
| 0015 | MTBLS36 | Waters ACQUITY UPLC HSS T3 | 72 | 25 | 2.72 | 1.47 | 1.25 | 2.90 | 1.60 | 1.30 | 0.45 | 3.68 | 2.33 | 1.35 | 15 |
| 0420 | ROASMI_60 | Waters SunFire C18 | 169 | 30 | 4.32 | 2.11 | 2.21 | 3.32 | 1.82 | 1.50 | 0.45 | 3.37 | 1.87 | 1.50 | 110 |
| 0421 | ROASMI_54 | Waters SunFire C18 | 84 | 25 | 3.14 | 1.72 | 1.41 | 3.05 | 1.64 | 1.42 | 0.46 | 3.51 | 2.07 | 1.44 | 10 |
| 0039 | CS21 | Merck LiChrospher RP-18 | 50 | 135 | 8.08 | 4.56 | 3.52 | 7.23 | 3.42 | 3.81 | 0.53 |  |  |  |  |
| 0405 | ROASMI_51 | Waters ACQUITY UPLC BEH C8 | 84 | 28 | 2.70 | 1.37 | 1.33 | 3.01 | 1.40 | 1.61 | 0.54 | 2.87 | 0.80 | 2.07 | 21 |
| 0411 | ROASMI_18 | Thermo Scientific Hypersil GOLD aQ | 50 | 10 | 1.38 | 0.47 | 0.91 | 1.03 | 0.38 | 0.65 | 0.63 | 1.04 | 0.41 | 0.62 | 11 |
| 0432 | ROASMI_65 | Merck Supelco Ascentis Express C18 | 216 | 25 | 8.83 | 2.55 | 6.29 | 9.63 | 2.01 | 7.62 | 0.79 | 10.35 | 1.19 | 9.16 | 30 |
| | **median** | | | | 1.91 | 1.36 | 0.70 | 1.72 | 1.35 | 0.65 | 0.34 | 2.01 | 1.24 | 0.50 | |

**HILIC** (spread anchors, minutes; order = O1, warp = I1 - O1, total = I1; purged columns blank when the dataset has < 20 compounds outside the training pool)

| id | dataset | column | n | grad (min) | k=15 total | k=15 order | k=15 warp | k=30 total | k=30 order | k=30 warp | k=30 warp share | purged k=30 total | purged k=30 order | purged k=30 warp | purged n_q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0231 | HPOS_ROI_V_3_1_1 | Waters ACQUITY UPLC BEH HILIC | 115 | 12.65 | 1.09 | 0.89 | 0.20 | 1.06 | 0.92 | 0.14 | 0.13 | 1.04 | 0.87 | 0.17 | 36 |
| 0183 | aHILIC Geneva DB maXisQTOF pos | Waters ACQUITY UPLC BEH Amide | 195 | 31 | 1.86 | 1.42 | 0.44 | 1.67 | 1.44 | 0.23 | 0.14 | 2.25 | 1.69 | 0.56 | 18 |
| 0044 | HILIC_BDD_2 | Merck SeQuant ZIC-HILIC | 453 | 25 | 3.65 | 2.75 | 0.90 | 3.40 | 2.76 | 0.65 | 0.19 | 3.77 | 2.69 | 1.08 | 186 |
| 0228 | AnnieTestFiehn | Waters ACQUITY UPLC BEH Amide | 993 | 17 | 2.36 | 1.75 | 0.61 | 2.33 | 1.74 | 0.59 | 0.25 | 2.46 | 1.66 | 0.80 | 744 |
| 0027 | KI_GIAR_zic_HILIC_ph2_7 | Merck SeQuant ZIC-HILIC | 165 | 41.5 | 1.80 | 1.40 | 0.40 | 1.83 | 1.36 | 0.47 | 0.25 | 1.87 | 1.51 | 0.36 | 55 |
| 0225 | AjsTestF | Waters ACQUITY UPLC BEH Amide | 996 | 17 | 2.48 | 1.74 | 0.74 | 2.36 | 1.74 | 0.62 | 0.26 | 2.53 | 1.64 | 0.89 | 744 |
| 0178 | HILIC19 neg | Merck SeQuant ZIC-HILIC | 85 | 40 | 1.47 | 0.99 | 0.48 | 1.49 | 1.05 | 0.44 | 0.29 | 2.58 | 2.72 | -0.14 | 6 |
| 0020 | Cao_HILIC | Merck SeQuant ZIC-pHILIC | 116 | 24 | 2.66 | 1.70 | 0.96 | 2.45 | 1.72 | 0.73 | 0.30 | 3.26 | 2.79 | 0.47 | 33 |
| 0185 | zHILIC Geneva DB VION QTOF neg | Merck SeQuant ZIC-pHILIC | 274 | 20 | 1.36 | 0.95 | 0.42 | 1.34 | 0.92 | 0.42 | 0.31 | 1.15 | 0.92 | 0.23 | 55 |
| 0111 | HILIC9 pos | Waters XBridge BEH Amide | 83 | 40 | 0.93 | 0.56 | 0.36 | 0.90 | 0.54 | 0.36 | 0.40 |  |  |  |  |
| 0007 | MTBLS87 | Merck SeQuant ZIC-pHILIC | 125 | 30 | 2.18 | 1.31 | 0.88 | 2.12 | 1.21 | 0.91 | 0.43 | 2.12 | 1.54 | 0.58 | 20 |
| | **median** | | | | 1.86 | 1.40 | 0.48 | 1.83 | 1.36 | 0.47 | 0.26 | 2.35 | 1.65 | 0.51 | |


**Reading the decomposition.** On the median RP-2step dataset at k = 30, total error is 1.72 min, order error is 1.35 min and warp error is 0.65 min. On purged compounds, total is 2.01 min, order 1.24 min and warp 0.50 min. Order error is roughly twice the warp error. **The remaining error lives mostly in the ordering, not the time axis.** The same holds for HILIC (order 1.36 min, warp 0.47 min) and for RP-proxy. It reverses only for PUB-2step, where the checkpoint has seen the datasets (order 0.35 min, warp 0.25 min). A few datasets are dominated by warp error, most clearly ROASMI_65 (I1 9.6 min vs O1 2.0 min on a 25 min gradient), where the non-negative-coefficient quadratic fits poorly. These cases drive the upper tail of H, not its median.

Retention order is *not* near-solved on systems the checkpoint never saw. Spearman between ROI and observed RT has a median of 0.83 on RP-2step (0.97 on the in-sample published targets, 0.71 for HILIC with the proxy). This is the only place Spearman is reported.

## Incumbent variants (I2 to I4)

None of the flexible maps beats the quadratic at k <= 30. Median G against I1, RP-2step purged with spread anchors: I2 isotonic -0.02 / -0.03, I3 linear -0.02 / -0.03, I4 monotone spline -0.00 / -0.02 at k = 15 / 30. This reproduces the 2-step authors' finding that flexible maps do no better than a quadratic at about 15 anchors. One exception, outside the decision range: at k = 50 (purged, spread) the monotone spline beats the quadratic (1.88 vs 2.13 min). The incumbent is unchanged for k <= 30.

## Step 5: is the time axis latent?

"Well documented" means t0 is present and at least 9 of the 10 tracked metadata fields are populated (`metadata_completeness.csv`, `step5_documentation_split.csv`). Purged, spread anchors:

| Group | k | H, well documented (n) | H, poorly documented (n) | Median Z1 / I1 per dataset | Median MAE I1 vs Z1 (min) |
|---|---|---|---|---|---|
| RP-2step | 15 | 0.32 (17) | 0.34 (7) | 1.15 | 1.98 vs 3.21 |
| RP-2step | 30 | 0.27 (17) | 0.28 (6) | 1.05 | 2.01 vs 3.05 |
| RP-proxy | 30 | 0.26 (19) | 0.28 (6) | 1.33 | 1.95 vs 3.27 |
| HILIC | 30 | 0.19 (5) | 0.28 (5) | 1.30 | 2.35 vs 3.01 |

- **RP: H does not shrink with documentation.** That is the pattern expected if the quantity that sets the time axis is latent and metadata cannot recover it. The poorly documented group is small (6 or 7 datasets), and no RP-2step dataset had recorded pH, so this is weak evidence.
- **HILIC: H is lower for the well-documented half** (0.19 vs 0.28, 5 vs 5 datasets). This suggests some time-axis information is observable there, but the sample is too small to lean on.
- **Z1 (metadata only, no anchors) is worse than the anchored incumbent in every group.** Its median MAE is 3.0 to 3.3 min, against 2.0 min for I1. The anchors carry information the recorded metadata does not. The per-dataset median ratio for RP-2step purged at k = 30 is close to 1 (1.05), so on roughly half the datasets the metadata model is competitive on novel compounds. On all compounds the ratio is 1.48.

## Step 6: order flips

A flip is a compound pair with |delta RT| >= 10 s whose observed order disagrees with the ROI order (`step6_order_flips.csv`). Recovery is measured on query-query pairs at k = 30, spread anchors, averaged over seeds and then taking the median over datasets.

| Group | Median flip fraction | Flipped pairs ordered correctly: I1 / T1 / T2 / R1 | Non-flipped pairs kept: I1 / T1 / T2 / R1 | Overall pair accuracy: I1 / T1 / T2 |
|---|---|---|---|---|
| RP-2step | 15.6% | 0.00 / 0.40 / 0.41 / 0.28 | 1.00 / 0.95 / 0.95 / 0.94 | 0.82 / 0.86 / 0.86 |
| RP-proxy | 15.0% | 0.00 / 0.38 / 0.38 / 0.22 | 1.00 / 0.95 / 0.95 / 0.96 | 0.84 / 0.85 / 0.85 |
| HILIC | 22.6% | 0.00 / 0.33 / 0.32 / 0.16 | 1.00 / 0.95 / 0.95 / 0.95 | 0.76 / 0.78 / 0.78 |
| PUB-2step (in-sample) | 4.8% | 0.00 / 0.21 / 0.20 / 0.23 | 1.00 / 0.98 / 0.98 / 0.98 | 0.95 / 0.94 / 0.94 |

**Yes, the anchor-using in-context arms recover some flips.** No monotone warp can. TabPFN puts about 40% of flipped RP pairs in the right order, at the cost of breaking about 5% of correctly ordered pairs, for a net gain of about 3.5 points of pairwise accuracy. This is where most of the in-context gain in G comes from. The in-context arms are partly re-ordering compounds using the support set, not learning a better time axis. That is the one mechanism the spec identified by which an in-context model could beat a pure warp, and it is real but small. On the in-sample published targets, where the ordering is already good, it does not pay.

## Published targets (indicative only)

I1 error as % of gradient length (MAE on non-anchor compounds, spread anchors), next to the published 2-step numbers. **PUB-2step is in-sample for the checkpoint** and is therefore optimistic. PUB-proxy is clean but uses a much weaker, condition-unaware ROI. Neither is a reproduction of the published benchmark.

| Dataset | Published | PUB-2step k=15 | PUB-2step k=30 | PUB-proxy k=15 | PUB-proxy k=30 | PUB-2step k=30, minutes |
|---|---|---|---|---|---|---|
| FEM long (0002) | 5.59 | 5.29 | 4.50 | 10.44 | 10.48 | 2.70 |
| IPB Halle (0003) | 1.78 | 2.05 | 1.66 | 6.58 | 5.99 | 0.33 |
| UniToyama Atlantis (0018) | 6.28 | 3.53 | 2.61 | 7.89 | 6.91 | 1.04 |
| Eawag XBridgeC18 (0019) | 4.68 | 4.00 | 2.98 | 4.73 | 3.96 | 0.90 |
| LIFE old (0054) | 3.12 | 3.63 | 3.42 | 7.82 | 7.93 | 0.21 |
| LIFE new (0055) | 5.59 | 6.31 | 6.28 | 14.74 | 13.45 | 0.44 |

The released pipeline lands in the published range, so the I1 implementation behaves as intended. On the 60 min FEM gradient, even this in-sample setting leaves about 2.7 min MAE, far from the "seconds to a minute" practitioners need. On unseen RP systems the incumbent's median is 8 to 10% of the gradient (1.7 to 2.0 min).

## Figures

- `figures/error_vs_k_purged.png` and `figures/error_vs_k_all.png`: median over datasets of MAE (minutes) versus k, one line per arm, one panel per held-out group, spread anchors. PUB-2step has no purged panel by construction.
- `figures/headroom_per_dataset.png`: H per dataset at k = 15 and k = 30 (all compounds), with k = 30 purged as diamonds. Dashed lines mark 0.15 and 0.25.

## Timing (4 CPU cores, no GPU)

| Component | Cost |
|---|---|
| 2-step checkpoint ROI inference | 234 s for 60,177 compounds over 331 RP datasets (about 4 ms per compound) |
| Proxy ROI training (LightGBM) | RP 26 s (45k rows), HILIC 2.4 s (3.5k rows), plus 5-fold out-of-fold refits |
| Z1 training (LightGBM) | 0.7 to 2.2 s per model |
| I1 quadratic (LAD + OLS, CBC solver) | 0.014 s per fit |
| I2 / I3 / I4 / O1 / R1 | 0.0001 to 0.003 s per fit |
| T1 / T2 TabPFN v2 fit + predict (1 thread) | 2.3 to 4.3 s per call; most of the 6,661 s wall-clock for the full grid |

Per-group, per-arm figures: `timing_per_fit_seconds.csv`, `timing_total_seconds.csv`, `timing_roi_2step.csv`, `timing_roi_proxy.csv`, `timing_z1.csv`.

## What this means for the investment decision

- **By the rule:** CONDITIONAL for RP, INDETERMINATE for HILIC, GO nowhere. Generic in-context inference (TabPFN v2, ridge residuals) captures 5 to 10% of the incumbent's error at the decision points, against a required 15%, and its CI includes zero at k = 30.
- **Beyond the rule:** the headroom the rule measures is about half warp and half in-sample oracle optimism (H falls from 0.27 to 0.14 at k = 30 with a cross-fitted oracle). The larger error component on unseen systems is retention *order* (Spearman 0.83, 16% of pairs flipped). A purpose-built model aimed at the time axis alone would chase the smaller half of the problem. If anything is worth building, it is one that uses the anchors to correct order on the target system, which is what the in-context arms did partially.

## Files

- `report.md`: this report.
- `prereg_interpretation.md`: operational choices frozen before running, plus amendment 1. Both committed before any error was computed.
- `results.csv`: every run. Columns: group, dataset, roi_source, scheme, k, seed, arm, eval_set, n_query, MAE, median AE, MAE as % of gradient.
- `decision_table.csv`: H, G per arm, CIs and verdict for every group, eval set, scheme and k.
- `error_decomposition.csv`, `step5_documentation_split.csv`, `step6_order_flips.csv`.
- `provenance.csv`: every dataset used, with RepoRT id, n compounds, column, mode, gradient length, t0 presence, metadata fields populated, pH used, role.
- `metadata_completeness.csv`, `overlap_summary.csv`, `ph_imputation_table.csv`, `2step_train_sets.json`.
- `code/`: survey, prep, ROI (checkpoint and proxy), arms, runner, analysis, provenance.

Provenance of external artifacts: RepoRT `420b550eca799be282a69c5cd7c2911a1593cbfa`; 2-step `41a423392bad5478d3cf053edd7f17a2ca8c347c`, checkpoint `models/2-step0525.pt` (trained from 2-step commit `341af06`, RepoRT version `94f43c1b` per the README); TabPFN v2 regressor `tabpfn-v2-regressor.ckpt` from Prior-Labs/TabPFN-v2-reg, run with the `tabpfn` 10.2.0 package and `ModelVersion.V2`.

## Appendix: decomposition for the other groups (minutes)
**PUB-2step** (spread anchors, minutes; order = O1, warp = I1 - O1, total = I1; purged columns blank when the dataset has < 20 compounds outside the training pool)

| id | dataset | column | n | grad (min) | k=15 total | k=15 order | k=15 warp | k=30 total | k=30 order | k=30 warp | k=30 warp share | purged k=30 total | purged k=30 order | purged k=30 warp | purged n_q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0019 | Eawag_XBridgeC18 | Waters XBridge C18 | 355 | 30 | 1.20 | 0.69 | 0.51 | 0.89 | 0.68 | 0.22 | 0.23 |  |  |  |  |
| 0054 | LIFE_old | Waters ACQUITY UPLC BEH C18 | 172 | 6 | 0.22 | 0.14 | 0.08 | 0.21 | 0.14 | 0.07 | 0.33 |  |  |  |  |
| 0002 | FEM_long | Waters ACQUITY UPLC HSS T3 | 346 | 60 | 3.17 | 1.76 | 1.42 | 2.70 | 1.75 | 0.95 | 0.35 |  |  |  |  |
| 0018 | UniToyama_Atlantis | Waters Atlantis T3 | 78 | 40 | 1.41 | 0.62 | 0.79 | 1.04 | 0.56 | 0.49 | 0.46 |  |  |  |  |
| 0003 | IPB_Halle | Waters ACQUITY UPLC HSS T3 | 80 | 20 | 0.41 | 0.14 | 0.27 | 0.33 | 0.12 | 0.21 | 0.64 |  |  |  |  |
| 0055 | LIFE_new | Waters ACQUITY UPLC HSS T3 | 162 | 7 | 0.44 | 0.15 | 0.29 | 0.44 | 0.15 | 0.29 | 0.66 |  |  |  |  |
| | **median** | | | | 0.82 | 0.38 | 0.40 | 0.67 | 0.35 | 0.25 | 0.40 | | | | |

**PUB-proxy** (spread anchors, minutes; order = O1, warp = I1 - O1, total = I1; purged columns blank when the dataset has < 20 compounds outside the training pool)

| id | dataset | column | n | grad (min) | k=15 total | k=15 order | k=15 warp | k=30 total | k=30 order | k=30 warp | k=30 warp share | purged k=30 total | purged k=30 order | purged k=30 warp | purged n_q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0002 | FEM_long | Waters ACQUITY UPLC HSS T3 | 346 | 60 | 6.27 | 5.36 | 0.91 | 6.29 | 5.36 | 0.92 | 0.15 | 6.02 | 5.52 | 0.50 | 57 |
| 0019 | Eawag_XBridgeC18 | Waters XBridge C18 | 355 | 30 | 1.42 | 0.92 | 0.49 | 1.19 | 0.92 | 0.27 | 0.22 | 1.45 | 1.34 | 0.11 | 45 |
| 0003 | IPB_Halle | Waters ACQUITY UPLC HSS T3 | 80 | 20 | 1.32 | 0.89 | 0.43 | 1.20 | 0.91 | 0.28 | 0.24 | 1.48 | 0.90 | 0.58 | 10 |
| 0054 | LIFE_old | Waters ACQUITY UPLC BEH C18 | 172 | 6 | 0.47 | 0.36 | 0.11 | 0.48 | 0.35 | 0.12 | 0.25 | 0.51 | 0.41 | 0.11 | 35 |
| 0018 | UniToyama_Atlantis | Waters Atlantis T3 | 78 | 40 | 3.16 | 1.89 | 1.26 | 2.76 | 1.67 | 1.09 | 0.39 | 1.55 | 1.46 | 0.09 | 15 |
| 0055 | LIFE_new | Waters ACQUITY UPLC HSS T3 | 162 | 7 | 1.03 | 0.52 | 0.51 | 0.94 | 0.51 | 0.43 | 0.45 | 1.21 | 0.63 | 0.58 | 14 |
| | **median** | | | | 1.37 | 0.91 | 0.50 | 1.19 | 0.92 | 0.36 | 0.24 | 1.47 | 1.12 | 0.30 | |

**RP-proxy** (spread anchors, minutes; order = O1, warp = I1 - O1, total = I1; purged columns blank when the dataset has < 20 compounds outside the training pool)

| id | dataset | column | n | grad (min) | k=15 total | k=15 order | k=15 warp | k=30 total | k=30 order | k=30 warp | k=30 warp share | purged k=30 total | purged k=30 order | purged k=30 warp | purged n_q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0382 | AA_LC1 | Thermo Scientific Acclaim RSLC 120 C18 | 2143 | 20.1 | 1.07 | 0.82 | 0.25 | 0.92 | 0.81 | 0.11 | 0.11 | 1.29 | 1.18 | 0.10 | 406 |
| 0400 | ROASMI_79 | Waters ACQUITY UPLC HSS C18 | 305 | 5.4 | 0.54 | 0.47 | 0.07 | 0.55 | 0.48 | 0.07 | 0.12 | 0.53 | 0.47 | 0.05 | 14 |
| 0383 | AA_LC2 | Thermo Scientific Acclaim RSLC 120 C18 | 1587 | 20.1 | 1.08 | 0.88 | 0.20 | 1.05 | 0.88 | 0.17 | 0.16 | 1.37 | 1.24 | 0.13 | 366 |
| 0262 | SPECIMEN_Pos | Waters ACQUITY UPLC BEH C18 | 883 | 30 | 1.50 | 1.30 | 0.20 | 1.61 | 1.30 | 0.31 | 0.19 | 1.89 | 1.79 | 0.10 | 239 |
| 0414 | ROASMI_67 | Agilent ZORBAX Eclipse Plus C18 | 225 | 28 | 1.88 | 1.48 | 0.40 | 1.87 | 1.48 | 0.39 | 0.21 | 2.26 | 1.92 | 0.33 | 8 |
| 0392 | UIUC_Carver | Waters ACQUITY Premier BEH C18 | 433 | 14.5 | 0.82 | 0.66 | 0.16 | 0.83 | 0.66 | 0.18 | 0.21 | 0.85 | 0.68 | 0.17 | 41 |
| 0047 | Acquity HSST3-RP-60mn | Waters ACQUITY UPLC HSS T3 | 159 | 60 | 4.95 | 3.66 | 1.29 | 4.64 | 3.61 | 1.03 | 0.22 | 7.10 | 5.31 | 1.78 | 7 |
| 0405 | ROASMI_51 | Waters ACQUITY UPLC BEH C8 | 84 | 28 | 2.42 | 1.79 | 0.63 | 2.32 | 1.80 | 0.52 | 0.22 | 2.99 | 2.29 | 0.70 | 12 |
| 0419 | ROASMI_59 | Waters SunFire C18 | 155 | 8 | 0.73 | 0.55 | 0.17 | 0.70 | 0.54 | 0.16 | 0.23 | 0.70 | 0.55 | 0.16 | 121 |
| 0427 | ROASMI_50 | GL Sciences Inertsil ODS-4 | 98 | 50 | 2.94 | 2.29 | 0.65 | 2.94 | 2.23 | 0.70 | 0.24 | 2.72 | 1.41 | 1.30 | 9 |
| 0263 | SPECIMEN_Neg | Waters ACQUITY UPLC BEH C18 | 725 | 30 | 2.11 | 1.75 | 0.36 | 2.34 | 1.74 | 0.60 | 0.24 | 2.06 | 1.79 | 0.28 | 184 |
| 0415 | ROASMI_83 | Agilent ZORBAX Eclipse Plus C18 | 650 | 22 | 1.57 | 1.16 | 0.40 | 1.53 | 1.14 | 0.38 | 0.25 | 1.73 | 1.37 | 0.36 | 274 |
| 0434 | ROASMI_71 | Waters XSelect HSS T3 XP | 192 | 25 | 2.02 | 1.55 | 0.47 | 2.05 | 1.51 | 0.54 | 0.25 | 2.95 | 2.65 | 0.30 | 11 |
| 0438 | Ion Pairing C18 HSS T3 | Waters ACQUITY Premier HSS T3 | 180 | 40 | 3.51 | 2.50 | 1.01 | 3.31 | 2.46 | 0.85 | 0.26 | 3.30 | 2.39 | 0.90 | 54 |
| 0404 | ROASMI_49 | Thermo Scientific Accucore C18 | 97 | 27.5 | 2.07 | 1.52 | 0.55 | 2.07 | 1.47 | 0.61 | 0.29 |  |  |  |  |
| 0428 | ROASMI_62 | Phenomenex Kinetex C18 | 180 | 25 | 0.89 | 0.51 | 0.38 | 0.74 | 0.52 | 0.22 | 0.30 | 0.70 | 0.47 | 0.23 | 94 |
| 0437 | ROASMI_86 | Agilent ZORBAX 300 SB-C18 | 1953 | 100 | 1.91 | 1.35 | 0.56 | 1.97 | 1.34 | 0.63 | 0.31 | 1.95 | 1.35 | 0.60 | 1722 |
| 0057 | CS5 | Agilent ZORBAX RRHD Eclipse Plus C18 | 82 | 12 | 1.46 | 0.93 | 0.54 | 1.31 | 0.90 | 0.41 | 0.31 | 1.30 | 0.85 | 0.45 | 19 |
| 0433 | ROASMI_66 | Merck Supelco Ascentis Express C18 | 119 | 25 | 2.17 | 1.40 | 0.77 | 1.81 | 1.25 | 0.57 | 0.31 | 2.31 | 1.87 | 0.44 | 8 |
| 0435 | ROASMI_77 | Waters XSelect HSS T3 XP | 178 | 25 | 6.90 | 4.90 | 2.00 | 6.88 | 4.66 | 2.22 | 0.32 | 7.07 | 4.23 | 2.84 | 22 |
| 0056 | CS4 | Agilent ZORBAX RRHD Eclipse Plus C18 | 89 | 12 | 1.11 | 0.68 | 0.44 | 1.07 | 0.72 | 0.35 | 0.33 | 0.87 | 0.63 | 0.23 | 24 |
| 0261 | Synergi_Neg | Phenomenex Synergi Polar-RP | 103 | 23 | 1.21 | 0.87 | 0.34 | 1.22 | 0.81 | 0.41 | 0.33 | 1.27 | 0.77 | 0.49 | 29 |
| 0318 | MSMS_pos_30_MeOH_Raptor Biphenyl | Waters CORTECS UPLC C18 | 407 | 10.01 | 0.43 | 0.26 | 0.18 | 0.38 | 0.25 | 0.13 | 0.34 |  |  |  |  |
| 0320 | MSMS_pos_40_MeOH_Raptor Biphenyl | Waters CORTECS UPLC C18 | 394 | 10.01 | 0.44 | 0.26 | 0.18 | 0.41 | 0.26 | 0.16 | 0.37 |  |  |  |  |
| 0338 | MSMS_pos_50_MeOH_CORTECS C18 | Restek Raptor Biphenyl | 391 | 10.01 | 0.58 | 0.32 | 0.26 | 0.52 | 0.31 | 0.21 | 0.41 |  |  |  |  |
| 0246 | MSMS_pos_40_ACN_XB-C18 | Phenomenex Kinetex XB-C18 | 429 | 10.01 | 0.49 | 0.26 | 0.23 | 0.45 | 0.26 | 0.19 | 0.41 |  |  |  |  |
| 0394 | ACQUITY_UPLC_HSSC18_QTOF_20min | Waters ACQUITY UPLC HSS T3 | 183 | 17 | 2.78 | 1.62 | 1.17 | 2.74 | 1.59 | 1.15 | 0.42 | 2.86 | 1.96 | 0.90 | 105 |
| 0248 | MSMS_pos_50_ACN_XB-C18 | Phenomenex Kinetex XB-C18 | 420 | 10.01 | 0.48 | 0.25 | 0.23 | 0.44 | 0.25 | 0.19 | 0.43 |  |  |  |  |
| 0336 | MSMS_pos_40_MeOH_CORTECS C18 | Restek Raptor Biphenyl | 394 | 10.01 | 0.62 | 0.33 | 0.30 | 0.58 | 0.32 | 0.26 | 0.45 |  |  |  |  |
| 0409 | ROASMI_46 | Waters Atlantis T3 | 91 | 29 | 1.26 | 0.75 | 0.51 | 1.33 | 0.65 | 0.68 | 0.49 |  |  |  |  |
| 0412 | ROASMI_61 | Thermo Scientific Hypersil GOLD aQ | 173 | 12 | 0.88 | 0.44 | 0.45 | 0.81 | 0.41 | 0.40 | 0.50 | 0.74 | 0.49 | 0.24 | 8 |
| 0420 | ROASMI_60 | Waters SunFire C18 | 169 | 30 | 4.20 | 2.18 | 2.01 | 4.07 | 1.85 | 2.22 | 0.54 | 4.42 | 1.93 | 2.49 | 110 |
| 0432 | ROASMI_65 | Merck Supelco Ascentis Express C18 | 216 | 25 | 8.35 | 2.54 | 5.80 | 8.73 | 2.00 | 6.73 | 0.77 | 8.88 | 1.35 | 7.53 | 27 |
| | **median** | | | | 1.46 | 0.93 | 0.40 | 1.33 | 0.90 | 0.40 | 0.31 | 1.95 | 1.37 | 0.36 | |
