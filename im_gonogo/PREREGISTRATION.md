# Pre-registered analysis choices (written before any evidence file was opened)

Timestamp: 2026-11-23, before MQ_output.zip finished downloading.

Decision rule: copied verbatim from the brief (NO-GO if A <= 2F; NO-GO if (A-O)/A < 0.25;
GO only if A > 2F and (A-O)/A >= 0.25 and Step 5 shows peptide-dependent residual).

Operational definitions fixed now:

- Study = one PRIDE project (one sub-folder of IM2Deep's Zenodo MQ_output.zip).
- PSM filter: MS/MS-identified rows only (Type contains "MSMS"; MaxQuant match-between-runs
  rows dropped), Reverse != "+", Potential contaminant != "+", PEP <= 0.01, 1/K0 > 0.
- Labels: rows with non-zero SILAC labeling state dropped; studies with TMT/iTRAQ reporter
  columns or dimethyl labels dropped (fixed label mods are invisible in "Modified sequence").
- Cysteine: MaxQuant omits fixed modifications from "Modified sequence", so a cysteine's
  alkylation state cannot be verified from the evidence table. All C-containing precursors
  are dropped unless the project's parameters file proves identical fixed modifications.
- Precursor identity = stripped sequence + canonicalized variable-mod string + charge.
  Allowed variable mods: Oxidation (M), Acetyl (Protein N-term). Any other mod -> dropped.
  Old and new MaxQuant notations ("(ox)", "(Oxidation (M))", "(ac)", "(Acetyl (Protein N-term))")
  are mapped to one canonical token.
- Per-run value: 1/K0 of the most intense MS/MS-identified evidence row of that precursor in that run.
- Per-study value: median of per-run values.
- F: within a study, for each (precursor, run) with >= 3 runs observing the precursor,
  |x_run - median(other runs)|; F_study = median of those; F pooled = median over all.
  Studies whose runs never share precursors are excluded from F.
- Consensus for precursor p and held-out study s = median over the other studies' values.
- D = median over (study, precursor) of |x_sp - consensus_sp| on shared precursors (>= 4 studies).
- A = same after per-study affine map x -> a + b x fit (OLS) on 20 random anchors, evaluated on
  non-anchor shared precursors; 5 seeds; median over seeds. A50 same with 50 anchors.
- O = same after per-study isotonic map fit on all shared precursors (in-sample, as the brief
  specifies). Sensitivity: 5-fold cross-fitted isotonic, and per-charge isotonic.
- All of the above repeated in CCS (Mason-Schamp, IM2Deep's im2ccs constants) and as relative
  deviation. Note: for a fixed precursor CCS = g(m/z, z) * (1/K0), so relative drift is identical
  in both spaces; only affine fitting across precursors differs.
- Multi-conformer flag: a precursor is flagged if, in any single run, two MS/MS-identified rows
  differ in 1/K0 by > 2 % relative, or if its study values split into two groups (each >= 2
  studies) separated by a gap > 2 % relative. Headline numbers reported with and without.
- Step 4 query set: all of the held-out study's precursors that are not anchors, capped at 3000
  randomly chosen per study and seed (fixed seed) for CPU cost; the cap is not tuned.
- O1 oracle: LightGBM on [Z1, features], 5-fold cross-fitted over all the study's precursors.
- Step 5: residual after affine (fit on all shared precursors of the study) regressed on
  charge, m/z, length, GRAVY hydrophobicity, 20 AA fractions; 5-fold CV R^2 for Ridge and LightGBM.
  Study-specificity test: within-study CV R^2 vs. leave-one-study-out pooled-model R^2.
