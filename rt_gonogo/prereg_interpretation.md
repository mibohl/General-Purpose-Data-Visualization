# Operational choices, frozen before any arm was evaluated

This file fixes every choice the go/no-go spec leaves open. It was written and
committed after inspecting only dataset metadata (sizes, columns, gradients,
which datasets the released 2-step checkpoint was trained on), and before any
retention time prediction or error was computed.

## ROI source

- Primary: released 2-step checkpoint `models/2-step0525.pt`, repo commit
  `41a423392bad5478d3cf053edd7f17a2ca8c347c` (boecker-lab/2-step). Its stored
  config lists the 171 RepoRT datasets it was trained on.
- Finding before any run: **all six published targets (0002, 0003, 0018, 0019,
  0054, 0055) are in that training set.** The checkpoint has seen their
  compounds and retention orders. They cannot be held-out systems for this
  checkpoint.
- Consequence: the RP headline uses RP datasets the checkpoint never saw. The
  six published targets are reported as a separate, flagged group, once with
  the checkpoint ROI (in-sample for ordering) and once with the proxy ROI
  (clean).
- Proxy ROI (spec fallback): LightGBM on 12 RDKit descriptors + ECFP4
  (radius 2, 1024 bits), target = within-dataset percentile rank of RT, pooled
  over all non-held-out datasets of the same mode. Used for all HILIC results
  (the checkpoint is RP-only) and as a clean second ROI for RP.
- When a column has no HSM/Tanaka entry, the 2-step code falls back to the
  descriptors of dataset 0045. This is the released tool's behaviour and is kept;
  the fallback is flagged per dataset.

## Data cleaning (per dataset)

- RepoRT commit `420b550eca799be282a69c5cd7c2911a1593cbfa`, isomeric
  `*_rtdata_isomeric_success.tsv` (as 2-step uses).
- Drop duplicate (SMILES, RT) rows; drop SMILES with more than one distinct RT
  (2-step `remove_doublets`).
- Void: drop compounds with RT < 2 * t0 (2-step `void_factor = 2`). t0 from
  `column.t0`; if missing or 0, estimated as 0.68 * pi * (id/2)^2 * L / F; if
  that is impossible, no void removal (flagged).
- Gradient length = last time point of `*_gradient.tsv`.
- SMRT (0186) excluded.

## Held-out selection (deterministic)

- RP-unseen (headline): RP, not in the checkpoint training set, gradient known,
  >= 80 compounds after cleaning, at most 2 datasets per recorded column name
  (largest first, ties by id). Cap raised to 3 only if this gives < 25.
- PUB (secondary): the six published targets.
- HILIC (secondary): HILIC, gradient known, >= 80 compounds after cleaning, at
  most 2 per column name (largest first). Cap raised only if < 10.
- Training pool for a held-out RP dataset = compounds of all RP datasets that
  are not held out (RP-unseen or PUB), plus the 171 checkpoint training sets.
  For HILIC = all non-held-out HILIC datasets.

## Anchors, arms, metrics

- k in {10, 15, 20, 30, 50}, seeds 0..4, schemes random and spread.
  Spread: stratified targets t_i = min + (i + u_i)/k * range with u_i ~ U(0,1),
  nearest unused compound per target.
- Anchors are drawn from all cleaned compounds; queries are all non-anchors.
- Purged evaluation = queries whose InChIKey first block is not in the
  training pool. Scaffold-purged = queries whose Bemis-Murcko scaffold is not in
  the pool. A dataset enters a purged analysis only if it has >= 20 compounds
  outside the pool (>= 10 query compounds remain at k = 50 in the worst case is
  not guaranteed, so runs with < 5 purged queries are dropped).
- Error = MAE in minutes on the query set; also 100 * MAE / gradient length.
  Median absolute error stored as a secondary column.
- I1 = 2-step `LADModel` exactly as in `predict.py` (LAD quadratic with
  non-negative coefficients, OLS refit on points with error <= 2 * median,
  OLS discarded if any coefficient is negative).
- I2 = isotonic (L2) on anchors, clipped extrapolation. I3 = OLS line.
  I4 = monotone piecewise-linear regression spline, knots at anchor ROI
  quantiles, n_segments = clip(k // 8, 2, 5), non-negative segment slopes,
  linear extrapolation.
- O1 = L1 isotonic regression of RT on ROI fitted on all cleaned compounds
  (MAE-optimal monotone map), scored on the same query set as the other arms.
  An L2 isotonic variant is stored for reference only.
- Z1 = LightGBM (500 trees, lr 0.05, 63 leaves, fixed) on [ROI, descriptors,
  metadata features], trained on the training-pool datasets, target = RT in
  minutes. No anchors.
- T1 = TabPFN v2 regressor (official v2 checkpoint, n_estimators = 8,
  random_state = seed) on [ROI, 12 descriptors, metadata features].
  T2 = T1 features + I1 prediction. Metadata features are constant within a
  dataset and therefore uninformative in context; included because the spec
  lists them.
- R1 = RidgeCV (alphas 10^-2..10^3, LOO on anchors only) on standardized
  [ROI, ROI^2, descriptors] predicting the I1 anchor residual; prediction =
  I1 + ridge.

## Decision quantities

- Per run: H_run = (MAE_I1 - MAE_O1) / MAE_I1; G_run(arm) =
  (MAE_I1 - MAE_arm) / MAE_I1. Per dataset: mean over the 5 seeds.
  Headline = median over held-out datasets.
- G = the arm among T1, T2, R1 with the highest median G in that cell. 95% CI
  = percentile bootstrap over datasets (10,000 resamples) of that arm's median.
- Decision cells: k = 15 and k = 30. Primary scheme = spread anchors
  (what a careful practitioner does); random = sensitivity.
- Headline = purged evaluation. If k = 15 and k = 30 disagree, the more
  conservative verdict is the headline (NO-GO > CONDITIONAL > GO).
- Verdicts: GO if H >= 0.25 and G >= 0.15 and CI excludes 0. NO-GO if
  H < 0.15. CONDITIONAL if H >= 0.25 and the G criterion fails.
  0.15 <= H < 0.25 is not covered by the rule; it is reported as
  "INDETERMINATE (rule gap)" with the numbers.

## Secondary analyses

- Step 5: well documented = t0 present in metadata and >= 9 of 10 tracked
  metadata fields populated. Compare median H between groups; Z1 vs I1.
- Step 6: pairs with |delta RT| >= 10 s and delta ROI != 0; a flip is a pair
  whose RT order disagrees with ROI order. Recovery measured at k = 30 (spread)
  on query-query pairs: fraction of flipped pairs each arm orders correctly,
  and fraction of non-flipped pairs it breaks.
- Spearman(ROI, RT) reported once, one line.
