For each held-out study (28 studies with ≥ 300 precursors), k ∈ {20, 50, 200, 1000} reference precursors drawn from that study (5 seeds), predict 1/K0 for up to 1,500 other precursors of the same study. Error = |prediction − study-level observed 1/K0|. Z1 = IM2Deep 2.0.2, default single-conformer model, CCS converted to 1/K0 with IM2Deep's `ccs2im` (the same constants MaxQuant uses). Features for R1/T1/O1: Z1, charge, m/z, mass, length, GRAVY, 20 amino-acid fractions, basic residues, C-terminal K/R, oxidation, N-terminal acetyl. R1 is fit on the residual (observed − Z1). T1 is TabPFN v2 (`tabpfn-v2-regressor.ckpt`, CPU, 4 estimators) on the same residual. O1 is a 5-fold cross-fitted LightGBM residual model on all of the study's precursors.

**Post-hoc control (not pre-registered), G1b:** a LightGBM residual model trained once on all *other* studies (target: observed − per-study affine-calibrated Z1), applied to the held-out study, followed by the same 2-parameter affine on the same k anchors. It is a one-time global retrain plus C2. It needs no in-context learning, and it tests whether any gain over C2 is study-specific. Anchors and queries are identical to the other arms (C2 reproduced to 0.0).

Median over the 26 studies that have results at every k (median over 5 seeds within study); in brackets, multiples of F = 0.0034:

{{T4_POOLED}}

![residual vs k](residual_vs_k.png)

Head-to-head against C2 on identical draws (all 28 studies):

{{T4_H2H}}

Readings:

- **IM2Deep as shipped is off by a per-study offset** (−0.042 to +0.043 Vs/cm²), so some calibration is mandatory. A one-parameter shift (C1) on 20 precursors removes almost all of what can be removed. C2 is no better than C1, and isotonic (C3) is worse at every k (it overfits).
- **Every arm, including the oracle, sits at 3.2 to 3.9 F**, above the published fine-tuned range of 0.007 to 0.010. What remains is IM2Deep's per-peptide error (and multi-conformer ambiguity), not study drift. The noise-limited cross-study residual from Step 3 (A = 0.0068) is already well below these prediction errors.
- **At the realistic k = 20 to 50, the learned in-context arms do not beat affine.** R1 ridge: +0.3 to +1.3 %. T1 TabPFN: +1.2 to +2.5 % (median relative improvement). They fall within anchor noise.
- **TabPFN helps only with many references.** It is +5.6 % at k = 200 and +10.6 % at k = 1000 (better than C2 in 98.5 % of draws at k = 1000). That is still 3.2 F.
- **Most of that gain is not study-specific.** The global retrain + affine control (G1b) is 7 to 8 % better than C2 at *every* k, including k = 20. It beats TabPFN at k ≤ 200, and TabPFN beats it by a median of only 2.6 % at k = 1000. So the study-specific, in-context-learnable part of the prediction error is at most ~3 %, and it only appears when 1,000 identified precursors from the target study are available. At that point per-study fine-tuning (current practice) is equally available.
- **O1 is not a strict ceiling.** The cross-fitted LightGBM "best within-study model" is beaten by TabPFN at k = 1000 (0.0111 vs 0.0115), so the true within-study ceiling is somewhat lower. It is still above 3 F.
- The relevant comparison in the brief (R1/T1 vs C2): a two-parameter affine map matches them at k ≤ 50, and a global model plus affine beats them. The direction is dead by the brief's own criterion.

Full per-study × k table: see the appendix at the end of this report and `step4_per_study_k.csv`.
