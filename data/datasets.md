# Dataset selection (fixed before any experiment was run)

Source: Zenodo record 5171792 (DeepLC paper, Bouwmeester et al. 2021), file
`train_eval.zip` (118 MB, md5 `ac1b54af23bdafb1c90c9698eea4ec52`), folder
`train_eval/datasets/`. Each CSV has `seq`, `modifications` (PEPREC format,
`pos|name|pos|name`), `tr` (observed retention time, units differ per set).
The 34 GB `figures.zip` was not downloaded.

## Selected: 16 datasets

### Easy (9): standard acidic reversed-phase C18, mostly unmodified or common modifications only

| Dataset | Peptides | Reason for label |
|---|---|---|
| hela_hf_psms_aligned | 161,195 | Standard RP C18 (HeLa, Q Exactive HF). Only CAM / Met-ox / N-term acetyl. |
| hela_lumos_2h_psms_aligned | 57,889 | Standard RP C18, 2 h gradient. Common mods only (3% non-CAM). |
| pancreas_psms_aligned | 50,296 | Standard RP C18 tissue proteome. Common mods only. |
| yeast_120min_psms_aligned | 27,501 | Standard RP C18, 2 h gradient. Common mods only. |
| yeast_60min_psms_aligned | 18,506 | Standard RP C18, 1 h gradient. Common mods only. |
| arabidopsis_psms_aligned | 15,568 | Standard RP C18 plant proteome. Common mods only (Met-ox 16%). |
| proteometools_library | 156,120 | Standard RP C18 synthetic peptide library (ProteomeTools). Common mods only. |
| dia_fixed_mods | 146,587 | Standard RP C18 DIA spectral library. CAM only. |
| unmod_fixed_mods | 14,266 | Standard RP C18; unmodified half of the phospho study whose modified half is `mod_fixed_mods` (direct control for that hard set). |

### Hard (7): non-standard chromatography or modification-rich

| Dataset | Peptides | Reason for label |
|---|---|---|
| ATLANTIS_SILICA_fixed_mods | 39,091 | Bare-silica column (HILIC-type / normal-phase selectivity): elution order roughly inverted vs. RP. |
| LUNA_SILICA_fixed_mods | 36,823 | Bare-silica column, as above, different vendor. |
| LUNA_HILIC_fixed_mods | 36,080 | HILIC column: hydrophilic-interaction selectivity, unlike RP training data. |
| SCX_fixed_mods | 30,471 | Strong cation exchange: retention driven mainly by charge (K/R/H, N-term), not hydrophobicity. |
| Xbridge_fixed_mods | 40,037 | RP C18 at high pH (pH 10). Borderline: still RP, but charge states of D/E/H/K and termini change, so selectivity departs from the acidic-RP training setup. |
| prosit_ptm_2020 | 5,694 | Standard RP, but 84% of peptides carry non-CAM PTMs (13 types: methyl, di/tri-methyl, formyl, propionyl, succinyl, malonyl, crotonyl, nitro, ...), most absent from the generic model's training data. |
| mod_fixed_mods | 3,413 | Standard RP, but 66% of peptides phosphorylated. Smallest set: n = 1,000 context still leaves 2,413 queries. |

## Excluded (and why)

- **PXD005573_mcp**: the training set of the generic DeepLC model used as
  feature (see below). Evaluating on it would be circular.
- hela_lumos_1h, plasma_lumos_1h, plasma_lumos_2h: same LC platforms as
  hela_lumos_2h with fewer peptides (plasma sets have only 3.5k to 4.3k
  peptides); dropped to keep the run on CPU within budget, not because of
  their results (no experiment had been run).

## Known data quirks (decided before running)

- The five 2D-LC sets (silica x2, HILIC, SCX, Xbridge) report first-dimension
  RTs that are strongly discretized (most values are integers, about 2.3k to
  4.3k distinct values for 30k to 40k peptides). This puts a floor on the
  achievable error for every method equally.
- Several "aligned" sets contain a few far outliers (negative RTs from
  alignment). Kept as is. The relative metrics divide by the RT range of
  the whole dataset (max - min), as stated in the plan.
- Duplicated (seq, modifications) pairs (only 18 in PXD005573_mcp, none in
  selected sets) would be averaged.

## Generic RT predictor (decided before running)

DeepLC 4.5.0 (current PyPI) ships a default multitask model with 6,543 LC
setup heads. **All 20 datasets in this Zenodo record, including every
dataset above, are named heads of that model** (checked via
`model.task_names`), so its shared trunk was trained on the very
peptide/RT pairs we would use as queries, and DeepLC's built-in
(multi-head ridge) calibration would simply select the target's own head.
That would leak query labels into every method.

Therefore the generic, uncalibrated prediction is the **legacy DeepLC
ensemble `full_hc_PXD005573_pub_*` (3 models, averaged as DeepLC 3.x does
when uncalibrated)**, trained on PXD005573 only. These weights ship inside
the deeplc 4.5 wheel (as onnx2torch-converted `.pt` files); features are
computed with DeepLC 3.1.13's own `feat_extractor.py` (vendored in
`rtpoc/legacy_deeplc/`). Inference only.

## Erratum (added after computing generic predictions, before any C1-C5 run)

`Xbridge_fixed_mods` is not high-pH reversed phase as guessed above. Its RTs
correlate 0.97 (Spearman, 14,692 shared peptides) with `LUNA_HILIC_fixed_mods`,
hydrophobic peptides elute first, and correlation with the RP generic
predictor is negative (-0.40), so it is almost certainly a HILIC column
(Waters XBridge Amide). The label stays **hard**; only the reason changes,
and the borderline caveat no longer applies.

Generic predictor vs. observed RT (Pearson): easy sets +0.947 to +0.988;
prosit_ptm_2020 +0.906; mod_fixed_mods +0.914; the five 2D-LC sets -0.23 to
-0.40 (inverted / charge-driven selectivity).
