# In-context RT prediction with TabPFN v2: PoC summary

**Verdict: promising.** All three parts of the decision rule pass, with margin.
The large hard-set gain over B2 comes mostly from setups where the generic
predictor is nearly useless, so the more telling comparison is with the
classical residual model B3. T2 beats B3 there too, by about 10 to 12%.

**Follow-up (see the last section): against DeepLC 4.5's own calibration on
setups published after its release, TabPFN with the legacy generic feature
loses by 2 to 3x.** The generic model matters more than the in-context
learner. With the DeepLC 4.5 head predictions as features (T1_d45), TabPFN
beats DeepLC 4.5's calibration by 17 to 22% on hard and 3 to 9% on easy
held-out setups (n >= 100), and by 27 to 50% at n = 25. DeepLC 4.5 fine-tuned on the
same context does not beat it either, at any n from 100 to 1,000
(`results/finetune/report.md`).

Setup in one paragraph: 16 datasets from the DeepLC Zenodo record (9 easy,
7 hard; list and reasons fixed beforehand in `data/datasets.md`). For each
dataset and seed, 2,000 peptides are held out as queries, and the contexts
are nested random subsets of the rest (n = 25 to 1,000, 5 seeds). Every
method is fit on the context only. The generic feature is an uncalibrated
DeepLC prediction (see "Deviations" for which model and why). Metrics are
relative MAE and relative delta-t95, both divided by the dataset's RT range.
Full tables: `results/tables.md`. Raw rows: `results/results.csv`.
Figures: `results/curves_easy.png`, `results/curves_hard.png`,
`results/curves_per_dataset.png`.

## Decision rule

> Promising if, on hard datasets at n = 100-250, T2 lowers relative MAE by
> >= 10% vs. B2 and beats B3, while being no more than 2% worse than B2 on
> easy datasets.

How it is evaluated: per dataset, relative MAE is averaged over 5 seeds, and
ratios are then aggregated as geometric means over datasets. "Best B3" is the
better of ridge and XGBoost **chosen per seed after seeing query error**. That
is an oracle choice that favors B3, so it makes the test harder for T2.

| Condition | n = 100 | n = 250 | Pass? |
|---|---|---|---|
| Hard: T2 / B2 (need <= 0.90) | **0.412** (-59%) | **0.324** (-68%) | yes |
| Hard: T2 / best B3 (need < 1) | **0.898** (-10%) | **0.877** (-12%) | yes |
| Hard: T2 beats best B3, share of dataset x seed pairs | 89% | 100% | |
| Easy: T2 / B2 (need <= 1.02) | **0.896** (-10%) | **0.849** (-15%) | yes |
| Easy: worst single dataset, T2 / B2 | 0.990 (arabidopsis) | 0.972 (arabidopsis) | yes |

**Hard-set breakdown.** The hard group mixes two kinds of difficulty:

| Hard subgroup | n | T2 / B2 | T2 / best B3 | T2 wins vs best B3 |
|---|---|---|---|---|
| 2D-LC first dimension (silica x2, HILIC x2, SCX) | 100 | 0.323 | 0.877 | 96% |
| | 250 | 0.245 | 0.858 | 100% |
| RP with many PTMs (prosit_ptm_2020, mod_fixed_mods) | 100 | 0.757 | 0.954 | 70% |
| | 250 | 0.649 | 0.925 | 100% |

- **2D-LC sets.** These are not reversed-phase. The generic RP predictor
  correlates negatively with observed RT (Pearson -0.23 to -0.40), so B1 and
  B2 can do little: about 0.11 to 0.14 relative MAE, close to predicting a
  constant. Any model that learns from composition wins big here. Against B3,
  which learns from the same features, T2 is 12 to 14% better.
- **PTM-rich RP sets.** This is the more realistic "hard" case: a good
  generic model with gaps. T2 still clears the 10% bar against B2 on both
  datasets (prosit -13% / -26%, mod -35% / -43% at n = 100 / 250). Against
  best B3 the margin is small at n = 100 (prosit -1%, mod -9%) and grows at
  n = 250 (-7%, -10%).

## C1 Headroom (context = all but 2,000 queries, 3 seeds)

There is real signal left after B2 everywhere, so the experiment was
continued.

| Datasets | Residual R^2, B3 ridge | Residual R^2, B3 XGBoost | MAE gain vs B2 (XGBoost) |
|---|---|---|---|
| Easy (9) | 0.07 to 0.36 | 0.12 to 0.50 | 8% to 33% |
| Hard, 2D-LC (5) | 0.86 to 0.97 | 0.93 to 0.98 | 74% to 87% |
| prosit_ptm_2020 | 0.44 | 0.46 | 35% |
| mod_fixed_mods | 0.53 | 0.62 | 47% |

On easy sets the generic predictor already explains most of the variance, but
a useful setup-specific residual remains, largest on dia, hela_lumos_2h and
proteometools (R^2 0.4 to 0.5).

## C2 Main curve

Geometric-mean ratio of relative MAE to B2 (from `curves_*.png`, right panel):

| | n = 25 | 50 | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|---|---|
| Easy, T2 | 0.987 | 0.948 | 0.896 | 0.849 | 0.806 | 0.780 |
| Easy, T1 | 0.950 | 0.923 | 0.877 | 0.840 | 0.799 | 0.776 |
| Easy, best of B3 (ridge / XGBoost) | 1.03 | 1.04 | 0.98 | 0.95 | 0.90 | 0.86 |
| Hard, T2 | 0.598 | 0.489 | 0.412 | 0.324 | 0.280 | 0.248 |
| Hard, T1 | 0.488 | 0.407 | 0.340 | 0.286 | 0.255 | 0.236 |
| Hard, B3 XGBoost | 0.683 | 0.559 | 0.468 | 0.377 | 0.334 | 0.303 |

- TabPFN is never worse than B2 on average at any n. On easy sets it starts
  paying off at about 50 to 100 peptides, where B3 still loses to B2 (B3
  XGBoost needs about 250 just to break even).
- **T1 (direct RT) is as good as T2 or better (residual) almost
  everywhere.** On easy sets the two are tied within 2%. On hard sets T1 is
  clearly better (0.34 vs 0.41 at n = 100). This is expected on 2D-LC: there
  the B2 offset is mostly noise, so subtracting it gives TabPFN a noisier
  target. T2 wins only on modified peptides in easy sets (below).
- Relative delta-t95 follows the same ordering (easy, n = 250: B2 0.148,
  T2 0.131; hard: B2 0.532, B3 XGBoost 0.214, T2 0.192, T1 0.173).
- DeepLC's own spline calibration (B2_deeplc) tracks our LOWESS from n >= 100
  and is unstable at n = 25.

**Modified vs. unmodified queries** (mean relative MAE, n = 100 / 250).
"Modified" means any modification other than carbamidomethyl. Easy sets have
few such queries (about 60 to 360 per dataset).

| | B2 | best B3 (mean) | T1 | T2 |
|---|---|---|---|---|
| Easy, modified | 0.079 / 0.078 | 0.074 / 0.071 | 0.072 / 0.067 | **0.067 / 0.064** |
| Easy, unmodified | 0.027 / 0.026 | 0.026 / 0.025 | 0.024 / 0.022 | 0.025 / 0.022 |
| Hard, modified (prosit + mod only) | 0.073 / 0.071 | 0.060 / 0.053 | 0.056 / 0.047 | 0.059 / 0.050 |
| Hard, unmodified | 0.109 / 0.107 | 0.049 / 0.038 | 0.034 / 0.029 | 0.042 / 0.032 |

Modified peptides are the main error source on easy setups (about 3x the
unmodified error). T2 cuts that error by 15 to 19% at n = 100 to 250, against 8 to 13% for
unmodified peptides; there T1 and T2 are tied.

## C3 Label noise (n = 250; noisy RTs uniform in the context range; scored on clean labels)

Mean relative MAE over datasets:

| | clean | 5% | 10% | 20% |
|---|---|---|---|---|
| Easy, B2 | 0.0280 | 0.0281 | 0.0282 | 0.0282 |
| Easy, B3 ridge | 0.0268 | 0.0294 | 0.0325 | 0.0394 |
| Easy, B3 XGBoost | 0.0284 | 0.0416 | 0.0526 | 0.0705 |
| Easy, T1 | **0.0239** | 0.0248 | 0.0254 | 0.0277 |
| Easy, T2 | 0.0242 | **0.0245** | **0.0251** | **0.0256** |
| Hard, B2 | 0.1088 | 0.1092 | 0.1090 | 0.1093 |
| Hard, B3 ridge | 0.0446 | 0.0529 | 0.0588 | 0.0682 |
| Hard, B3 XGBoost | 0.0401 | 0.0514 | 0.0591 | 0.0726 |
| Hard, T1 | **0.0313** | **0.0339** | **0.0364** | **0.0414** |
| Hard, T2 | 0.0349 | 0.0374 | 0.0390 | 0.0431 |

T2 is fairly robust: +6% error at 20% noise on easy sets, +23% on hard sets.
B3 degrades much faster (XGBoost up to 2.5x on easy, where it ends far worse
than plain B2). B2's robust LOWESS barely moves. At 20% false IDs, T2 still
beats B2 on easy sets (-9%) and stays far ahead of B3 on hard sets.

T1 was run later on the identical corruptions. It degrades more than T2:
+16% on easy and +32% on hard sets at 20% noise. It stays ahead of T2 on hard
sets at every noise level. On easy sets T2 is better from 5% noise on,
because its LOWESS anchor is robust.

## C4 Uncertainty (90% intervals)

TabPFN uses its 5% and 95% predictive quantiles; for T2 they are shifted by
the B2 prediction. B2 uses a constant window of +/- the 90th percentile of the
cross-fitted absolute context residuals.

| | n | B2 coverage / width | T1 coverage / width | T2 coverage / width |
|---|---|---|---|---|
| Easy | 100 | 0.892 / 0.112 | 0.888 / 0.100 | 0.890 / 0.103 |
| Easy | 250 | 0.897 / 0.110 | 0.890 / 0.092 | 0.888 / 0.093 |
| Hard | 100 | 0.884 / 0.435 | 0.903 / 0.160 | 0.897 / 0.189 |
| Hard | 250 | 0.890 / 0.436 | 0.909 / 0.132 | 0.910 / 0.151 |

(Width is the mean interval width divided by RT range.) On easy sets TabPFN
intervals are about 1 point under nominal coverage (0.886 to 0.890) and 8 to
24% narrower than the B2 window (n = 100 to 1,000). They are also adaptive: their width varies
per peptide, while the B2 window is constant. On hard sets they meet nominal
coverage while being 2.3 to 2.9x (T2) or 2.7 to 3.3x (T1) narrower. At
n = 25 to 50, coverage is 0.86 to 0.90 for all methods.

## C5 Runtime (CPU)

Wall-clock seconds per fit + predict of 2,000 queries. TabPFN v2 runs with
n_estimators = 4 on 1 torch thread; 4 such jobs shared 4 cores. Figures are
the median (max) over datasets x seeds.

| n | 25 | 100 | 250 | 1000 |
|---|---|---|---|---|
| B2 (LOWESS incl. span CV) | <0.1 | 0.1 | 0.2 (0.4) | 0.8 (1.1) |
| B3 XGBoost (incl. B2 + cross-fit) | 0.1 | 0.3 | 0.9 (1.4) | 3.7 (4.4) |
| T1 | 20 (31) | 25 (39) | 36 (50) | 73 (93) |
| T2 (incl. B2 + cross-fit) | 21 (30) | 25 (38) | 37 (48) | 76 (94) |

That is about 3 to 4 minutes per TabPFN method per dataset per seed for
all 6 context sizes. Most of it is a fixed overhead at small n and grows with context and
query count. With 4 threads, a single fit was about 2.4x faster in a spot
check. For a full run (e.g. 50k queries) cost scales roughly linearly in
queries, so batching or a GPU would be needed for routine use. B2 is about 100 to
200x cheaper.

## Deviations from the plan and caveats

1. **Generic predictor is the legacy DeepLC ensemble, not DeepLC 4.5's
   default model.** All 20 Zenodo datasets are named heads of the 4.5
   multitask model, so using it would leak query labels into every method.
   The legacy `full_hc_PXD005573_pub` ensemble ships inside the 4.5 wheel. It
   was trained on PXD005573 only, which is excluded from the evaluation (see
   `data/datasets.md`). A stronger modern generic model would shrink the
   residual headroom on easy sets, so the easy-set gains here are likely an
   upper bound.
2. **TabPFN settings:** v2 regressor (ungated weights; newer default weights
   are license-gated), CPU, n_estimators = 4 (default 8) to fit the CPU
   budget, `ignore_pretraining_limits=True`. N- and C-terminal residues are
   passed as categorical features. The point prediction is the predictive
   median.
3. **delta-t95 definition:** 2 x the 95th percentile of |error|, divided by
   RT range. This is the symmetric-window convention; check it before
   comparing with numbers from the DeepLC transfer learning paper.
4. **B3 is trained on cross-fitted B2 residuals,** the same targets as T2,
   so B3 and T2 differ only in the learner.
5. **Data budget:** data and weights are 160 MB (Zenodo `train_eval.zip`
   118 MB, TabPFN weights 44 MB). However, `pip install tabpfn deeplc` pulled
   a CUDA build of torch plus NVIDIA libraries, which on its own exceeds the
   ~2 GB guideline. It is unused at runtime (CPU only).
6. **Xbridge label:** I first guessed high-pH RP, but its RTs correlate 0.97
   with LUNA_HILIC, so it is HILIC (erratum in `data/datasets.md`, written
   before any check ran). It stays hard.
7. **Hard sets are not independent.** The five 2D-LC sets share one study and
   largely the same peptides, and their first-dimension RTs are discretized.
   In effect there are three distinct hard situations: non-RP selectivity,
   phospho-rich RP, and multi-PTM RP. That is narrow evidence.
8. **The modification features are coarse** (counts per modification type,
   no positions), and TabPFN sees no sequence order. Gains on modified
   peptides come from learning per-setup shifts for modification counts.
9. C3 runs B2, B3 and T2, as the plan says; T1 was added in the follow-up.
10. Four smaller easy sets (hela_lumos_1h, plasma x2) were left out for
    runtime only; none were run.

## Follow-up: the suggested next steps

The first version of this summary ended with four next steps. All four were run
on the same CPU-only machine (4 cores; no GPU available). Tables are in
`results/tables.md` (sections "C3", "TabPFN n_estimators", "Per-dataset ratios",
"Held-out"), the held-out figure is `results/curves_heldout.png`, and timing
rows are in `results/timing.csv`.

### 1. Held-out comparison with DeepLC 4.5's own calibration

**Data.** There are 8 LC setups from PRIDE projects published 2026-09-21 to
2026-11-19, after DeepLC 4.5.0 was released on 2026-09-11. None of their PXD
accessions occurs among the 958 PXDs in DeepLC 4.5's training task names.
Each set is one run with RT taken as the per-run median, so there is no
alignment. The labels were fixed before any run (`data/datasets.md`, addendum;
details in `data/heldout/heldout.md`):

- 3 easy: tryptic, mostly unmodified.
- 5 hard: crotonyl-enriched, TMTpro phospho, two HLA-I sets, and an
  endogenous peptidome.

All 8 are reversed phase, because no new non-RP dataset could be found. Same
protocol as C2: nested contexts n = 25 to 1,000, 5 seeds, up to 2,000 queries,
every fit on the context only.

**New methods.**

- **B4**: the DeepLC 4.5 default model (6,543 setup heads) with its own
  `MultiHeadRidgeCalibration` fit on the context. This is what
  `predict_and_calibrate(..., psm_list_reference=context)` does. Its 90% window
  is +/- the 90th percentile of 5-fold cross-fitted context residuals.
- **T1_d45**: T1 plus the raw predictions of the 8 DeepLC 4.5 heads that B4
  ranks highest on the context. The ranking uses context labels only.

Geometric-mean ratio of relative MAE to B4 (below 1 = better than B4):

| | n=25 | 50 | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|---|---|
| Easy, T1 (legacy generic) | 1.88 | 2.34 | 3.10 | 3.22 | 3.12 | 3.13 |
| Easy, T1_d45 | **0.73** | **0.73** | **0.95** | **0.97** | **0.92** | **0.91** |
| Hard, T1 (legacy generic) | 1.38 | 2.48 | 2.41 | 2.38 | 2.34 | 2.22 |
| Hard, T1_d45 | **0.50** | **0.79** | **0.78** | **0.83** | **0.83** | **0.81** |

Mean relative MAE at n = 250 was as follows:

- Easy: B4 0.0044, T1_d45 0.0042, T1 0.0140, B2 0.0193.
- Hard: B4 0.0150, T1_d45 0.0122, T1 0.0358, B2 0.0574.

What this shows:

- **DeepLC 4.5 + its calibration is the baseline to beat, and the original
  TabPFN pipeline does not beat it.** With the legacy generic feature, every
  method is 2 to 7x worse than B4 once n >= 100. That includes T1, T2 and B3,
  and holds on both easy and hard sets. The main study's gains were relative
  to calibrating an old model. The newer generic model gives a much larger
  gain than switching calibration to in-context learning.
- **TabPFN on top of DeepLC 4.5 helps.** T1_d45 is better than B4 in every
  group and at every n. On hard sets the gain is 17 to 22% for n >= 50.
  Per dataset, the 95% CI of T1_d45/B4 lies below 1:
  - n = 250: 4 of 5 hard sets (PXD061731 is at 0.86 [0.73, 1.01]);
  - n = 100: 3 of 5.
  At n = 100 and 250, T1_d45 wins in 80 to 100% of seeds on every hard set. On easy sets the gain
  is 3 to 9% for n >= 100: clear on PXD062983 (0.88 to 0.89, CI below 1),
  and a tie on PXD076643 and PXD079751 (CI includes 1).
- **Small contexts: B4 is fragile, T1_d45 is not.** At n = 25, B4 has single
  seeds with relative MAE up to 0.35 (PXD085246) and a hard-set mean of
  0.047, against 0.020 for T1_d45. At n = 50, B4 recovers.
- **Where the gain is.** On hard sets T1_d45 helps unmodified peptides most:
  at n = 250, 0.0116 vs 0.0158 for B4, which covers the non-tryptic HLA and
  peptidome sets. On modified peptides the gain is 0.0183 vs 0.0205. On easy
  sets modified peptides are a tie.
- **Intervals.** T1_d45 coverage is 0.89 to 0.91 for n >= 100 on easy sets
  and 0.886 to 0.890 on hard sets, so 1 point under nominal on hard sets. At
  n = 25 it undercovers (0.85 to 0.86). On hard sets its intervals are
  20 to 35% narrower than B4's constant window (n >= 50); B4 covers 0.89
  to 0.90. On easy sets the widths are about equal.
- **Cost.** B4 takes 0.1 to 0.6 s per fit. T1_d45 takes 28 to 100 s for
  2,000 queries (1 thread, 4 jobs on 4 cores), so 100 to 300x more.

DeepLC's transfer learning (fine-tuning the network on the context) was not
run. The plan forbids neural-network training, and DeepLC's own documentation
recommends calibration over fine-tuning for contexts below about 500 peptides,
which is the range studied here. So B4 is the relevant baseline in this
regime, but it is not the strongest possible one for n = 500 to 1,000.

**Update: fine-tuning was run later, see `results/finetune/report.md`.**
DeepLC 4.5 fine-tuned on the context (its shipped API, settings fixed on
Part 1) does not beat T1_d45 at any k from 100 to 1,000: 21 to 27% worse on
hard and 1 to 12% worse on easy setups. Fine-tuned on each setup's largest
context (11,000 to 30,000 peptides), it beats T1_d45 at k = 1,000 on 2 of 3
easy setups but on none of the 5 hard ones. Adding the B4 prediction as an
extra TabPFN input lowers T1_d45's error by 5 to 6% on easy setups (k >= 250)
and leaves hard setups unchanged.

### 2. T1 as default, and a predictor CLI

T1 beats T2 in the main study on hard sets (per-dataset CIs below) and under
label noise on hard sets (C3 above). The held-out results add that the best
variant is T1 with DeepLC 4.5 head features. `rtpoc/predict.py` implements
it, with T1_d45 as the default and T1 / T2 as options. It takes a context CSV
(`seq, modifications, tr`) and a query CSV, and writes `rt_pred` and a 90%
interval. On PXD084949 (250 context peptides, 500 queries) it gave relative
MAE 0.017 and 89% coverage.

### 3. Per-dataset confidence intervals (main study, 5 seeds)

Ratio of relative MAE with a 95% t-interval over seeds on the log ratio
(full table in `results/tables.md`):

- T1 vs best B3 (per-seed min of ridge and XGBoost), n = 250: the CI lies
  below 1 on 15 of 16 datasets. The exception is arabidopsis
  (0.97 [0.94, 1.01]). On the non-RP hard sets the ratio is 0.64 to 0.79.
- T2 vs best B3 is smaller: 0.81 to 0.95 on hard sets at n = 250. At
  n = 100, prosit_ptm_2020 is not significant (0.995 [0.946, 1.047]).
- The two new independent hard groups, PTM-rich RP and HLA/peptidome RP, now
  come from the held-out sets in section 1. Non-RP remains covered by one
  study only (the 5 HILIC/SCX/silica sets from the same lab).

### 4. n_estimators = 8 and speed for large query sets

**n_estimators 8 vs 4** (all 16 datasets, seeds 0 to 2, n = 100 and 250):
the geometric-mean ratio of relative MAE is 0.995 to 1.000 in every group.
The best single dataset gains 0.8% (Xbridge, T1). Runtime doubles. Keep 4.

**Timing** (`rtpoc/timing.py`, hela_hf, 4 torch threads on an idle machine,
queries predicted in blocks of 10,000; pruned grid to fit the budget):

| context n | queries | n_estimators | ms / query | total |
|---|---|---|---|---|
| 250 | 2,000 | 4 | 7.3 | 15 s |
| 250 | 20,000 | 4 | 5.8 | 1.9 min |
| 250 | 100,000 | 4 | 6.9 | 11.6 min |
| 1,000 | 2,000 | 4 | 15.6 | 31 s |
| 1,000 | 20,000 | 4 | 9.8 | 3.3 min |
| 1,000 | 100,000 | 4 | 9.3 | 15.5 min |
| 250 | 20,000 | 8 | 15.0 | 5.0 min |
| 1,000 | 20,000 | 8 | 21.2 | 7.1 min |

- Cost is linear in queries, at about 6 to 10 ms per query on 4 CPU threads.
  A 100k-peptide library takes 12 to 16 minutes per setup. That is usable
  for offline library building, not for interactive use.
- `fit_with_cache` (context cache reused across query blocks) gave 9.3
  ms/query vs 9.8 at n = 1,000 / 20k queries, with the same accuracy. That
  is a 5% saving, not worth the extra 9 s fit and memory.
- n_estimators = 8 costs 2.2 to 2.6x. On this single split it lowered
  relative MAE by 3 to 4%. The 16-dataset NE8 test above (< 0.5%) is the
  better estimate.
- GPU: none is available in this environment (no CUDA device), so GPU
  speed was not measured.

### Revised recommendation

1. For a new LC setup, use **DeepLC 4.5 + TabPFN on its head predictions
   (T1_d45, `rtpoc/predict.py`)**. It was better than DeepLC 4.5's own
   calibration in both groups at every n, by about 20% on hard setups.
   At n = 250 it had lower error on 7 of 8 held-out datasets; the eighth is
   a tie. It is robust at n = 25 to 50, where the
   calibration can fail.
2. If runtime matters more than the last 5 to 20%, DeepLC 4.5's
   calibration (B4) is 100 to 300x cheaper and already 2 to 3x better than
   anything built on the legacy model.
3. Open questions:
   - non-RP held-out setups (none were found);
   - DeepLC fine-tuning at n >= 500, excluded here by the no-training rule;
   - adding B4's calibrated prediction itself as a T1_d45 feature;
   - GPU throughput.

Data budget: Zenodo 118 MB, TabPFN weights 43 MB (the DeepLC models ship
with the pip package), held-out downloads
352 MB kept (517 MB on disk after extraction in `data/heldout/raw/`, ignored
by git). That is under 2 GB in total.
