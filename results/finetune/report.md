# DeepLC 4.5 fine-tuning vs TabPFN stacking on the 8 held-out setups

**Fine-tuning does not beat TabPFN stacking at any context size tested, so
there is no crossover in k = 100 to 1,000.** DeepLC 4.5 fine-tuned on the k
context peptides (FT) is worse than TabPFN on the DeepLC 4.5 heads (T1+heads,
the old winner) at every k in both groups. On the 5 hard setups FT is 21 to 27%
worse (geometric mean over setups), and the 95% CI over seeds lies above 1 on
19 of 20 (setup, k) cells. On the 3 easy setups the gap shrinks with k: 11.5%
at k = 100, 3.6% at 250, 4.8% at 500 and 1.4% at 1,000. At k = 1,000 FT is
significantly worse on 1 easy setup and indistinguishable on the other 2. On
easy data the curves are converging: fine-tuned on each setup's largest
available context (11,000 to 30,000 peptides), FT beats T1+heads at k = 1,000 on 2 of 3 easy
setups, so the easy-setup crossover lies between 1,000 and about 10,000 context
peptides (TabPFN itself was not run on contexts that large). On hard setups T1+heads at k = 1,000 still beats FT on the whole pool
on all 5 setups. Against DeepLC 4.5's default calibration (B4), FT
only pays off on easy setups from k = 500 (3% better at 500, 8% at 1,000); on
hard setups FT and B4 are tied on average at every k. The hybrid, TabPFN
stacked on the fine-tuned prediction (T1_d45_FT), is the best method overall:
5 to 7% better than T1+heads on easy setups at k >= 250 (CI below 1 on all 3
setups at k = 500 and 1,000), and identical to T1+heads on hard setups. Feeding
TabPFN the B4 calibrated prediction instead (T1_d45_B4, no fine-tuning) gives
nearly the same gain as the hybrid at a fraction of the cost.

![error vs k](curves_finetune.png)

## Main table: relative MAE, % of RT range

Mean over setups of the 5-seed mean. Lower is better. Per-setup mean and
standard deviation over the 5 seeds are in `tables.md` and `results.csv`.

| group | k | B4 calibration | FT fine-tuned | T1+heads | hybrid T1_d45_FT | T1_d45_B4 |
|---|---|---|---|---|---|---|
| easy | 100 | 0.504 | 0.534 | 0.479 | **0.468** | 0.472 |
| easy | 250 | 0.436 | 0.437 | 0.422 | 0.400 | **0.396** |
| easy | 500 | 0.426 | 0.412 | 0.393 | **0.372** | **0.372** |
| easy | 1000 | 0.413 | 0.379 | 0.375 | **0.351** | 0.356 |
| hard | 100 | 1.723 | 1.739 | **1.354** | 1.359 | 1.373 |
| hard | 250 | 1.498 | 1.491 | 1.219 | 1.211 | **1.207** |
| hard | 500 | 1.393 | 1.421 | 1.125 | 1.126 | **1.122** |
| hard | 1000 | 1.351 | 1.353 | 1.072 | 1.070 | **1.062** |

Geometric-mean ratio over setups (below 1 = the first method is better), and
the number of setups whose 95% CI over seeds (t-interval on the log ratio) lies
fully below / above 1:

| group | k | FT / B4 | FT / T1+heads | hybrid / FT | hybrid / T1+heads | T1_d45_B4 / T1+heads |
|---|---|---|---|---|---|---|
| easy (3) | 100 | 1.063 (0 / 0) | 1.115 (0 / 2) | 0.876 (3 / 0) | 0.977 (1 / 0) | 0.985 (0 / 0) |
| easy (3) | 250 | 1.004 (0 / 0) | 1.036 (0 / 0) | 0.915 (1 / 0) | 0.948 (1 / 0) | 0.936 (2 / 0) |
| easy (3) | 500 | 0.970 (0 / 0) | 1.048 (0 / 1) | 0.902 (3 / 0) | 0.945 (3 / 0) | 0.941 (3 / 0) |
| easy (3) | 1000 | 0.920 (2 / 0) | 1.014 (0 / 1) | 0.920 (3 / 0) | 0.933 (3 / 0) | 0.946 (2 / 0) |
| hard (5) | 100 | 1.007 (0 / 0) | 1.272 (0 / 4) | 0.788 (4 / 0) | 1.002 (0 / 0) | 1.012 (0 / 0) |
| hard (5) | 250 | 0.998 (1 / 0) | 1.210 (0 / 5) | 0.823 (5 / 0) | 0.996 (0 / 0) | 0.993 (0 / 0) |
| hard (5) | 500 | 1.007 (1 / 1) | 1.220 (0 / 5) | 0.818 (5 / 0) | 0.999 (0 / 0) | 0.997 (1 / 0) |
| hard (5) | 1000 | 0.987 (2 / 1) | 1.225 (0 / 5) | 0.814 (5 / 0) | 0.997 (0 / 0) | 0.992 (0 / 0) |

## Findings

1. **FT vs B4 (does fine-tuning beat DeepLC's own calibration?).** Only with
   enough data and only on easy setups. At k = 100 FT is 6% worse on easy
   setups and tied on hard ones; DeepLC itself warns below 500 reference
   peptides (`MIN_FINETUNE_REFERENCE = 500`; here the training part is 85% of
   k). At k = 1,000 FT is 8% better on easy setups (CI below 1 on PXD062983
   and PXD076643). On hard setups the picture is mixed: FT is better on
   PXD050641 (5 to 9% from k = 250) and PXD085246 at k = 1,000, but 13% worse
   on PXD084949 at k = 500 and 1,000 (CI above 1). On average over hard setups
   FT and B4 are within 1.3% at every k.
2. **FT vs T1+heads (the question asked).** T1+heads wins at every k in both
   groups. No setup has a CI for FT / T1+heads below 1 at any k. The hard-group
   gap does not close with k (1.27 at 100, 1.23 at 1,000). On easy setups it
   narrows to 1.4% at k = 1,000, where only PXD062983 is still significant.
3. **The hybrid (TabPFN on top of the fine-tuned model).** Since FT does not
   win, the hybrid is mainly a check of whether fine-tuning adds information
   that TabPFN cannot get from the frozen heads. On easy setups it does: the
   hybrid is 5 to 7% better than T1+heads at k >= 250 and the best method in
   the study. On hard setups it adds nothing (ratio 0.996 to 1.002, no CI away
   from 1). The hybrid is 8 to 21% better than FT itself everywhere, with the
   CI below 1 on 29 of 32 (setup, k) cells, so in-context calibration still
   adds a lot on top of a fine-tuned DeepLC.
4. **B4 prediction as an extra TabPFN feature (open item from the old
   study).** It helps on easy setups (6% at k = 250 to 1,000, CI below 1 on 2
   or 3 of 3 setups from k = 250) and does nothing on hard setups (within 1%,
   1 of 20 cells significant). It gets about the same gain as the hybrid on
   easy setups without any fine-tuning: T1_d45_B4 costs 34 to 94 s per run
   against 92 to 813 s for the hybrid. **Recommendation: make T1_d45_B4 the
   new default of `predict.py`** (not yet done here; the change is one extra
   column).
5. **Modified vs unmodified queries.** FT is relatively better on modified
   peptides on hard setups (1.88 vs 1.95 for B4 at k = 1,000), but still worse
   than T1+heads (1.75). On unmodified peptides on hard setups FT is worse than
   B4 at every k. See `tables.md`.
6. **Side observation, not part of the protocol.** The shipped
   `finetune_and_predict` path recalibrates the fine-tuned column on the
   context (`calibrate` with MultiHeadRidgeCalibration). Skipping that step
   (FT_raw, the fine-tuned head's raw output) is 2 to 5% better than FT in
   every cell (geometric mean). FT_raw still loses to T1+heads on hard setups
   by 17 to 21% at every k; on easy setups it is 2% better at k = 1,000
   (mostly PXD076643). This variant was not chosen beforehand, so treat it as
   a hint for DeepLC users, not as a headline number.
7. **Intervals.** Adding the FT or B4 column leaves TabPFN's 90% interval
   coverage (0.88 to 0.91) and width unchanged.

## Headroom: DeepLC 4.5 with no target labels vs fine-tuned on the largest context

Relative MAE in % of the RT range on the seed-0 queries. "No labels" is
`deeplc.predict` with the default model and no reference peptides: DeepLC then
reports its PXD005573_mcp setup, a 200-minute gradient. "Fine-tuned, largest
context" is fine-tuned on the whole non-query pool (85% train, 15% validation),
except PXD079751 (see below). Headroom = no labels minus fine-tuned. The k =
1,000 columns are 5-seed means from the grid, for reference.

| setup | group | no labels | fine-tuned, largest context | headroom (points) | context used | epochs (cap) | B4 on same context | FT k=1000 | T1+heads k=1000 |
|---|---|---|---|---|---|---|---|---|---|
| PXD062983 | easy | 30.91 | 0.31 | 30.60 | 24,941 | 27 (27) | 0.35 | 0.34 | 0.32 |
| PXD076643 | easy | 40.36 | 0.37 | 39.99 | 10,944 | 59 (61) | 0.41 | 0.40 | 0.41 |
| PXD079751 | easy | 19.57 | 0.38 | 19.18 | 30,000 of 90,685 | 22 (22) | 0.39 (all 90,685) | 0.40 | 0.40 |
| PXD050641 | hard | 36.92 | 1.51 | 35.41 | 12,275 | 54 (54) | 1.66 | 1.57 | 1.24 |
| PXD061731 | hard | 22.21 | 0.89 | 21.32 | 3,021 | 112 (150) | 2.03 | 0.88 | 0.78 |
| PXD074962 | hard | 25.82 | 1.20 | 24.62 | 2,006 | 88 (150) | 1.45 | 1.16 | 1.02 |
| PXD084949 | hard | 52.87 | 2.00 | 50.87 | 1,575 | 52 (150) | 1.98 | 2.14 | 1.41 |
| PXD085246 | hard | 19.93 | 1.00 | 18.93 | 3,800 | 73 (150) | 1.02 | 1.04 | 0.89 |

- The headroom is 19 to 51 points of the RT range. Almost all of it is the
  time axis: without labels DeepLC predicts on its own 200-minute gradient. The
  order of peptides is already good without labels (Pearson 0.96 to 0.995), so
  even a handful of labeled peptides recovers most of it, as the k = 100 rows
  of the main table show.
- **Easy setups: with their whole pool (11k to 30k peptides), fine-tuning beats
  T1+heads at k = 1,000** on 2 of 3 setups (PXD062983 0.31 vs 0.32,
  PXD076643 0.37 vs 0.41) and ties on the third. So the crossover on easy data
  sits somewhere between 1,000 and about 10,000 context peptides. TabPFN was
  not run at those sizes.
- **Hard setups: T1+heads at k = 1,000 is still better than fine-tuning on the
  whole pool** on all 5 setups (by 11 to 30%). Fine-tuning gains little from
  more data on these (PXD061731: 0.88 at k = 1,000, 0.89 with 3,021).
- PXD079751: fine-tuning on its 90,685-peptide pool was killed twice by the
  16 GB memory limit (about 14 GB resident, DeepLC's encoded training data),
  so it was fine-tuned on the first 30,000 pool peptides, the largest context
  that fits. The B4 column for this setup uses the whole pool.
- The three largest contexts had their epochs capped (27, 54, 22) to keep one
  fine-tune under 10 minutes on 4 threads; PXD076643 early-stopped before its
  cap. These caps were set by a fixed seconds-per-peptide-epoch budget before
  the runs, not tuned.
- Side finding: B4 on the whole pool of PXD061731 (3,021 peptides) is worse
  than B4 at k = 1,000 (2.03 vs 0.90, Pearson 0.61), so DeepLC's default
  calibration can fail with more data. Fine-tuning on the same pool did not.

## Timing

**4 threads, idle machine, seed 0, one fine-tune per (setup, k)**,
seconds for `deeplc.finetune` alone (calibrate and predict add under 2 s):

| setup | k = 100 | 250 | 500 | 1,000 |
|---|---|---|---|---|
| PXD062983 | 9 | 30 | 45 | 141 |
| PXD076643 | 8 | 29 | 51 | 85 |
| PXD079751 | 15 | 16 | 42 | 78 |
| PXD050641 | 10 | 29 | 66 | 78 |
| PXD061731 | 9 | 13 | 42 | 136 |
| PXD074962 | 6 | 12 | 71 | 156 |
| PXD084949 | 6 | 37 | 49 | 68 |
| PXD085246 | 5 | 26 | 43 | 84 |
| **median** | **9** | **27** | **47** | **84** |

Median epochs run: 96, 136, 110 and 100 at k = 100, 250, 500 and 1,000.

No fine-tune in the grid came near 10 minutes (longest 156 s on 4 threads),
so no epochs were reduced there. The headroom fine-tunes on 2,006 to 30,000
peptides took 64 to 587 s on 4 threads; the three largest ran with an epoch
cap, listed in the headroom table.

**Grid runs** (1 thread per job, 4 jobs on 4 cores), median seconds per
(setup, k, seed):

| | k = 100 | 250 | 500 | 1,000 |
|---|---|---|---|---|
| B4 calibration | 0.3 | 0.4 | 0.4 | 0.5 |
| FT fine-tune + calibrate + predict | 20 | 67 | 150 | 309 |
| T1_d45_B4 (TabPFN, 2,000 queries) | 34 | 45 | 58 | 94 |
| Hybrid T1_d45_FT (full fine-tune + 3 cross-fit fine-tunes + TabPFN) | 92 | 211 | 464 | 813 |

The whole study (grid, cross-fitting, timing, headroom) took about 7 hours
of wall-clock on the 4 cores.

## Protocol

- **Data, splits, seeds, metric:** reused unchanged from the old study
  (`rtpoc/run.py` `split`, contexts are the first k of the pool, up to 2,000
  queries, relative MAE divided by the setup's full RT range). B4 and T1+heads
  numbers are the stored ones from `results/parts/H_<setup>.csv`. B4 was also
  recomputed here (B4_check) and matches the stored values to the 4th digit on
  every cell, which confirms the splits are identical.
- **Fine-tuning API:** DeepLC 4.5.0's shipped transfer-learning path.
  `deeplc.finetune_and_predict` is `deeplc.finetune` followed by
  `deeplc.calibrate` on the reference; it does not accept a validation set, so
  `rtpoc/ft_d45.py` calls the same two public functions in the same order and
  passes our own 15% validation split to `finetune(psm_list_validation=...)`.
  There is no training loop of our own. For the default model DeepLC adds a new
  low-rank setup head (66 trainable parameters) and keeps the encoder and the
  pretrained heads frozen; it anchors the head by least squares before
  training, early-stops on validation L1 and restores the best weights.
- **Validation split:** 15% of the context (at least 2 peptides), drawn with a
  seed-dependent permutation; fine-tuning sees the other 85%. The final
  calibration step uses the whole context. Query peptides are never touched.
- **Hyperparameters, fixed once on Part 1 datasets** (yeast_120min_psms_aligned
  and prosit_ptm_2020, which are not among the held-out setups), k = 100 and
  1,000, seeds 0 and 1, geometric-mean relative MAE:

  | config | lr | max epochs | patience | geo-mean rel. MAE |
  |---|---|---|---|---|
  | A (DeepLC defaults) | 0.05 | 50 | 10 | 0.02045 |
  | B | 0.05 | 150 | 20 | 0.02009 |
  | C | 0.2 | 100 | 15 | 0.02188 |
  | **D (chosen)** | **0.01** | **150** | **20** | **0.01988** |

  Batch size and everything else are DeepLC defaults. The same config is used
  for every held-out setup and k.
- **Stacking without leakage:** TabPFN never sees an in-sample prediction. For
  the hybrid, the context column is 3-fold cross-fitted (each fold fine-tuned
  on the other 2/3 of the context, with its own 15% validation split) and the
  query column comes from the full-context fine-tune. For T1_d45_B4 the context
  column is 5-fold cross-fitted B4. 3 folds rather than 5 for the hybrid was a
  compute choice.
- **5 seeds** per (setup, k); `torch.manual_seed(seed)` is set before each
  fine-tune.
- **Grid runs** used 1 torch thread per job and 4 jobs on the 4 cores; the
  separate timing pass used 4 threads on an otherwise idle machine.

## Caveats and negative results

- The hybrid's gain on easy setups is real but small in absolute terms
  (0.024 percentage points of the RT range at k = 1,000, 1 to 2 s on these
  59- to 127-min gradients).
- Only 3 easy and 5 hard setups; the per-setup CIs are over 5 seeds of the
  same setup, not over setups.
- Config D hit its 150-epoch cap in 32% of grid fine-tunes (45% at k = 500
  and 1,000). A check on the Part 1 tuning datasets (k = 500 and 1,000, seeds
  0 and 1, `epochcap_*.csv`) reran D with a 500-epoch cap: only 2 of 8 runs
  went past 150 epochs and error changed by at most 0.3%. So the cap did not
  hold FT back in any way that matters here. The held-out grid was not rerun.
- The default model's fine-tuning only trains a new 66-parameter head. A
  deeper fine-tune (unfreezing the encoder) is not what DeepLC ships for this
  model and was not tried.
- The epoch counter in `rtpoc/ft_d45.py` counts one extra validation pass
  (DeepLC runs a post-training sanity check); the aggregate script subtracts
  it, so the epoch numbers in the tables are correct.

## Files

- `results.csv`: one row per (setup, k, seed, method, query subset), with
  `rel_mae_pct` = relative MAE in % of RT range.
- `tables.md`: all tables, including per-setup mean ± SD and CIs.
- `curves_finetune.png`: the plot above.
- `timing_4threads.csv`, `parts/HR_*.csv`: timing pass and headroom runs
  (epochs in these raw files include the off-by-one noted above).
- `epochcap_*.csv`: epoch-cap check on Part 1.
- `tune_*.csv`: hyperparameter selection on Part 1.
- Code: `rtpoc/ft_d45.py`, `rtpoc/ft_tune.py`, `rtpoc/run_ft.py`,
  `rtpoc/ft_timing.py`, `rtpoc/ft_headroom.py`, `rtpoc/ft_epochcap.py`,
  `rtpoc/ft_aggregate.py`, `rtpoc/run_ft_all.sh`, `rtpoc/run_ft_after.sh`,
  `rtpoc/run_ft_epochcap.sh`, `rtpoc/run_hr_big.sh`.
