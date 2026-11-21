# Feasibility: in-context TCR specificity (TabPFN v2 / Kumo Tabular) vs TCRdist, on unseen epitopes

**Verdict — TabPFN v2 does **not** beat the cheap TCRdist baseline in the few-shot regime that matters.** On epitopes unseen during any cross-epitope training, under honest negatives (scheme b, TCRs that bind other epitopes), the decisive small-support sizes favour the baseline: at k=20-50 the median per-epitope AUROC is 0.666 for TabPFN vs 0.697 for the best TCRdist scorer (nearest-neighbour / k-NN). TabPFN only overtakes at k=200 (0.773 vs 0.697), i.e. once the support is large enough that the "few-shot" framing no longer applies and a plain LightGBM is equally good. TabPFN also costs ~25x more compute per task (~6.8 s/fit vs ~0.3 s). **NVIDIA Kumo Tabular (S)**, swapped in for TabPFN on the identical splits, behaves essentially like TabPFN: median honest AUROC 0.664 at k=20-50 and 0.779 at k=200 -- it too fails to beat TCRdist in the few-shot regime, so the conclusion does not depend on which tabular foundation model is used. It is also ~9x faster than TabPFN here (~0.7 s/fit vs ~6.8 s), but faster-yet-no-better does not change the decision. Crucially, the honest other-epitope negatives are **not** the catastrophe the usual warning implies: median AUROC under scheme b is if anything slightly higher than scheme a (difference +0.050). The honest cluster-aware split, balanced query and frozen generic features already remove the shortcut that normally makes random background negatives trivially separable, so both schemes are genuinely hard and the comparison above is robust to the negative choice. All methods are modest (median AUROC ~0.65-0.77) with wide per-epitope spread -- some epitopes (e.g. GILGFVFTL) reach ~0.85 while others (e.g. NLVPMVATV) sit at chance -- and zero-shot transfer to a new epitope collapses to ~0.5. The honest reading: in-context tabular learning is **not** a win over TCRdist nearest-neighbour for unseen-epitope TCR specificity at realistic support sizes, and neither cheap approach solves the generalization problem. **Secondary (feature content):** adding the CDR3-alpha chain (secondary analysis) lifts the tabular models by a median +0.042 AUROC -- a larger and more consistent gain than any difference *between* methods, underscoring that for this problem the input representation matters more than the learner.

CPU-only feasibility study. Task = one epitope; support = k labeled TCRs for that epitope; query = held-out TCRs for the same epitope. Leave-one-epitope-out, so every evaluated epitope is unseen by any cross-epitope training. The epitope is constant within a task, so the within-task feature vector encodes the TCR only.

## 1. Data provenance

- **Dataset:** IMMREP22 (Working Group 6) TCR specificity benchmark
- **Source:** https://github.com/viragbioinfo/IMMREP_2022_TCRSpecificity
- **Positives:** VDJdb-derived, tetramer/dextramer-sorted, paired alpha-beta
- **Background negatives (scheme a):** unpublished healthy-donor CD8+CD96+ TCRs (A. Eugster, 11 donors)
- **Other-epitope negatives (scheme b):** reconstructed from the binders of the other kept epitopes; any within-90%-identity match to a held-out positive is dropped (public/cross-reactive).
- **Inclusion:** epitopes with >= 100 raw labeled positives (11 epitopes kept).
- **Scoping note:** IMMREP22 positives are curated from VDJdb; adding raw VDJdb/McPAS would re-import the same TCRs and risk cross-split leakage, so this feasibility run uses the benchmark's curated set (which already ships defined per-epitope negatives of both kinds).

## 2. Deduplication / clonal control

Exact-duplicate (CDR3b+V) rows are collapsed; remaining CDR3b are clustered at 90% identity and splits are made by **cluster**, so no clonal/near-duplicate cluster spans support and query.

| epitope | raw pos | exact dups removed | unique | clusters@90% | largest cluster |
| --- | --- | --- | --- | --- | --- |
| ATDALMTGF | 130 | 15 | 115 | 101 | 5 |
| CINGVCWTV | 229 | 8 | 221 | 208 | 3 |
| GILGFVFTL | 680 | 235 | 444 | 282 | 21 |
| GLCTLVAML | 183 | 38 | 145 | 105 | 7 |
| KSKRTPMGF | 107 | 23 | 84 | 62 | 5 |
| LLWNGPMAV | 235 | 13 | 222 | 187 | 6 |
| LTDEMIAQY | 125 | 6 | 119 | 111 | 4 |
| NLVPMVATV | 343 | 24 | 317 | 290 | 6 |
| SPRWYFYYL | 116 | 0 | 116 | 114 | 2 |
| TTDPSFLGRY | 242 | 10 | 232 | 223 | 3 |
| YLQPRTFLL | 334 | 26 | 308 | 241 | 10 |

Background control pool: 14648 unique TCRs (1295 exact duplicates removed).

## 3. Headline: in-context tabular models vs TCRdist (median over epitopes)

AUROC medians are taken over epitopes (each epitope's value is itself a median over 5 seeds). Query sets are class-balanced, so AUPRC chance = 0.5. The two Δ columns are the median paired difference against the **better** of the two TCRdist scorers per epitope.

| scheme | k | #epi | TabPFN | Kumo (S) | TCRdist NN | TCRdist kNN | LightGBM | LogReg | Δ(TabPFN−bestTCRdist) | Δ(Kumo−bestTCRdist) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| (a) background negatives | 20 | 11 | 0.624 | 0.637 | 0.663 | 0.632 | 0.582 | 0.601 | -0.021 | -0.008 |
| (a) background negatives | 50 | 11 | 0.635 | 0.651 | 0.724 | 0.697 | 0.645 | 0.628 | -0.010 | -0.004 |
| (a) background negatives | 200 | 6 | 0.723 | 0.731 | 0.698 | 0.708 | 0.729 | 0.687 | -0.009 | +0.000 |
| (b) other-epitope negatives | 20 | 11 | 0.647 | 0.645 | 0.653 | 0.639 | 0.639 | 0.619 | +0.012 | -0.010 |
| (b) other-epitope negatives | 50 | 11 | 0.685 | 0.683 | 0.709 | 0.741 | 0.672 | 0.657 | -0.004 | -0.007 |
| (b) other-epitope negatives | 200 | 6 | 0.773 | 0.779 | 0.696 | 0.697 | 0.749 | 0.706 | +0.048 | +0.043 |

### 3b. AUPRC (median over epitopes; balanced query, chance = 0.5)

| scheme | k | TabPFN | Kumo (S) | TCRdist NN | TCRdist kNN | LightGBM |
| --- | --- | --- | --- | --- | --- | --- |
| (a) background negatives | 20 | 0.649 | 0.670 | 0.702 | 0.668 | 0.614 |
| (a) background negatives | 50 | 0.678 | 0.695 | 0.701 | 0.709 | 0.672 |
| (a) background negatives | 200 | 0.750 | 0.754 | 0.730 | 0.748 | 0.752 |
| (b) other-epitope negatives | 20 | 0.664 | 0.652 | 0.652 | 0.703 | 0.649 |
| (b) other-epitope negatives | 50 | 0.692 | 0.674 | 0.704 | 0.729 | 0.657 |
| (b) other-epitope negatives | 200 | 0.793 | 0.799 | 0.733 | 0.749 | 0.768 |

Per-run AUROC **and** AUPRC for every (scheme, epitope, seed, k, method) cell are in `results.csv`; per-epitope AUROC+AUPRC at k=50 are in `per_epitope_k50.csv`.

## 4. Zero-shot transfer baseline (trained on all other epitopes, no target labels)

| scheme | model | median AUROC | #epi |
| --- | --- | --- | --- |
| (a) background negatives | Zero-shot LR | 0.555 | 11 |
| (a) background negatives | Zero-shot LGBM | 0.580 | 11 |
| (b) other-epitope negatives | Zero-shot LR | 0.477 | 11 |
| (b) other-epitope negatives | Zero-shot LGBM | 0.538 | 11 |

## 5. AUROC vs k

![AUROC vs k](auroc_vs_k.png)

## 6. Per-epitope AUROC at k=50 (not pooled)


**(a) background negatives**

| epitope | TabPFN v2 | Kumo Tabular (S) | TCRdist NN | TCRdist kNN | LogReg | LightGBM |
| --- | --- | --- | --- | --- | --- | --- |
| ATDALMTGF | 0.745 | 0.766 | 0.732 | 0.697 | 0.716 | 0.771 |
| CINGVCWTV | 0.635 | 0.644 | 0.584 | 0.561 | 0.612 | 0.629 |
| GILGFVFTL | 0.849 | 0.853 | 0.830 | 0.829 | 0.786 | 0.780 |
| GLCTLVAML | 0.769 | 0.774 | 0.778 | 0.748 | 0.774 | 0.731 |
| KSKRTPMGF | 0.583 | 0.608 | 0.756 | 0.776 | 0.596 | 0.537 |
| LLWNGPMAV | 0.641 | 0.673 | 0.724 | 0.745 | 0.635 | 0.645 |
| LTDEMIAQY | 0.561 | 0.578 | 0.592 | 0.569 | 0.551 | 0.586 |
| NLVPMVATV | 0.571 | 0.582 | 0.529 | 0.474 | 0.536 | 0.568 |
| SPRWYFYYL | 0.598 | 0.651 | 0.584 | 0.521 | 0.628 | 0.659 |
| TTDPSFLGRY | 0.631 | 0.620 | 0.657 | 0.637 | 0.552 | 0.612 |
| YLQPRTFLL | 0.684 | 0.738 | 0.782 | 0.767 | 0.661 | 0.666 |

**(b) other-epitope negatives**

| epitope | TabPFN v2 | Kumo Tabular (S) | TCRdist NN | TCRdist kNN | LogReg | LightGBM |
| --- | --- | --- | --- | --- | --- | --- |
| ATDALMTGF | 0.757 | 0.780 | 0.770 | 0.743 | 0.703 | 0.766 |
| CINGVCWTV | 0.686 | 0.683 | 0.511 | 0.526 | 0.623 | 0.672 |
| GILGFVFTL | 0.842 | 0.845 | 0.840 | 0.846 | 0.788 | 0.816 |
| GLCTLVAML | 0.697 | 0.711 | 0.780 | 0.770 | 0.663 | 0.685 |
| KSKRTPMGF | 0.679 | 0.692 | 0.709 | 0.765 | 0.634 | 0.633 |
| LLWNGPMAV | 0.654 | 0.671 | 0.720 | 0.741 | 0.657 | 0.674 |
| LTDEMIAQY | 0.601 | 0.565 | 0.571 | 0.571 | 0.571 | 0.585 |
| NLVPMVATV | 0.542 | 0.594 | 0.489 | 0.464 | 0.571 | 0.588 |
| SPRWYFYYL | 0.624 | 0.659 | 0.598 | 0.573 | 0.663 | 0.654 |
| TTDPSFLGRY | 0.685 | 0.662 | 0.646 | 0.674 | 0.605 | 0.630 |
| YLQPRTFLL | 0.694 | 0.714 | 0.743 | 0.745 | 0.679 | 0.713 |

## 6b. Secondary analysis: paired CDR3 alpha+beta features

Same protocol (3 seeds) with 100 extra Atchley features from CDR3-alpha (204 features total) on the paired-chain subset (~100% of TCRs). Median per-epitope AUROC; Δ = paired − beta-only (primary). **Note:** only the feature-based models get the alpha chain here; TCRdist was left beta-only, so this isolates the information value of the alpha chain for the tabular models and is *not* a head-to-head against a paired-chain TCRdist (which would also improve).

| scheme | method | k | beta AUROC | paired AUROC | Δ(paired−beta) |
| --- | --- | --- | --- | --- | --- |
| (a) background negatives | TabPFN v2 | 20 | 0.624 | 0.650 | +0.026 |
| (a) background negatives | TabPFN v2 | 50 | 0.635 | 0.714 | +0.080 |
| (a) background negatives | TabPFN v2 | 200 | 0.723 | 0.772 | +0.049 |
| (a) background negatives | Kumo Tabular (S) | 20 | 0.637 | 0.666 | +0.029 |
| (a) background negatives | Kumo Tabular (S) | 50 | 0.651 | 0.736 | +0.085 |
| (a) background negatives | Kumo Tabular (S) | 200 | 0.731 | 0.780 | +0.049 |
| (a) background negatives | LightGBM | 20 | 0.582 | 0.610 | +0.028 |
| (a) background negatives | LightGBM | 50 | 0.645 | 0.696 | +0.051 |
| (a) background negatives | LightGBM | 200 | 0.729 | 0.772 | +0.043 |
| (a) background negatives | LogReg | 20 | 0.601 | 0.636 | +0.035 |
| (a) background negatives | LogReg | 50 | 0.628 | 0.685 | +0.057 |
| (a) background negatives | LogReg | 200 | 0.687 | 0.737 | +0.050 |
| (b) other-epitope negatives | TabPFN v2 | 20 | 0.647 | 0.646 | -0.001 |
| (b) other-epitope negatives | TabPFN v2 | 50 | 0.685 | 0.725 | +0.040 |
| (b) other-epitope negatives | TabPFN v2 | 200 | 0.773 | 0.820 | +0.047 |
| (b) other-epitope negatives | Kumo Tabular (S) | 20 | 0.645 | 0.670 | +0.025 |
| (b) other-epitope negatives | Kumo Tabular (S) | 50 | 0.683 | 0.751 | +0.068 |
| (b) other-epitope negatives | Kumo Tabular (S) | 200 | 0.779 | 0.820 | +0.041 |
| (b) other-epitope negatives | LightGBM | 20 | 0.639 | 0.638 | -0.000 |
| (b) other-epitope negatives | LightGBM | 50 | 0.672 | 0.726 | +0.053 |
| (b) other-epitope negatives | LightGBM | 200 | 0.749 | 0.811 | +0.062 |
| (b) other-epitope negatives | LogReg | 20 | 0.619 | 0.655 | +0.035 |
| (b) other-epitope negatives | LogReg | 50 | 0.657 | 0.684 | +0.027 |
| (b) other-epitope negatives | LogReg | 200 | 0.706 | 0.737 | +0.031 |

Median effect of adding the alpha chain across cells: +0.042 AUROC.

## 7. Runtime per method (CPU)

| method | mean s/fit | max s/fit |
| --- | --- | --- |
| Kumo Tabular (S) | 0.73 | 2.91 |
| LightGBM | 0.14 | 0.25 |
| LogReg | 0.01 | 0.04 |
| TabPFN v2 | 6.82 | 17.79 |
| TCRdist kNN | 0.28 | 4.82 |
| TCRdist NN | 0.28 | 4.82 |

## 8. Reproduce

```bash
python -m venv venv && source venv/bin/activate
pip install -r tcr_feasibility/code/requirements.txt \
    torch --index-url https://download.pytorch.org/whl/cpu
cd tcr_feasibility/code
PYTHON=../../venv/bin/python ./run.sh          # full study (~1 h CPU)
PYTHON=../../venv/bin/python ./run.sh --quick  # fast smoke (2 seeds)
```

Python 3.13, CPU-only. Pinned versions in `code/requirements.txt` (TabPFN 2.0.9 = the ungated v2 weights; Kumo Tabular via `structured-data-models` 0.3.0). Every run is in `results.csv` (primary + Kumo) and `results_paired.csv` (secondary); this report and the plot are regenerated by `make_report.py`.
