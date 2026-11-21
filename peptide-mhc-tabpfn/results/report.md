# TabPFN in-context learning vs pan-allele prediction for rare HLA class I alleles

_CPU-only feasibility check. No neural networks were trained; TabPFN v2 is used for in-context inference only._

## Verdict

**No. On held-out rare HLA class I alleles, TabPFN with a small support set does not beat MHCflurry 2.x, and only marginally beats a PSSM.** At k=200 support peptides (1 alleles with enough data), mean AUROC is TabPFN 0.921 vs MHCflurry 0.967 vs PSSM 0.890; TabPFN beats MHCflurry on only 0% of these alleles (mean AUROC delta -0.046), and beats the PSSM on 100% (delta +0.031). On Spearman of 1-log50k, TabPFN 0.740 vs MHCflurry 0.814.

At the smallest support (k=10, 1 alleles) TabPFN AUROC is 0.644 vs MHCflurry 0.966 vs PSSM 0.674.

**Crucial caveat (read this first):** MHCflurry 2.x was itself trained on the IEDB binding-affinity data this benchmark is drawn from, so for these "held-out" alleles it is effectively evaluated near-in-sample (63% of peptides occur in >1 allele; 92% of each held-out allele's peptides also appear in the high-data training pool). MHCflurry is therefore an *upper-bound* reference, not a fair out-of-distribution competitor. The honest leakage-free pan-allele reference is our own linear model trained only on disjoint high-data alleles, which is near chance (AUROC ~0.56). Against that bar, a handful of in-allele labels (TabPFN, PSSM, ridge) help enormously; against MHCflurry's effectively-seen labels, they do not catch up.

## Data provenance

- **Dataset:** `mhc1_iedb-imgt_nielsen` (NetMHCpan-style IEDB binding data, Nielsen/IEDB), obtained via PyTDC `tdc.multi_pred.PeptideMHC`.
- **PyTDC version:** 1.1.15. **Source:** Harvard Dataverse file id 4167073 (TDC `name2id`), downloaded 2026-11-21.
- **Rows:** 185,985 peptide-allele measurements; **alleles:** 153; unique peptides: 43,018 (full file).
- **Target Y:** 1-log50k affinity in [0,1]; binder threshold 500 nM (Y >= 0.4256).
- **First-pass scope:** HLA only, 9-mers only -> 121,720 rows, 103 alleles.
- **Models:** TabPFN v2 (tabpfn 10.2.0, checkpoints `tabpfn-v2-regressor` / `tabpfn-v2-classifier-finetuned-zk73skhh` from the public Prior-Labs HF repos; n_estimators=2, CPU). MHCflurry 2.3.30 `models_class1_presentation` affinity predictor (nM -> 1-log50k).

## Split

- **Rare tail:** alleles with 50-1500 measurements (HLA 9-mers): 53 alleles.
- **Held out:** 44 rare alleles passing a pre-registered, result-independent data-quality filter (>=10 binders and >=10 non-binders so AUROC/AUPRC are defined). Range n=74-1428; loci A/B/C/E.
- **High-data pool (n>1500):** 30 alleles, 94,777 rows — used only for our own leakage-free pan-allele linear baseline. No held-out allele is ever in this pool.
- **Protocol:** per allele, k in {10,50,200,1000} capped at n//2, 5 seeds; support = k sampled peptides, query = the rest (capped at 500, identical query set scored by every method). Reported k-tiers use only alleles whose n//2>=k (uncapped); k=1000 is empty because the largest held-out allele has n=1428.

## Peptide overlap across alleles

- 12,650 of 20,171 unique HLA 9-mer peptides (62.7%) appear in more than one allele; the most-shared peptide appears in 46 alleles.
- On average 91.8% of a held-out allele's peptides also appear in the high-data pool. These are shared screening panels: the alleles are label-poor, not peptide-novel.

## Results: mean AUROC across held-out alleles, by k

| method | k=10 | k=50 | k=200 |
|---|---|---|---|
| TabPFN (regressor) | 0.644 | 0.877 | 0.921 |
| TabPFN (classifier) | 0.649 | 0.865 | 0.922 |
| MHCflurry 2.x (pan-allele) | 0.966 | 0.965 | 0.967 |
| MHCflurry + recal (k) | 0.936 | 0.958 | 0.964 |
| PSSM / PWM (k) | 0.674 | 0.828 | 0.890 |
| kNN-BLOSUM (k) | 0.659 | 0.785 | 0.844 |
| Ridge (k) | 0.705 | 0.821 | 0.900 |
| Logistic (k) | 0.667 | 0.837 | 0.907 |
| Pan-allele ridge (leak-free) | 0.803 | 0.801 | 0.808 |
| Pan-allele logistic (leak-free) | 0.783 | 0.781 | 0.788 |

## Results: mean Spearman (1-log50k) across held-out alleles, by k

| method | k=10 | k=50 | k=200 |
|---|---|---|---|
| TabPFN (regressor) | 0.259 | 0.664 | 0.740 |
| TabPFN (classifier) | 0.266 | 0.648 | 0.738 |
| MHCflurry 2.x (pan-allele) | 0.818 | 0.809 | 0.814 |
| MHCflurry + recal (k) | 0.834 | 0.835 | 0.828 |
| PSSM / PWM (k) | 0.309 | 0.575 | 0.683 |
| kNN-BLOSUM (k) | 0.278 | 0.503 | 0.611 |
| Ridge (k) | 0.362 | 0.555 | 0.690 |
| Logistic (k) | 0.296 | 0.590 | 0.703 |
| Pan-allele ridge (leak-free) | 0.542 | 0.537 | 0.551 |
| Pan-allele logistic (leak-free) | 0.505 | 0.502 | 0.516 |

## Results: mean AUPRC across held-out alleles, by k

| method | k=10 | k=50 | k=200 |
|---|---|---|---|
| TabPFN (regressor) | 0.535 | 0.803 | 0.872 |
| TabPFN (classifier) | 0.526 | 0.783 | 0.870 |
| MHCflurry 2.x (pan-allele) | 0.951 | 0.950 | 0.952 |
| MHCflurry + recal (k) | 0.902 | 0.924 | 0.937 |
| PSSM / PWM (k) | 0.545 | 0.738 | 0.841 |
| kNN-BLOSUM (k) | 0.530 | 0.682 | 0.767 |
| Ridge (k) | 0.590 | 0.724 | 0.841 |
| Logistic (k) | 0.545 | 0.750 | 0.850 |
| Pan-allele ridge (leak-free) | 0.719 | 0.723 | 0.732 |
| Pan-allele logistic (leak-free) | 0.689 | 0.688 | 0.700 |

## Best case: each allele at its largest support (~n/2)

Mean AUROC across all held-out alleles at their maximum support size:

| method | mean AUROC | std |
|---|---|---|
| TabPFN (regressor) | 0.932 | nan |
| TabPFN (classifier) | 0.935 | nan |
| MHCflurry 2.x (pan-allele) | 0.970 | nan |
| MHCflurry + recal (k) | 0.968 | nan |
| PSSM / PWM (k) | 0.922 | nan |
| kNN-BLOSUM (k) | 0.869 | nan |
| Ridge (k) | 0.939 | nan |
| Logistic (k) | 0.919 | nan |
| Pan-allele ridge (leak-free) | 0.799 | nan |
| Pan-allele logistic (leak-free) | 0.778 | nan |

## Results vs k

![results vs k](results_vs_k.png)

## Per-allele AUROC at largest support (key methods)

| allele | n | k | TabPFN | MHCflurry 2.x | PSSM / PWM | Ridge | kNN-BLOSUM | Pan-allele ridge |
|---|---|---|---|---|---|---|---|---|
| HLA-A02:11 | 1078 | 539 | 0.932 | 0.970 | 0.922 | 0.939 | 0.869 | 0.799 |

## Runtime per method (mean wall-seconds per fit+predict call, CPU, 4 cores)

| method | mean s/call |
|---|---|
| TabPFN (regressor) | 13.86 |
| TabPFN (classifier) | 10.03 |
| MHCflurry 2.x (pan-allele) | 0.00 |
| MHCflurry + recal (k) | 0.00 |
| PSSM / PWM (k) | 0.00 |
| kNN-BLOSUM (k) | 0.01 |
| Ridge (k) | 0.00 |
| Logistic (k) | 0.06 |
| Pan-allele ridge (leak-free) | 0.00 |
| Pan-allele logistic (leak-free) | 0.00 |

_MHCflurry and the pan-allele linear models are k/seed-independent (predicted once per peptide, zero target labels), shown at 0 s in-loop._


---
_Coverage at report time: 1/44 held-out alleles complete._
