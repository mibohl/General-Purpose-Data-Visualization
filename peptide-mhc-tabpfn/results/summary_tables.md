## AUROC by k

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

## AUPRC by k

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

## Spearman by k

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
