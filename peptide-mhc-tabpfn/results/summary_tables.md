## AUROC by k

| method | k=10 | k=50 | k=200 |
|---|---|---|---|
| TabPFN (regressor) | 0.629 | 0.788 | 0.880 |
| TabPFN (classifier) | 0.636 | 0.791 | 0.905 |
| Kumo Tabular (regressor) | 0.645 | 0.822 | 0.914 |
| Kumo Tabular (classifier) | 0.631 | 0.790 | 0.911 |
| MHCflurry 2.x (pan-allele) | 0.902 | 0.922 | 0.954 |
| MHCflurry + recal (k) | 0.863 | 0.909 | 0.943 |
| PSSM / PWM (k) | 0.708 | 0.804 | 0.892 |
| kNN-BLOSUM (k) | 0.641 | 0.741 | 0.804 |
| Ridge (k) | 0.658 | 0.775 | 0.858 |
| Logistic (k) | 0.648 | 0.784 | 0.887 |
| Pan-allele ridge (leak-free) | 0.559 | 0.566 | 0.612 |
| Pan-allele logistic (leak-free) | 0.553 | 0.560 | 0.600 |

## AUPRC by k

| method | k=10 | k=50 | k=200 |
|---|---|---|---|
| TabPFN (regressor) | 0.442 | 0.625 | 0.704 |
| TabPFN (classifier) | 0.450 | 0.612 | 0.699 |
| Kumo Tabular (regressor) | 0.453 | 0.649 | 0.757 |
| Kumo Tabular (classifier) | 0.451 | 0.616 | 0.716 |
| MHCflurry 2.x (pan-allele) | 0.820 | 0.839 | 0.853 |
| MHCflurry + recal (k) | 0.727 | 0.794 | 0.814 |
| PSSM / PWM (k) | 0.545 | 0.606 | 0.671 |
| kNN-BLOSUM (k) | 0.450 | 0.535 | 0.541 |
| Ridge (k) | 0.465 | 0.574 | 0.625 |
| Logistic (k) | 0.452 | 0.591 | 0.664 |
| Pan-allele ridge (leak-free) | 0.381 | 0.373 | 0.312 |
| Pan-allele logistic (leak-free) | 0.376 | 0.368 | 0.303 |

## Spearman by k

| method | k=10 | k=50 | k=200 |
|---|---|---|---|
| TabPFN (regressor) | 0.211 | 0.482 | 0.606 |
| TabPFN (classifier) | 0.220 | 0.433 | 0.567 |
| Kumo Tabular (regressor) | 0.235 | 0.506 | 0.628 |
| Kumo Tabular (classifier) | 0.214 | 0.427 | 0.573 |
| MHCflurry 2.x (pan-allele) | 0.616 | 0.635 | 0.640 |
| MHCflurry + recal (k) | 0.637 | 0.666 | 0.679 |
| PSSM / PWM (k) | 0.327 | 0.447 | 0.552 |
| kNN-BLOSUM (k) | 0.222 | 0.367 | 0.448 |
| Ridge (k) | 0.249 | 0.423 | 0.523 |
| Logistic (k) | 0.238 | 0.416 | 0.535 |
| Pan-allele ridge (leak-free) | 0.097 | 0.109 | 0.174 |
| Pan-allele logistic (leak-free) | 0.085 | 0.098 | 0.160 |
