## Spearman on seed-0 queries

|                       |   B_knn zero-shot |   B_lgbm zero-shot |   DeepLC 4.5 zero-shot |   B4 k=100 |   FT k=100 |   T1_d45 k=100 |   B4 k=1000 |   FT k=1000 |   T1_d45 k=1000 |   best zero-shot |   best at k=1000 |   Spearman gap (k=1000 best - zero-shot best) |
|:----------------------|------------------:|-------------------:|-----------------------:|-----------:|-----------:|---------------:|------------:|------------:|----------------:|-----------------:|-----------------:|----------------------------------------------:|
| ('easy', 'PXD079751') |            0.9989 |             0.9985 |                 0.9983 |     0.9997 |     0.9997 |         0.9997 |      0.9998 |      0.9998 |          0.9997 |           0.9989 |           0.9998 |                                        0.0008 |
| ('easy', 'PXD062983') |            0.9977 |             0.9947 |                 0.9970 |     0.9993 |     0.9993 |         0.9994 |      0.9996 |      0.9995 |          0.9996 |           0.9977 |           0.9996 |                                        0.0019 |
| ('easy', 'PXD076643') |            0.9936 |             0.9929 |                 0.9931 |     0.9948 |     0.9948 |         0.9947 |      0.9947 |      0.9948 |          0.9948 |           0.9936 |           0.9948 |                                        0.0012 |
| ('hard', 'PXD050641') |            0.9784 |             0.9663 |                 0.9790 |     0.9911 |     0.9910 |         0.9936 |      0.9923 |      0.9927 |          0.9940 |           0.9790 |           0.9940 |                                        0.0150 |
| ('hard', 'PXD085246') |          nan      |           nan      |                 0.9890 |     0.9919 |     0.9918 |         0.9919 |      0.9919 |      0.9919 |          0.9921 |           0.9890 |           0.9921 |                                        0.0031 |
| ('hard', 'PXD061731') |            0.9848 |             0.9856 |                 0.9870 |     0.9923 |     0.9924 |         0.9920 |      0.9928 |      0.9927 |          0.9931 |           0.9870 |           0.9931 |                                        0.0061 |
| ('hard', 'PXD074962') |            0.9875 |             0.9823 |                 0.9842 |     0.9932 |     0.9929 |         0.9938 |      0.9910 |      0.9929 |          0.9924 |           0.9875 |           0.9929 |                                        0.0054 |
| ('hard', 'PXD084949') |            0.9884 |             0.9899 |                 0.9895 |     0.9914 |     0.9912 |         0.9922 |      0.9903 |      0.9896 |          0.9927 |           0.9899 |           0.9927 |                                        0.0028 |

Group medians:

|                                             |   easy |   hard |
|:--------------------------------------------|-------:|-------:|
| B_knn zero-shot                             | 0.9977 | 0.9861 |
| B_lgbm zero-shot                            | 0.9947 | 0.9840 |
| DeepLC 4.5 zero-shot                        | 0.9970 | 0.9870 |
| B4 k=100                                    | 0.9993 | 0.9919 |
| FT k=100                                    | 0.9993 | 0.9918 |
| T1_d45 k=100                                | 0.9994 | 0.9922 |
| B4 k=1000                                   | 0.9996 | 0.9919 |
| FT k=1000                                   | 0.9995 | 0.9927 |
| T1_d45 k=1000                               | 0.9996 | 0.9927 |
| best zero-shot                              | 0.9977 | 0.9875 |
| best at k=1000                              | 0.9996 | 0.9929 |
| Spearman gap (k=1000 best - zero-shot best) | 0.0012 | 0.0054 |

## Easy: compute ratio to the cheapest adequate arm (core-seconds; adequate = within 5% of the best error)

| arm                               |   100 |   250 |   500 |   1000 |
|:----------------------------------|------:|------:|------:|-------:|
| B4                                |  0.01 |  0.01 |  0.01 |   0.00 |
| B_knn+ridge (model trained once)  |  0.00 |  0.00 |  0.00 |   0.00 |
| B_knn+ridge (per new setup)       |  0.01 |  0.01 |  0.01 |   0.00 |
| B_lgbm+ridge (model trained once) |  0.00 |  0.00 |  0.00 |   0.00 |
| B_lgbm+ridge (per new setup)      |  1.96 |  1.42 |  1.10 |   0.70 |
| FT                                |  0.93 |  1.94 |  2.57 |   2.94 |
| T1_d45                            |  1.00 |  0.97 |  0.95 |   0.94 |
| T1_d45_B4                         |  1.02 |  1.00 |  1.00 |   1.00 |
| T1_d45_FT                         |  2.20 |  3.66 |  6.27 |   6.98 |

Adequate arms per k: k=100: T1_d45, T1_d45_B4, T1_d45_FT; k=250: T1_d45_B4, T1_d45_FT; k=500: T1_d45_B4, T1_d45_FT; k=1000: T1_d45_B4, T1_d45_FT

## Hard: compute ratio to the cheapest adequate arm (core-seconds; adequate = within 5% of the best error)

| arm                               |   100 |   250 |   500 |   1000 |
|:----------------------------------|------:|------:|------:|-------:|
| B4                                |  0.01 |  0.01 |  0.01 |   0.01 |
| B_knn+ridge (model trained once)  |  0.00 |  0.00 |  0.00 |   0.00 |
| B_knn+ridge (per new setup)       |  0.01 |  0.01 |  0.01 |   0.00 |
| B_lgbm+ridge (model trained once) |  0.00 |  0.00 |  0.00 |   0.00 |
| B_lgbm+ridge (per new setup)      |  1.96 |  1.47 |  1.16 |   0.75 |
| FT                                |  0.93 |  2.01 |  2.70 |   3.14 |
| T1_d45                            |  1.00 |  1.00 |  1.00 |   1.00 |
| T1_d45_B4                         |  1.02 |  1.03 |  1.05 |   1.07 |
| T1_d45_FT                         |  2.20 |  3.79 |  6.58 |   7.45 |

Adequate arms per k: k=100: T1_d45, T1_d45_B4, T1_d45_FT; k=250: T1_d45, T1_d45_B4, T1_d45_FT; k=500: T1_d45, T1_d45_B4, T1_d45_FT; k=1000: T1_d45, T1_d45_B4, T1_d45_FT

## Seconds per fit + 2,000-query prediction at k = 250 (median over setups and seeds)

| arm                               |   seconds |   threads |   core_s | note                                                                       |
|:----------------------------------|----------:|----------:|---------:|:---------------------------------------------------------------------------|
| B4                                |       0.4 |         1 |      0.4 | DeepLC 4.5 calibration fit + predict (head matrix precomputed)             |
| FT                                |      28.0 |         4 |    112.0 | finetune + calibrate + predict, idle machine                               |
| T1_d45                            |      55.7 |         1 |     55.7 | TabPFN fit + predict (4 estimators) + B4 head ranking                      |
| T1_d45_B4                         |      57.6 |         1 |     57.6 | as T1_d45 plus 5-fold cross-fitted B4                                      |
| T1_d45_FT                         |     211.0 |         1 |    211.0 | full fine-tune + 3 cross-fit fine-tunes + TabPFN (grid, 4 jobs on 4 cores) |
| B_knn+ridge (per new setup)       |       0.5 |         1 |      0.5 | zero-shot model 1 s per new setup + anchored ridge                         |
| B_knn+ridge (model trained once)  |       0.0 |         1 |      0.0 | anchored ridge only                                                        |
| B_lgbm+ridge (per new setup)      |      20.4 |         4 |     81.8 | zero-shot model 20 s per new setup + anchored ridge                        |
| B_lgbm+ridge (model trained once) |       0.0 |         1 |      0.0 | anchored ridge only                                                        |

All k:

| arm                               |   100 |   250 |   500 |   1000 |
|:----------------------------------|------:|------:|------:|-------:|
| B4                                |   0.3 |   0.4 |   0.5 |    0.6 |
| B_knn+ridge (model trained once)  |   0.0 |   0.0 |   0.0 |    0.0 |
| B_knn+ridge (per new setup)       |   0.5 |   0.5 |   0.5 |    0.5 |
| B_lgbm+ridge (model trained once) |   0.0 |   0.0 |   0.0 |    0.0 |
| B_lgbm+ridge (per new setup)      |  20.4 |  20.4 |  20.4 |   20.4 |
| FT                                |   9.7 |  28.0 |  47.7 |   85.5 |
| T1_d45                            |  41.6 |  55.7 |  70.6 |  109.0 |
| T1_d45_B4                         |  42.5 |  57.6 |  74.1 |  116.4 |
| T1_d45_FT                         |  91.5 | 211.0 | 464.2 |  812.5 |
