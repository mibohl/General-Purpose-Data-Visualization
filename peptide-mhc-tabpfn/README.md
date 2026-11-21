# TabPFN in-context learning vs pan-allele prediction for rare HLA class I alleles

CPU-only feasibility study. **Question:** on held-out rare HLA class I alleles
(few measurements), does TabPFN v2 with *k* measured support peptides match or
beat (a) MHCflurry 2.x pan-allele predictions and (b) a PSSM built from the same
*k* peptides?

No neural network is trained. TabPFN v2 is used for in-context inference only;
all other learned baselines are ridge / logistic / PSSM / kNN.

## Run

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu   # CPU torch
./run.sh
```

`run.sh` is idempotent and the experiment is resumable (it appends to
`results/results.csv` and skips alleles already present). Outputs:

- `results/report.md` — provenance, verdict, tables, plot.
- `results/results.csv` — every run (allele x k x seed x method).
- `results/results_vs_k.png` — metric vs k.

## Design

- **Data:** TDC `mhc1_iedb-imgt_nielsen` (NetMHCpan-style IEDB, Nielsen), 185,985
  rows / 153 alleles. First pass: HLA, 9-mers (121,720 rows / 103 alleles).
- **Task = one allele.** Allele is fixed per task, so the MHC pseudosequence is
  **not** a feature (keeps feature count small). Peptides encoded as BLOSUM62 row
  per residue -> 180 features for 9-mers.
- **Held out:** 44 rare-tail HLA alleles (50-1500 measurements) with >=10 of each
  class. **Pan-allele pool:** 30 high-data alleles (n>1500); no held-out allele
  is ever in it.
- **Per allele x k in {10,50,200,1000} (capped at n/2) x 5 seeds:** TabPFN
  regressor (1-log50k) + classifier (500 nM binder), PSSM log-odds, kNN-BLOSUM,
  ridge, logistic, MHCflurry 2.x affinity, MHCflurry recalibrated on the k
  support peptides, and a leakage-free pan-allele ridge/logistic reference.
- **Metrics:** AUROC, AUPRC, Spearman on 1-log50k; mean +/- std across alleles.

## Files

```
src/encoders.py        BLOSUM62 + physchem peptide encoders
src/data.py            load/filter/splits/overlap diagnostics
src/baselines.py       PSSM, kNN, ridge, logistic, recalibration, metrics
src/precompute.py      MHCflurry + pan-allele linear predictions (k-independent)
src/run_experiment.py  main sweep, incremental + resumable
src/make_report.py     aggregation, plot, report.md
```
