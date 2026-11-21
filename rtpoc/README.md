# In-context RT prediction PoC (TabPFN v2 vs. calibration)

Code for the research plan. Results and conclusions: `results/SUMMARY.md`.
Dataset choice and labels: `data/datasets.md`.

## Reproduce (CPU only, about 5 h on 4 cores)

```bash
pip install deeplc==4.5.0 tabpfn==10.2.0 onnx2torch pandas numpy scikit-learn statsmodels \
    xgboost matplotlib pyarrow tabulate

# 1. Data: train_eval.zip (118 MB) from the DeepLC Zenodo record 5171792
mkdir -p data/raw data/zenodo
curl -L -C - -o data/raw/train_eval.zip \
  "https://zenodo.org/records/5171792/files/train_eval.zip?download=1"
md5sum data/raw/train_eval.zip   # ac1b54af23bdafb1c90c9698eea4ec52
unzip -q data/raw/train_eval.zip -d data/zenodo

# 2. TabPFN v2 regressor weights (ungated). tabpfn downloads them to
#    ~/.cache/tabpfn on first use of ModelVersion.V2; copy them next to the data.
mkdir -p data/models
python -c "from tabpfn import TabPFNRegressor; from tabpfn.constants import ModelVersion; \
import numpy as np; m=TabPFNRegressor.create_default_for_version(ModelVersion.V2, device='cpu'); \
m.fit(np.random.rand(20,3), np.random.rand(20))"
cp ~/.cache/tabpfn/tabpfn-v2-regressor.ckpt data/models/

# 3. Uncalibrated generic DeepLC prediction per peptide -> data/processed/<ds>.parquet
python rtpoc/generic_pred.py

# 4. Checks C1-C3 (C4/C5 are computed inside C2): results/parts/*.csv
mkdir -p results/logs && ./rtpoc/run_all.sh

# 5. results/results.csv, tables.md, curves_*.png
python rtpoc/aggregate.py
```

## Files

| File | Role |
|---|---|
| `datasets.py` | the 16 datasets and their easy / hard label |
| `generic_pred.py` | legacy DeepLC `full_hc_PXD005573_pub` ensemble, uncalibrated, inference only |
| `legacy_deeplc/` | DeepLC 3.1.13 feature extractor (vendored, needed by the legacy models) |
| `methods.py` | features, B1, B2 (LOWESS, CV-chosen span), B2 DeepLC spline, B3 ridge / XGBoost, TabPFN |
| `run.py` | C1, C2 (+ C4 intervals, C5 timings), C3 for one dataset |
| `aggregate.py` | tables and figures |

Every fit, scaler, LOWESS span choice and hyperparameter search uses context
peptides only; query peptides are only used to compute metrics.

## Follow-up runs (held-out DeepLC 4.5 comparison, noise, n_estimators, timing)

```bash
# 6. Held-out PRIDE sets (published after DeepLC 4.5): download the files listed
#    in data/heldout/heldout.md into data/heldout/raw/<PXD>/, then
python3 -I rtpoc/heldout_prepare.py --check          # -> data/heldout/<PXD>.csv
python rtpoc/heldout_generic.py PXD079751 PXD062983 PXD076643 PXD050641 \
    PXD085246 PXD061731 PXD074962 PXD084949           # legacy + DeepLC 4.5 predictions
./rtpoc/run_heldout_all.sh                            # results/parts/H_<PXD>.csv

# 7. T1 under label noise (C3T1) and n_estimators=8 (NE8) on the 16 Zenodo sets
./rtpoc/run_extra.sh                                  # results/parts/{C3T1,NE8}_<ds>.csv

# 8. Speed for 2k to 100k queries (idle machine)
./rtpoc/run_timing.sh                                 # pruned grid -> results/timing.csv

python rtpoc/aggregate.py
```

| File | Role |
|---|---|
| `heldout_prepare.py` | PRIDE identification tables (MaxQuant, FragPipe, DIA-NN, PEAKS, PD) to `seq,modifications,tr` |
| `heldout_generic.py` | legacy DeepLC prediction and the DeepLC 4.5 multi-head prediction matrix (stored as its low-rank factors) |
| `run_heldout.py`, `run_heldout_all.sh` | C2 protocol plus B4 (DeepLC 4.5 `MultiHeadRidgeCalibration` on the context) and T1_d45 |
| `run_extra.sh` | C3T1 and NE8 checks (via `run.py --checks C3T1 NE8`) |
| `timing.py`, `run_timing.sh` | wall-clock and accuracy for large query sets (pruned grid) |
| `predict.py` | command-line predictor for a new setup (default T1_d45) |

### Using the predictor

```bash
python rtpoc/predict.py --context ctx.csv --query query.csv --out pred.csv [--method T1_d45|T1|T2]
```

`ctx.csv` has `seq,modifications,tr` (a few hundred identified peptides from the
new setup, any RT unit); `query.csv` has `seq,modifications`. The output adds
`rt_pred` and a 90% interval `rt_lo90`, `rt_hi90`. Modifications outside the
built-in vocabulary can be added as count features with `--mods TMTpro ...`.
