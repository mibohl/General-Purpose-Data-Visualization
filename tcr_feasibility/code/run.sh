#!/usr/bin/env bash
# One-command runner for the TCR in-context-vs-TCRdist feasibility study.
#
#   PYTHON=/path/to/venv/bin/python ./run.sh           # full study
#   PYTHON=/path/to/venv/bin/python ./run.sh --quick   # 2 seeds, skip paired
#
# Env deps: tabpfn==2.0.9, structured-data-models (Kumo Tabular), tcrdist3,
# scikit-learn, lightgbm, rapidfuzz, pyarrow, matplotlib, torch (CPU).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
PY="${PYTHON:-python}"
RES="$HERE/../results"

QUICK=""
[ "$1" = "--quick" ] && QUICK="--seeds 2"

echo "[1/4] prepare data (IMMREP22 load + dedup + 90% clustering)"
$PY prepare_data.py

echo "[2/4] primary: leave-one-epitope-out, all methods incl. TabPFN v2 + Kumo (S)"
# default --methods = all (tabpfn kumo tcrdist logreg lightgbm zeroshot) -> results.csv
$PY run_experiment.py --fresh $QUICK

echo "[3/4] secondary: paired CDR3 alpha+beta features (feature-based methods, 3 seeds)"
if [ -z "$QUICK" ]; then
  $PY run_experiment.py --feature-mode paired \
      --methods tabpfn kumo lightgbm logreg --seeds 3 \
      --out "$RES/results_paired.csv" --fresh
fi

echo "[4/4] build report + plot"
$PY make_report.py
echo "done -> $RES/report.md , $RES/results.csv , $RES/auroc_vs_k.png"
