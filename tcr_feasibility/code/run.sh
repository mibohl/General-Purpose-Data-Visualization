#!/usr/bin/env bash
# One-command runner for the TCR TabPFN-vs-TCRdist feasibility study.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
PY="${PYTHON:-python}"
echo "[1/3] preparing data (IMMREP22 load + dedup + clustering)"
$PY prepare_data.py
echo "[2/3] running leave-one-epitope-out protocol (this is the long step)"
$PY run_experiment.py "$@"
echo "[3/3] building report + plot"
$PY make_report.py
echo "done. See ../results/report.md and ../results/results.csv"
