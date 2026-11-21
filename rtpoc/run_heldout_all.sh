#!/usr/bin/env bash
# Held-out DeepLC 4.5 comparison: 8 datasets, 4 single-threaded workers.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p results/logs
cat <<LIST | xargs -P 4 -L 1 sh -c 'python3 rtpoc/run_heldout.py "$0" --label "$1" --mods TMTpro Cysteinyl > results/logs/H_$0.log 2>&1'
PXD079751 easy
PXD062983 easy
PXD050641 hard
PXD076643 easy
PXD085246 hard
PXD061731 hard
PXD074962 hard
PXD084949 hard
LIST
