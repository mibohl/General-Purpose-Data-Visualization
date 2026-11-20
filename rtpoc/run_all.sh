#!/usr/bin/env bash
# Runs C1-C3 for all datasets, 4 single-threaded workers in parallel.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -c "import sys; sys.path.insert(0,'rtpoc'); from datasets import ALL; print('\n'.join(ALL))" |
  xargs -P 4 -I{} sh -c 'python3 rtpoc/run.py {} --threads 1 > results/logs/{}.log 2>&1'
