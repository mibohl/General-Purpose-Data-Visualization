#!/usr/bin/env bash
# Follow-up checks: T1 under label noise (C3T1) and TabPFN n_estimators=8 (NE8).
# 3 single-threaded workers, leaving one core free.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p results/logs
python3 -c "import sys; sys.path.insert(0,'rtpoc'); from datasets import ALL; print('\n'.join(ALL))" |
  xargs -P 3 -I{} sh -c 'python3 rtpoc/run.py {} --checks C3T1 NE8 --threads 1 > results/logs/extra_{}.log 2>&1'
