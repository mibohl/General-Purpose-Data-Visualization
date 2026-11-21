#!/usr/bin/env bash
# Follow-up checks: T1 under label noise (C3T1) and TabPFN n_estimators=8 (NE8).
# Single-threaded workers; a new one starts whenever fewer than MAXJ rtpoc
# jobs are running (so it can share the machine with other runs).
# Datasets whose NE8 part already exists are skipped.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MAXJ=${MAXJ:-4}
mkdir -p results/logs
for ds in $(python3 -c "import sys; sys.path.insert(0,'rtpoc'); from datasets import ALL; print('\n'.join(ALL))"); do
  [ -f results/parts/NE8_$ds.csv ] && continue
  while [ "$(pgrep -fc '^python3 rtpoc/run')" -ge "$MAXJ" ]; do sleep 20; done
  python3 rtpoc/run.py "$ds" --checks C3T1 NE8 --threads 1 > results/logs/extra_$ds.log 2>&1 &
  sleep 5
done
wait
