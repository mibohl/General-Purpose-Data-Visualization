#!/usr/bin/env bash
# Interval rerun for all 8 held-out setups, 4 jobs x 1 thread.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
printf "%s\n" "PXD079751 easy" "PXD062983 easy" "PXD076643 easy" "PXD050641 hard" \
  "PXD085246 hard" "PXD061731 hard" "PXD074962 hard" "PXD084949 hard" |
  xargs -P 4 -L 1 sh -c 'python3 rtpoc/interval_rerun.py $0 --label $1 > results/logs/iv_$0.log 2>&1'
echo IV DONE >> results/logs/iv_all.log
