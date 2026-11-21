#!/usr/bin/env bash
# Timing grid (pruned to fit the compute budget), run after all other jobs.
cd "$(dirname "$0")/.."
while pgrep -f '^python3 rtpoc/run' >/dev/null || pgrep -f run_extra.sh >/dev/null; do sleep 30; done
python3 rtpoc/timing.py --threads 4 --ns 250 1000 --n_est 4 --modes fit_preprocessors --queries 2000 20000 100000
python3 rtpoc/timing.py --threads 4 --ns 250 1000 --n_est 8 --modes fit_preprocessors --queries 20000
python3 rtpoc/timing.py --threads 4 --ns 1000 --n_est 4 --modes fit_with_cache --queries 20000
echo TIMING DONE
