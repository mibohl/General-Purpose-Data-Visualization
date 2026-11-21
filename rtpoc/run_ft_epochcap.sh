#!/usr/bin/env bash
# Runs the epoch-cap sensitivity check once the headroom runs have finished.
cd "$(dirname "$0")/.."
until grep -q "AFTER DONE" results/logs/ft_after.log 2>/dev/null; do sleep 30; done
for n in yeast_120min_psms_aligned prosit_ptm_2020; do
  python3 rtpoc/ft_epochcap.py $n --threads 4 > results/logs/ft_epochcap_$n.log 2>&1
done
echo EPOCHCAP DONE >> results/logs/ft_after.log
