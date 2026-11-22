#!/usr/bin/env bash
# Rerun Competitor B on an idle machine (timings, saved zero-shot predictions) after the interval rerun.
cd "$(dirname "$0")/.."
until grep -q "IV DONE" results/logs/iv_all.log 2>/dev/null; do sleep 60; done
python3 rtpoc/meta_arms.py > results/logs/meta_arms_rerun.log 2>&1
echo META DONE >> results/logs/iv_all.log
