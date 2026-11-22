#!/usr/bin/env bash
# Headroom for PXD079751 in two processes (one process ran out of memory), after the epoch-cap check.
cd "$(dirname "$0")/.."
until grep -q "EPOCHCAP DONE" results/logs/ft_after.log 2>/dev/null; do sleep 30; done
A="--lr 0.01 --epochs 150 --patience 20 --threads 4"
python3 rtpoc/ft_headroom.py PXD079751 $A --parts zs ft > results/logs/ft_hr_PXD079751.log 2>&1
python3 rtpoc/ft_headroom.py PXD079751 $A --parts b4 >> results/logs/ft_hr_PXD079751.log 2>&1
echo HRBIG DONE >> results/logs/ft_after.log
