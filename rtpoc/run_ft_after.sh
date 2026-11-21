#!/usr/bin/env bash
# After the grid: 4-thread timing pass, then headroom runs (one at a time, 4 threads).
cd "$(dirname "$0")/.."
while pgrep -f run_ft_all.sh >/dev/null; do sleep 30; done
python3 rtpoc/ft_timing.py > results/logs/ft_timing.log 2>&1
for n in PXD079751 PXD062983 PXD076643 PXD050641 PXD085246 PXD061731 PXD074962 PXD084949; do
  python3 rtpoc/ft_headroom.py $n --lr 0.01 --epochs 150 --patience 20 --threads 4 > results/logs/ft_hr_$n.log 2>&1
done
echo AFTER DONE
