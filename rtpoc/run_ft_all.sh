#!/usr/bin/env bash
# Fine-tuning study on the 8 held-out setups: phase main (FT, B4 check, T1_d45_B4),
# then phase stack (T1_d45_FT). 4 single-threaded workers. Hyperparameters fixed
# once by ft_tune.py on Part 1 datasets (config D).
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
CFG="--lr 0.01 --epochs 150 --patience 20"
mkdir -p results/logs
LIST="PXD079751 easy
PXD062983 easy
PXD076643 easy
PXD050641 hard
PXD085246 hard
PXD061731 hard
PXD074962 hard
PXD084949 hard"
for PH in main stack; do
  echo "$LIST" | xargs -P 4 -L 1 sh -c "python3 rtpoc/run_ft.py \"\$0\" --label \"\$1\" --phase $PH $CFG > results/logs/ft_${PH}_\$0.log 2>&1"
done
echo ALL DONE
