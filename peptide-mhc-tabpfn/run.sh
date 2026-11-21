#!/usr/bin/env bash
# One-command runner for the TabPFN-vs-pan-allele feasibility study (CPU only).
# Idempotent: steps whose outputs already exist are skipped.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data cache results

REG_CKPT="cache/tabpfn-v2-regressor.ckpt"
CLF_CKPT="cache/tabpfn-v2-classifier-finetuned-zk73skhh.ckpt"
DATA="data/mhc1_iedb-imgt_nielsen.tab"

echo "== 1/6 install deps =="
pip install -q -r requirements.txt

echo "== 2/6 dataset (TDC mhc1_iedb-imgt_nielsen, Harvard Dataverse 4167073) =="
if [ ! -f "$DATA" ]; then
  curl -sSL -o "$DATA" "https://dataverse.harvard.edu/api/access/datafile/4167073"
fi
echo "   rows: $(($(wc -l < "$DATA") - 1))"

echo "== 3/6 TabPFN v2 checkpoints (public Prior-Labs HF repos) =="
[ -f "$REG_CKPT" ] || curl -sSL -o "$REG_CKPT" \
  "https://huggingface.co/Prior-Labs/TabPFN-v2-reg/resolve/main/tabpfn-v2-regressor.ckpt"
[ -f "$CLF_CKPT" ] || curl -sSL -o "$CLF_CKPT" \
  "https://huggingface.co/Prior-Labs/TabPFN-v2-clf/resolve/main/tabpfn-v2-classifier-finetuned-zk73skhh.ckpt"

echo "== 4/6 MHCflurry models =="
python3 -c "from mhcflurry import Class1AffinityPredictor as P; P.load()" 2>/dev/null \
  || mhcflurry-downloads fetch models_class1_presentation

echo "== 5/6 precompute k-independent predictions (MHCflurry + pan-allele linear) =="
[ -f cache/precomputed.csv ] || python3 src/precompute.py

echo "== 6/6 run experiment (resumable) + build report =="
python3 src/run_experiment.py
python3 src/make_report.py

echo "Done. See results/report.md, results/results.csv, results/results_vs_k.png"
