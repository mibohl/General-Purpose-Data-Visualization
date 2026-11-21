"""Fix DeepLC fine-tuning hyperparameters ONCE on Part 1 datasets (not held-out).

Grid of 4 configs (A = DeepLC 4.5's shipped defaults for its low-rank task head:
lr 0.05, 50 epochs, patience 10). Contexts k = 100 and 1000, seeds 0-1, split
and 15% validation exactly as in the main run. Score = geometric-mean relative
MAE on these Part 1 datasets' query peptides (none are held-out data).
Output: results/finetune/tune_<dataset>.csv
Usage: python rtpoc/ft_tune.py <dataset> [--threads 2]
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import run as R, ft_d45 as F  # noqa: E402

CONFIGS = {"A": dict(lr=0.05, epochs=50, patience=10),
           "B": dict(lr=0.05, epochs=150, patience=20),
           "C": dict(lr=0.2, epochs=100, patience=15),
           "D": dict(lr=0.01, epochs=150, patience=20)}

ap = argparse.ArgumentParser(); ap.add_argument("name"); ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args()
ROOT = os.path.dirname(HERE); out = os.path.join(ROOT, "results", "finetune"); os.makedirs(out, exist_ok=True)
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", a.name + ".parquet"))
df = df.reset_index(drop=True)
tr = df.tr.values.astype(float); rng = tr.max() - tr.min()
rows = []
for seed in (0, 1):
    q, pool = R.split(len(df), seed)
    for k in (100, 1000):
        for c, cfg in CONFIGS.items():
            pred, raw, tft, tall, ep = F.finetune_predict(df, pool[:k], q, seed, cfg, a.threads)
            r = dict(dataset=a.name, seed=seed, k=k, config=c, **cfg, rel_mae=np.abs(pred - tr[q]).mean() / rng,
                     rel_mae_raw=np.abs(raw - tr[q]).mean() / rng, ft_s=tft, epochs_run=ep)
            rows.append(r); print(r, flush=True)
            pd.DataFrame(rows).to_csv(os.path.join(out, f"tune_{a.name}.csv"), index=False)
