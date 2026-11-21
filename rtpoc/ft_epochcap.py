"""Sensitivity check for the epoch cap, on Part 1 datasets only (not held-out).

Config D (lr 0.01, patience 20) hit its 150-epoch cap in about a third of the
held-out fine-tunes. This reruns D with the cap at 150 and at 500 on the same
Part 1 datasets used to choose D, k = 500 and 1000, seeds 0-1, and reports
whether the longer schedule lowers query error. Not used for any headline number.
Output: results/finetune/epochcap_<dataset>.csv
Usage: python rtpoc/ft_epochcap.py <dataset> [--threads 4]
"""
import argparse, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import run as R, ft_d45 as F  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("name"); ap.add_argument("--threads", type=int, default=4)
a = ap.parse_args()
ROOT = os.path.dirname(HERE)
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", a.name + ".parquet")).reset_index(drop=True)
tr = df.tr.values.astype(float); rng = tr.max() - tr.min()
rows = []
for seed in (0, 1):
    q, pool = R.split(len(df), seed)
    for k in (500, 1000):
        for cap in (150, 500):
            cfg = dict(lr=0.01, epochs=cap, patience=20)
            pred, raw, tft, tall, ep = F.finetune_predict(df, pool[:k], q, seed, cfg, a.threads)
            r = dict(dataset=a.name, seed=seed, k=k, epochs_cap=cap, rel_mae=np.abs(pred - tr[q]).mean() / rng,
                     rel_mae_raw=np.abs(raw - tr[q]).mean() / rng, ft_s=tft,
                     epochs_run=ep - 1)  # ft_d45 counts DeepLC's post-training check as an epoch
            rows.append(r); print(r, flush=True)
            pd.DataFrame(rows).to_csv(os.path.join(ROOT, "results", "finetune", f"epochcap_{a.name}.csv"),
                                      index=False)
