"""Wall-clock of one DeepLC 4.5 fine-tune on 4 CPU threads (idle machine).

Seed 0 split, k in {100, 250, 500, 1000}, every held-out setup, fixed config.
Output: results/finetune/timing_4threads.csv
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ft_d45 as FT, run as R  # noqa: E402,E401
ROOT = os.path.dirname(HERE)
CFG = dict(lr=0.01, epochs=150, patience=20)
NAMES = ["PXD079751", "PXD062983", "PXD076643", "PXD050641", "PXD085246", "PXD061731",
         "PXD074962", "PXD084949"]
rows = []
for name in NAMES:
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{name}.parquet"))
    tr = df.tr.values; rng = tr.max() - tr.min()
    q, pool = R.split(len(df), 0)
    for k in (100, 250, 500, 1000):
        pred, raw, t_ft, t_all, ep = FT.finetune_predict(df, pool[:k], q, 0, CFG, threads=4)
        r = dict(dataset=name, n=k, threads=4, ft_train_s=round(t_ft, 1), total_s=round(t_all, 1),
                 epochs_run=ep, rel_mae=np.abs(pred - tr[q]).mean() / rng)
        rows.append(r); print(r, flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(ROOT, "results", "finetune", "timing_4threads.csv"),
                                  index=False)
