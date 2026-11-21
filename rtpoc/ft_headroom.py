"""Headroom per held-out setup: DeepLC 4.5 with no target labels vs fine-tuned
on the largest available context.

- zero_shot: deeplc.predict with the default model and no reference (DeepLC
  reports its PXD005573_mcp setup, a 200-min gradient), on the seed's queries.
- B4_full: DeepLC 4.5 calibration fit on the whole pool (all non-query peptides).
- FT_full: fine-tuned on the whole pool (85/15 split), then calibrated on it.
  Epochs are capped so one fine-tune stays under about 9 minutes on 4 threads;
  the cap and the measured time are written to the output.
Seed 0 only for FT_full (cost); zero-shot and B4_full for seed 0 as well.
Output: results/finetune/parts/HR_<name>.csv
Usage: python rtpoc/ft_headroom.py <name> --lr .. --epochs .. --patience .. [--threads 4]
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import deeplc  # noqa: E402
import ft_d45 as FT  # noqa: E402
import run as R  # noqa: E402
from heldout_generic import load_d45  # noqa: E402
from deeplc.calibration import MultiHeadRidgeCalibration  # noqa: E402

SEC_PER_PEPTIDE_EPOCH = 0.0008   # measured: k=1000, 4 threads, about 0.8 s per epoch
BUDGET_S = 540

ap = argparse.ArgumentParser(); ap.add_argument("name")
ap.add_argument("--lr", type=float, required=True); ap.add_argument("--epochs", type=int, required=True)
ap.add_argument("--patience", type=int, required=True); ap.add_argument("--threads", type=int, default=4)
a = ap.parse_args()
import torch; torch.set_num_threads(a.threads)  # noqa: E401,E702
ROOT = os.path.dirname(HERE)
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{a.name}.parquet"))
tr = df.tr.values.astype(float); rng = tr.max() - tr.min()
seed = 0
q, pool = R.split(len(df), seed)
rows = []
t = time.time()
zs = np.asarray(deeplc.predict(FT.to_psms(df, q, False)), dtype=float)
rows.append(dict(dataset=a.name, method="zero_shot", n=0, rel_mae=np.abs(zs - tr[q]).mean() / rng,
                 pearson=np.corrcoef(zs, tr[q])[0, 1], seconds=time.time() - t))
D = load_d45(a.name)
t = time.time()
yc = tr[pool]
cal = MultiHeadRidgeCalibration(); cal.fit(target=yc, source=D[pool])
pb = np.asarray(cal.transform(D[q]), dtype=float)
rows.append(dict(dataset=a.name, method="B4_full", n=len(pool), rel_mae=np.abs(pb - tr[q]).mean() / rng,
                 pearson=np.corrcoef(pb, tr[q])[0, 1], seconds=time.time() - t))
cap = int(min(a.epochs, max(3, BUDGET_S // (SEC_PER_PEPTIDE_EPOCH * len(pool)))))
cfg = dict(lr=a.lr, epochs=cap, patience=a.patience)
pred, raw, t_ft, t_all, ep = FT.finetune_predict(df, pool, q, seed, cfg, a.threads)
rows.append(dict(dataset=a.name, method="FT_full", n=len(pool), rel_mae=np.abs(pred - tr[q]).mean() / rng,
                 pearson=np.corrcoef(pred, tr[q])[0, 1], seconds=t_all, ft_train_s=t_ft,
                 epochs_cap=cap, epochs_cfg=a.epochs, epochs_run=ep))
out = pd.DataFrame(rows); out["threads"] = a.threads
os.makedirs(os.path.join(ROOT, "results", "finetune", "parts"), exist_ok=True)
out.to_csv(os.path.join(ROOT, "results", "finetune", "parts", f"HR_{a.name}.csv"), index=False)
print(out.to_string(), flush=True)
