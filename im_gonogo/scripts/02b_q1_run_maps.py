"""Question 1: per run, is reported 1/K0 a deterministic function of the raw TIMS scan index? Also MaxQuant's
'Calibrated 1/K0' (cross-run alignment) shift relative to '1/K0'.
Usage: python -I 02b_q1_run_maps.py <work_dir>   (reads evidence_all.parquet, writes q1_run_index_map.csv)
"""
import sys
import numpy as np, pandas as pd
W = sys.argv[1]
df = pd.read_parquet(f'{W}/evidence_all.parquet', columns=['study', 'Raw file', 'Type', 'Ion mobility index', '1/K0', 'Calibrated 1/K0'])
df = df[df['1/K0'].notna() & (df['Type'] == 'TIMS-MULTI-MSMS')]
rows = []
for (s, r), g in df.groupby(['study', 'Raw file']):
    if len(g) < 200: continue
    x = g['Ion mobility index'].values.astype(float); y = g['1/K0'].values
    c1 = np.polyfit(x, y, 1); c3 = np.polyfit(x, y, 3); res = y - np.polyval(c3, x)
    cal = g['Calibrated 1/K0'].values - y
    rows.append(dict(study=s, run=r, n=len(g), slope=c1[0], k0_at_idx500=np.polyval(c1, 500), idx_min=x.min(), idx_max=x.max(),
                     lin_mad=np.median(np.abs(y - np.polyval(c1, x))), cubic_mad=np.median(np.abs(res)), cubic_p99=np.quantile(np.abs(res), .99),
                     mq_align_shift_med=np.median(cal), mq_align_shift_mad=np.median(np.abs(cal - np.median(cal)))))
pd.DataFrame(rows).to_csv(f'{W}/q1_run_index_map.csv', index=False)
