"""Predict ROIs with the released 2-step checkpoint for every cleaned RP dataset. Runs in the py3.12 chemprop env."""
import sys, os, io, time, contextlib, logging, warnings
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
sys.path.insert(0, '/home/user/data/2-step')
os.chdir('/home/user/data/2-step')
import pandas as pd, numpy as np, torch
torch.set_num_threads(4)
from predict import load_model, load_metadata
from predict_rois import predict_rois

R = '/home/user/data/RepoRT/processed_data'
OUT = '/home/user/work/out'
TMP = '/home/user/work/tmp_roi'; os.makedirs(TMP, exist_ok=True)
comp = pd.read_parquet(f'{OUT}/compounds.parquet')
meta = pd.read_csv(f'{OUT}/datasets_meta.csv', dtype={'id': str}).set_index('id')
rp = [d for d in meta.index if meta.loc[d, 'mode'] == 'RP']
model = load_model('models/2-step0525.pt', all_in_one=True)
model.eval()
res, timing = [], []
for i, d in enumerate(rp):
    sub = comp[comp.dataset == d]
    if len(sub) < 5:
        continue
    fn = f'{TMP}/{d}.tsv'
    sub[['smiles', 'rt']].rename(columns={'smiles': 'smiles.std'}).to_csv(fn, sep='\t', index=False)
    md = load_metadata(f'{R}/{d}/{d}_metadata.yaml')
    if bool(meta.loc[d, 'pH_imputed']):
        md['eluent.A.pH'] = float(meta.loc[d, 'pH_2step_input'])
        md['eluent.B.pH'] = float(meta.loc[d, 'pH_2step_input'])
    t = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = predict_rois(model, fn, md, '/home/user/data/RepoRT', 512, False)
    except Exception as e:
        print(d, 'FAILED', repr(e)[:200], flush=True)
        continue
    dt = time.time() - t
    m = out[['smiles', 'roi']].drop_duplicates('smiles')
    s = sub.merge(m, on='smiles', how='left')
    res.append(s[['dataset', 'cid', 'roi']])
    timing.append(dict(dataset=d, n=len(sub), seconds=dt))
    if i % 20 == 0:
        print(i, d, len(sub), f'{dt:.1f}s', 'nan', s.roi.isna().sum(), flush=True)
pd.concat(res).to_parquet(f'{OUT}/roi_2step.parquet')
pd.DataFrame(timing).to_csv(f'{OUT}/timing_roi_2step.csv', index=False)
print('done', len(res))
