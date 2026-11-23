"""Z1: IM2Deep (as shipped, default single-conformer model, no calibration) CCS predictions -> 1/K0.
Usage: python -I 05_im2deep_predict.py <work_dir>
"""
import sys, re
import numpy as np, pandas as pd
from psm_utils import PSM, PSMList
from im2deep import predict
from im2deep.utils import ccs2im
W = sys.argv[1]
sv = pd.read_parquet(f'{W}/study_values.parquet', columns=['precursor', 'modseq', 'charge', 'mz']).drop_duplicates('precursor')
def proforma(modseq, z):
    s = modseq
    nterm = ''
    if s.startswith('[ac]'):
        nterm = '[Acetyl]-'; s = s[4:]
    s = s.replace('[ox]', '[Oxidation]')
    return f'{nterm}{s}/{int(z)}'
sv['proforma'] = [proforma(m, z) for m, z in zip(sv.modseq, sv.charge)]
preds = []
B = 20000
for i in range(0, len(sv), B):
    chunk = sv.iloc[i:i + B]
    pl = PSMList(psm_list=[PSM(peptidoform=p, spectrum_id=str(j)) for j, p in enumerate(chunk.proforma)])
    preds.append(np.asarray(predict(pl)))
    print(i + len(chunk), '/', len(sv), flush=True)
sv['ccs_pred'] = np.concatenate(preds)
sv['k0_pred'] = ccs2im(sv['ccs_pred'].values, sv['mz'].values, sv['charge'].values)
sv[['precursor', 'proforma', 'ccs_pred', 'k0_pred']].to_parquet(f'{W}/im2deep_pred.parquet')
print(sv.head())
