"""Parse MaxQuant evidence tables from IM2Deep's Zenodo MQ_output.zip into one filtered parquet.
Usage: python -I 01_parse.py <MQ_output.zip> <out_dir>
"""
import sys, re, zipfile, hashlib, json
import numpy as np, pandas as pd

zpath, out = sys.argv[1], sys.argv[2]
z = zipfile.ZipFile(zpath)

def study_of(fn):
    b = fn.split('/')[-1]
    if 'PXD019086' in b:
        return 'PXD019086'
    if b.startswith('evidence_') and not re.search(r'\d{6}', b):
        return 'SYN21PTM'          # synthetic 21-PTM peptide library (AnWi), one lab/project
    b2 = b.replace('PXXD', 'PXD')
    m = re.search(r'(\d{5,6})', b2)
    return 'PXD' + m.group(1).zfill(6)

WANT = ['Sequence', 'Modified sequence', 'Modifications', 'Charge', 'm/z', 'Mass', 'Raw file', 'Type',
        'Reverse', 'Potential contaminant', 'PEP', 'Score', 'Intensity', '1/K0', 'Calibrated 1/K0',
        'CCS', 'Calibrated CCS', 'Ion mobility index', 'K0', 'Calibrated K0', 'Retention time',
        'Labeling state', 'Experiment', 'id']
frames, log, seen_hash = [], [], {}
for info in z.infolist():
    if info.file_size == 0:
        continue
    fn = info.filename.split('/')[-1]
    # detect exact duplicate files by content hash
    h = hashlib.md5()
    with z.open(info) as f:
        for chunk in iter(lambda: f.read(1 << 24), b''):
            h.update(chunk)
    h = h.hexdigest()
    if h in seen_hash:
        log.append(dict(file=fn, action='skip_duplicate_of', other=seen_hash[h])); continue
    seen_hash[h] = fn
    with z.open(info) as f:
        hdr = f.readline().decode().rstrip('\r\n').split('\t')
    cols = [c for c in WANT if c in hdr]
    labelcols = [c for c in hdr if c.lower().startswith('reporter intensity') or 'ratio h/l' in c.lower()]
    with z.open(info) as f:
        df = pd.read_csv(f, sep='\t', usecols=cols, low_memory=False)
    n0 = len(df)
    df['study'] = study_of(fn); df['file'] = fn
    df['has_reporter'] = len(labelcols) > 0
    log.append(dict(file=fn, study=df['study'].iat[0], rows=n0, md5=h, n_reporter_or_silac_cols=len(labelcols),
                    has_1k0='1/K0' in cols, has_cal1k0='Calibrated 1/K0' in cols, has_ccs='CCS' in cols,
                    n_raw=df['Raw file'].nunique()))
    frames.append(df)
    print(fn, n0, flush=True)
df = pd.concat(frames, ignore_index=True)
for c in ['Reverse', 'Potential contaminant', 'Modifications', 'Modified sequence', 'Experiment']:
    if c in df: df[c] = df[c].astype('string')
df.to_parquet(f'{out}/evidence_all.parquet')
pd.DataFrame(log).to_csv(f'{out}/parse_log.csv', index=False)
