import pandas as pd, numpy as np, json, os, glob
R='/home/user/data/RepoRT/processed_data'
train2s=set(json.load(open('/home/user/work/2step_train_sets.json')))
hsm=pd.read_csv('/home/user/data/RepoRT/resources/hsm_database/hsm_database.tsv',sep='\t')
tan=pd.read_csv('/home/user/data/RepoRT/resources/tanaka_database/tanaka_database.tsv',sep='\t')
rows=[]
for d in sorted(os.listdir(R)):
    p=f'{R}/{d}/{d}'
    if not os.path.exists(p+'_rtdata_canonical_success.tsv'): continue
    info=pd.read_csv(p+'_info.tsv',sep='\t',dtype=str).iloc[0]
    md=pd.read_csv(p+'_metadata.tsv',sep='\t').iloc[0]
    rt=pd.read_csv(p+'_rtdata_canonical_success.tsv',sep='\t')
    try: gr=pd.read_csv(p+'_gradient.tsv',sep='\t')
    except Exception: gr=pd.DataFrame()
    glen=float(gr['t [min]'].max()) if len(gr) and 't [min]' in gr else np.nan
    def has(f):
        v=md.get(f,np.nan); return pd.notna(v) and v!=0 and v!=''
    solv=any(has(f'eluent.A.{s}') for s in ['h2o','meoh','acn','iproh']) and any(has(f'eluent.B.{s}') for s in ['h2o','meoh','acn','iproh','acetone'])
    fields={'column.name':has('column.name'),'column.length':has('column.length'),'column.id':has('column.id'),
            'column.particle.size':has('column.particle.size'),'column.temperature':has('column.temperature'),
            'column.flowrate':has('column.flowrate'),'column.t0':has('column.t0'),'eluent.composition':solv,
            'eluent.pH':has('eluent.A.pH') or has('eluent.B.pH'),'gradient':len(gr)>1}
    cname=md.get('column.name',None)
    rows.append(dict(id=d,name=info.get('name'),mode=info.get('method.type'),n_entries=len(rt),
        n_compounds=rt['inchikey.std'].str[:14].nunique(),column=cname,usp=md.get('column.usp.code'),
        gradient_min=glen,t0=md.get('column.t0'),t0_present=has('column.t0'),
        n_meta_fields=sum(fields.values()),meta_fields_total=len(fields),
        hsm=cname in set(hsm['name_normalized']) if 'name_normalized' in hsm else cname in set(hsm.iloc[:,0]),
        tanaka=cname in set(tan['name_normalized']) if 'name_normalized' in tan else cname in set(tan.iloc[:,0]),
        in_2step_train=d in train2s, rt_max=rt.rt.max(), **{'f_'+k:v for k,v in fields.items()}))
df=pd.DataFrame(rows)
df.to_csv('/home/user/work/out/survey_all.csv',index=False)
print(df['mode'].value_counts(dropna=False))
print(pd.crosstab(df['mode'],df.in_2step_train))
print(df[['n_compounds','gradient_min','n_meta_fields']].describe())
print(df.t0_present.value_counts())
print(hsm.columns.tolist()[:8], tan.columns.tolist()[:8])
