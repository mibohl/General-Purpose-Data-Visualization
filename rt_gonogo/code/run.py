"""Main experiment runner. One 'unit' = (group, dataset, roi_source). Results per unit go to out/runs/*.parquet."""
import os, sys, time, json, zlib, argparse, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
OUT = '/home/user/work/out'
RUNS = f'{OUT}/runs'; os.makedirs(RUNS, exist_ok=True)
KS = [10, 15, 20, 30, 50]
SEEDS = range(5)
SCHEMES = ['random', 'spread']
EPS_FLIP = 10 / 60  # 10 s, as in 2-step
DESC = ['MolLogP', 'TPSA', 'MolWt', 'NumHDonors', 'NumHAcceptors', 'NumRotatableBonds', 'NumAromaticRings',
        'FractionCSP3', 'HeavyAtomCount', 'FormalCharge', 'NumHeteroatoms', 'MolMR']
META = ['col_length', 'col_id', 'col_particle.size', 'col_temperature', 'col_flowrate', 't0_eff', 'pH_2step_input',
        'B_acn', 'B_meoh', 'B_iproh', 'A_formic', 'A_acetic', 'A_nh4ac', 'A_nh4form', 'A_trifluoroacetic',
        'A_nh4bicarb', 'A_nh4oh', 'hsm_H', 'hsm_S*', 'hsm_A', 'hsm_B', 'hsm_C (pH 2.8)', 'hsm_C (pH 7.0)',
        'tan_kPB', 'tan_αCH2', 'tan_αT/O', 'tan_αC/P', 'tan_αB/P', 'tan_αB/P.1', 'grad_len', 'grad_org_start',
        'grad_org_max', 'grad_hold', 'grad_t50', 'grad_t90', 'grad_slope']


def load():
    comp = pd.read_parquet(f'{OUT}/compounds.parquet')
    meta = pd.read_csv(f'{OUT}/datasets_meta.csv', dtype={'id': str}).set_index('id')
    D = pd.read_parquet(f'{OUT}/descriptors.parquet')
    comp = comp.join(D[DESC], on='smiles')
    r2 = pd.read_parquet(f'{OUT}/roi_2step.parquet').rename(columns={'roi': 'roi_2step'})
    rp = pd.read_parquet(f'{OUT}/roi_proxy.parquet').rename(columns={'roi': 'roi_proxy'})
    comp = comp.merge(r2, on=['dataset', 'cid'], how='left').merge(rp, on=['dataset', 'cid'], how='left')
    train2s = set(json.load(open('/home/user/work/2step_train_sets.json')))
    return comp, meta, train2s


def units(meta):
    U = []
    for d in meta.index[meta.g_rp2step]: U.append(('RP-2step', d, '2step'))
    for d in meta.index[meta.g_rpproxy]: U.append(('RP-proxy', d, 'proxy'))
    for d in meta.index[meta.g_pub]:
        U.append(('PUB-2step', d, '2step')); U.append(('PUB-proxy', d, 'proxy'))
    for d in meta.index[meta.g_hilic]: U.append(('HILIC', d, 'proxy'))
    return U


def pools(comp, meta, train2s):
    """Training pools (datasets) per (mode, roi_source): used for Z1 training and for purging."""
    held = set(meta.index[meta.heldout_any])
    ok2 = comp.groupby('dataset').roi_2step.apply(lambda s: s.notna().all())
    rp = set(meta.index[meta['mode'] == 'RP']); hi = set(meta.index[meta['mode'] == 'HILIC'])
    z1 = {('RP', '2step'): sorted(d for d in rp - held if ok2.get(d, False)),
          ('RP', 'proxy'): sorted(rp - held),
          ('HILIC', 'proxy'): sorted(hi - held)}
    purge = {('RP', '2step'): set(z1[('RP', '2step')]) | train2s,   # checkpoint training sets + Z1 pool
             ('RP', 'proxy'): set(z1[('RP', 'proxy')]),
             ('HILIC', 'proxy'): set(z1[('HILIC', 'proxy')])}
    return z1, purge


def train_z1(comp, meta, z1):
    import lightgbm as lgb
    models, tt = {}, {}
    for key, ds in z1.items():
        src = key[1]
        sub = comp[comp.dataset.isin(ds)].dropna(subset=[f'roi_{src}'])
        X = np.column_stack([sub[f'roi_{src}'].values, sub[DESC].values, meta.loc[sub.dataset, META].values])
        t = time.time()
        models[key] = lgb.LGBMRegressor(n_estimators=500, learning_rate=0.05, num_leaves=63, verbose=-1,
                                        n_jobs=4).fit(X, sub.rt.values)
        tt[key] = time.time() - t
        print('Z1', key, 'rows', len(sub), 'datasets', sub.dataset.nunique(), f'{tt[key]:.1f}s', flush=True)
    return models, tt


def anchor_sets(rt, d, scheme, k, seed):
    h = zlib.crc32(f'{d}|{scheme}|{k}|{seed}'.encode())
    rng = np.random.default_rng(h)
    n = len(rt)
    if scheme == 'random':
        return np.sort(rng.choice(n, k, replace=False))
    lo, hi = rt.min(), rt.max()
    targets = lo + (np.arange(k) + rng.uniform(size=k)) / k * (hi - lo)
    used = np.zeros(n, bool); idx = []
    for t in targets:
        cand = np.where(~used)[0]
        j = cand[np.argmin(np.abs(rt[cand] - t) + 1e-12 * rng.uniform(size=len(cand)))]
        used[j] = True; idx.append(j)
    return np.sort(np.array(idx))


def pair_stats(rt, roi, pred, q):
    """Over query-query pairs with |drt| >= 10 s and droi != 0: concordance on flipped and non-flipped pairs."""
    r, o, p = rt[q], roi[q], pred[q]
    iu = np.triu_indices(len(r), 1)
    drt = r[iu[0]] - r[iu[1]]; dro = o[iu[0]] - o[iu[1]]; dp = p[iu[0]] - p[iu[1]]
    m = (np.abs(drt) >= EPS_FLIP) & (dro != 0)
    drt, dro, dp = drt[m], dro[m], dp[m]
    flip = np.sign(drt) != np.sign(dro)
    conc = np.where(dp == 0, 0.5, (np.sign(dp) == np.sign(drt)).astype(float))
    return dict(n_pairs=int(m.sum()), n_flipped=int(flip.sum()),
                conc_flipped=float(conc[flip].mean()) if flip.any() else np.nan,
                conc_nonflipped=float(conc[~flip].mean()) if (~flip).any() else np.nan,
                conc_all=float(conc.mean()) if len(conc) else np.nan)


def run_unit(args):
    group, d, src, sub, mrow, z1pred, pool_ik, pool_scaf, tabpfn_ks = args
    import arms
    roi = sub[f'roi_{src}'].values.astype(float)
    rt = sub.rt.values.astype(float)
    Xd = sub[DESC].values.astype(float)
    Xm = np.tile(mrow[META].values.astype(float), (len(sub), 1))
    glen = float(mrow.grad_len)
    in_pool = sub.ik14.isin(pool_ik).values
    in_scaf = sub.scaffold.isin(pool_scaf).values & (sub.scaffold.values != '')
    n = len(sub)
    t = time.time(); o1 = arms.isotonic_l1(roi, rt); t_o1 = time.time() - t
    o1_l2 = arms.isotonic_l2(roi, rt)
    o1_cv = arms.isotonic_l1_crossfit(roi, rt)
    rows, prow, timing = [], [], {}
    def add_t(a, dt): timing[a] = timing.get(a, 0.0) + dt
    for scheme in SCHEMES:
        for k in KS:
            if k > n - 10 or (k == 50 and n < 80):
                continue
            for seed in SEEDS:
                a = anchor_sets(rt, d, scheme, k, seed)
                q = np.ones(n, bool); q[a] = False
                preds = {}
                t = time.time(); preds['I1'] = arms.i1_quadratic(roi[a], rt[a], roi); add_t('I1', time.time() - t)
                t = time.time(); preds['I2'] = arms.i2_isotonic(roi[a], rt[a], roi); add_t('I2', time.time() - t)
                t = time.time(); preds['I3'] = arms.i3_linear(roi[a], rt[a], roi); add_t('I3', time.time() - t)
                t = time.time(); preds['I4'] = arms.i4_monotone_spline(roi[a], rt[a], roi); add_t('I4', time.time() - t)
                preds['O1'] = o1; preds['O1_L2'] = o1_l2; preds['O1_CV'] = o1_cv; preds['Z1'] = z1pred
                t = time.time()
                preds['R1'] = arms.r1_ridge_residual(roi[a], rt[a], Xd[a], preds['I1'][a], roi, Xd, preds['I1'])
                add_t('R1', time.time() - t)
                if k in tabpfn_ks:
                    F = np.column_stack([roi, Xd, Xm])
                    t = time.time(); preds['T1'] = arms.tabpfn_arm(F[a], rt[a], F, seed); add_t('T1', time.time() - t)
                    F2 = np.column_stack([F, preds['I1']])
                    t = time.time(); preds['T2'] = arms.tabpfn_arm(F2[a], rt[a], F2, seed); add_t('T2', time.time() - t)
                for ev, mask in [('all', q), ('purged', q & ~in_pool), ('scaffold', q & ~in_scaf)]:
                    nq = int(mask.sum())
                    for arm, p in preds.items():
                        if nq == 0:
                            mae = medae = np.nan
                        else:
                            e = np.abs(p[mask] - rt[mask]); mae, medae = float(e.mean()), float(np.median(e))
                        rows.append(dict(group=group, dataset=d, roi_source=src, scheme=scheme, k=k, seed=seed,
                                         arm=arm, eval_set=ev, n_query=nq, mae_min=mae, medae_min=medae,
                                         mae_pct_grad=100 * mae / glen if nq else np.nan, grad_len=glen))
                if k == 30:
                    for arm, p in preds.items():
                        prow.append(dict(group=group, dataset=d, roi_source=src, scheme=scheme, seed=seed, arm=arm,
                                         **pair_stats(rt, roi, p, q)))
    timing['O1'] = t_o1
    tim = pd.DataFrame([dict(group=group, dataset=d, roi_source=src, arm=a, seconds=s, n=n) for a, s in timing.items()])
    tag = f'{group}_{d}'
    pd.DataFrame(rows).to_parquet(f'{RUNS}/{tag}.parquet')
    pd.DataFrame(prow).to_parquet(f'{RUNS}/{tag}_pairs.parquet')
    tim.to_parquet(f'{RUNS}/{tag}_timing.parquet')
    return tag, n, sum(timing.values())


def init_worker():
    import torch; torch.set_num_threads(1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--no_tabpfn', action='store_true')
    args = ap.parse_args()
    comp, meta, train2s = load()
    z1, purge = pools(comp, meta, train2s)
    json.dump({f'{k[0]}|{k[1]}': v for k, v in z1.items()}, open(f'{OUT}/z1_pools.json', 'w'))
    models, z1t = train_z1(comp, meta, z1)
    pd.DataFrame([dict(model=f'Z1 {k[0]} {k[1]}', train_seconds=v) for k, v in z1t.items()]).to_csv(f'{OUT}/timing_z1.csv', index=False)
    jobs = []
    for group, d, src in units(meta):
        if args.only and group not in args.only: continue
        tag = f'{group}_{d}'
        if os.path.exists(f'{RUNS}/{tag}_timing.parquet'): continue
        mode = meta.loc[d, 'mode']
        sub = comp[comp.dataset == d].dropna(subset=[f'roi_{src}']).reset_index(drop=True)
        X = np.column_stack([sub[f'roi_{src}'].values, sub[DESC].values,
                             np.tile(meta.loc[d, META].values.astype(float), (len(sub), 1))])
        z1pred = models[(mode, src)].predict(X)
        pd_ = purge[(mode, src)] - {d} if group != 'PUB-2step' else purge[(mode, src)]
        pool_c = comp[comp.dataset.isin(pd_)]
        ks = [] if args.no_tabpfn else KS
        jobs.append((group, d, src, sub, meta.loc[d], z1pred, set(pool_c.ik14), set(pool_c.scaffold) - {''}, ks))
    jobs.sort(key=lambda j: -len(j[3]))
    print('units to run', len(jobs), flush=True)
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(args.workers, initializer=init_worker) as pool:
        for i, (tag, n, s) in enumerate(pool.imap_unordered(run_unit, jobs)):
            print(f'[{i+1}/{len(jobs)}] {tag} n={n} {s:.0f}s  elapsed {time.time()-t0:.0f}s', flush=True)
