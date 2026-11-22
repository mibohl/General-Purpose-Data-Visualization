"""Competitor B: explicit LC-metadata conditioning, zero-shot and with k anchors.

Training collection: the PRIDE projects that DeepLC 4.5 was trained on (its
6,543 setup heads, 5,553 of them named by a PXD accession). Each project's
LC metadata comes from its PRIDE record (rtpoc/meta_parse.py). A project is a
training setup if its metadata is "usable" and states one gradient length.
Its RT behaviour is the median over its heads of the head predictions (DeepLC
4.5 fitted these heads to the project's observed RTs, in minutes). The
observed-RT datasets of this repo cannot be pooled instead: the Zenodo sets
are in seconds, aligned units or iRT, and only 7 other held-out setups are in
minutes.

Held-out metadata: data/metadata/curated_lc_metadata.csv (PRIDE read by hand
plus publications), i.e. the best metadata a user could give.

Zero-shot arms (no labels of the held-out setup):
- B_knn: the K = 25 training projects nearest in metadata space (z-scored
  log gradient, log column length, log inner diameter, particle size, log
  flow, temperature, additive, and enrichment tags; mean squared difference
  over the fields both sides report; gradient length required on both
  sides); prediction = per-peptide median of their project predictions.
- B_lgbm: LightGBM fit on (peptide, project) rows: peptide features (as for
  TabPFN), the legacy generic DeepLC prediction and DeepLC 4.5's default
  (PXD005573_mcp) head prediction, plus the metadata fields; target = project
  prediction. Peptide panel: 300 random peptides from each of the other 7
  held-out setups (sequences only, no labels), so no held-out label or
  sequence of the target enters training.
All settings (K, distance, LightGBM parameters, panel) were fixed before
running and are the same for every setup.

Anchored arms at k in {100, 250, 500, 1000}, same splits and seeds:
- <arm>+ridge: least-squares line from the zero-shot prediction to the
  context RTs, then ridge (alpha = 10 on standardized features) on the
  residual with peptide features and the zero-shot prediction as inputs.
- <arm>+iso: isotonic regression from the zero-shot prediction to the
  context RTs (clipped outside the context range).
Output: results/finetune/parts/MB_<name>.csv and results/finetune/meta_zero_shot.csv
Usage: python rtpoc/meta_arms.py
"""
from __future__ import annotations

import os
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import methods as M  # noqa: E402
import run as R  # noqa: E402
from heldout_generic import load_d45  # noqa: E402
from meta_parse import enrichment  # noqa: E402

ROOT = os.path.dirname(HERE)
MD = os.path.join(ROOT, "data", "metadata")
OUT = os.path.join(ROOT, "results", "finetune", "parts")
HELDOUT = {"PXD079751": "easy", "PXD062983": "easy", "PXD076643": "easy", "PXD050641": "hard",
           "PXD085246": "hard", "PXD061731": "hard", "PXD074962": "hard", "PXD084949": "hard"}
KS = [100, 250, 500, 1000]
K_NN = 25
PANEL_PER_SETUP = 300
TAGS = ["TMT", "phospho-enriched", "HLA/MHC peptides", "peptidome", "acetyl-enriched", "diGly", "glyco"]
NUMF = ["log_grad", "log_len", "log_id", "particle_um", "log_flow", "temperature_c", "additive_fa"]
LGB = dict(objective="l1", learning_rate=0.05, num_leaves=63, min_data_in_leaf=50,
           feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, n_estimators=600,
           num_threads=4, verbose=-1)


def meta_frame(df):
    out = pd.DataFrame(index=df.index)
    for src, dst in (("gradient_min", "log_grad"), ("column_length_cm", "log_len"),
                     ("inner_diameter_um", "log_id"), ("flow_nl_min", "log_flow")):
        out[dst] = np.log(pd.to_numeric(df[src], errors="coerce"))
    out["particle_um"] = pd.to_numeric(df["particle_um"], errors="coerce")
    out["temperature_c"] = pd.to_numeric(df["temperature_c"], errors="coerce")
    add = df["additive"].astype(str)
    out["additive_fa"] = np.where(add.str.contains("FA|formic", case=False, regex=True)
                                  & ~add.str.contains("TFA"), 1.0,
                                  np.where(add.str.contains("TFA"), 0.0, np.nan))
    enr = df["enrichment"].fillna("").astype(str)
    for t in TAGS:
        out["tag_" + t] = enr.str.contains(re.escape(t), regex=True).astype(float)
    return out


def training_projects():
    md = pd.read_csv(os.path.join(MD, "pride_lc_metadata.csv"))
    md = md[md.usable & ~md.multi_gradient.fillna(False).astype(bool)].set_index("pxd")
    tn = [t for t in open(os.path.join(MD, "deeplc45_task_names.txt")).read().split("\n") if t]
    heads = {}
    for i, t in enumerate(tn):
        m = re.findall(r"PXD\d{6}", t)
        if m and m[0] in md.index:
            heads.setdefault(m[0], []).append(i)
    md = md.loc[sorted(heads)]
    mcp = tn.index("PXD005573_mcp")
    return md, heads, mcp


def heldout_meta():
    c = pd.read_csv(os.path.join(MD, "curated_lc_metadata.csv"))
    c = c[c.pxd.isin(HELDOUT)].drop_duplicates("pxd").set_index("pxd")
    # tags with the same regexes as for the training projects
    c["enrichment"] = [enrichment(str(e) + " " + str(n), [], []) for e, n in
                       zip(c["enrichment"].fillna(""), c.get("notes", pd.Series("", index=c.index)).fillna(""))]
    return c


def project_preds(D, idx, md, heads, keep=None):
    """Per-project median head prediction on peptides idx -> (len(idx), n_projects)."""
    P, A = [], []
    for c in range(0, len(idx), 4000):
        Ac = np.asarray(D[idx[c:c + 4000]], dtype=np.float32)
        P.append(np.stack([np.median(Ac[:, heads[p]], axis=1) for p in md.index], axis=1))
        A.append(Ac[:, :0] if keep is None else Ac[:, keep])
    return np.vstack(P), np.vstack(A)


def knn_predict(mx, mtrain, P):
    z = mtrain[NUMF + ["tag_" + t for t in TAGS]]
    mu, sd = z.mean(), z.std().replace(0, 1).fillna(1)
    zt = (z - mu) / sd
    zx = (mx[zt.columns] - mu) / sd
    diff = (zt - zx.values) ** 2
    d = diff.mean(axis=1, skipna=True)
    d[zt["log_grad"].isna() | np.isnan(zx["log_grad"])] = np.inf
    nn = np.argsort(d.values)[:K_NN]
    return np.median(P[:, nn], axis=1), list(mtrain.index[nn]), float(np.sort(d.values)[K_NN - 1])


def anchored(zs_c, zs_q, Fc, Fq, y_c):
    from sklearn.isotonic import IsotonicRegression
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    A = np.vstack([zs_c, np.ones_like(zs_c)]).T
    a, b = np.linalg.lstsq(A, y_c, rcond=None)[0]
    lin_c, lin_q = a * zs_c + b, a * zs_q + b
    Xc = np.hstack([Fc, zs_c[:, None]])
    Xq = np.hstack([Fq, zs_q[:, None]])
    sc = StandardScaler().fit(Xc)
    r = Ridge(alpha=10.0).fit(sc.transform(Xc), y_c - lin_c)
    ridge = lin_q + r.predict(sc.transform(Xq))
    iso = IsotonicRegression(out_of_bounds="clip").fit(zs_c, y_c).predict(zs_q)
    return ridge, iso


def main():
    import lightgbm as lgb
    os.makedirs(OUT, exist_ok=True)
    md, heads, mcp = training_projects()
    mtrain = meta_frame(md)
    hm = heldout_meta()
    mheld = meta_frame(hm)
    print(f"training projects {len(md)}, heads {sum(len(v) for v in heads.values())}", flush=True)
    data = {}
    for name in HELDOUT:
        df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{name}.parquet"))
        F = M.peptide_features(df, extra_mods=["TMTpro", "Cysteinyl"])
        data[name] = (df, F, load_d45(name))
    zrows = []
    for name, grp in HELDOUT.items():
        df, F, D = data[name]
        tr = df["tr"].values.astype(float)
        rng = tr.max() - tr.min()
        mx = mheld.loc[name]
        # zero-shot predictions on ALL peptides of the setup (labels unused)
        t0 = time.time()
        used = np.unique(np.concatenate([np.concatenate(R.split(len(df), s_)[:1] + (R.split(len(df), s_)[1][:max(KS)],))
                                         for s_ in range(5)]))
        P_u, A_u = project_preds(D, used, md, heads, keep=[mcp])
        P = np.full((len(df), P_u.shape[1]), np.nan, dtype=np.float32)
        P[used] = P_u
        zk, nn, dk = knn_predict(mx, mtrain, P_u)
        zs_knn = np.full(len(df), np.nan)
        zs_knn[used] = zk
        t_knn = time.time() - t0
        # pooled LightGBM on the panel from the other 7 setups
        t0 = time.time()
        Xs, ys = [], []
        for other in HELDOUT:
            if other == name:
                continue
            odf, oF, oD = data[other]
            pidx = np.random.default_rng(0).choice(len(odf), min(PANEL_PER_SETUP, len(odf)), replace=False)
            oP, oA = project_preds(oD, pidx, md, heads, keep=[mcp])
            pep = np.hstack([oF.values[pidx].astype(np.float32), oA])
            nproj = oP.shape[1]
            Xs.append(np.hstack([np.repeat(pep, nproj, axis=0),
                                 np.tile(mtrain.values.astype(np.float32), (len(pidx), 1))]))
            ys.append(oP.reshape(-1))
        X = np.vstack(Xs)
        y = np.concatenate(ys)
        model = lgb.LGBMRegressor(**LGB).fit(X, y)
        pep_x = np.hstack([F.values[used].astype(np.float32), A_u])
        Xx = np.hstack([pep_x, np.tile(mx[mtrain.columns].values.astype(np.float32), (len(used), 1))])
        zs_lgb = np.full(len(df), np.nan)
        zs_lgb[used] = model.predict(Xx)
        t_lgb = time.time() - t0
        zs_dlc = np.full(len(df), np.nan)
        zs_dlc[used] = A_u[:, 0]
        q0, _ = R.split(len(df), 0)
        for meth, z, t in (("B_knn", zs_knn, t_knn), ("B_lgbm", zs_lgb, t_lgb), ("DeepLC45_zero_shot", zs_dlc, 0.0)):
            zrows.append(dict(dataset=name, group=grp, method=meth, rel_mae=np.abs(z[q0] - tr[q0]).mean() / rng,
                              pearson=np.corrcoef(z[q0], tr[q0])[0, 1], seconds=t,
                              knn_dist=dk if meth == "B_knn" else np.nan,
                              nn_projects=";".join(nn[:10]) if meth == "B_knn" else "",
                              n_train_rows=len(y) if meth == "B_lgbm" else np.nan))
        print(f"{name}: zero-shot rel MAE knn {zrows[-3]['rel_mae']:.4f} lgbm {zrows[-2]['rel_mae']:.4f} "
              f"deeplc {zrows[-1]['rel_mae']:.4f} (rows {len(y)}, {t_lgb:.0f}s)", flush=True)
        # anchored arms
        rows = []
        for seed in range(5):
            q, pool = R.split(len(df), seed)
            for k in KS:
                ctx = pool[:k]
                base = dict(check="MB", dataset=name, group=grp, n=k, seed=seed)
                for arm, z in (("B_knn", zs_knn), ("B_lgbm", zs_lgb)):
                    t0 = time.time()
                    ridge, iso = anchored(z[ctx], z[q], F.values[ctx].astype(float),
                                          F.values[q].astype(float), tr[ctx])
                    dt = time.time() - t0
                    for suffix, p in (("+ridge", ridge), ("+iso", iso)):
                        rows.append(dict(base, method=arm + suffix, subset="all",
                                         rel_mae=np.abs(p - tr[q]).mean() / rng, runtime_s=dt))
                    rows.append(dict(base, method=arm + "_zs", subset="all",
                                     rel_mae=np.abs(z[q] - tr[q]).mean() / rng, runtime_s=0.0))
        pd.DataFrame(rows).to_csv(os.path.join(OUT, f"MB_{name}.csv"), index=False)
    pd.DataFrame(zrows).to_csv(os.path.join(ROOT, "results", "finetune", "meta_zero_shot.csv"), index=False)


if __name__ == "__main__":
    main()
