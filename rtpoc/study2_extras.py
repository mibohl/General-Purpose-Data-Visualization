"""Spearman headroom and timing / compute-ratio tables (results/metadata_study/).

Spearman (seed-0 queries; the small-molecule RepoRT result was stated in
Spearman): DeepLC 4.5 zero-shot, metadata kNN and LightGBM zero-shot, and at
k = 100 and 1,000: B4, FT, TabPFN on heads (T1_d45). Gap = best k = 1,000
method minus best zero-shot.

Timing: wall-clock per (setup, k, seed) fit-and-predict for 2,000 queries,
median over setups and seeds, with the thread setting it was measured at,
and core-seconds (seconds x threads). "Adequate" (rule fixed before
computing): geometric-mean error ratio to B4 within 5% of the best arm of
that group at that k. Compute ratio = core-seconds / core-seconds of the
cheapest adequate arm.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run as R  # noqa: E402
from heldout_generic import load_d45  # noqa: E402
from run_ft import b4  # noqa: E402

ROOT = os.path.dirname(HERE)
FTD = os.path.join(ROOT, "results", "finetune")
OUT = os.path.join(ROOT, "results", "metadata_study")
HELDOUT = {"PXD079751": "easy", "PXD062983": "easy", "PXD076643": "easy", "PXD050641": "hard",
           "PXD085246": "hard", "PXD061731": "hard", "PXD074962": "hard", "PXD084949": "hard"}


def spearman(a, b):
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def spearman_table():
    rows = []
    for name, g in HELDOUT.items():
        df = pd.read_parquet(os.path.join(ROOT, "data", "processed", f"heldout_{name}.parquet"))
        tr = df.tr.values.astype(float)
        q, pool = R.split(len(df), 0)
        r = dict(dataset=name, group=g)
        mp = os.path.join(FTD, "metapreds", f"{name}.npz")
        D = load_d45(name)
        if os.path.exists(mp):
            z = np.load(mp)
            pos = {v: i for i, v in enumerate(z["used"])}
            ii = [pos[v] for v in q]
            for m in ("B_knn", "B_lgbm"):
                r[m + " zero-shot"] = spearman(z[m][ii], tr[q])
        tn = [t for t in open(os.path.join(ROOT, "data", "metadata", "deeplc45_task_names.txt")).read().split("\n") if t]
        r["DeepLC 4.5 zero-shot"] = spearman(np.asarray(D[q])[:, tn.index("PXD005573_mcp")], tr[q])
        for k in (100, 1000):
            ctx = pool[:k]
            lo, hi = tr[ctx].min(), tr[ctx].max()
            yc = (tr[ctx] - lo) / (hi - lo)
            pq, _, _, _ = b4(D, ctx, q, yc, 0)
            r[f"B4 k={k}"] = spearman(pq, tr[q])
            iv = np.load(os.path.join(FTD, "ivpreds", f"{name}_s0_k{k}.npz"))
            r[f"FT k={k}"] = spearman(iv["FT"], tr[q])
            r[f"T1_d45 k={k}"] = spearman(iv["T1_d45"][0], tr[q])
        rows.append(r)
    t = pd.DataFrame(rows).set_index(["group", "dataset"])
    zs = [c for c in t if "zero-shot" in c]
    t["best zero-shot"] = t[zs].max(axis=1)
    t["best at k=1000"] = t[[c for c in t if c.endswith("k=1000")]].max(axis=1)
    t["Spearman gap (k=1000 best - zero-shot best)"] = t["best at k=1000"] - t["best zero-shot"]
    return t


def timing_table():
    h = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(ROOT, "results", "parts", "H_*.csv"))])
    h = h[(h.subset == "all") & h.dataset.isin(HELDOUT) & h.n.isin([100, 250, 500, 1000])]
    iv = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "IV_*.csv"))])
    st = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "ST_*.csv"))])
    st = st[st.subset == "all"]
    t4 = pd.read_csv(os.path.join(FTD, "timing_4threads.csv"))
    mb = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "MB_*.csv"))])
    zs = pd.read_csv(os.path.join(FTD, "meta_zero_shot.csv"))
    b4_s = h[h.method == "B4"].groupby("n").runtime_s.median()
    rows = []

    def add(arm, ser, threads, note):
        for n, v in ser.items():
            rows.append(dict(arm=arm, n=n, seconds=v, threads=threads, core_s=v * threads, note=note))

    add("B4", b4_s, 1, "DeepLC 4.5 calibration fit + predict (head matrix precomputed)")
    add("FT", t4.groupby("n").total_s.median(), 4, "finetune + calibrate + predict, idle machine")
    add("T1_d45", iv[iv.method == "T1_d45"].groupby("n").seconds.median() + b4_s, 1,
        "TabPFN fit + predict (4 estimators) + B4 head ranking")
    add("T1_d45_B4", iv[iv.method == "T1_d45_B4"].groupby("n").seconds.median() + 6 * b4_s, 1,
        "as T1_d45 plus 5-fold cross-fitted B4")
    add("T1_d45_FT", st[st.method == "T1_d45_FT"].groupby("n").runtime_s.median(), 1,
        "full fine-tune + 3 cross-fit fine-tunes + TabPFN (grid, 4 jobs on 4 cores)")
    for arm in ("B_knn", "B_lgbm"):
        fit = zs[zs.method == arm].seconds.median()
        anc = mb[mb.method == arm + "+ridge"].groupby("n").runtime_s.median()
        add(arm + "+ridge (per new setup)", anc + fit, 4 if arm == "B_lgbm" else 1,
            f"zero-shot model {fit:.0f} s per new setup + anchored ridge")
        add(arm + "+ridge (model trained once)", anc, 1, "anchored ridge only")
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    sp = spearman_table()
    sp.to_csv(os.path.join(OUT, "spearman_headroom.csv"))
    out = ["## Spearman on seed-0 queries\n", sp.to_markdown(floatfmt=".4f"), "\nGroup medians:\n",
           sp.groupby(level="group").median().T.to_markdown(floatfmt=".4f")]
    tt = timing_table()
    hh = pd.read_csv(os.path.join(OUT, "ratio_to_b4_ci.csv"))
    gm = hh.groupby(["group", "method", "n"]).ratio.apply(lambda v: float(np.exp(np.log(v).mean()))).reset_index()
    gm = pd.concat([gm, pd.DataFrame([dict(group=g, method="B4", n=n, ratio=1.0)
                                      for g in ("easy", "hard") for n in (100, 250, 500, 1000)])])
    arm_of = {"B_knn+ridge (per new setup)": "B_knn+ridge", "B_knn+ridge (model trained once)": "B_knn+ridge",
              "B_lgbm+ridge (per new setup)": "B_lgbm+ridge", "B_lgbm+ridge (model trained once)": "B_lgbm+ridge"}
    tt["method"] = tt.arm.map(lambda a: arm_of.get(a, a))
    res = []
    for g in ("easy", "hard"):
        for n in (100, 250, 500, 1000):
            e = gm[(gm.group == g) & (gm.n == n)].set_index("method").ratio
            best = e.min()
            x = tt[tt.n == n].copy()
            x["err_ratio_to_B4"] = x.method.map(e)
            x["adequate"] = x.err_ratio_to_B4 <= best * 1.05
            cheap = x[x.adequate].core_s.min()
            x["compute_ratio_to_cheapest_adequate"] = x.core_s / cheap
            x["group"] = g
            res.append(x)
    res = pd.concat(res)
    res.to_csv(os.path.join(OUT, "timing.csv"), index=False)
    for g in ("easy", "hard"):
        x = res[res.group == g]
        piv = x.pivot_table(index="arm", columns="n", values="compute_ratio_to_cheapest_adequate")
        sec = x[x.n == 250].set_index("arm")[["seconds", "threads", "core_s", "note"]]
        out.append(f"\n## {g.capitalize()}: compute ratio to the cheapest adequate arm (core-seconds; "
                   "adequate = within 5% of the best error)\n")
        out.append(piv.to_markdown(floatfmt=".2f"))
        ad = x[x.adequate].groupby("n").apply(lambda d: ", ".join(sorted(d.arm)), include_groups=False)
        out.append("\nAdequate arms per k: " + "; ".join(f"k={k}: {v}" for k, v in ad.items()))
    out.append("\n## Seconds per fit + 2,000-query prediction at k = 250 (median over setups and seeds)\n")
    out.append(res[(res.n == 250) & (res.group == "easy")].set_index("arm")[["seconds", "threads", "core_s", "note"]]
               .to_markdown(floatfmt=".1f"))
    out.append("\nAll k:\n")
    out.append(res[res.group == "easy"].pivot_table(index="arm", columns="n", values="seconds").to_markdown(floatfmt=".1f"))
    open(os.path.join(OUT, "extras.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
