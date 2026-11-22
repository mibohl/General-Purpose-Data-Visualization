"""Aggregate the fine-tuning + metadata-conditioning study (results/metadata_study/).

Arms (all on the same 8 held-out setups, splits, seeds and k):
  B4            DeepLC 4.5 default calibration           results/parts/H_*.csv
  T1_d45        TabPFN on DeepLC 4.5 heads (old winner)  results/parts/H_*.csv
  FT            DeepLC 4.5 fine-tuned on the context     results/finetune/parts/FT_*.csv
  T1_d45_FT     TabPFN on heads + fine-tuned prediction  results/finetune/parts/ST_*.csv
  T1_d45_B4     TabPFN on heads + B4 prediction          results/finetune/parts/FT_*.csv
  B_knn / B_lgbm (+ridge, +iso, _zs)  metadata conditioning  results/finetune/parts/MB_*.csv
Intervals: results/finetune/parts/IV_*.csv. Headroom: results/finetune/parts/HR_*.csv
and results/finetune/meta_zero_shot.csv. Timing: results/finetune/timing_4threads.csv
plus runtimes logged by each arm.
Easy (3) and hard (5) setups are never pooled. Competitor B is untestable on
PXD085246 (no recorded gradient), so hard-group comparisons that involve B
arms use the 4 testable hard setups and say so.
"""
from __future__ import annotations

import glob
import os

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FTD = os.path.join(ROOT, "results", "finetune")
OUT = os.path.join(ROOT, "results", "metadata_study")
KS = [100, 250, 500, 1000]
ARMS = ["B4", "FT", "T1_d45", "T1_d45_FT", "T1_d45_B4", "B_knn+ridge", "B_knn+iso",
        "B_lgbm+ridge", "B_lgbm+iso", "B_knn_zs", "B_lgbm_zs"]
NAME = {"B4": "DeepLC 4.5 calibration (B4)", "FT": "DeepLC 4.5 fine-tuned (FT)",
        "T1_d45": "TabPFN on DeepLC 4.5 heads", "T1_d45_FT": "TabPFN on heads + fine-tuned (hybrid)",
        "T1_d45_B4": "TabPFN on heads + B4 prediction",
        "B_knn+ridge": "Metadata kNN + ridge residual", "B_knn+iso": "Metadata kNN + isotonic",
        "B_lgbm+ridge": "Metadata LightGBM + ridge residual", "B_lgbm+iso": "Metadata LightGBM + isotonic",
        "B_knn_zs": "Metadata kNN, zero-shot", "B_lgbm_zs": "Metadata LightGBM, zero-shot"}
# validated categorical palette (dataviz reference instance, fixed slot order)
# 5 hues pass the validator (all pairs); the two metadata arms share the green hue
# and differ by dash and marker. T1_d45_B4 is in the tables, not the plot (it
# overlaps the hybrid).
COLOR = {"B4": "#2a78d6", "FT": "#eda100", "T1_d45": "#1baf7a", "T1_d45_FT": "#e34948",
         "B_knn+ridge": "#008300", "B_lgbm+ridge": "#008300"}
DASH = {"T1_d45_FT": "--", "B_lgbm+ridge": "--"}
MARK = {"B4": "v", "FT": "s", "T1_d45": "o", "T1_d45_FT": "D", "T1_d45_B4": "^",
        "B_knn+ridge": "P", "B_lgbm+ridge": "X"}
UNTESTABLE_B = {"PXD085246"}


def load():
    h = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(ROOT, "results", "parts", "H_*.csv"))])
    h = h[h.method.isin(["B4", "T1_d45"]) & h.n.isin(KS) & (h.subset == "all")]
    ft = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "FT_*.csv"))]
                   + [pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "ST_*.csv"))])
    ft = ft[ft.method.isin(["FT", "T1_d45_FT", "T1_d45_B4"]) & (ft.subset == "all")]
    mb = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "MB_*.csv"))])
    h = h[h.dataset.isin(ft.dataset.unique())]
    cols = ["dataset", "group", "n", "seed", "method", "rel_mae", "runtime_s"]
    d = pd.concat([h[cols + ["coverage90", "width90"]], ft[cols + ["coverage90", "width90"]],
                   mb[cols]], ignore_index=True)
    d["rel_mae_pct"] = 100 * d.rel_mae
    return d


def ratio_ci(ps, num, den):
    rows = []
    for (g, ds, n), x in ps.groupby(level=["group", "dataset", "n"]):
        if num not in x or den not in x:
            continue
        x = x[[num, den]].dropna()
        if len(x) < 2:
            continue
        lr = np.log(x[num] / x[den])
        k = len(lr)
        se = lr.std(ddof=1) / np.sqrt(k)
        t = stats.t.ppf(0.975, k - 1)
        rows.append(dict(group=g, dataset=ds, n=n, method=num, ratio=np.exp(lr.mean()),
                         lo=np.exp(lr.mean() - t * se), hi=np.exp(lr.mean() + t * se),
                         wins=float((x[num] < x[den]).mean())))
    return pd.DataFrame(rows)


def gmean(v):
    return float(np.exp(np.log(v).mean()))


def fig(ps, path):
    meth = [m for m in COLOR if m in ps.columns]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    pdm = ps.groupby(level=["group", "dataset", "n"]).mean()
    for ax, g in zip(axes, ("easy", "hard")):
        x = pdm.xs(g, level="group")
        x = x[~x.index.get_level_values("dataset").isin(UNTESTABLE_B)]
        agg = x.groupby(level="n").mean()
        for m in meth:
            ax.plot(agg.index, agg[m], color=COLOR[m], marker=MARK[m], ms=7, lw=2.0,
                    ls=DASH.get(m, "-"), label=NAME[m], markeredgecolor="white", markeredgewidth=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(KS)
        ax.set_xticklabels([str(k) for k in KS])
        ax.minorticks_off()
        yt = {"easy": [0.35, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2], "hard": [1.0, 1.2, 1.5, 2.0, 2.5, 3.0]}[g]
        ax.set_yticks(yt)
        ax.set_yticklabels([f"{v:g}" for v in yt])
        ax.set_xlabel("context size k (identified peptides)")
        ax.set_ylabel("relative MAE, % of RT range (log scale)")
        nds = x.index.get_level_values("dataset").nunique()
        extra = "; PXD085246 left out: no LC metadata" if g == "hard" else ""
        ax.set_title(f"{g.capitalize()} held-out setups ({nds}{extra})", fontsize=10, loc="left")
        ax.grid(True, axis="y", color="#e4e3df", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle("All arms vs context size, mean over setups of the 5-seed mean", fontsize=11, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    d = load()
    d.to_csv(os.path.join(OUT, "results.csv"), index=False)
    ps = d.pivot_table(index=["group", "dataset", "n", "seed"], columns="method", values="rel_mae_pct")
    meth = [m for m in ARMS if m in ps.columns]
    out = []
    pdm = ps.groupby(level=["group", "dataset", "n"]).mean()
    sd = ps.groupby(level=["group", "dataset", "n"]).std()
    for g in ("easy", "hard"):
        for subset, lab in ((False, "all setups"), (True, "B-testable setups only")):
            x = pdm.xs(g, level="group", drop_level=False)
            if subset:
                if g == "easy":
                    continue
                x = x[~x.index.get_level_values("dataset").isin(UNTESTABLE_B)]
            m = x.groupby(level="n").mean()[meth]
            out.append(f"## {g.capitalize()}, {lab} ({x.index.get_level_values('dataset').nunique()}): "
                       "mean relative MAE, % of RT range\n")
            out.append(m.T.to_markdown(floatfmt=".3f"))
            out.append("")
    tab = pdm[meth].round(3).astype(str) + " ± " + sd[meth].round(3).astype(str)
    out.append("## Per setup: relative MAE % (mean ± SD over 5 seeds)\n")
    out.append(tab.to_markdown())
    # ratios to B4
    cis = pd.concat([ratio_ci(ps, m, "B4") for m in meth if m != "B4"], ignore_index=True)
    cis.to_csv(os.path.join(OUT, "ratio_to_b4_ci.csv"), index=False)
    s = cis.groupby(["group", "method", "n"]).agg(
        gmean_ratio=("ratio", gmean), ci_below_1=("hi", lambda v: int((v < 1).sum())),
        ci_above_1=("lo", lambda v: int((v > 1).sum())), setups=("ratio", "size")).reset_index()
    for g in ("easy", "hard"):
        t = s[s.group == g].copy()
        t["cell"] = t.apply(lambda r: f"{r.gmean_ratio:.3f} ({r.ci_below_1}/{r.ci_above_1} of {r.setups})", axis=1)
        out.append(f"\n## {g.capitalize()}: ratio to B4 (geometric mean over setups), in brackets: "
                   "setups with 95% CI below / above 1\n")
        out.append(t.pivot(index="method", columns="n", values="cell").reindex(
            [m for m in meth if m != "B4"]).to_markdown())
    out.append("\n## Per setup ratio to B4 with 95% CI (t-interval on the log ratio over 5 seeds)\n")
    cis["cell"] = cis.apply(lambda r: f"{r.ratio:.3f} [{r.lo:.3f}, {r.hi:.3f}]", axis=1)
    out.append(cis.pivot_table(index=["group", "dataset", "method"], columns="n", values="cell",
                               aggfunc="first").to_markdown())
    # head-to-head ratios of interest
    pairs = [("FT", "T1_d45"), ("T1_d45_FT", "FT"), ("T1_d45_FT", "T1_d45"), ("T1_d45_B4", "T1_d45"),
             ("B_knn+ridge", "T1_d45"), ("B_lgbm+ridge", "T1_d45"), ("B_knn+ridge", "FT"),
             ("B_lgbm+ridge", "FT")]
    hh = []
    for a, b in pairs:
        c = ratio_ci(ps, a, b)
        if c.empty:
            continue
        for (g, n), x in c.groupby(["group", "n"]):
            hh.append(dict(pair=f"{a} / {b}", group=g, n=n, gmean=gmean(x.ratio),
                           below=int((x.hi < 1).sum()), above=int((x.lo > 1).sum()), setups=len(x)))
    hh = pd.DataFrame(hh)
    hh["cell"] = hh.apply(lambda r: f"{r.gmean:.3f} ({r.below}/{r.above} of {r.setups})", axis=1)
    out.append("\n## Head-to-head geometric-mean ratios (CI below / above 1)\n")
    out.append(hh.pivot_table(index=["group", "pair"], columns="n", values="cell", aggfunc="first").to_markdown())
    hh.to_csv(os.path.join(OUT, "head_to_head.csv"), index=False)
    # intervals
    ivf = glob.glob(os.path.join(FTD, "parts", "IV_*.csv"))
    if ivf:
        iv = pd.concat([pd.read_csv(f) for f in ivf])
        b4 = d[(d.method == "B4")][["dataset", "group", "n", "seed", "coverage90", "width90", "rel_mae"]]
        b4 = b4.assign(method="B4 window", width90_median=b4.width90)
        iv = pd.concat([iv, b4], ignore_index=True)
        t = iv.groupby(["group", "method", "n"]).agg(coverage=("coverage90", "mean"),
                                                    median_width_pct=("width90_median", lambda v: 100 * v.mean()),
                                                    runs=("seed", "size")).reset_index()
        out.append("\n## 90% intervals: empirical coverage and median width (% of RT range), mean over "
                   "setups and seeds. FT = split-conformal on its 15% validation part; FT_spread = the "
                   "central 90% range of FT's query residuals (an oracle reference, not a method); B4 "
                   "window = constant cross-fitted window from the old study\n")
        t["cell"] = t.apply(lambda r: f"{r.coverage:.3f} / {r.median_width_pct:.2f}", axis=1)
        out.append(t.pivot_table(index=["group", "method"], columns="n", values="cell", aggfunc="first").to_markdown())
        if "max_abs_diff_vs_main" in iv:
            out.append(f"\nFT rerun vs the saved predictions: max |difference| = "
                       f"{iv.max_abs_diff_vs_main.max():.2e} of the RT range.")
        chk = iv[iv.method.isin(["T1_d45", "T1_d45_B4"])].merge(
            d[d.method.isin(["T1_d45", "T1_d45_B4"])][["dataset", "n", "seed", "method", "rel_mae"]],
            on=["dataset", "n", "seed", "method"], suffixes=("", "_stored"))
        if len(chk):
            out.append(f"TabPFN rerun vs stored rel. MAE: max |difference| = "
                       f"{(chk.rel_mae - chk.rel_mae_stored).abs().max():.2e} over {len(chk)} runs.")
    # headroom
    hr = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(FTD, "parts", "HR_*.csv"))])
    zs = pd.read_csv(os.path.join(FTD, "meta_zero_shot.csv"))
    w = hr.pivot_table(index="dataset", columns="method", values="rel_mae") * 100
    z = zs.pivot_table(index="dataset", columns="method", values="rel_mae") * 100
    w = w.join(z[[c for c in ("B_knn", "B_lgbm") if c in z]])
    w["metadata zero-shot (best)"] = w[["B_knn", "B_lgbm"]].min(axis=1)
    w["best zero-shot"] = w[["zero_shot", "metadata zero-shot (best)"]].min(axis=1)
    w["gap = best zero-shot - FT largest"] = w["best zero-shot"] - w["FT_full"]
    w["FT context"] = hr[hr.method == "FT_full"].set_index("dataset").n
    grp = d.drop_duplicates("dataset").set_index("dataset").group
    w.insert(0, "group", grp)
    w = w.rename(columns={"zero_shot": "DeepLC 4.5 zero-shot", "FT_full": "FT largest context",
                          "B_knn": "metadata kNN zero-shot", "B_lgbm": "metadata LightGBM zero-shot"})
    keep = ["group", "DeepLC 4.5 zero-shot", "metadata kNN zero-shot", "metadata LightGBM zero-shot",
            "best zero-shot", "FT largest context", "FT context", "gap = best zero-shot - FT largest"]
    w = w[keep].sort_values(["group", "gap = best zero-shot - FT largest"])
    w.to_csv(os.path.join(OUT, "headroom.csv"))
    out.append("\n## Headroom (relative MAE %, seed-0 queries)\n")
    out.append(w.to_markdown(floatfmt=".2f"))
    med = w.groupby("group")[["DeepLC 4.5 zero-shot", "best zero-shot", "FT largest context",
                              "gap = best zero-shot - FT largest"]].median()
    out.append("\nGroup medians:\n")
    out.append(med.to_markdown(floatfmt=".2f"))
    open(os.path.join(OUT, "tables.md"), "w").write("\n".join(out) + "\n")
    fig(ps, os.path.join(OUT, "curves_all_arms.png"))
    print("\n".join(out[:8]))


if __name__ == "__main__":
    main()
