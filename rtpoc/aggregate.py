"""Collect results/parts/*.csv into results/results.csv, tables, and figures."""
from __future__ import annotations

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from datasets import EASY, HARD  # noqa: E402

ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
COLS = ["dataset", "group", "n", "seed", "method", "rel_mae", "rel_dt95", "coverage90",
        "width90", "runtime_s"]
METHODS = ["B1", "B2", "B2_deeplc", "B3_ridge", "B3_xgb", "T1", "T2"]
# Reference categorical palette (dataviz skill), fixed slot order; validated adjacent pairs.
COLORS = dict(zip(METHODS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                            "#008300", "#6250d6"]))
MARKERS = dict(zip(METHODS, ["o", "s", "D", "^", "v", "P", "X"]))
LABELS = {"B1": "B1 linear", "B2": "B2 LOWESS", "B2_deeplc": "B2 DeepLC spline",
          "B3_ridge": "B3 ridge", "B3_xgb": "B3 XGBoost", "T1": "T1 TabPFN direct",
          "T2": "T2 TabPFN residual"}


def load():
    parts = {c: pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{RES}/parts/{c}_*.csv"))],
                          ignore_index=True)
             for c in ("C1", "C2", "C3") if glob.glob(f"{RES}/parts/{c}_*.csv")}
    return parts


def fmt(df, digits=4):
    return df.to_markdown(floatfmt=f".{digits}f")


def curves(c2, group, path):
    d = c2[(c2.group == group) & (c2.subset == "all")]
    # per dataset: mean over seeds; then across datasets
    per_ds = d.groupby(["dataset", "method", "n"]).rel_mae.mean().reset_index()
    b2 = per_ds[per_ds.method == "B2"].set_index(["dataset", "n"]).rel_mae
    per_ds["ratio"] = per_ds.rel_mae.values / b2.loc[list(zip(per_ds.dataset, per_ds.n))].values
    agg = per_ds.groupby(["method", "n"]).agg(
        rel_mae=("rel_mae", "mean"), ratio=("ratio", lambda x: np.exp(np.mean(np.log(x)))))
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#8a8986", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ns = sorted(d.n.unique())
    for ax, col, ylab in ((axes[0], "rel_mae", "relative MAE (mean over datasets)"),
                          (axes[1], "ratio", "relative MAE / B2 (geometric mean)")):
        for m in METHODS:
            if m not in agg.index.get_level_values(0):
                continue
            s = agg.loc[m][col]
            lw = 2.4 if m in ("B2", "T2") else 1.6
            ax.plot(s.index, s.values, color=COLORS[m], marker=MARKERS[m], ms=6, lw=lw,
                    label=LABELS[m], markeredgecolor="white", markeredgewidth=0.8)
        ax.set_xscale("log")
        ax.set_xticks(ns)
        ax.set_xticklabels([str(n) for n in ns])
        ax.set_xlabel("context size n (identified peptides)")
        ax.set_ylabel(ylab)
        ax.grid(True, axis="y", color="#e4e3df", lw=0.8)
        if col == "ratio":
            ax.axhline(1.0, color="#8a8986", lw=1, ls="--")
            ax.axhline(0.9, color="#c3c2b7", lw=1, ls=":")
            ax.text(ns[-1], 0.9, " -10%", va="center", fontsize=8, color="#52514e")
    axes[0].legend(frameon=False, fontsize=8.5)
    names = EASY if group == "easy" else HARD
    fig.suptitle(f"{group.capitalize()} datasets (k={len(names)}): RT error vs. context size, "
                 f"2,000 queries, 5 seeds", fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return agg


def per_dataset_grid(c2, path):
    d = c2[c2.subset == "all"].groupby(["dataset", "method", "n"]).rel_mae.mean().reset_index()
    names = EASY + HARD
    fig, axes = plt.subplots(4, 4, figsize=(14, 12), sharex=True)
    for ax, name in zip(axes.flat, names):
        x = d[d.dataset == name]
        for m in METHODS:
            s = x[x.method == m]
            ax.plot(s.n, s.rel_mae, color=COLORS[m], marker=MARKERS[m], ms=4,
                    lw=2 if m in ("B2", "T2") else 1.2, label=LABELS[m])
        ax.set_xscale("log")
        ax.set_title(("[hard] " if name in HARD else "[easy] ") + name.replace("_psms_aligned", "")
                     .replace("_fixed_mods", ""), fontsize=9)
        ax.grid(True, axis="y", color="#e4e3df", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[-1]:
        ax.set_xticks([25, 100, 250, 1000])
        ax.set_xticklabels(["25", "100", "250", "1000"])
        ax.set_xlabel("n")
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle("Relative MAE vs. context size, per dataset (mean of 5 seeds)", x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main():
    P = load()
    out = []
    c2 = P.get("C2")
    c3 = P.get("C3")
    c1 = P.get("C1")
    frames = []
    for k, df in P.items():
        df = df.copy()
        df["check"] = k
        frames.append(df)
    allr = pd.concat(frames, ignore_index=True)
    extra = [c for c in allr.columns if c not in COLS]
    allr = allr[COLS + extra]
    allr.to_csv(os.path.join(RES, "results.csv"), index=False)

    if c1 is not None:
        t = c1.groupby(["group", "dataset", "subset", "method"]).agg(
            resid_r2=("resid_r2", "mean"), rel_mae=("rel_mae", "mean"),
            b2_rel_mae=("b2_rel_mae", "mean"), n_ctx=("n", "mean")).reset_index()
        t["mae_gain_vs_b2"] = 1 - t.rel_mae / t.b2_rel_mae
        out.append("## C1 headroom (context = all but 2,000; mean of 3 seeds)\n")
        piv = t[t.subset == "all"].pivot_table(index=["group", "dataset"], columns="method",
                                              values=["resid_r2", "mae_gain_vs_b2"])
        out.append(fmt(piv, 3))
        piv2 = t[t.subset != "all"].pivot_table(index=["group", "dataset", "subset"],
                                               columns="method", values="resid_r2")
        out.append("\nResidual R^2 by modified / unmodified queries:\n")
        out.append(fmt(piv2, 3))

    if c2 is not None:
        for g in ("easy", "hard"):
            agg = curves(c2, g, os.path.join(RES, f"curves_{g}.png"))
            out.append(f"\n## C2 {g}: mean relative MAE over datasets\n")
            out.append(fmt(agg.rel_mae.unstack("n").reindex(METHODS)))
            out.append(f"\n## C2 {g}: geometric-mean ratio of relative MAE to B2\n")
            out.append(fmt(agg.ratio.unstack("n").reindex(METHODS), 3))
            dt = c2[(c2.group == g) & (c2.subset == "all")].groupby(
                ["method", "n"]).rel_dt95.mean().unstack("n").reindex(METHODS)
            out.append(f"\n## C2 {g}: mean relative delta-t95 (2 x q95|err| / range)\n")
            out.append(fmt(dt))
        per_dataset_grid(c2, os.path.join(RES, "curves_per_dataset.png"))
        # per-dataset T2 vs B2 / B3 at n in {100,250}
        d = c2[(c2.subset == "all") & c2.n.isin([100, 250])]
        pd_ = d.groupby(["group", "dataset", "n", "method"]).rel_mae.mean().unstack("method")
        pd_["T2_vs_B2"] = pd_["T2"] / pd_["B2"] - 1
        pd_["T1_vs_B2"] = pd_["T1"] / pd_["B2"] - 1
        pd_["bestB3"] = pd_[["B3_ridge", "B3_xgb"]].min(axis=1)
        pd_["T2_vs_bestB3"] = pd_["T2"] / pd_["bestB3"] - 1
        out.append("\n## C2 per dataset at n=100, 250 (mean of 5 seeds)\n")
        out.append(fmt(pd_[["B2", "B2_deeplc", "B3_ridge", "B3_xgb", "T1", "T2", "T2_vs_B2",
                            "T1_vs_B2", "T2_vs_bestB3"]]))
        # modified vs unmodified
        sub = c2[c2.subset != "all"].groupby(["group", "subset", "method", "n"]).rel_mae.mean() \
            .unstack("n").reindex(METHODS, level="method")
        out.append("\n## C2 relative MAE by modified / unmodified queries (mean over datasets)\n")
        out.append(fmt(sub))
        # C4 intervals
        iv = c2[(c2.subset == "all") & c2.method.isin(["B2", "T1", "T2"])].groupby(
            ["group", "method", "n"]).agg(coverage90=("coverage90", "mean"),
                                          width90=("width90", "mean")).unstack("n")
        out.append("\n## C4 90% intervals: coverage and mean width / range\n")
        out.append(fmt(iv, 3))
        # C5 runtime
        rt = c2[(c2.subset == "all") & c2.method.isin(["B2", "B3_xgb", "T1", "T2"])].groupby(
            ["method", "n"]).runtime_s.agg(["median", "max"]).unstack("n")
        out.append("\n## C5 wall-clock seconds per fit+predict (2,000 queries; 1 CPU thread, "
                   "4 jobs sharing 4 cores; median and max over datasets x seeds)\n")
        out.append(fmt(rt, 1))
        rtd = c2[(c2.subset == "all") & c2.method.isin(["T1", "T2"])].groupby(
            ["dataset", "method"]).runtime_s.sum().unstack("method") / 5
        out.append("\n## C5 TabPFN seconds per dataset per seed (all 6 n summed)\n")
        out.append(fmt(rtd, 0))

    if c3 is not None:
        d = c3[c3.subset == "all"]
        clean = c2[(c2.subset == "all") & (c2.n == 250) & c2.method.isin(["B2", "B3_ridge",
                                                                             "B3_xgb", "T2"])]
        dd = pd.concat([clean.assign(noise=0.0), d])
        per = dd.groupby(["group", "dataset", "method", "noise"]).rel_mae.mean().reset_index()
        t3 = per.groupby(["group", "method", "noise"]).rel_mae.mean().unstack("noise")
        out.append("\n## C3 label noise at n=250: mean relative MAE over datasets\n")
        out.append(fmt(t3))

    with open(os.path.join(RES, "tables.md"), "w") as f:
        f.write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
