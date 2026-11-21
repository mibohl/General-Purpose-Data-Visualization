"""Aggregate the DeepLC fine-tuning study.

Inputs: results/finetune/parts/{FT,ST,HR}_<name>.csv, results/parts/H_<name>.csv
(B4 and T1_d45 from the earlier held-out run, same splits and seeds), and
results/finetune/timing_4threads.csv.
Outputs (results/finetune/): results.csv, tables.md, curves_finetune.png
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
OUT = os.path.join(ROOT, "results", "finetune")
KS = [100, 250, 500, 1000]
ORDER = ["B4", "FT", "T1_d45", "T1_d45_FT", "T1_d45_B4", "FT_raw", "B4_check"]
LABEL = {"B4": "DeepLC 4.5 calibration (B4)", "FT": "DeepLC 4.5 fine-tuned (FT)",
         "T1_d45": "TabPFN on DeepLC 4.5 heads (T1+heads)",
         "T1_d45_FT": "TabPFN on heads + fine-tuned DeepLC (hybrid)",
         "T1_d45_B4": "TabPFN on heads + B4 prediction"}
COLOR = {"B4": "#2a78d6", "FT": "#eb6834", "T1_d45": "#1baf7a", "T1_d45_FT": "#6250d6",
         "T1_d45_B4": "#8b939b"}
MARK = {"B4": "v", "FT": "s", "T1_d45": "o", "T1_d45_FT": "D", "T1_d45_B4": "^"}


def fmt(df, nd=2):
    return df.to_markdown(floatfmt=f".{nd}f")


def load():
    parts = [pd.read_csv(f) for f in glob.glob(os.path.join(OUT, "parts", "FT_*.csv"))]
    parts += [pd.read_csv(f) for f in glob.glob(os.path.join(OUT, "parts", "ST_*.csv"))]
    ft = pd.concat(parts, ignore_index=True)
    h = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(ROOT, "results", "parts",
                                                                   "H_*.csv"))])
    h = h[h.method.isin(["B4", "T1_d45"]) & h.n.isin(KS) & h.dataset.isin(ft.dataset.unique())]
    h = h.assign(check="H_reused")
    allr = pd.concat([h, ft], ignore_index=True)
    allr["rel_mae_pct"] = 100 * allr["rel_mae"]
    # ft_d45 counts _validate_epoch calls minus one; deeplc.finetune also calls it once
    # after training (_warn_if_fit_collapsed), so the logged count is one too high.
    allr["ft_epochs"] = allr["ft_epochs"] - 1
    return allr


def per_seed(d, subset="all"):
    return d[d.subset == subset].pivot_table(index=["group", "dataset", "n", "seed"],
                                             columns="method", values="rel_mae_pct")


def ratio_ci(ps, num, den):
    rows = []
    for (g, ds, n), x in ps.groupby(level=["group", "dataset", "n"]):
        x = x[[num, den]].dropna()
        if len(x) < 2:
            continue
        lr = np.log(x[num] / x[den])
        k = len(lr)
        se = lr.std(ddof=1) / np.sqrt(k)
        t = stats.t.ppf(0.975, k - 1)
        rows.append(dict(group=g, dataset=ds, n=n, ratio=np.exp(lr.mean()),
                         lo=np.exp(lr.mean() - t * se), hi=np.exp(lr.mean() + t * se),
                         wins=float((x[num] < x[den]).mean()), seeds=k))
    return pd.DataFrame(rows)


def fig(ps, path):
    meth = [m for m in ["B4", "FT", "T1_d45", "T1_d45_FT"] if m in ps.columns]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0))
    pdm = ps.groupby(level=["group", "dataset", "n"]).mean()
    for ax, g in zip(axes, ("easy", "hard")):
        x = pdm.xs(g, level="group")
        agg = x.groupby(level="n").mean()
        for m in meth:
            ax.plot(agg.index, agg[m], color=COLOR[m], marker=MARK[m], ms=7,
                    lw=2.2, ls="--" if m == "T1_d45_FT" else "-", label=LABEL[m],
                    markeredgecolor="white", markeredgewidth=0.8)
        ax.set_xscale("log")
        ax.set_xticks(KS)
        ax.set_xticklabels([str(k) for k in KS])
        ax.minorticks_off()
        ax.set_xlabel("context size k (identified peptides)")
        ax.set_ylabel("relative MAE, % of RT range (mean over setups)")
        nds = x.index.get_level_values("dataset").nunique()
        ax.set_title(f"{g.capitalize()} held-out setups ({nds} setups)", fontsize=10, loc="left")
        ax.grid(True, axis="y", color="#e4e3df", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=9)
    axes[1].annotate("hybrid and T1+heads overlap", xy=(500, 1.13), xytext=(1000, 1.25), ha="right",
                     fontsize=8, color="#555", arrowprops=dict(arrowstyle="-", color="#999"))
    fig.suptitle("DeepLC 4.5 fine-tuning vs calibration vs TabPFN stacking, 5 seeds",
                 fontsize=11, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    d = load()
    d.to_csv(os.path.join(OUT, "results.csv"), index=False)
    ps = per_seed(d)
    out = []
    meth = [m for m in ORDER if m in ps.columns]
    # group means
    pdm = ps.groupby(level=["group", "dataset", "n"]).mean()
    grp = pdm.groupby(level=["group", "n"]).mean()[meth]
    out.append("## Mean relative MAE, % of RT range (mean over setups of the 5-seed mean)\n")
    out.append(fmt(grp, 3))
    # per setup mean +- sd
    sd = ps.groupby(level=["group", "dataset", "n"]).std()
    tab = pdm[meth].round(3).astype(str) + " ± " + sd[meth].round(3).astype(str)
    out.append("\n## Per setup: relative MAE % (mean ± SD over 5 seeds)\n")
    out.append(tab.to_markdown())
    # geometric-mean ratios per group
    out.append("\n## Geometric-mean ratio over setups (below 1 = numerator better)\n")
    pairs = [("FT", "B4"), ("FT", "T1_d45"), ("FT_raw", "T1_d45"), ("FT_raw", "FT"),
             ("T1_d45_FT", "FT"), ("T1_d45_FT", "T1_d45"),
             ("T1_d45_B4", "T1_d45")]
    gr = {}
    cis = {}
    for a, b in pairs:
        if a in ps.columns and b in ps.columns:
            c = ratio_ci(ps, a, b)
            if c.empty:
                continue
            cis[(a, b)] = c
            gr[f"{a} / {b}"] = c.groupby(["group", "n"]).ratio.apply(
                lambda r: np.exp(np.log(r).mean()))
    out.append(fmt(pd.DataFrame(gr), 3))
    for (a, b), c in cis.items():
        c = c.copy()
        c["ratio [95% CI]"] = c.apply(lambda r: f"{r.ratio:.3f} [{r.lo:.3f}, {r.hi:.3f}]", axis=1)
        c["CI<1"] = c.hi < 1
        c["CI>1"] = c.lo > 1
        out.append(f"\n## Per setup: {a} / {b} with 95% CI over seeds (t on log ratio), "
                   f"share of seeds where {a} wins\n")
        out.append(c.set_index(["group", "dataset", "n"])[["ratio [95% CI]", "wins", "CI<1", "CI>1"]]
                   .to_markdown(floatfmt=".2f"))
        s = c.groupby(["group", "n"]).agg(better=("hi", lambda v: int((v < 1).sum())),
                                          worse=("lo", lambda v: int((v > 1).sum())),
                                          setups=("ratio", "size"))
        out.append(f"\nCount of setups with CI fully below 1 (better) / above 1 (worse):\n")
        out.append(s.to_markdown())
    # modified / unmodified
    for sub in ("modified", "unmodified"):
        p2 = per_seed(d, sub)
        if len(p2):
            m2 = [m for m in meth if m in p2.columns]
            g2 = p2.groupby(level=["group", "dataset", "n"]).mean().groupby(level=["group", "n"]).mean()
            out.append(f"\n## {sub.capitalize()} queries: mean relative MAE %\n")
            out.append(fmt(g2[m2], 3))
    # intervals
    iv = d[(d.subset == "all") & d.coverage90.notna() & d.method.isin(["T1_d45", "T1_d45_FT",
                                                                         "T1_d45_B4"])]
    if len(iv):
        t = iv.groupby(["group", "method", "n"])[["coverage90", "width90"]].mean().unstack("n")
        out.append("\n## TabPFN 90% intervals: coverage and width / range\n")
        out.append(fmt(t, 3))
    # timing (grid runs: 1 thread, 4 jobs on 4 cores)
    a = d[(d.subset == "all") & d.method.isin(["FT", "T1_d45_FT", "T1_d45_B4", "B4_check"])]
    tt = a.groupby(["method", "n"]).runtime_s.median().unstack("n")
    if "ft_train_s" in a:
        tt.loc["FT training only"] = a[a.method == "FT"].groupby("n").ft_train_s.median()
        tt.loc["FT epochs run (median)"] = a[a.method == "FT"].groupby("n").ft_epochs.median()
        tt.loc["FT epochs run (max)"] = a[a.method == "FT"].groupby("n").ft_epochs.max()
    out.append("\n## Grid run seconds per (setup, k, seed), median (1 thread, 4 jobs on 4 "
               "cores; T1_d45_FT includes the full fine-tune and its 3 cross-fit fine-tunes)\n")
    out.append(fmt(tt, 1))
    tp = os.path.join(OUT, "timing_4threads.csv")
    if os.path.exists(tp):
        t4 = pd.read_csv(tp)
        t4["epochs_run"] = t4["epochs_run"] - 1  # same off-by-one as above
        out.append("\n## Fine-tune wall-clock on 4 threads, idle machine (seed 0)\n")
        out.append(fmt(t4.groupby("n")[["ft_train_s", "total_s", "epochs_run"]]
                       .agg(["median", "max"]), 1))
    # headroom
    hr = glob.glob(os.path.join(OUT, "parts", "HR_*.csv"))
    if hr:
        h = pd.concat([pd.read_csv(f) for f in hr])
        h["rel_mae_pct"] = 100 * h.rel_mae
        h["epochs_run"] = h["epochs_run"] - 1  # same off-by-one as above
        w = h.pivot_table(index="dataset", columns="method", values="rel_mae_pct")
        k1 = pdm.xs(1000, level="n").droplevel("group")
        w["FT_k1000"] = k1["FT"]
        w["B4_k1000"] = k1["B4"]
        w["T1_d45_k1000"] = k1["T1_d45"]
        w["headroom (zero_shot - FT_full)"] = w["zero_shot"] - w["FT_full"]
        info = h[h.method == "FT_full"].set_index("dataset")[["n", "epochs_cap", "epochs_run",
                                                              "ft_train_s"]]
        w = w.join(info)
        grp_map = d.drop_duplicates("dataset").set_index("dataset").group
        w.insert(0, "group", grp_map)
        out.append("\n## Headroom: relative MAE % on seed-0 queries\n")
        out.append(w.sort_values("group").to_markdown(floatfmt=".2f"))
    open(os.path.join(OUT, "tables.md"), "w").write("\n".join(out) + "\n")
    fig(ps, os.path.join(OUT, "curves_finetune.png"))
    print("\n".join(out[:3]))


if __name__ == "__main__":
    main()
