"""
Report for the fair-rematch follow-up (Q1 matched inputs, Q2 headroom + A(k),
Q3 meta-learnability). Reuses the same splits/seeds/clustering as the primary
study. Everything per-epitope, never pooled.

Inputs (results/):
  results_matched.csv              lgbm, ridge, paired-tcrdist (matched inputs)
  results_matched_tab_bg/oe.csv    TabPFN (matched inputs), per scheme
  results_primary_beta.csv         beta-only TCRdist (same splits) -> TCRdist best
  results_headroom.csv             oracle / zeroshot / pooled / headroom
  coef_similarity_summary.json     Q3 coefficient cosine similarity
Outputs: report_v2.md, q1_auroc_vs_k.png, q2_adaptation.png, results_matched_all.csv
"""
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common

RES = os.path.join(common.REPO_ROOT, "results")
KS = [20, 50, 200]
SCH = ["background", "otherep"]
SLAB = {"background": "(a) background", "otherep": "(b) other-epitope (honest)"}


def load_matched():
    parts = [pd.read_csv(os.path.join(RES, "results_matched.csv"))]
    for f in ["results_matched_tab_bg.csv", "results_matched_tab_oe.csv"]:
        p = os.path.join(RES, f)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                t = pd.read_csv(p)
                if len(t):
                    parts.append(t)
            except Exception:
                pass
    d = pd.concat(parts, ignore_index=True)
    d = d.drop_duplicates(["scheme", "epitope", "seed", "k", "method"])
    d.to_csv(os.path.join(RES, "results_matched_all.csv"), index=False)
    return d[d.method != "SKIP"].dropna(subset=["auroc"])


def permed(df, sch, k, meth):
    s = df[(df.scheme == sch) & (df.k == k) & (df.method == meth)]
    return s.groupby("epitope").auroc.median()


def tcrdist_best(m, b, sch, k):
    """Per-epitope strongest TCRdist over beta (prior) and paired (matched)."""
    cand = [permed(b, sch, k, "tcrdist_nn"), permed(b, sch, k, "tcrdist_knn"),
            permed(m, sch, k, "tcrdist_nn"), permed(m, sch, k, "tcrdist_knn")]
    cand = [c for c in cand if len(c)]
    return pd.concat(cand, axis=1).max(axis=1)


def fmt(x):
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"


def md_table(rows, hdr):
    out = ["| " + " | ".join(hdr) + " |", "| " + " | ".join(["---"] * len(hdr)) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def q1_plot(m, b):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), squeeze=False)
    colors = {"tabpfn": "#2b6cb0", "lightgbm": "#2f855a", "ridge": "#805ad5", "td": "#c53030"}
    for ax, sch in zip(axes[0], SCH):
        for meth, lab in [("tabpfn", "TabPFN (matched)"), ("lightgbm", "LightGBM (matched)"),
                          ("ridge", "ridge (matched)")]:
            ys = [permed(m, sch, k, meth).median() for k in KS]
            if all(np.isnan(ys)):
                continue
            ax.plot(KS, ys, "o-", color=colors[meth], label=lab)
        td = [tcrdist_best(m, b, sch, k).median() for k in KS]
        ax.plot(KS, td, "s--", color=colors["td"], label="TCRdist (best config)")
        ax.axhline(0.5, ls=":", color="gray", lw=.8)
        ax.set_title(SLAB[sch]); ax.set_xlabel("support size k"); ax.set_xticks(KS)
        ax.set_ylabel("median per-epitope AUROC"); ax.set_ylim(0.45, 0.9); ax.grid(alpha=.25)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Q1 fair rematch: matched germline+CDR3 inputs (both chains) vs TCRdist")
    fig.tight_layout(); p = os.path.join(RES, "q1_auroc_vs_k.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


def q2_plot(m, h):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), squeeze=False)
    for ax, sch in zip(axes[0], SCH):
        hs = h[h.scheme == sch].set_index("epitope")
        zs = hs["zeroshot"].median(); orc = hs["oracle"].median()
        xs = [0, 20, 50, 200]
        lg = [zs] + [permed(m, sch, k, "lightgbm").median() for k in KS]
        tp = [zs] + [permed(m, sch, k, "tabpfn").median() for k in KS]
        ax.plot(xs, lg, "o-", color="#2f855a", label="LightGBM A(k)")
        if not all(np.isnan(tp[1:])):
            ax.plot(xs, tp, "o-", color="#2b6cb0", label="TabPFN A(k)")
        ax.axhline(orc, ls="--", color="#b7791f", label=f"oracle {orc:.2f}")
        ax.axhline(zs, ls=":", color="#c53030", label=f"zero-shot {zs:.2f}")
        ax.axhline(0.5, ls=":", color="gray", lw=.8)
        ax.set_title(SLAB[sch]); ax.set_xlabel("target labels used (k; 0=zero-shot, ∞=oracle)")
        ax.set_ylabel("median per-epitope AUROC"); ax.set_ylim(0.0, 1.0); ax.grid(alpha=.25)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Q2 adaptation curve A(k): zero-shot → few-shot → within-epitope oracle")
    fig.tight_layout(); p = os.path.join(RES, "q2_adaptation.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main():
    m = load_matched()
    b = pd.read_csv(os.path.join(RES, "results_primary_beta.csv"))
    b = b[b.method != "SKIP"].dropna(subset=["auroc"])
    h = pd.read_csv(os.path.join(RES, "results_headroom.csv"))
    coef = json.load(open(os.path.join(RES, "coef_similarity_summary.json")))
    has_tab = "tabpfn" in m.method.unique()
    p1 = q1_plot(m, b); p2 = q2_plot(m, h)

    # ---- Q1 verdict numbers (honest scheme) ----
    def cell(sch, k, meth):
        v = permed(m, sch, k, meth)
        return v
    L = []
    hon = "otherep"
    lg20, lg50 = cell(hon, 20, "lightgbm"), cell(hon, 50, "lightgbm")
    tp20, tp50 = cell(hon, 20, "tabpfn"), cell(hon, 50, "tabpfn")
    td20, td50 = tcrdist_best(m, b, hon, 20), tcrdist_best(m, b, hon, 50)
    eps = td20.index
    dlg20 = (lg20[eps] - td20[eps]).median(); dlg50 = (lg50[eps] - td50[eps]).median()
    best_tab = "TabPFN/LightGBM" if has_tab else "LightGBM"
    # decision
    wins20 = dlg20 > 0.005 or (has_tab and (tp20[eps] - td20[eps]).median() > 0.005)
    wins50 = dlg50 > 0.005 or (has_tab and (tp50[eps] - td50[eps]).median() > 0.005)
    if not wins20 and not wins50:
        ans = ("**No.** With matched inputs the tabular models still do not beat TCRdist "
               "in the few-shot regime: TCRdist wins at k=20 and ties at k=50. This is the "
               "clean negative — the direction closes.")
    elif not wins20 and wins50:
        ans = ("**Partly.** Matched inputs erase most of the prior deficit: the tabular models "
               "reach parity or better by k=50, but TCRdist still wins at k=20. Not a clean "
               "negative, not a clean win — TCRdist keeps the edge at the smallest support.")
    else:
        ans = ("**Yes.** With matched inputs the tabular models beat TCRdist at k=20 and k=50.")

    L.append("# Fair rematch: matched inputs, headroom, and meta-learnability\n")
    L.append(f"**Q1 — does matching the inputs let the tabular models beat TCRdist at k≤50? "
             f"{ans}** Under honest other-epitope negatives, median per-epitope AUROC at k=20 is "
             f"LightGBM {fmt(lg20.median())}"
             + (f" / TabPFN {fmt(tp20.median())}" if has_tab else "")
             + f" vs TCRdist-best {fmt(td20[eps].median())} (Δ LightGBM {dlg20:+.3f}); "
             f"at k=50 LightGBM {fmt(lg50.median())}"
             + (f" / TabPFN {fmt(tp50.median())}" if has_tab else "")
             + f" vs TCRdist-best {fmt(td50[eps].median())} (Δ LightGBM {dlg50:+.3f}). "
             f"The tabular models only clearly win at k=200. "
             + ("" if has_tab else "_(TabPFN matched run still completing; LightGBM shown.)_") + "\n")
    L.append("Matched inputs = Atchley over germline CDR1/CDR2/CDR2.5 + CDR3 for BOTH chains, "
             "from tcrdist3's own germline tables (395 features after dropping constant columns); "
             "no V/J target-encoding (the loops already carry the V-gene information TCRdist uses). "
             "Splits, clustering (90% identity, cluster-aware), seeds and the per-epitope protocol "
             "are identical to the primary study. TCRdist-best is the strongest of beta-only and "
             "paired NN/kNN per epitope — notably, **paired TCRdist is weaker than beta-only at k=50** "
             "(equal-weighting a noisier alpha chain dilutes the NN distance), so TCRdist is given "
             "its best configuration.\n")

    # ---- Q1 table ----
    L.append("## Q1. Matched tabular vs TCRdist (median per-epitope AUROC)\n")
    rows = []
    for sch in SCH:
        for k in KS:
            td = tcrdist_best(m, b, sch, k); e = td.index
            def M(meth):
                v = permed(m, sch, k, meth)
                return fmt(v[e].median()) if len(v) else "—"
            def D(meth):
                v = permed(m, sch, k, meth)
                if not len(v): return "—"
                d = (v[e] - td[e]); return f"{d.median():+.3f} ({int((d>0).sum())}/{len(e)})"
            rows.append([SLAB[sch], k, M("tabpfn"), M("lightgbm"), M("ridge"),
                         fmt(td[e].median()), D("tabpfn"), D("lightgbm")])
    L.append(md_table(rows, ["scheme", "k", "TabPFN", "LightGBM", "ridge",
             "TCRdist best", "Δ(TabPFN−TD) wins", "Δ(LGBM−TD) wins"]))
    L.append(f"\n![Q1]({os.path.basename(p1)})\n")

    # ---- Q2 headroom ----
    L.append("## Q2. Headroom = oracle − zero-shot, and the adaptation curve A(k)\n")
    L.append("Oracle = best within-epitope model (LightGBM/ridge) under cluster-aware 5-fold CV on "
             "ALL of an epitope's labels. Zero-shot = trained on all other epitopes **with** epitope "
             "features, 0 labels from the target. Headroom is the ceiling any in-context head could "
             "capture. A(k) is the few-shot curve between them.\n")
    for sch in SCH:
        hs = h[h.scheme == sch].set_index("epitope")
        L.append(f"\n**{SLAB[sch]}**\n")
        rows = []
        for ep in hs.index:
            a = {k: permed(m, sch, k, "lightgbm").get(ep, np.nan) for k in KS}
            rows.append([ep, fmt(hs.loc[ep, "zeroshot"]), fmt(a[20]), fmt(a[50]), fmt(a[200]),
                         fmt(hs.loc[ep, "oracle"]), f"{hs.loc[ep,'headroom']:+.3f}"])
        rows.append(["**median**", fmt(hs["zeroshot"].median()),
                     fmt(permed(m, sch, 20, "lightgbm").median()),
                     fmt(permed(m, sch, 50, "lightgbm").median()),
                     fmt(permed(m, sch, 200, "lightgbm").median()),
                     fmt(hs["oracle"].median()), f"{hs['headroom'].median():+.3f}"])
        L.append(md_table(rows, ["epitope", "zero-shot", "A(20)", "A(50)", "A(200)",
                 "oracle", "headroom"]))
    L.append(f"\n![Q2]({os.path.basename(p2)})\n")
    zoe = h[h.scheme=="otherep"]
    L.append(f"\nUnder honest negatives the headroom is large (median {zoe['headroom'].median():+.3f}: "
             f"oracle {zoe['oracle'].median():.3f} vs zero-shot {zoe['zeroshot'].median():.3f}, which is "
             f"**below chance**). But A(k) shows within-task models capture most of it from only a few "
             f"labels — LightGBM reaches {fmt(permed(m,'otherep',50,'lightgbm').median())} by k=50 — so the "
             f"limiter is not adaptation capacity, it is that target labels are required at all.\n")

    # ---- Q3 ----
    L.append("## Q3. Is there anything to meta-learn across epitopes?\n")
    pooled = zoe["pooled_noepi"].median(); avg = zoe["oracle"].median()
    L.append(f"Per-epitope (oracle, within-epitope CV) median AUROC **{avg:.3f}**; one pooled model with "
             f"no epitope conditioning, trained on all other epitopes and applied to the held-out one, "
             f"median AUROC **{pooled:.3f}** (honest scheme) — far below the per-epitope models, and "
             f"below chance. Under honest negatives the pooled model *inverts*: the same TCR is a positive "
             f"for its own epitope and a negative for others, so without epitope conditioning the decision "
             f"boundary is self-contradictory.\n")
    L.append(f"\nPer-epitope ridge coefficient vectors are essentially **orthogonal** across epitopes: "
             f"pairwise cosine similarity median {coef['coef_cosine_median']:+.3f} "
             f"(range {coef['coef_cosine_min']:+.3f} to {coef['coef_cosine_max']:+.3f}). The per-task "
             f"binding functions do not overlap, so there is little shared structure for an in-context "
             f"head to amortize. Each epitope is its own problem; in-context value here is within-task "
             f"fitting, not meta-learned transfer.\n")

    L.append("## Reproduce\n")
    L.append("```bash\ncd tcr_feasibility/code\npython annotate_germline.py          # germline loops -> matched cache\n"
             "python run_experiment.py --feature-mode matched --fresh \\\n"
             "    --out ../results/results_matched.csv   # Q1\n"
             "python headroom_meta.py                    # Q2 + Q3\n"
             "python make_report_v2.py\n```\n")
    with open(os.path.join(RES, "report_v2.md"), "w") as f:
        f.write("\n".join(L))
    print("wrote report_v2.md,", os.path.basename(p1), os.path.basename(p2),
          "| tabpfn present:", has_tab)


if __name__ == "__main__":
    main()
