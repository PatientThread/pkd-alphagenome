#!/usr/bin/env python3
"""
19_figures_ejhg.py

The three figures for the EJHG submission, rebuilt after review (amendment A7).
Figure numbering follows order of appearance. The former Figure 2 is replaced by
Table 4, which states the assay endpoint for each case rather than implying that
a predicted RNA change and a measured protein output are comparable.

  Figure 1  discrimination within the splicing class
  Figure 2  what drives the AVI score, by variant class
  Figure 3  rank within the locus, SpliceAI against AVI

Colour is assigned by identity and validated for colour-vision deficiency; every
figure also ships a CSV of its values. Outputs PNG (review), PDF and EPS.

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "publication" / "figures_ejhg"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "analysis"))
from importlib import import_module  # noqa: E402
M = import_module("18_stratified_analysis")

POS_C, NEG_C = "#2a78d6", "#eb6834"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK3,
                     "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.titlesize": 9, "axes.labelsize": 8, "figure.dpi": 300})


def save(fig, name, table):
    for ext in ("png", "pdf", "eps"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    table.to_csv(FIG / f"{name}_values.csv", index=False)
    plt.close(fig)
    print("wrote", name)


def load():
    d = pd.read_csv(RES / "per_variant.tsv", sep="\t")
    d["assay_class"] = [M.assay_class(str(a), str(m)) for a, m in zip(d.assay, d.mechanism)]
    return d


def fig1(d, strat):
    sp = d[(d.assay_class == "splicing") & d.label.isin(["POS", "NEG"])]
    ncs = sp[sp.consequence_class.isin(["noncoding", "synonymous"])]
    panels = [("spliceai_max", "SpliceAI", "delta score", sp, "SpliceAI"),
              (M.SPLICE_COL, "AlphaGenome splicing composite", "composite score", sp,
               "AlphaGenome splicing composite"),
              ("avi", "AVI, non-coding and synonymous only", "raw AVI", ncs, None)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
    rows = []
    rng = np.random.default_rng(7)
    for ax, (col, title, ylab, data, key) in zip(axes, panels):
        for i, (lab, c, nm) in enumerate((("POS", POS_C, "splicing shown"), ("NEG", NEG_C, "no effect detected"))):
            v = pd.to_numeric(data[data.label == lab][col], errors="coerce").dropna().to_numpy()
            ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v, s=9, color=c, alpha=0.75,
                       linewidths=0.4, edgecolors="white", zorder=3, label=nm)
            if len(v):
                ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=1.4, zorder=4)
            rows += [{"panel": title, "label": lab, "value": float(x)} for x in v]
        if key:
            a = strat["primary_splicing_class"]["auc"][key]
            sub = f"AUC {a['auc']:.2f} ({a['ci95'][0]:.2f}-{a['ci95'][1]:.2f})"
        else:
            a = json.load(open(RES / "stratified.json"))["composite_scores_noncoding_synonymous_splicing_class"]["AVI"]
            sub = f"AUC {a['auc']:.2f} ({a['ci95'][0]:.2f}-{a['ci95'][1]:.2f})"
        ax.set_title(f"{title}\n{sub}", color=INK, loc="left", fontsize=8)
        ax.set_ylabel(ylab, color=INK2)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["positive", "negative"])
        ax.set_xlim(-0.5, 1.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color="#e8e7e3", lw=0.5, zorder=0); ax.set_axisbelow(True)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, ncol=2, loc="upper right", bbox_to_anchor=(1.0, 1.06),
               labelcolor=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "figure1_discrimination", pd.DataFrame(rows))


def fig2(df: pd.DataFrame) -> None:
    att = pd.read_csv(RES / "therapeutic_sites_attribution.tsv", sep="\t")
    groups = ["splicing", "expression_chromatin", "protein", "conservation"]
    # validated set (blue, aqua, yellow, violet): CVD dE 9.1, normal-vision dE 22.9, both PASS.
    # Orange is not used here: beside yellow it fails the normal-vision floor in a stacked bar.
    cols = {"splicing": POS_C, "expression_chromatin": "#1baf7a", "protein": "#eda100", "conservation": "#4a3aa7"}
    sets = {"miR-17 seed site\n(18 possible changes)": att[att.site == "miR17_seed"][groups].mean(),
            "uORF start codons\n(18 possible changes)": att[att.site == "uORF_ATG"][groups].mean()}
    cod = df[(df.label == "POS") & df.consequence_class.isin(["missense", "nonsense"]) & df.avi.notna()]
    nc = df[(df.label == "POS") & (df.consequence_class == "noncoding") & df.avi.notna()]
    for name, sub in (("Non-coding positives\n(splice variants)", nc), ("Coding positives\n(minigene studies)", cod)):
        sets[name] = pd.Series({g: sub[f"avi_attr_{g}"].mean() for g in groups})
    fig, ax = plt.subplots(figsize=(7.4, 2.6))
    y = np.arange(len(sets))[::-1]
    left = np.zeros(len(sets))
    for g in groups:
        vals = np.array([sets[k][g] for k in sets])
        ax.barh(y, vals, left=left, height=0.55, color=cols[g], label=g.replace("_", " and "),
                edgecolor="white", linewidth=1.2, zorder=3)
        for yy, v, l in zip(y, vals, left):
            if v > 0.08:
                ax.text(l + v / 2, yy, f"{v:.2f}", ha="center", va="center", color="white", fontsize=7)
        left += vals
    ax.set_yticks(y); ax.set_yticklabels(list(sets), fontsize=7)
    ax.set_xlabel("mean contribution to the AVI score", color=INK2)
    ax.set_title("What drives AVI, by variant class", color=INK, loc="left")
    ax.legend(frameon=False, fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.28), labelcolor=INK2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#e8e7e3", lw=0.5, zorder=0); ax.set_axisbelow(True)
    save(fig, "figure2_avi_attribution", pd.DataFrame(sets).T.reset_index().rename(columns={"index": "variant_class"}))


def fig3(df: pd.DataFrame) -> None:
    if not ((RES / "spliceai_locus_PKD1.tsv.gz").exists() or (RES / "spliceai_locusbg_PKD1.tsv.gz").exists()):
        print("fig4 skipped: gene-wide SpliceAI not finished")
        return
    # Percentile compresses badly here: most of PKD1 scores near zero with SpliceAI, so small score
    # differences swing the percentile. Plot instead the quantity a laboratory acts on, the share of
    # the gene that outranks the variant, on a log scale. Lower and further left is better.
    # Splicing evidence class only: a translation-reporter variant has no place in a ranking
    # against a background of splice predictions (amendment A8.5).
    sub = df[(df.gene == "PKD1") & df.label.isin(["POS", "NEG"]) & df.is_snv
             & (df.assay_class == "splicing")
             & df.consequence_class.isin(["noncoding", "synonymous"])
             & df.avi_locus_pct.notna() & df.spliceai_locus_pct.notna()].copy()
    # Each panel is censored at the resolution of ITS OWN background: one allele in 39,144 for the
    # whole locus (0.0026%) and one in 28,401 for the non-coding background (0.0035%). A shared floor
    # would imply a resolution one of the two backgrounds does not have.
    sbg = pd.read_csv(RES / "spliceai_locusbg_PKD1.tsv.gz", sep="\t")
    floor = round(100 / len(sbg), 4)
    sub["x"] = (100 - sub.spliceai_locus_pct).clip(lower=floor)
    sub["y"] = (100 - sub.avi_locus_pct).clip(lower=floor)

    # Panel b: the same variants ranked against a NON-CODING background, matched on the proportion
    # of positions flagged. This is the comparison that is fair to a whole-genome pathogenicity
    # score asked a non-coding question, and it reverses the apparent result of panel a.
    ncb_path = RES / "noncoding_background_ranks.tsv"
    ncb = pd.read_csv(ncb_path, sep="\t") if ncb_path.exists() else pd.DataFrame()
    floor_b = floor
    if len(ncb):
        ncb = ncb[ncb.gene == "PKD1"].copy()
        rm = json.loads((RES / "ranking_manifest.json").read_text())
        n_b = rm["backgrounds"]["PKD1"]["non-coding positions only"]["spliceai"]["alleles_scored"]
        floor_b = round(100 / n_b, 4)
        ncb["x"] = (100 - ncb.spliceai_pct_nc_bg).clip(lower=floor_b)
        ncb["y"] = (100 - ncb.avi_pct_nc_bg).clip(lower=floor_b)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
    panels = [(axes[0], sub, "a", "all possible PKD1 variants", floor),
              (axes[1], ncb, "b", "non-coding PKD1 positions only", floor_b)]
    lim = (min(floor, floor_b) * 0.7, 120)
    for ax, data, tag, bg, fl in panels:
        if not len(data):
            ax.set_visible(False)
            continue
        for lab, c, nm in (("POS", POS_C, "tested positive"), ("NEG", NEG_C, "tested negative")):
            d_ = data[data.label == lab]
            ax.scatter(d_.x, d_.y, s=26, color=c, alpha=0.8, edgecolors="white", linewidths=0.5,
                       label=nm, zorder=3)
        ax.plot(lim, lim, color=INK3, lw=0.6, ls=(0, (4, 3)), zorder=1)
        for v in (1, 0.1):
            ax.axvline(v, color="#e8e7e3", lw=0.8, zorder=0)
            ax.axhline(v, color="#e8e7e3", lw=0.8, zorder=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(*lim); ax.set_ylim(*lim)
        ax.invert_xaxis(); ax.invert_yaxis()
        ax.set_xlabel("percent scoring higher, SpliceAI", color=INK2, fontsize=8)
        ax.set_ylabel("percent scoring higher, AVI", color=INK2, fontsize=8)
        ax.set_title(f"({tag}) background: {bg}", color=INK, loc="left", fontsize=8.5)
        ax.annotate(f"censored at {fl}%", xy=(0.03, 0.03), xycoords="axes fraction",
                    color=INK3, fontsize=6)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=7)
    axes[0].legend(frameon=False, fontsize=7, loc="upper left", labelcolor=INK2)
    fig.suptitle("How far down the list each variant sits (better is towards the top right)",
                 color=INK, x=0.02, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.text(0.02, 0.005, "SpliceAI background is an evenly spaced sample of the locus (amendment "
             "A5). Each panel is censored at the resolution of its own background, one allele in "
             f"{len(sbg):,} ({floor}%) in (a) and one in {n_b:,} ({floor_b}%) in (b). Both panels "
             "match the two tools on the PROPORTION of their own background flagged, never the raw "
             "count; the unit throughout is the alternate allele.",
             color=INK3, fontsize=6.5, ha="left")
    save(fig, "figure3_locus_rank",
         sub[["gene", "hgvs_c", "label", "spliceai_locus_pct", "avi_locus_pct", "x", "y"]]
         .rename(columns={"x": "pct_of_gene_above_spliceai", "y": "pct_of_gene_above_avi"}))




def main() -> None:
    d = load()
    strat = json.load(open(RES / "stratified.json"))
    fig1(d, strat)
    fig2(d)
    fig3(d)


if __name__ == "__main__":
    main()
