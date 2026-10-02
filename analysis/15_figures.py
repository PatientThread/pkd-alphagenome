#!/usr/bin/env python3
"""
15_figures.py

Manuscript figures. Static (print), so no hover layer; every figure ships a CSV
"table view" beside it so no value is available only as colour or position.

Design rules followed: categorical colour assigned by identity in fixed order and
validated for colour-vision deficiency (blue #2a78d6 / orange #eb6834, checked with
the dataviz validator: CVD dE 9.2+, normal-vision dE 27.6, both PASS); one axis per
panel, never two scales; thin marks; recessive axes; legend plus direct labels;
text in ink colours, never the series colour.

Figures
  fig1_scores_by_label      what each tool scores the tested-positive and
                            tested-negative variants (one panel per tool)
  fig2_blind_spots          measured effect versus AlphaGenome's prediction for the
                            five variants that act through translation or miRNA
  fig3_avi_attribution      what drives AVI at the drug-target sites and at coding
                            positives
  fig4_locus_percentile     (only if the gene-wide SpliceAI run has finished)
                            AVI versus SpliceAI on the same within-gene scale

Output: publication/figures/*.png, *.pdf and *_table.csv

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "publication" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

POS_C, NEG_C = "#2a78d6", "#eb6834"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK3,
                     "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.titlesize": 9, "axes.labelsize": 8, "figure.dpi": 200})


def save(fig, name: str, table: pd.DataFrame) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    table.to_csv(FIG / f"{name}_table.csv", index=False)
    plt.close(fig)
    print("wrote", name)


def strip(ax, vals_pos, vals_neg, title, ylabel, logy=False):
    rng = np.random.default_rng(7)
    for i, (v, c, lab) in enumerate(((vals_pos, POS_C, "tested positive"), (vals_neg, NEG_C, "tested negative"))):
        v = pd.to_numeric(pd.Series(v), errors="coerce").dropna().to_numpy()
        x = i + rng.uniform(-0.16, 0.16, len(v))
        ax.scatter(x, v, s=9, color=c, alpha=0.75, linewidths=0.4, edgecolors="white", zorder=3, label=lab)
        if len(v):
            ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=1.4, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["positive", "negative"])
    ax.set_xlim(-0.5, 1.5); ax.set_title(title, color=INK, loc="left")
    ax.set_ylabel(ylabel, color=INK2)
    if logy:
        ax.set_yscale("symlog", linthresh=0.01)
        ax.set_ylim(bottom=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#e8e7e3", lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def fig1(df: pd.DataFrame, S: dict) -> None:
    # AVI is a composite pathogenicity score, so per amendment A1 its panel uses the variant classes
    # where a high score cannot come from a protein consequence.
    panels = [("spliceai_max", "SpliceAI", "delta score", False, None),
              ("splice_site_usage_maxabs_all", "AlphaGenome splice-site usage", "max |score|", True, None),
              ("splice_junctions_maxabs_all", "AlphaGenome splice junctions", "max |score|", True, None),
              ("avi", "AVI (Atlas)\nnon-coding and synonymous only", "raw AVI", False, "nc_syn")]
    auc = S.get("A3_auc", {})
    fig, axes = plt.subplots(1, 4, figsize=(9.2, 2.9))
    rows = []
    for ax, (col, title, ylab, logy, subset) in zip(axes, panels):
        d = df[df.consequence_class.isin(["noncoding", "synonymous"])] if subset else df
        p = d[d.label == "POS"][col]; n = d[d.label == "NEG"][col]
        strip(ax, p, n, title, ylab, logy)
        a = (auc.get("composite_noncoding_synonymous", {}).get(col) if subset
             else auc.get("splicing_specific_all_variants", {}).get(col)) or {}
        if a.get("auc") is not None:
            ax.set_title(f"{title}\nAUC {a['auc']:.2f} ({a['ci95'][0]:.2f}-{a['ci95'][1]:.2f}), "
                         f"n={a['n_pos']} vs {a['n_neg']}", color=INK, loc="left", fontsize=8)
        for lab, v in (("POS", p), ("NEG", n)):
            for x in pd.to_numeric(v, errors="coerce").dropna():
                rows.append({"tool": title, "label": lab, "value": x})
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, ncol=2, loc="upper right",
               bbox_to_anchor=(1.0, 1.04), labelcolor=INK2)
    fig.suptitle("What each tool scores variants that were tested in the laboratory",
                 x=0.005, ha="left", color=INK, fontsize=10)
    fig.text(0.005, -0.06, "Each dot is one variant; black bar is the median. AVI is a composite "
             "pathogenicity score, so its panel is restricted to non-coding and synonymous variants, where a "
             "high score cannot come from a protein consequence.",
             color=INK3, fontsize=7, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, "fig1_scores_by_label", pd.DataFrame(rows))


def fig2(df: pd.DataFrame, TS: dict) -> None:
    # Only three of these five have a percentage in their source paper. The other two report a
    # direction only, so the measured effect is shown as text, never as an invented bar length.
    items = [("c.-69dupG\n5'UTR uORF extension", "87% less translation",
              float(df.loc[df.hgvs_c == "c.-69dupG", "rna_seq_maxabs_all"].iloc[0])),
             ("c.-87A>T\nuORF1 start loss", "9% more translation",
              abs(TS["Q4_uORF1_c.-87A>T"]["rna_seq"]["mean_all"])),
             ("c.-20A>T\nuORF2 start loss", "34% more translation",
              abs(TS["Q4_uORF2_c.-20A>T"]["rna_seq"]["mean_all"])),
             ("miR-17 seed\n6-base edit", "more PKD1 mRNA and protein (mouse; not given as %)",
              abs(TS["Q3_seed_6nt_edit"]["rna_seq"]["mean_all"])),
             ("rs3874648\nintron 30 modifier", "less full-length PKD1 (not given as %)",
              float(df.loc[df.hgvs_c == "c.10051-239G>A", "rna_seq_maxabs_all"].iloc[0]))]
    fig = plt.figure(figsize=(9.6, 3.0))
    ax = fig.add_axes((0.17, 0.24, 0.30, 0.60))
    y = np.arange(len(items))[::-1]
    ax.barh(y, [i[2] for i in items], color=NEG_C, height=0.5, zorder=3)
    for yy, it in zip(y, items):
        ax.text(it[2] + 0.0015, yy, f"{it[2]:.3f}", va="center", color=INK2, fontsize=7)
        ax.text(1.03, yy, it[1], va="center", color=INK2, fontsize=7.5, transform=ax.get_yaxis_transform(),
                clip_on=False)
    ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items], fontsize=7.5)
    ax.set_xlabel("AlphaGenome predicted change in PKD1 RNA (log2, absolute)", color=INK2)
    ax.set_xlim(0, 0.13)
    ax.text(1.03, len(items) - 0.35, "measured in the laboratory", transform=ax.get_yaxis_transform(),
            color=INK, fontsize=8, va="center", clip_on=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#e8e7e3", lw=0.5, zorder=0); ax.set_axisbelow(True)
    fig.text(0.02, 0.04, "For scale, a doubling of RNA would be 1.0 on this axis. Two of the five sources "
             "report a direction only, so their measured effect is given as text.",
             color=INK3, fontsize=7, ha="left")
    fig.text(0.02, 0.93, "Five variants that act after the message is made, and are invisible to the model",
             ha="left", color=INK, fontsize=10)
    save(fig, "fig2_blind_spots", pd.DataFrame(
        [{"variant": i[0].replace("\n", " "), "measured_effect": i[1],
          "alphagenome_rna_abs": i[2]} for i in items]))


def fig3(df: pd.DataFrame) -> None:
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
    save(fig, "fig3_avi_attribution", pd.DataFrame(sets).T.reset_index().rename(columns={"index": "variant_class"}))


def fig4(df: pd.DataFrame) -> None:
    if not ((RES / "spliceai_locus_PKD1.tsv.gz").exists() or (RES / "spliceai_locusbg_PKD1.tsv.gz").exists()):
        print("fig4 skipped: gene-wide SpliceAI not finished")
        return
    # Percentile compresses badly here: most of PKD1 scores near zero with SpliceAI, so small score
    # differences swing the percentile. Plot instead the quantity a laboratory acts on, the share of
    # the gene that outranks the variant, on a log scale. Lower and further left is better.
    sub = df[(df.gene == "PKD1") & df.label.isin(["POS", "NEG"]) & df.is_snv
             & df.consequence_class.isin(["noncoding", "synonymous"])
             & df.avi_locus_pct.notna() & df.spliceai_locus_pct.notna()].copy()
    floor = 0.001  # 1 in 100,000 of the locus; below this the sampled background cannot resolve
    sub["x"] = (100 - sub.spliceai_locus_pct).clip(lower=floor)
    sub["y"] = (100 - sub.avi_locus_pct).clip(lower=floor)
    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    for lab, c, nm in (("POS", POS_C, "tested positive"), ("NEG", NEG_C, "tested negative")):
        s = sub[sub.label == lab]
        ax.scatter(s.x, s.y, s=26, color=c, alpha=0.8, edgecolors="white", linewidths=0.5,
                   label=nm, zorder=3)
    lim = (floor * 0.7, 120)
    ax.plot(lim, lim, color=INK3, lw=0.6, ls=(0, (4, 3)), zorder=1)
    for v in (1, 0.1):
        ax.axvline(v, color="#e8e7e3", lw=0.8, zorder=0); ax.axhline(v, color="#e8e7e3", lw=0.8, zorder=0)
    ax.annotate("top 1%", xy=(1, 0.0016), color=INK3, fontsize=6.5, ha="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.invert_xaxis(); ax.invert_yaxis()
    ax.set_xlabel("percent of PKD1 scoring higher, SpliceAI", color=INK2)
    ax.set_ylabel("percent of PKD1 scoring higher, AVI", color=INK2)
    ax.set_title("How far down the gene-wide list each variant sits\n(better is towards the top right)",
                 color=INK, loc="left", fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc="upper left", labelcolor=INK2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_axisbelow(True)
    fig.text(0.02, -0.02, "SpliceAI background is an evenly spaced sample of the locus (amendment A5); "
             "values are clipped at 0.001%.", color=INK3, fontsize=6.5, ha="left")
    save(fig, "fig4_locus_rank",
         sub[["gene", "hgvs_c", "label", "spliceai_locus_pct", "avi_locus_pct", "x", "y"]]
         .rename(columns={"x": "pct_of_gene_above_spliceai", "y": "pct_of_gene_above_avi"}))


def main() -> None:
    df = pd.read_csv(RES / "per_variant.tsv", sep="\t")
    S = json.loads((RES / "summary.json").read_text())
    TS = json.loads((RES / "therapeutic_sites.json").read_text())
    fig1(df, S); fig2(df, TS); fig3(df); fig4(df)


if __name__ == "__main__":
    main()
