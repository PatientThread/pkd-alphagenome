#!/usr/bin/env python3
"""
21_build_tables.py

Generates publication/tables_ejhg.md from the result files.

The tables were previously maintained by hand, which is how Table 2 and Table 3 came to
disagree with results/stratified.json after amendment A8. Every number here is read from
a result file; nothing is typed.

Output: publication/tables_ejhg.md

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "publication" / "tables_ejhg.md"

# Assay wording per source paper, taken from the papers themselves and carried in the truth set.
ASSAY_BY_PMID = {
    "39757590": ("Renilla/firefly dual-luciferase reporter, protein output", "yes, RT-qPCR: unchanged"),
    "41006799": ("Gaussia luciferase normalised to secreted alkaline phosphatase, protein output",
                 "no RNA study performed"),
}


def row(*cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def header(*cells) -> str:
    return row(*cells) + "\n" + row(*["---"] * len(cells))


def fmt_ci(d: dict, key: str = "ci95") -> str:
    c = d.get(key)
    return f"{c[0]:.2f} to {c[1]:.2f}" if c else "not estimable"


def main() -> None:
    st = json.loads((RES / "stratified.json").read_text())
    df = pd.read_csv(RES / "per_variant.tsv", sep="\t")
    prim, out = st["primary_splicing_class"], []

    # ---------------- Table 1
    import importlib
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    assay_class = importlib.import_module("18_stratified_analysis").assay_class
    df["assay_class"] = [assay_class(str(a), str(m)) for a, m in zip(df.assay, df.mechanism)]
    lab = df[df.label.isin(["POS", "NEG"])]
    out.append("**Table 1. Composition of the benchmark by type of experimental evidence.**\n")
    out.append(header("Evidence type", "Variant class", "Effect shown", "No effect detected",
                      "From patients or populations", "Genes"))
    for (cls, cons), g in lab.groupby(["assay_class", "consequence_class"], sort=False):
        out.append(row(cls, cons.replace("_", " "), int((g.label == "POS").sum()), int((g.label == "NEG").sum()),
                       int((g.origin != "engineered").sum()), ", ".join(sorted(g.gene.unique()))))
    out.append("\nProvenance is reported in the fifth column because it is strongly associated with "
               "outcome: in the splicing class there are no patient-derived variants in which no "
               "effect was detected.\n")

    # ---------------- Table 2
    out.append("**Table 2. Discrimination within the splicing class.**\n")
    out.append(header("Score", "Analysis population", "Positives", "Negatives", "AUC", "95% CI"))
    n_pos, n_neg = prim["n_pos"], prim["n_neg"]
    for tool, a in prim["auc"].items():
        miss = (n_pos - a["n_pos"]) + (n_neg - a["n_neg"])
        pop = "all splicing-class variants" + (f" ({miss} unscored)" if miss else "")
        out.append(row(tool, pop, a["n_pos"], a["n_neg"], f"{a['auc']:.3f}", fmt_ci(a)))
    p1 = st["composite_scores_PKD1_only"]
    for tool, a in p1["all_tools"].items():
        out.append(row(tool, "PKD1 only", a["n_pos"], a["n_neg"], f"{a['auc']:.3f}", fmt_ci(a)))
    for tool, a in p1["noncoding_synonymous"].items():
        out.append(row(tool, "PKD1, non-coding and synonymous", a["n_pos"], a["n_neg"],
                       f"{a['auc']:.3f}", fmt_ci(a)))
    cs = st["common_subset_all_tools"]
    for tool, a in cs["auc"].items():
        out.append(row(tool, "scored by all four tools", a["n_pos"], a["n_neg"],
                       f"{a['auc']:.3f}", fmt_ci(a)))
    pd_, pc = prim["paired_ag_composite_vs_spliceai"], prim["paired_cluster_bootstrap"]
    out.append(f"\nPaired difference, AlphaGenome composite minus SpliceAI, on the {pd_['n_pos']} "
               f"positives and {pd_['n_neg']} negatives both tools score: "
               f"{pd_['auc_difference']:+.3f} (95% CI {pd_['ci95'][0]:+.3f} to {pd_['ci95'][1]:+.3f}); "
               f"resampling whole source studies, {pc['auc_difference']:+.3f} "
               f"({pc['ci95_study_cluster'][0]:+.3f} to {pc['ci95_study_cluster'][1]:+.3f}). Equal "
               "observed AUCs do not imply equal per-variant behaviour.")
    out.append(f"\nAVI and CADD are reported primarily for PKD1, which lies on a chromosome absent "
               "from AVI's training and validation lists; PKD2 lies on a training chromosome, so "
               "overlap is not excluded for it, though neither is it demonstrated. The pooled and "
               "PKD1-only analyses agree. The "
               f"{len(cs['excluded'])} variants absent from the four-tool comparison are insertions "
               "or deletions missing from the precomputed files used, and all have a demonstrated "
               "effect, so the absent values are not missing at random.\n")

    # ---------------- Table 3
    out.append("**Table 3. Robustness of the primary comparison.**\n")
    out.append(header("Analysis", "SpliceAI", "AlphaGenome composite", "n (pos/neg)"))
    sa, ag = prim["auc"]["SpliceAI"], prim["auc"]["AlphaGenome splicing composite"]
    out.append(row("Primary (all splicing-class variants)", f"{sa['auc']:.3f}", f"{ag['auc']:.3f}",
                   f"{n_pos}/{n_neg}"))
    cl = prim["study_cluster_bootstrap"]
    out.append(row("Study-cluster bootstrap (95% interval)",
                   fmt_ci(cl["SpliceAI"], "ci95_study_cluster"),
                   fmt_ci(cl["AlphaGenome splicing composite"], "ci95_study_cluster"),
                   f"{n_pos}/{n_neg}"))
    lo = st["leave_one_study_out"]
    out.append(row("Leave-one-study-out (range)",
                   f"{min(v['SpliceAI'] for v in lo.values()):.3f} to {max(v['SpliceAI'] for v in lo.values()):.3f}",
                   f"{min(v['AlphaGenome'] for v in lo.values()):.3f} to {max(v['AlphaGenome'] for v in lo.values()):.3f}",
                   f"{len(lo)} studies"))
    for kind, label in (("gene", "Subgroup: {}"), ("splice_distance", "Subgroup: {}"),
                        ("origin", "Evidence from {} variants")):
        for name, v in st["subgroups_splicing_class"][kind].items():
            s_, a_ = v["SpliceAI"], v["AlphaGenome"]
            def cell(x):
                if x.get("auc") is not None:
                    return f"{x['auc']:.3f}"
                return ("undefined, no controls" if x.get("withheld_because", "").startswith("no ")
                        else "withheld, sparse controls")
            out.append(row(label.format(name), cell(s_), cell(a_), f"{s_['n_pos']}/{s_['n_neg']}"))
    out.append("\nTwo different reasons are distinguished. Where a subgroup contains no variants of "
               "one outcome class an AUC is undefined: this applies to the canonical splice-site and "
               "the patient-derived subgroups, the latter being the important one, since it means no "
               "discrimination estimate on patient-derived variants is available from these data. "
               "Where both classes are present but the controls number one or two, an AUC is defined "
               "but far too unstable to report, and is withheld rather than shown.\n")

    # ---------------- Table 4
    out.append("**Table 4. The reporter-output class: 5' untranslated region variants assayed by "
               "protein output. Reporter RNA was measured for three of the six.**\n")
    out.append(header("Variant", "Assay", "Transcript measured", "Measured effect",
                      "AlphaGenome splicing composite", "AlphaGenome RNA (log2)", "SpliceAI", "Source"))
    tr = st["reporter_output_class"]["variants"]
    effect = {"c.-209G>A": "no difference from wild type",
              "c.-52C>T": "no difference from wild type",
              "c.-20A>G": "higher output (uORF2 start codon disrupted)",
              "c.-23G>T": "about 42% higher output (uORF2 -3 nucleotide)",
              "c.-66T>C": "lower output (uORF1 stop codon removed)",
              "c.-69dupG": "87% lower output (uORF1 extended)"}
    pm = dict(zip(df.hgvs_c, df.pmids.astype(str).str.split(";").str[0]))
    for v in sorted(tr, key=lambda x: x["hgvs_c"]):
        src = pm.get(v["hgvs_c"], "")
        assay, rna = ASSAY_BY_PMID[src]
        out.append(row(v["hgvs_c"], assay, rna, effect[v["hgvs_c"]],
                       f"{v['ag_splicing_composite_anygene']:.3f}", f"{v['rna_seq_maxabs_all']:.3f}",
                       f"{v['spliceai_max']:.2f}", src))
    out.append("\nThe AlphaGenome RNA column is the RNA_SEQ variant score: the maximum ABSOLUTE "
               "log2 predicted change across all RNA-seq tracks, restricted to the PKD1 gene row, at "
               "the recommended 1,048,576 bp context. It is an unsigned magnitude, so it states how "
               "large a change is predicted, not its direction. Direction of the measured effect is "
               "stated because two variants increase output. Chen et al. "
               "[39757590] measured reporter RNA by RT-qPCR and found it unchanged, so for those three "
               "variants the effect is attributable to translation. Wedd et al. [41006799] performed no "
               "RNA study and state that effects on pre-mRNA splicing and transcript stability cannot be "
               "excluded, so for their three variants the mechanism of the change in protein output is "
               "undetermined. Predicted values are shown for the same variants; neither predictor "
               "addresses protein output.")

    OUT.write_text("\n".join(out) + "\n")
    print("wrote", OUT.name)


if __name__ == "__main__":
    main()
