#!/usr/bin/env python3
"""
18_stratified_analysis.py

Rebuild of the primary analysis after external review (amendment A7).

The reviewers' central objection was that one reference standard was built from
three different kinds of experiment and then used for a single splicing AUC. A
variant shown to reduce translation in a reporter is not evidence about splicing,
and a normal translation reporter is not a splicing-negative control.

This script therefore assigns every labelled variant to the assay class that
produced its evidence, and analyses each class separately:

  splicing    RT-PCR, cDNA sequencing, RNA-seq or minigene read out as splicing
  abundance   transcript-level measurement (allele-specific expression)
  reporter output  5'UTR luciferase reporter (protein output; RNA measured in one study of two)

PRIMARY: splicing class only, AlphaGenome's documented splicing composite
    max(splice_sites) + max(splice_site_usage) + max(splice_junctions)/5
(aggregated across all tracks and genes, per alphagenomedocs.com) against
SpliceAI, on identical variants, as a PAIRED difference with a bootstrap
interval. Composite scores (AVI, CADD) are reported on the non-coding and
synonymous subset only, as before.

Also produced, because the negatives are concentrated in a few papers:
  leave-one-study-out AUCs, a study-cluster bootstrap, and subgroup results by
  gene, by distance from the splice site, and by patient RNA versus minigene.

Coverage is reported per tool, so a variant a tool cannot score is counted as
unscored and never as a low score.

Output: results/stratified.json, results/coverage_matrix.csv

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
SPLICE_COL = "ag_splicing_composite_anygene"


def assay_class(assay: str, mechanism: str) -> str:
    """Class of experiment that produced the evidence.

    Corrected 21 September 2026 after review: rs3874648 was previously placed in the abundance
    class because its mechanism string begins "expression modifier", but Zhang et al. demonstrated
    an atypical splice form by RT-PCR of patient leukocytes and activation of a cryptic acceptor.
    Splicing evidence therefore takes precedence over the downstream abundance consequence.
    """
    a, m = assay.lower(), mechanism.lower()
    if "luciferase" in a:
        return "reporter output"
    splice_words = ("splice", "splicing", "exon", "intron", "pseudoexon", "cryptic", "branchpoint")
    if any(w in m for w in splice_words) or any(w in a for w in ("minigene", "rt-pcr", "rna-seq", "cdna")):
        return "splicing"
    if "allele-specific expression" in a or "expression modifier" in m:
        return "abundance"
    return "splicing"


def auc(x: np.ndarray, y: np.ndarray) -> float:
    return float(((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum()) / (len(x) * len(y)))


def auc_ci(p: np.ndarray, n: np.ndarray, n_boot=2000, seed=20260920) -> dict:
    if len(p) < 3 or len(n) < 3:
        # An AUC is undefined without both outcome classes, but merely unstable when both are
        # present and one is tiny. The reviewer asked that these be distinguished, and they are
        # different statements: one is arithmetic, the other is a judgement about precision.
        reason = ("no variants of one outcome class" if len(p) == 0 or len(n) == 0
                  else "not reported because of sparse controls")
        return {"n_pos": int(len(p)), "n_neg": int(len(n)), "auc": None, "withheld_because": reason}
    rng = np.random.default_rng(seed)
    b = [auc(rng.choice(p, len(p)), rng.choice(n, len(n))) for _ in range(n_boot)]
    return {"n_pos": int(len(p)), "n_neg": int(len(n)), "auc": round(auc(p, n), 3),
            "ci95": [round(float(np.percentile(b, 2.5)), 3), round(float(np.percentile(b, 97.5)), 3)]}


def cluster_auc_ci(df: pd.DataFrame, col: str, n_boot=2000, seed=20260920) -> dict:
    """Bootstrap resampling whole SOURCE STUDIES, not variants, so that concentration of the
    controls in a few papers is reflected in the interval."""
    d = df[[col, "label", "first_pmid"]].dropna()
    studies = d.first_pmid.unique()
    rng = np.random.default_rng(seed)
    obs_p = d[d.label == "POS"][col].to_numpy(float)
    obs_n = d[d.label == "NEG"][col].to_numpy(float)
    if len(obs_p) < 3 or len(obs_n) < 3:
        return {"auc": None}
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(studies, len(studies), replace=True)
        s = pd.concat([d[d.first_pmid == q] for q in pick])
        p, n = s[s.label == "POS"][col].to_numpy(float), s[s.label == "NEG"][col].to_numpy(float)
        if len(p) and len(n):
            boots.append(auc(p, n))
    return {"auc": round(auc(obs_p, obs_n), 3),
            "ci95_study_cluster": [round(float(np.percentile(boots, 2.5)), 3),
                                   round(float(np.percentile(boots, 97.5)), 3)],
            "n_studies": int(len(studies))}


def paired_diff(df: pd.DataFrame, a: str, b: str, n_boot=2000, seed=20260920) -> dict:
    d = df[[a, b, "label"]].dropna()
    ap, an = d[d.label == "POS"][a].to_numpy(float), d[d.label == "NEG"][a].to_numpy(float)
    bp, bn = d[d.label == "POS"][b].to_numpy(float), d[d.label == "NEG"][b].to_numpy(float)
    if len(ap) < 3 or len(an) < 3:
        return {"auc_difference": None}
    rng = np.random.default_rng(seed)
    obs = auc(ap, an) - auc(bp, bn)
    boots = []
    for _ in range(n_boot):
        ip, iq = rng.integers(0, len(ap), len(ap)), rng.integers(0, len(an), len(an))
        boots.append(auc(ap[ip], an[iq]) - auc(bp[ip], bn[iq]))
    return {"n_pos": int(len(ap)), "n_neg": int(len(an)),
            "auc_a": round(auc(ap, an), 3), "auc_b": round(auc(bp, bn), 3),
            "auc_difference": round(obs, 3),
            "ci95": [round(float(np.percentile(boots, 2.5)), 3), round(float(np.percentile(boots, 97.5)), 3)]}


def splice_distance(h: str) -> str:
    """Canonical (+/-1,2), near-splice (within 20 nt of a junction), or other/exonic."""
    m = re.search(r"[+-](\d+)", h.split("c.")[-1])
    if m:
        d = int(m.group(1))
        return "canonical" if d <= 2 else ("near-splice (3-20 nt)" if d <= 20 else "deep intronic (>20 nt)")
    return "exonic or UTR"


def main() -> None:
    df = pd.read_csv(RES / "per_variant.tsv", sep="\t")
    df["assay_class"] = [assay_class(str(a), str(m)) for a, m in zip(df.assay, df.mechanism)]
    df["first_pmid"] = df.pmids.astype(str).str.split(";").str[0]
    df["splice_distance"] = df.hgvs_c.map(splice_distance)
    lab = df[df.label.isin(["POS", "NEG"])].copy()

    out: dict = {"assay_class_counts": lab.groupby(["assay_class", "label"]).size().unstack(fill_value=0).to_dict()}

    # coverage per tool, per assay class
    tools = {"SpliceAI": "spliceai_max", "AlphaGenome splicing composite": SPLICE_COL,
             "AVI": "avi", "CADD": "cadd_raw"}
    cov = []
    for cls, g in lab.groupby("assay_class"):
        for tool, col in tools.items():
            cov.append({"assay class": cls, "tool": tool, "variants": len(g),
                        "scored": int(g[col].notna().sum()), "unscored": int(g[col].isna().sum())})
    pd.DataFrame(cov).to_csv(RES / "coverage_matrix.csv", index=False)
    out["coverage"] = cov

    # PRIMARY: splicing class
    sp = lab[lab.assay_class == "splicing"]
    out["primary_splicing_class"] = {
        "definition": "variants whose evidence is a splicing assay (patient RNA or minigene)",
        "n_pos": int((sp.label == "POS").sum()), "n_neg": int((sp.label == "NEG").sum()),
        "auc": {t: auc_ci(sp[sp.label == "POS"][c].dropna().to_numpy(float),
                          sp[sp.label == "NEG"][c].dropna().to_numpy(float))
                for t, c in tools.items()},
        "paired_ag_composite_vs_spliceai": paired_diff(sp, SPLICE_COL, "spliceai_max"),
        "study_cluster_bootstrap": {t: cluster_auc_ci(sp, c) for t, c in
                                    (("SpliceAI", "spliceai_max"), ("AlphaGenome splicing composite", SPLICE_COL))}}

    # leave one study out
    loo = {}
    for pmid in sorted(sp.first_pmid.unique()):
        s = sp[sp.first_pmid != pmid]
        if (s.label == "POS").sum() >= 3 and (s.label == "NEG").sum() >= 3:
            loo[pmid] = {t: auc_ci(s[s.label == "POS"][c].dropna().to_numpy(float),
                                   s[s.label == "NEG"][c].dropna().to_numpy(float))["auc"]
                         for t, c in (("SpliceAI", "spliceai_max"), ("AlphaGenome", SPLICE_COL))}
    out["leave_one_study_out"] = loo

    # subgroups
    subs = {}
    for key, col in (("gene", "gene"), ("splice_distance", "splice_distance"), ("origin", "origin")):
        subs[key] = {}
        for v, g in sp.groupby(col):
            subs[key][str(v)] = {t: auc_ci(g[g.label == "POS"][c].dropna().to_numpy(float),
                                           g[g.label == "NEG"][c].dropna().to_numpy(float))
                                 for t, c in (("SpliceAI", "spliceai_max"), ("AlphaGenome", SPLICE_COL))}
    out["subgroups_splicing_class"] = subs

    # composite scores on non-coding and synonymous splicing-class variants
    ncs = sp[sp.consequence_class.isin(["noncoding", "synonymous"])]
    out["composite_scores_noncoding_synonymous_splicing_class"] = {
        t: auc_ci(ncs[ncs.label == "POS"][c].dropna().to_numpy(float),
                  ncs[ncs.label == "NEG"][c].dropna().to_numpy(float))
        for t, c in (("AVI", "avi"), ("CADD", "cadd_raw"))}

    # other classes, described not tested
    for cls in ("abundance", "reporter output"):
        g = lab[lab.assay_class == cls]
        out[f"{cls.replace(chr(32), chr(95))}_class"] = {
            "n": int(len(g)),
            "variants": g[["gene", "hgvs_c", "label", "mechanism", "spliceai_max", SPLICE_COL,
                           "splice_sites_maxabs_anygene", "rna_seq_maxabs_all", "avi_locus_pct"]]
            .round(3).to_dict("records")}

    # Common analysis set: every tool scored (reviewer request). Score availability is
    # outcome-dependent here, so the restricted comparison is reported beside the full one.
    cols = [c for c in tools.values()]
    common = sp.dropna(subset=cols)
    out["common_subset_all_tools"] = {
        "n_pos": int((common.label == "POS").sum()), "n_neg": int((common.label == "NEG").sum()),
        "excluded": sp[sp[cols].isna().any(axis=1)][["gene", "hgvs_c", "label"]].to_dict("records"),
        "auc": {t: auc_ci(common[common.label == "POS"][c].dropna().to_numpy(float),
                          common[common.label == "NEG"][c].dropna().to_numpy(float))
                for t, c in tools.items()},
        "paired_ag_vs_spliceai": paired_diff(common, SPLICE_COL, "spliceai_max")}

    # Paired difference under the SAME study draws for both tools.
    def cluster_paired(dfi, a, b, n_boot=2000, seed=20260921):
        d = dfi[[a, b, "label", "first_pmid"]].dropna()
        studies = d.first_pmid.unique()
        rng = np.random.default_rng(seed)
        def diff(s):
            p, n = s[s.label == "POS"], s[s.label == "NEG"]
            if len(p) < 1 or len(n) < 1:
                return None
            return auc(p[a].to_numpy(float), n[a].to_numpy(float)) - auc(p[b].to_numpy(float), n[b].to_numpy(float))
        obs = diff(d)
        boots = [x for x in (diff(pd.concat([d[d.first_pmid == q] for q in rng.choice(studies, len(studies), True)]))
                             for _ in range(n_boot)) if x is not None]
        return {"auc_difference": round(obs, 3),
                "ci95_study_cluster": [round(float(np.percentile(boots, 2.5)), 3),
                                       round(float(np.percentile(boots, 97.5)), 3)],
                "n_studies": int(len(studies)), "n_resamples": len(boots), "seed": 20260921}
    out["primary_splicing_class"]["paired_cluster_bootstrap"] = cluster_paired(sp, SPLICE_COL, "spliceai_max")

    # Threshold performance ON THE SPLICING SET, with full confusion matrices.
    def wilson(k, n):
        if n == 0:
            return [None, None]
        z, ph = 1.959964, k / n
        d = 1 + z * z / n
        c = (ph + z * z / (2 * n)) / d
        h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
        return [round(100 * (c - h), 1), round(100 * (c + h), 1)]
    thr_tbl = {}
    for t in (0.2, 0.5):
        g = sp.dropna(subset=["spliceai_max"])
        tp = int(((g.label == "POS") & (g.spliceai_max >= t)).sum()); fn = int((g.label == "POS").sum()) - tp
        fp = int(((g.label == "NEG") & (g.spliceai_max >= t)).sum()); tn = int((g.label == "NEG").sum()) - fp
        thr_tbl[f"SpliceAI >= {t}"] = {"TP": tp, "FN": fn, "FP": fp, "TN": tn,
                                       "sensitivity_%": round(100 * tp / (tp + fn), 1), "sens_CI": wilson(tp, tp + fn),
                                       "flagged_among_no_effect_%": round(100 * fp / (fp + tn), 1),
                                       "flag_CI": wilson(fp, fp + tn)}
    out["thresholds_splicing_class"] = thr_tbl

    # Provenance by outcome, which the reviewers asked to be made explicit.
    out["provenance_by_outcome"] = sp.groupby(["origin", "label"]).size().unstack(fill_value=0).to_dict()

    g1 = sp[sp.gene == "PKD1"]
    out["composite_scores_PKD1_only"] = {
        "rationale": "PKD1 lies on chromosome 16, which is absent from the AVI training and "
                     "validation chromosome lists; PKD2 lies on chromosome 4, which is a training "
                     "chromosome, so overlap cannot be excluded for PKD2 (it is not demonstrated "
                     "either: no position-level overlap check was possible).",
        "all_tools": {t: auc_ci(g1[g1.label == "POS"][c].dropna().to_numpy(float),
                                g1[g1.label == "NEG"][c].dropna().to_numpy(float))
                      for t, c in tools.items()},
        "noncoding_synonymous": {
            t: auc_ci(g1[(g1.label == "POS") & g1.consequence_class.isin(["noncoding", "synonymous"])][c].dropna().to_numpy(float),
                      g1[(g1.label == "NEG") & g1.consequence_class.isin(["noncoding", "synonymous"])][c].dropna().to_numpy(float))
            for t, c in (("AVI", "avi"), ("CADD", "cadd_raw"))}}

    # Resampling and sampling details, requested so the intervals can be reproduced exactly.
    dup = lab[lab.pmids.astype(str).str.contains(";")][["gene", "hgvs_c", "pmids"]]
    out["resampling_details"] = {
        "seeds": {"auc_ci": 20260920, "cluster_auc_ci": 20260920, "paired_diff": 20260920,
                  "paired_cluster_bootstrap": 20260921},
        "n_resamples": 2000,
        "interval_method": "percentile bootstrap, 2.5th and 97.5th centiles",
        "variant_level_bootstrap": "positives and negatives resampled independently with "
                                   "replacement, preserving the observed class sizes",
        "cluster_bootstrap": "whole source studies resampled with replacement; a resample "
                             "containing only one outcome class yields no AUC and is discarded, "
                             "and the number retained is reported alongside the interval",
        "paired_comparison": "both tools evaluated on the SAME resampled indices, which is why the "
                             "paired interval can be narrow while the marginal intervals are wide",
        "variants_with_two_sources": {
            "n": int(len(dup)),
            "handling": "assigned to the FIRST listed PMID for clustering, so each variant belongs "
                        "to exactly one cluster and no variant is counted twice",
            "variants": dup.to_dict("records")},
        "ties": "AUC counts ties as one half, so identical scores neither help nor harm a tool"}

    (RES / "stratified.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({k: out[k] for k in ("assay_class_counts", "primary_splicing_class")}, indent=1, default=str)[:2600])


if __name__ == "__main__":
    main()
