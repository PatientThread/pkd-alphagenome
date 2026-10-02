#!/usr/bin/env python3
"""
12_analyse.py

Implements docs/ANALYSIS_PLAN.md (with amendments A1 and A2) and nothing else.
Each section is labelled with the plan item it answers. No threshold, metric
or subset appears here that is not in the plan; where the plan left a choice
open, the choice is stated in the output.

Inputs: data/curated/{truthset_verified,annotations,gnomad_overlap}.tsv,
        results/{atlas_truthset,ondemand_scores,spliceai_scores,cadd_truthset}.tsv,
        results/{atlas_locus,cadd_locus}_<gene>.tsv.gz
Outputs: results/per_variant.tsv, results/summary.json, results/summary.md

AlphaGenome Output is non-commercial; see LEGALLY_BINDING_TERMS_OF_USE.txt.

Author: Christopher Lawrence
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CUR = ROOT / "data" / "curated"
RES = ROOT / "results"
FEATURES_SPLICE = "fi_MERGED_SPLICING"
FEATURES_EXPR = ["fi_MAX_ABS_RNA_SEQ", "fi_MAX_ABS_CAGE", "fi_MAX_ABS_PROCAP", "fi_MAX_ABS_POLYADENYLATION",
                 "fi_MAX_ABS_DNASE", "fi_MAX_ABS_ATAC", "fi_MAX_ABS_CHIP_TF", "fi_MAX_ABS_CHIP_HISTONE",
                 "fi_MAX_ABS_CONTACT_MAPS"]
FEATURES_PROTEIN = ["fi_PROTEIN_TERMINATION", "fi_ALPHAMISSENSE", "fi_START_LOST", "fi_STOP_LOST"]
FEATURES_CONS = ["fi_CACTUS_241_WAY", "fi_PHASTCONS_470_WAY"]


def key(chrom, pos, ref, alt) -> str:
    return f"chr{chrom}:{pos}:{ref}>{alt}"


def locus_pct(values: pd.Series, background: np.ndarray) -> pd.Series:
    """Percentile of each value within the locus background (mid-rank for ties), 100 = highest."""
    bg = np.sort(background[~np.isnan(background)])
    def one(x):
        if pd.isna(x):
            return np.nan
        lo = np.searchsorted(bg, x, "left"); hi = np.searchsorted(bg, x, "right")
        return 100.0 * (lo + 0.5 * (hi - lo)) / len(bg)
    return values.map(one)


def load() -> pd.DataFrame:
    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"]
    # one row per unique variant; keep every supporting source
    grp = ts.groupby(["gene", "chrom", "pos", "ref", "alt"], as_index=False).agg(
        entry_id=("entry_id", "first"), entry_ids=("entry_id", ";".join), hgvs_c=("hgvs_c", "first"),
        label=("label", lambda s: s.iloc[0] if s.nunique() == 1 else "DISCORDANT"),  # A3 rule
        labels_all=("label", ";".join), mechanism=("mechanism", "first"), assay=("assay", " | ".join),
        origin=("origin", lambda s: "patient" if "patient" in set(s) else s.iloc[0]),
        pmids=("pmid", ";".join))
    ann = pd.read_csv(CUR / "annotations.tsv", sep="\t", dtype=str).fillna("")
    gno = pd.read_csv(CUR / "gnomad_overlap.tsv", sep="\t", dtype=str).fillna("")
    df = grp.merge(ann[["entry_id", "protein", "consequence_class", "clinvar_present", "clinvar_classification"]],
                   on="entry_id", how="left")
    df = df.merge(gno[["entry_id", "gnomad_genome_af", "gnomad_genome_overlap", "avi_leakage_flag"]],
                  on="entry_id", how="left")
    df["key"] = [key(*t) for t in df[["chrom", "pos", "ref", "alt"]].itertuples(index=False)]
    df["is_snv"] = (df.ref.str.len() == 1) & (df.alt.str.len() == 1)

    od = pd.read_csv(RES / "ondemand_scores.tsv", sep="\t")
    od["key"] = [key(*t) for t in od[["chrom", "pos", "ref", "alt"]].astype(str).itertuples(index=False)]
    keep = [c for c in od.columns if c.startswith(("splice_", "rna_seq_")) and not c.endswith("n_kidney_tracks")]
    df = df.merge(od[["key"] + keep].drop_duplicates("key"), on="key", how="left")

    # Documented AlphaGenome splicing composite (alphagenomedocs.com FAQ, accessed 20 Sep 2026):
    #   max(splice_sites) + max(splice_site_usage) + max(splice_junctions)/5
    # with each component the maximum absolute score across all tracks AND genes. We report that
    # (suffix _anygene, the recommendation) and the gene-restricted variant used elsewhere here.
    for suff in ("anygene", "all"):
        cols = [f"splice_sites_maxabs_{suff}", f"splice_site_usage_maxabs_{suff}",
                f"splice_junctions_maxabs_{suff}"]
        if all(c in df.columns for c in cols):
            df[f"ag_splicing_composite_{suff}"] = (df[cols[0]] + df[cols[1]] + df[cols[2]] / 5)

    sa = pd.read_csv(RES / "spliceai_scores.tsv", sep="\t")
    sa["key"] = [key(*t) for t in sa[["chrom", "pos", "ref", "alt"]].astype(str).itertuples(index=False)]
    df = df.merge(sa[["key", "spliceai_max", "spliceai_ag", "spliceai_al", "spliceai_dg", "spliceai_dl",
                      "spliceai_mnv_extension"]].drop_duplicates("key"), on="key", how="left")

    cd = pd.read_csv(RES / "cadd_truthset.tsv", sep="\t", dtype=str)
    cd["key"] = [key(*t) for t in cd[["chrom", "pos", "ref", "alt"]].itertuples(index=False)]
    cd["cadd_phred"] = pd.to_numeric(cd.cadd_phred, errors="coerce")
    cd["cadd_raw"] = pd.to_numeric(cd.cadd_raw, errors="coerce")
    df = df.merge(cd[["key", "cadd_phred", "cadd_raw", "cadd_source"]].drop_duplicates("key"), on="key", how="left")

    at_path = RES / "atlas_truthset.tsv"
    if at_path.exists():
        at = pd.read_csv(at_path, sep="\t")
        at = at.rename(columns={"variant": "key"})
        df = df.merge(at, on="key", how="left")
    else:
        df["avi"] = np.nan

    # locus percentiles (plan section 5)
    df["avi_locus_pct"] = np.nan
    df["cadd_locus_pct"] = np.nan
    for gene in ("PKD1", "PKD2"):
        m = df.gene == gene
        ap = RES / f"atlas_locus_{gene}.tsv.gz"
        if ap.exists():
            bg = pd.read_csv(ap, sep="\t")["avi"].to_numpy(float)
            df.loc[m, "avi_locus_pct"] = locus_pct(df.loc[m, "avi"], bg)
            df.loc[m, "avi_locus_n"] = len(bg)
        sp = RES / f"spliceai_locus_{gene}.tsv.gz"
        if not sp.exists():
            sp = RES / f"spliceai_locusbg_{gene}.tsv.gz"   # evenly spaced sample (amendment A5)
        if sp.exists():
            bg = pd.read_csv(sp, sep="\t")["delta_max"].to_numpy(float)
            df.loc[m, "spliceai_locus_pct"] = locus_pct(df.loc[m, "spliceai_max"], bg)
            df.loc[m, "spliceai_locus_n"] = len(bg)
        cp = RES / f"cadd_locus_{gene}.tsv.gz"
        if cp.exists():
            bg = pd.read_csv(cp, sep="\t")["raw"].to_numpy(float)
            df.loc[m, "cadd_locus_pct"] = locus_pct(df.loc[m, "cadd_raw"], bg)
            df.loc[m, "cadd_locus_n"] = len(bg)

    # AVI attribution: which feature group drives each score (amendment A1 item 3)
    if FEATURES_SPLICE in df.columns:
        groups = {"splicing": [FEATURES_SPLICE], "expression_chromatin": FEATURES_EXPR,
                  "protein": FEATURES_PROTEIN, "conservation": FEATURES_CONS}
        for g, cols in groups.items():
            cols = [c for c in cols if c in df.columns]
            df[f"avi_attr_{g}"] = df[cols].sum(axis=1, min_count=1)
        gcols = [f"avi_attr_{g}" for g in groups]
        has = df[gcols].notna().any(axis=1)
        df["avi_top_attr_group"] = pd.Series([None] * len(df), dtype="object")
        df.loc[has, "avi_top_attr_group"] = (df.loc[has, gcols].idxmax(axis=1)
                                             .str.replace("avi_attr_", "", regex=False))
    return df


def auc_boot(pos_vals, neg_vals, n_boot: int = 2000, seed: int = 20260919) -> dict:
    """Mann-Whitney AUC with a stratified bootstrap 95% interval (amendment A3)."""
    a = np.asarray(pd.to_numeric(pd.Series(pos_vals), errors="coerce").dropna(), float)
    b = np.asarray(pd.to_numeric(pd.Series(neg_vals), errors="coerce").dropna(), float)
    if len(a) < 3 or len(b) < 3:
        return {"n_pos": int(len(a)), "n_neg": int(len(b)), "auc": None}
    def auc(x, y):
        return float(((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum()) / (len(x) * len(y)))
    rng = np.random.default_rng(seed)
    boots = [auc(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(n_boot)]
    return {"n_pos": int(len(a)), "n_neg": int(len(b)), "auc": round(auc(a, b), 3),
            "ci95": [round(float(np.percentile(boots, 2.5)), 3), round(float(np.percentile(boots, 97.5)), 3)]}


def describe(sub: pd.DataFrame, col: str) -> dict:
    x = pd.to_numeric(sub[col], errors="coerce").dropna()
    return {"n_scored": int(len(x)), "n_missing": int(len(sub) - len(x)),
            "median": None if x.empty else round(float(x.median()), 3)}


def main() -> None:
    df = load()
    df.to_csv(RES / "per_variant.tsv", sep="\t", index=False)

    # Amendment A8.5. The locus-rank and matched-burden analyses put a variant's score against a
    # background of splice predictions, so they answer a splicing question and are restricted to the
    # splicing evidence class. Until the Chen 5'UTR reporter variants were added on 20 September 2026
    # this restriction held by accident; it is now explicit, so that translation-reporter evidence
    # never enters a splicing statistic.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    assay_class = importlib.import_module("18_stratified_analysis").assay_class
    df["assay_class"] = [assay_class(str(a), str(m)) for a, m in zip(df.assay, df.mechanism)]
    df = df[df.assay_class == "splicing"].copy()
    S: dict = {"notes": [
        "Locus percentile: 100 = highest-scoring of all possible SNVs in gene span + 5 kb upstream.",
        "AlphaGenome on-demand splicing outputs are reported separately; no composite was pre-specified.",
        "AVI and CADD are composite pathogenicity scores; per amendment A1 their primary evaluation "
        "uses non-coding and synonymous variants only.",
        "Per amendment A2, AVI is reported primarily on PKD1 (chromosome 16, held out of AVI training "
        "and validation).",
        "locus_burden: SpliceAI locus percentile is pre-specified (plan 5). The equal-burden comparison "
        "(AVI threshold chosen to flag as many locus SNVs as SpliceAI >= 0.5) is EXPLORATORY, added "
        "19 Sep 2026 before the gene-wide SpliceAI results existed, and is not a pre-specified endpoint."]}

    pos, neg = df[df.label == "POS"], df[df.label == "NEG"]
    nc_syn = df.consequence_class.isin(["noncoding", "synonymous"])

    # ---------------- PRIMARY ENDPOINT (plan 5, A1, A2)
    prim = df[(df.label == "POS") & (df.gene == "PKD1") & df.is_snv & nc_syn]
    P = {"definition": "PKD1 POS SNVs, non-coding or synonymous; AVI percentile within the PKD1 locus",
         "n": int(len(prim))}
    if prim.avi_locus_pct.notna().any():
        p = prim.avi_locus_pct.dropna()
        P.update(n_scored=int(len(p)), median_pct=round(float(p.median()), 2),
                 top1pct=int((p >= 99).sum()), top0_1pct=int((p >= 99.9).sum()),
                 per_variant=prim[["hgvs_c", "avi", "avi_locus_pct", "avi_top_attr_group"]]
                 .round(3).to_dict("records"))
    c = prim.cadd_locus_pct.dropna()
    P["cadd_same_variants"] = {"n_scored": int(len(c)), "median_pct": round(float(c.median()), 2) if len(c) else None,
                               "top1pct": int((c >= 99).sum()), "top0_1pct": int((c >= 99.9).sum()),
                               "phred_ge_20": int((prim.cadd_phred >= 20).sum())}
    s = prim.spliceai_max.dropna()
    P["spliceai_same_variants"] = {"n_scored": int(len(s)), "ge_0_2": int((s >= 0.2).sum()),
                                   "ge_0_5": int((s >= 0.5).sum()), "median": round(float(s.median()), 3) if len(s) else None}
    S["primary"] = P

    # ---------------- secondary: every POS / NEG, by tool
    # ---------------- plan 5 secondary: SpliceAI on the same locus-percentile scale, and the
    # false-alarm burden each tool would put on a laboratory across the whole gene
    burden = {}
    for gene in ("PKD1", "PKD2"):
        ap = RES / f"atlas_locus_{gene}.tsv.gz"
        sp = RES / f"spliceai_locus_{gene}.tsv.gz"
        sampled = False
        if not sp.exists():
            sp = RES / f"spliceai_locusbg_{gene}.tsv.gz"
            sampled = True
        if not sp.exists():
            continue
        sb = pd.read_csv(sp, sep="\t")
        g = {"spliceai_background": "evenly spaced sample (A5)" if sampled else "every position",
             "n_snvs": int(len(sb)), "in_spliceai_annotation": int(sb.in_annotation.sum()),
             "spliceai_ge_0_2": int((sb.delta_max.round(2) >= 0.2).sum()),
             "spliceai_ge_0_5": int((sb.delta_max.round(2) >= 0.5).sum())}
        g["spliceai_ge_0_5_pct_of_locus"] = round(100 * g["spliceai_ge_0_5"] / len(sb), 3)
        g["spliceai_ge_0_2_pct_of_locus"] = round(100 * g["spliceai_ge_0_2"] / len(sb), 3)
        if ap.exists():
            ab = pd.read_csv(ap, sep="\t")["avi"].to_numpy(float)
            # AVI threshold that flags the same PROPORTION of the locus as SpliceAI >= 0.5.
            # Matching raw counts is wrong: the SpliceAI background is an evenly spaced sample
            # (amendment A5) whereas the AVI background is complete, so the same count is a four
            # times stricter threshold for AVI. This error made AVI recover 1 of 25 instead of 5.
            frac = g["spliceai_ge_0_5"] / g["n_snvs"]
            kk = int(round(frac * len(ab)))
            thr = float(np.sort(ab)[::-1][kk - 1]) if kk else None
            g["avi_matched_burden_pct_of_locus"] = round(100 * kk / len(ab), 3)
            g["avi_threshold_matching_spliceai_0_5_burden"] = thr
            gpos = df[(df.gene == gene) & (df.label == "POS") & df.is_snv & df.consequence_class.isin(["noncoding", "synonymous"])]
            g["pos_nc_syn_snvs"] = int(len(gpos))
            g["pos_caught_spliceai_0_5"] = int((gpos.spliceai_max >= 0.5).sum())
            g["pos_caught_avi_same_burden"] = int((gpos.avi >= thr).sum()) if thr is not None else None
            gneg = df[(df.gene == gene) & (df.label == "NEG") & df.is_snv & df.consequence_class.isin(["noncoding", "synonymous"])]
            g["neg_nc_syn_snvs"] = int(len(gneg))
            g["neg_flagged_spliceai_0_5"] = int((gneg.spliceai_max >= 0.5).sum())
            g["neg_flagged_avi_same_burden"] = int((gneg.avi >= thr).sum()) if thr is not None else None
        prim_g = df[(df.gene == gene) & (df.label == "POS") & df.is_snv
                    & df.consequence_class.isin(["noncoding", "synonymous"])]
        if "spliceai_locus_pct" in df:
            g["pos_nc_syn_spliceai_locus_pct_median"] = describe(prim_g, "spliceai_locus_pct")["median"]
            g["pos_nc_syn_avi_locus_pct_median"] = describe(prim_g, "avi_locus_pct")["median"]
        burden[gene] = g
    S["locus_burden"] = burden

    # ---------------- EXPLORATORY (added 20 Sep 2026, after seeing the first burden result, and
    # recorded as such): the same comparison with the background restricted to NON-CODING positions.
    # AVI is a whole-genome pathogenicity score, so ranking a non-coding candidate against a
    # background dominated by protein-altering variants understates it. A laboratory hunting a
    # non-coding cause would compare against non-coding variants. Exons are the canonical Ensembl
    # transcript (data/raw/exons_grch38.json), padded by 2 bp to exclude canonical splice sites.
    exons = json.loads((ROOT / "data" / "raw" / "exons_grch38.json").read_text())
    ncb, ncb_rows = {}, []
    for gene in ("PKD1", "PKD2"):
        ap = RES / f"atlas_locus_{gene}.tsv.gz"
        sp = RES / f"spliceai_locus_{gene}.tsv.gz"
        if not sp.exists():
            sp = RES / f"spliceai_locusbg_{gene}.tsv.gz"
        if not (ap.exists() and sp.exists()):
            continue
        iv = [(a - 2, b + 2) for a, b in exons[gene]["exons"]]
        def is_nc(pos: np.ndarray) -> np.ndarray:
            m = np.ones(len(pos), bool)
            for a, b in iv:
                m &= ~((pos >= a) & (pos <= b))
            return m
        ab = pd.read_csv(ap, sep="\t")
        ab["pos"] = ab.variant.str.split(":").str[1].astype(int)
        ab = ab[is_nc(ab.pos.to_numpy())]
        sb = pd.read_csv(sp, sep="\t")
        sb = sb[is_nc(sb.pos.to_numpy())]
        sub = df[(df.gene == gene) & df.is_snv & (df.consequence_class == "noncoding")]
        pos_s, neg_s = sub[sub.label == "POS"], sub[sub.label == "NEG"]
        k = int((sb.delta_max.round(2) >= 0.5).sum())
        frac = k / len(sb)                      # proportion, not count (see note above)
        kk = int(round(frac * len(ab)))
        thr = float(np.sort(ab.avi.to_numpy())[::-1][kk - 1]) if kk else None
        ncb[gene] = {
            "noncoding_background_avi_snvs": int(len(ab)),
            "noncoding_background_spliceai_snvs": int(len(sb)),
            "spliceai_ge_0_5_pct_of_noncoding": round(100 * k / len(sb), 3),
            "avi_matched_burden_pct_of_noncoding": round(100 * kk / len(ab), 3),
            "avi_threshold_matching_that_burden": thr,
            "n_pos_noncoding": int(len(pos_s)), "n_neg_noncoding": int(len(neg_s)),
            "pos_caught_spliceai_0_5": int((pos_s.spliceai_max >= 0.5).sum()),
            "pos_caught_avi_same_burden": int((pos_s.avi >= thr).sum()) if thr else None,
            "neg_flagged_spliceai_0_5": int((neg_s.spliceai_max >= 0.5).sum()),
            "neg_flagged_avi_same_burden": int((neg_s.avi >= thr).sum()) if thr else None,
            "pos_median_avi_pct_noncoding_bg": round(float(locus_pct(pos_s.avi, ab.avi.to_numpy()).median()), 2),
            "pos_median_spliceai_pct_noncoding_bg": round(float(locus_pct(pos_s.spliceai_max, sb.delta_max.to_numpy()).median()), 2)}
        # Per-variant ranks against the same non-coding background, so Figure 3b plots exactly the
        # comparison the matched-burden statistic above is computed from.
        ncb_rows.append(sub.assign(
            spliceai_pct_nc_bg=locus_pct(sub.spliceai_max, sb.delta_max.to_numpy()).values,
            avi_pct_nc_bg=locus_pct(sub.avi, ab.avi.to_numpy()).values)
            [["gene", "hgvs_c", "label", "spliceai_pct_nc_bg", "avi_pct_nc_bg"]])
    S["locus_burden_noncoding_background"] = ncb
    if ncb_rows:
        pd.concat(ncb_rows).to_csv(RES / "noncoding_background_ranks.tsv", sep="\t", index=False)

    S["pos_vs_neg"] = {}
    for col in ["avi_locus_pct", "spliceai_locus_pct", "cadd_locus_pct", "cadd_phred", "spliceai_max",
                "splice_sites_maxabs_all", "splice_site_usage_maxabs_all", "splice_junctions_maxabs_all",
                "splice_site_usage_maxabs_kidney", "splice_junctions_maxabs_kidney", "rna_seq_maxabs_all"]:
        if col in df:
            S["pos_vs_neg"][col] = {"POS": describe(pos, col), "NEG": describe(neg, col),
                                    "NEG_values": neg[["hgvs_c", "gene", "consequence_class", col]]
                                    .round(3).to_dict("records")}
    S["spliceai_thresholds_all_pos"] = {"n": int(pos.spliceai_max.notna().sum()),
                                        "ge_0_2": int((pos.spliceai_max >= 0.2).sum()),
                                        "ge_0_5": int((pos.spliceai_max >= 0.5).sum())}
    S["spliceai_thresholds_all_neg"] = {"n": int(neg.spliceai_max.notna().sum()),
                                        "ge_0_2": int((neg.spliceai_max >= 0.2).sum()),
                                        "ge_0_5": int((neg.spliceai_max >= 0.5).sum())}

    # ---------------- H1 splice-acting PKD1 non-coding positives
    h1 = df[(df.label == "POS") & (df.consequence_class == "noncoding") & (df.gene == "PKD1")
            & df.is_snv & ~df.hgvs_c.str.startswith("c.-")]
    S["H1"] = {"definition": "PKD1 non-coding POS SNVs excluding 5'UTR (splice-acting); "
                             "top 1% of locus by AVI with splicing as the largest attribution group",
               "n": int(len(h1))}
    if h1.avi_locus_pct.notna().any():
        S["H1"].update(top1pct=int((h1.avi_locus_pct >= 99).sum()),
                       splicing_dominant=int((h1.avi_top_attr_group == "splicing").sum()),
                       both=int(((h1.avi_locus_pct >= 99) & (h1.avi_top_attr_group == "splicing")).sum()))

    # ---------------- H2 c.-69dupG
    def row(h):
        r = df[df.hgvs_c == h]
        return {} if r.empty else r.iloc[0].replace({np.nan: None}).to_dict()
    cols_h = ["hgvs_c", "label", "avi", "avi_locus_pct", "cadd_phred", "cadd_locus_pct", "spliceai_max",
              "splice_sites_maxabs_all", "splice_site_usage_maxabs_all", "splice_junctions_maxabs_all",
              "rna_seq_maxabs_all", "rna_seq_signed_at_max_all", "rna_seq_maxabs_kidney", "rna_seq_signed_at_max_kidney"]
    S["H2_c.-69dupG"] = {k: row("c.-69dupG").get(k) for k in cols_h}
    S["H2_context_5UTR_negatives"] = [{k: row(h).get(k) for k in cols_h} for h in ("c.-52C>T", "c.-209G>A")]

    # ---------------- H3 pseudoexon haplotype (X1 from the side files)
    od = pd.read_csv(RES / "ondemand_scores.tsv", sep="\t")
    sa = pd.read_csv(RES / "spliceai_scores.tsv", sep="\t")
    cd = pd.read_csv(RES / "cadd_truthset.tsv", sep="\t")
    h3 = []
    for eid in ("E015", "X1", "E016"):
        o = od[od.entry_id == eid].iloc[0]; a = sa[sa.entry_id == eid].iloc[0]; d = cd[cd.entry_id == eid].iloc[0]
        h3.append({"entry": eid, "variant": o.hgvs_c,
                   "splice_sites": round(o.splice_sites_maxabs_all, 4),
                   "splice_site_usage": round(o.splice_site_usage_maxabs_all, 4),
                   "splice_junctions": round(o.splice_junctions_maxabs_all, 4),
                   "spliceai_max": a.spliceai_max, "cadd_phred": d.cadd_phred})
    S["H3_pseudoexon"] = h3

    # ---------------- H4 allele resolution
    S["H4_allele_resolution"] = [{k: row(h).get(k) for k in cols_h}
                                 for h in ("c.7866C>A", "c.7866C>T", "c.11257C>A", "c.11257C>G", "c.11257C>T")]

    # ---------------- H5 exonic minigene positives versus near-splice positives
    bmc = df[(df.label == "POS") & df.pmids.str.contains("37468838")]
    near = df[(df.label == "POS") & (df.consequence_class == "noncoding") & ~df.hgvs_c.str.startswith("c.-")]
    S["H5"] = {name: {c: describe(sub, c) for c in ["spliceai_max", "splice_site_usage_maxabs_all",
                                                    "avi_locus_pct", "cadd_locus_pct"]}
               for name, sub in (("BMC_exonic_POS", bmc), ("near_splice_noncoding_POS", near))}

    # ---------------- H6 common modifier
    S["H6_rs3874648"] = {k: row("c.10051-239G>A").get(k) for k in cols_h + ["gnomad_genome_af"]}

    # ---------------- A2 PKD2 split by leakage flag
    pk2 = df[(df.gene == "PKD2") & (df.label == "POS")]
    S["A2_PKD2"] = {f: {"n": int((pk2.avi_leakage_flag == f).sum()),
                        "avi_locus_pct_median": describe(pk2[pk2.avi_leakage_flag == f], "avi_locus_pct")["median"]}
                    for f in ("True", "False")}

    # ---------------- kidney vs all tracks (secondary)
    S["kidney_vs_all"] = {c: {"POS_median_all": describe(pos, f"{c}_maxabs_all")["median"],
                              "POS_median_kidney": describe(pos, f"{c}_maxabs_kidney")["median"]}
                          for c in ("splice_site_usage", "splice_junctions", "rna_seq")}

    # ---------------- ClinVar stratification (plan 6)
    if df.avi_locus_pct.notna().any():
        S["clinvar_strata_PKD1_POS"] = {
            s: describe(df[(df.gene == "PKD1") & (df.label == "POS") & (df.clinvar_present == ("True" if s == "in_clinvar" else "False"))],
                        "avi_locus_pct") for s in ("in_clinvar", "not_in_clinvar")}

    # ---------------- A3: AUC with bootstrap interval, on the classes each tool may be judged on
    lab = df[df.label.isin(["POS", "NEG"])]
    ncs = lab[lab.consequence_class.isin(["noncoding", "synonymous"])]
    A = {"splicing_specific_all_variants": {}, "composite_noncoding_synonymous": {}, "composite_all_variants_secondary": {}}
    for col in ["spliceai_max", "ag_splicing_composite_anygene", "ag_splicing_composite_all",
                "splice_sites_maxabs_all", "splice_site_usage_maxabs_all",
                "splice_junctions_maxabs_all", "splice_site_usage_maxabs_kidney", "splice_junctions_maxabs_kidney"]:
        if col in lab:
            A["splicing_specific_all_variants"][col] = auc_boot(lab[lab.label == "POS"][col], lab[lab.label == "NEG"][col])
    for col in ["avi", "cadd_raw"]:
        A["composite_noncoding_synonymous"][col] = auc_boot(ncs[ncs.label == "POS"][col], ncs[ncs.label == "NEG"][col])
        A["composite_all_variants_secondary"][col] = auc_boot(lab[lab.label == "POS"][col], lab[lab.label == "NEG"][col])
    # the same, restricted to the held-out gene
    p1 = lab[lab.gene == "PKD1"]
    A["PKD1_only_splicing_specific"] = {c: auc_boot(p1[p1.label == "POS"][c], p1[p1.label == "NEG"][c])
                                        for c in ["spliceai_max", "splice_sites_maxabs_all", "splice_junctions_maxabs_all"]}
    # sensitivity: c.11257C>T kept as POS (its original label) instead of DISCORDANT
    d2 = df.copy(); d2.loc[d2.hgvs_c == "c.11257C>T", "label"] = "POS"
    l2 = d2[d2.label.isin(["POS", "NEG"])]
    A["sensitivity_c11257CT_as_POS"] = {c: auc_boot(l2[l2.label == "POS"][c], l2[l2.label == "NEG"][c])
                                        for c in ["spliceai_max", "splice_sites_maxabs_all", "splice_junctions_maxabs_all"]}
    A["discordant"] = df[df.label == "DISCORDANT"][["hgvs_c", "labels_all"]].to_dict("records")
    A["counts"] = df.label.value_counts().to_dict()
    # Paired difference in AUC (same variants, resampled together). Reported with an interval
    # instead of inferring equivalence from overlapping marginal intervals.
    def paired_auc_diff(col_a, col_b, n_boot=2000, seed=20260920):
        sub = lab[[col_a, col_b, "label"]].dropna()
        a_p, a_n = sub[sub.label == "POS"][col_a].to_numpy(float), sub[sub.label == "NEG"][col_a].to_numpy(float)
        b_p, b_n = sub[sub.label == "POS"][col_b].to_numpy(float), sub[sub.label == "NEG"][col_b].to_numpy(float)
        if min(len(a_p), len(a_n)) < 3:
            return None
        def auc(x, y):
            return float(((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum()) / (len(x) * len(y)))
        rng = np.random.default_rng(seed)
        obs = auc(a_p, a_n) - auc(b_p, b_n)
        boots = []
        for _ in range(n_boot):
            ip = rng.integers(0, len(a_p), len(a_p)); ineg = rng.integers(0, len(a_n), len(a_n))
            boots.append(auc(a_p[ip], a_n[ineg]) - auc(b_p[ip], b_n[ineg]))
        return {"n_pos": int(len(a_p)), "n_neg": int(len(a_n)), "auc_difference": round(obs, 3),
                "ci95": [round(float(np.percentile(boots, 2.5)), 3), round(float(np.percentile(boots, 97.5)), 3)],
                "note": "positive favours the first score; interval crossing zero means the data do "
                        "not resolve which is better, which is not the same as equivalence"}
    A["paired_differences"] = {
        "ag_composite_anygene_minus_spliceai": paired_auc_diff("ag_splicing_composite_anygene", "spliceai_max"),
        "ag_usage_minus_spliceai": paired_auc_diff("splice_site_usage_maxabs_all", "spliceai_max")}
    S["A3_auc"] = A

    # ---------------- validation: local SpliceAI versus the SpliceAI scores Xie et al. published
    s2p = ROOT / "data" / "papers" / "34257392_DatasetS2_transcribed.tsv"
    if s2p.exists():
        s2 = pd.read_csv(s2p, sep="\t", comment="#")
        mm = s2.merge(df[["hgvs_c", "spliceai_max"]], left_on="location", right_on="hgvs_c", how="left")
        d = (mm.spliceai_max - mm.spliceai_published).abs()
        S["validation_spliceai_vs_published"] = {
            "n": int(mm.spliceai_max.notna().sum()), "exact_to_2dp": int((d < 0.005).sum()),
            "within_0_02": int((d <= 0.02).sum()), "max_abs_diff": round(float(d.max()), 3),
            "discrepancies": mm[d > 0.02][["location", "spliceai_published", "spliceai_max"]].to_dict("records")}

    (RES / "summary.json").write_text(json.dumps(S, indent=2, default=str))
    print(json.dumps(S, indent=1, default=str)[:12000])


if __name__ == "__main__":
    main()
