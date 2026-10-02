#!/usr/bin/env python3
"""
23_ranking_manifest.py

Documents the two ranking comparisons completely, because the manuscript previously
conflated two different operations:

  (i)  a DENOMINATOR CORRECTION, counts -> proportions, which was a bug and applies to BOTH
       backgrounds. It is history, recorded in the amendment log, and every number now
       reported already has it applied.
  (ii) a CHANGE OF CANDIDATE BACKGROUND, whole locus -> non-coding positions only, which is
       not a correction at all. It changes the question being asked, and it is the reason the
       two comparisons disagree.

For each background and tool this writes the number of scored alternate alleles, the
threshold, the achieved flag proportion, the eligible benchmark variants by identifier, and
the recovery counts. Units are alternate alleles throughout: both backgrounds carry the three
possible substitutions at each position.

Output: results/ranking_manifest.json

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
RES = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))


def exon_intervals(gene: str) -> list[tuple[int, int]]:
    ex = json.loads((ROOT / "data" / "raw" / "exons_grch38.json").read_text())
    return [(a - 2, b + 2) for a, b in ex[gene]["exons"]]


def main() -> None:
    assay_class = importlib.import_module("18_stratified_analysis").assay_class
    df = pd.read_csv(RES / "per_variant.tsv", sep="\t")
    df["assay_class"] = [assay_class(str(a), str(m)) for a, m in zip(df.assay, df.mechanism)]
    df = df[df.assay_class == "splicing"]

    out: dict = {"unit": "alternate allele (three substitutions per position in both backgrounds)",
                 "matching_rule": "identical for both backgrounds: the AVI threshold is set so that "
                                  "AVI flags the same PROPORTION of its background as SpliceAI >= 0.5 "
                                  "flags of its own background",
                 "history": "Until 20 September 2026 the two tools were matched on raw COUNTS. The "
                            "SpliceAI background is a 1-in-4 sample of positions and the AVI "
                            "background is complete, so equal counts set a fourfold stricter AVI "
                            "threshold. Amendment A7 corrected this. Every number below is "
                            "post-correction; the comparison between the two backgrounds is NOT a "
                            "consequence of that bug.",
                 "backgrounds": {}}

    for gene in ("PKD1", "PKD2"):
        ap, sp = RES / f"atlas_locus_{gene}.tsv.gz", RES / f"spliceai_locusbg_{gene}.tsv.gz"
        if not (ap.exists() and sp.exists()):
            continue
        ab = pd.read_csv(ap, sep="\t")
        ab["pos"] = ab.variant.str.split(":").str[1].astype(int)
        sb = pd.read_csv(sp, sep="\t")
        iv = exon_intervals(gene)

        def nc(pos: np.ndarray) -> np.ndarray:
            m = np.ones(len(pos), bool)
            for a, b in iv:
                m &= ~((pos >= a) & (pos <= b))
            return m

        g = df[(df.gene == gene) & df.is_snv]
        for name, abx, sbx, elig in (
                ("whole locus", ab, sb, g[g.consequence_class.isin(["noncoding", "synonymous"])]),
                ("non-coding positions only", ab[nc(ab.pos.to_numpy())], sb[nc(sb.pos.to_numpy())],
                 g[g.consequence_class == "noncoding"])):
            k = int((sbx.delta_max.round(2) >= 0.5).sum())
            frac = k / len(sbx)
            kk = int(round(frac * len(abx)))
            thr = float(np.sort(abx.avi.to_numpy())[::-1][kk - 1]) if kk else None
            p, n = elig[elig.label == "POS"], elig[elig.label == "NEG"]
            out["backgrounds"].setdefault(gene, {})[name] = {
                "spliceai": {"alleles_scored": int(len(sbx)), "positions": int(sbx.pos.nunique()),
                             "threshold": 0.5, "flagged": k,
                             "flag_proportion_pct": round(100 * frac, 4),
                             "resolution_one_allele_pct": round(100 / len(sbx), 4)},
                "avi": {"alleles_scored": int(len(abx)), "positions": int(abx.pos.nunique()),
                        "threshold": round(thr, 6) if thr else None, "flagged": kk,
                        "flag_proportion_pct": round(100 * kk / len(abx), 4)},
                "eligible_benchmark_variants": {
                    "n_pos": int(len(p)), "n_neg": int(len(n)),
                    "pos_ids": sorted(p.hgvs_c.tolist()), "neg_ids": sorted(n.hgvs_c.tolist())},
                "recovery": {
                    "spliceai_pos": int((p.spliceai_max >= 0.5).sum()),
                    "avi_pos": int((p.avi >= thr).sum()) if thr else None,
                    "spliceai_neg_flagged": int((n.spliceai_max >= 0.5).sum()),
                    "avi_neg_flagged": int((n.avi >= thr).sum()) if thr else None}}

    # Why the two backgrounds disagree: what the whole-locus top slice is made of.
    ab = pd.read_csv(RES / "atlas_locus_PKD1.tsv.gz", sep="\t")
    ab["pos"] = ab.variant.str.split(":").str[1].astype(int)
    iv = exon_intervals("PKD1")
    is_nc = np.ones(len(ab), bool)
    for a, b in iv:
        is_nc &= ~((ab.pos.to_numpy() >= a) & (ab.pos.to_numpy() <= b))
    thr = out["backgrounds"]["PKD1"]["whole locus"]["avi"]["threshold"]
    top = ab[ab.avi >= thr]
    out["why_the_backgrounds_disagree"] = {
        "explanation": "AVI is a whole-genome pathogenicity score and ranks amino-acid-changing "
                       "variants highly. Against a whole-locus background its top slice is filled "
                       "by coding positions, which leaves few places for non-coding positives. "
                       "Restricting the background to non-coding positions removes that competition. "
                       "The disagreement between the two comparisons is therefore a property of the "
                       "question asked, not an artefact.",
        "avi_top_slice_alleles": int(len(top)),
        "percent_of_avi_top_slice_that_is_non_coding": round(100 * float(is_nc[ab.avi.to_numpy() >= thr].mean()), 1),
        "percent_of_whole_locus_that_is_non_coding": round(100 * float(is_nc.mean()), 1)}

    (RES / "ranking_manifest.json").write_text(json.dumps(out, indent=1))
    b = out["backgrounds"]["PKD1"]
    for name, v in b.items():
        r, e = v["recovery"], v["eligible_benchmark_variants"]
        print(f"{name:28s} SpliceAI {v['spliceai']['flag_proportion_pct']:.3f}% of "
              f"{v['spliceai']['alleles_scored']:,} alleles -> {r['spliceai_pos']}/{e['n_pos']}; "
              f"AVI thr {v['avi']['threshold']:.3f} {v['avi']['flag_proportion_pct']:.3f}% of "
              f"{v['avi']['alleles_scored']:,} -> {r['avi_pos']}/{e['n_pos']}; "
              f"negatives flagged {r['spliceai_neg_flagged']}/{e['n_neg']} and "
              f"{r['avi_neg_flagged']}/{e['n_neg']}")
    print(json.dumps(out["why_the_backgrounds_disagree"], indent=1)[-320:])


if __name__ == "__main__":
    main()
