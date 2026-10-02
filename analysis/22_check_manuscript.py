#!/usr/bin/env python3
"""
22_check_manuscript.py

Checks that the numbers asserted in the manuscript still match the result files.

Twice now a figure quoted in the text has survived a reanalysis that changed it: the
count-versus-proportion burden correction, and the reclassification of rs3874648 under
amendment A8. This script fails loudly rather than relying on a careful reading.

Each check names the claim, the value read from results/, and where in the manuscript the
claim is made. Exit status is non-zero if any check fails.

Usage: python3 22_check_manuscript.py [publication/EJHG_manuscript_v2.1.md]

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
DEFAULT = ROOT / "publication" / "EJHG_manuscript_v2.1.md"


def main() -> None:
    ms = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    # The tables are part of the submission, so a number stated only in a table still counts as
    # reported. Banned phrases are checked against both for the same reason.
    companions = [ms.parent / "tables_ejhg.md"]
    stem = ms.stem.split("_manuscript")[0]
    companions += sorted(ms.parent.glob(f"{stem}_supplementary*.md"))
    text = ms.read_text()
    for c in companions:
        if c.exists():
            text += "\n" + c.read_text()
    st = json.loads((RES / "stratified.json").read_text())
    sm = json.loads((RES / "summary.json").read_text())
    prim = st["primary_splicing_class"]
    sa, ag = prim["auc"]["SpliceAI"], prim["auc"]["AlphaGenome splicing composite"]
    cs, th = st["common_subset_all_tools"], st["thresholds_splicing_class"]
    prov = st["provenance_by_outcome"]
    nb = sm["locus_burden_noncoding_background"]["PKD1"]
    wl = sm["locus_burden"]["PKD1"]
    rm = json.loads((RES / "ranking_manifest.json").read_text())
    p1 = st["composite_scores_PKD1_only"]
    wlb, ncb = rm["backgrounds"]["PKD1"]["whole locus"], rm["backgrounds"]["PKD1"]["non-coding positions only"]
    why = rm["why_the_backgrounds_disagree"]

    checks: list[tuple[str, str]] = [
        ("splicing class size", f"{prim['n_pos']} positives and {prim['n_neg']}"),
        ("SpliceAI AUC", f"{sa['auc']:.3f}"),
        ("AlphaGenome AUC", f"{ag['auc']:.3f}"),
        ("paired difference", f"{prim['paired_ag_composite_vs_spliceai']['auc_difference']:.3f}"),
        ("AVI AUC", f"{prim['auc']['AVI']['auc']:.3f}"),
        ("CADD AUC", f"{prim['auc']['CADD']['auc']:.3f}"),
        ("common subset n", f"{cs['n_pos'] + cs['n_neg']} variants scored by all four tools"),
        ("sensitivity at 0.2", f"{th['SpliceAI >= 0.2']['sensitivity_%']:.1f}"),
        ("flag rate at 0.2", f"{th['SpliceAI >= 0.2']['flagged_among_no_effect_%']:.1f}"),
        ("TP at 0.2", f"{th['SpliceAI >= 0.2']['TP']} of {prim['n_pos']} positives"),
        ("patient-derived positives", f"{prov['POS']['patient']} were carried by patients"),
        ("engineered negatives", f"{prov['NEG']['engineered']} were engineered"),
        ("matched-burden recovery", f"{wl['pos_caught_spliceai_0_5']} of {wl['pos_nc_syn_snvs']}"),
        ("non-coding burden", f"{nb['spliceai_ge_0_5_pct_of_noncoding']:.2f}%"),
        ("non-coding recovery SpliceAI", f"{nb['pos_caught_spliceai_0_5']} of {nb['n_pos_noncoding']}"),
        ("non-coding recovery AVI", f"{nb['pos_caught_avi_same_burden']} of {nb['n_pos_noncoding']}"),
        # These appear as prose in one journal format and as table rows in another, so only the
        # value is checked: that the value is stated correctly is the assertion, and the words
        # around it are the journal's house style.
        ("whole-locus SpliceAI alleles", f"{wlb['spliceai']['alleles_scored']:,}"),
        ("whole-locus AVI alleles", f"{wlb['avi']['alleles_scored']:,}"),
        ("whole-locus AVI threshold", f"{wlb['avi']['threshold']:.3f}"),
        ("non-coding SpliceAI alleles", f"{ncb['spliceai']['alleles_scored']:,}"),
        ("non-coding AVI threshold", f"{ncb['avi']['threshold']:.3f}"),
        ("AVI top slice size", f"{why['avi_top_slice_alleles']:,} alleles above"),
        ("AVI top slice non-coding", f"only {why['percent_of_avi_top_slice_that_is_non_coding']}% are non-coding"),
        ("locus non-coding share", f"{why['percent_of_whole_locus_that_is_non_coding']}% of the locus"),
        ("PKD1-only AVI", f"AVI gives {p1['all_tools']['AVI']['auc']:.3f}"),
        ("PKD1-only AVI n", f"on {p1['all_tools']['AVI']['n_pos']} positives"),
        # Stated in Table 2 rather than in prose, so matched as a table row.
        ("PKD1-only CADD", f"| CADD | PKD1 only | {p1['all_tools']['CADD']['n_pos']} | "
                           f"{p1['all_tools']['CADD']['n_neg']} | {p1['all_tools']['CADD']['auc']:.3f} |"),
        ("PKD1 nc/syn AVI", f"| AVI | PKD1, non-coding and synonymous | "
                            f"{p1['noncoding_synonymous']['AVI']['n_pos']} | "
                            f"{p1['noncoding_synonymous']['AVI']['n_neg']} | "
                            f"{p1['noncoding_synonymous']['AVI']['auc']:.3f} |"),
        ("splicing-class engineered", f"contains 92 engineered variants"),
    ]

    fails = 0
    for name, needle in checks:
        # Numbers are written with ordinary spacing in the text; collapse whitespace before matching
        # so that a line break inside a sentence does not produce a spurious failure.
        flat = re.sub(r"\s+", " ", text)
        ok = re.sub(r"\s+", " ", needle) in flat
        print(f"{'OK  ' if ok else 'FAIL'}  {name:30s} expects: {needle}")
        fails += not ok

    banned = {
        "transcript level shown to be unchanged": "withdrawn: Wedd performed no RNA studies",
        "Following external review": "this was an AI-assisted critique, not peer review",
        "rather than as evidence in itself": "ClinGen allows PP3/BP4; do not overstate",
        "dominates unresolved ADPKD": "not established by this study",
        "no accuracy reason to change": "clinical interchangeability not supported",
        "performed identically": "equal observed AUC does not mean identical per-variant behaviour",
        "was seen in training": "chromosome eligibility is not demonstrated exposure",
        "not independent of the model's training data": "overstates what the split shows",
        "indistinguishable from wild type": 'use "no effect detected in the reported assay"',
        "two tools cannot score insertions": "the precomputed files lack them; the tools can",
        "the translation class": "renamed reporter-output class; RNA unresolved for three variants",
        "Figure legends begin on a new page": "instruction to the journal, not manuscript text",
        "most often sought": "unsupported frequency claim",
        "other_coding": "render as 'other coding'",
    }
    flat_text = re.sub(r"\s+", " ", text)
    for phrase, why in banned.items():
        if re.sub(r"\s+", " ", phrase) in flat_text:
            print(f"FAIL  banned phrase present: {phrase!r} ({why})")
            fails += 1

    # Citation order: first appearance of each reference number must be 1, 2, 3, ... The numbering
    # script numbered PMID and named citations in separate passes until 21 September 2026, which put
    # reference 8 in the text before reference 7.
    ref = Path(str(ms).replace(".md", ".referenced.md"))
    if ref.exists():
        body = ref.read_text()
        body = body[:body.index("## References")] if "## References" in body else body
        seen: list[int] = []
        for m in re.finditer(r"\[(\d+(?:,\d+)*)\]", body):
            for n in (int(x) for x in m.group(1).split(",")):
                if n not in seen:
                    seen.append(n)
        expected = list(range(1, len(seen) + 1))
        ok = seen == expected
        print(f"{'OK  ' if ok else 'FAIL'}  {'citation order':30s} first appearances: {seen}")
        fails += not ok

        figs = [int(m.group(1)) for m in re.finditer(r"Figure (\d)", body)]
        first = sorted(set(figs), key=figs.index)
        ok = first == sorted(first)
        print(f"{'OK  ' if ok else 'FAIL'}  {'figure order':30s} first appearances: {first}")
        fails += not ok

    print(f"\n{len(checks)} numeric checks, {fails} failure(s)")
    raise SystemExit(1 if fails else 0)


if __name__ == "__main__":
    main()
