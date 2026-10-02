#!/usr/bin/env python3
"""
16_references.py

Builds the manuscript reference list from authoritative records, not from memory.

Every reference with a PMID is fetched from NCBI E-utilities (esummary) and
formatted in Vancouver style from the returned fields (authors, title, journal
abbreviation, year, volume, issue, pages, DOI). Nothing is typed by hand except
the two items that have no PubMed record (the AlphaGenome Atlas preprint and the
software resources), whose details are taken from the stored PDF and stated as
such.

It then rewrites the manuscript, replacing [PMID nnnnn] and [NAMED-KEY] markers
with numbered citations in order of first appearance, and appends the numbered
list.

Usage:  python3 16_references.py ../publication/DRAFT_manuscript_v1.0.md

Output: <manuscript>.referenced.md and publication/references.md

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
PUB = ROOT / "publication"
CACHE = ROOT / "data" / "raw" / "esummary_cache.json"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Items without a PubMed record. Details taken from the stored source file named in `source`.
MANUAL = {
    "ATLAS": {
        "text": ("AlphaGenome Atlas team; Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: "
                 "in silico mutagenesis of the entire human genome improves prioritization and "
                 "interpretation of non-coding variants. Google DeepMind; 2026. Available from: "
                 "https://deepmind.google/science/alphagenome/atlas"),
        "source": "docs/sources/alphagenome_atlas_preprint.pdf"},
    "ATLASDOCS": {
        "text": ("Google DeepMind. AlphaGenome documentation: recommended splicing score and variant "
                 "scoring definitions. Available from: https://www.alphagenomedocs.com/faqs.html "
                 "(accessed 20 September 2026)."),
        "source": "quoted formula: max(splice_sites) + max(splice_site_usage) + max(splice_junctions)/5"},
    "VARIANTVALIDATOR_WEB": {
        "text": ("VariantValidator [Internet]. Manchester: University of Manchester. Available from: "
                 "https://variantvalidator.org (accessed 19 September 2026)."),
        "source": "used via its REST API in analysis/05_verify.py"},
}

# PMIDs cited in the Methods/Introduction/Discussion that are not truth-set sources.
# Each is resolved from NCBI; the expectation recorded here is checked against the fetched title.
NAMED = {
    "ALPHAGENOME": (41606153, "advancing regulatory variant effect prediction"),
    "SPLICEAI": (30661751, "splicing from primary sequence with deep learning"),
    "CADD17": (38183205, "CADD"),
    "GNOMAD4": (38057664, "76,156 human genomes"),  # gnomAD v4 flagship
    "VARIANTVALIDATOR": (28967166, "variantvalidator"),
    "CLINGEN": (37352859, "impact on splicing"),
}


def esummary(pmids: list[int]) -> dict:
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    todo = [str(p) for p in pmids if str(p) not in cache]
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        r = requests.post(f"{EU}/esummary.fcgi", data={"db": "pubmed", "id": ",".join(chunk),
                                                      "retmode": "json"}, timeout=120)
        r.raise_for_status()
        res = r.json()["result"]
        for p in chunk:
            if p in res:
                cache[p] = res[p]
        time.sleep(0.4)
    CACHE.write_text(json.dumps(cache))
    return cache


def vancouver(rec: dict) -> str:
    auths = [a["name"] for a in rec.get("authors", []) if a.get("authtype") == "Author"]
    if len(auths) > 6:
        alist = ", ".join(auths[:6]) + ", et al"
    else:
        alist = ", ".join(auths)
    title = rec.get("title", "").rstrip(".")
    journal = rec.get("source", "")
    year = (rec.get("pubdate") or "")[:4]
    vol, issue, pages = rec.get("volume", ""), rec.get("issue", ""), rec.get("pages", "")
    doi = next((a["value"] for a in rec.get("articleids", []) if a["idtype"] == "doi"), "")
    bits = f"{alist}. {title}. {journal}. {year}"
    if vol:
        bits += f";{vol}"
        if issue:
            bits += f"({issue})"
        if pages:
            bits += f":{pages}"
        elif doi and "/" in doi:
            # Journals that number articles rather than pages (Oxford titles such as Nucleic Acids
            # Research) return an empty pages field from esummary. The article identifier is the
            # final segment of the DOI, which is what such journals expect in the citation.
            bits += f":{doi.rsplit('/', 1)[1]}"
    bits += "."
    if doi:
        bits += f" doi:{doi}."
    return bits


def main() -> None:
    ms = Path(sys.argv[1])
    text = ms.read_text()
    used_pmids = [int(p) for p in re.findall(r"\[PMID\s+([0-9;\s]+)\]", text) for p in re.findall(r"\d+", p)]
    named_pmids = [v[0] for v in NAMED.values()]
    recs = esummary(sorted(set(used_pmids + named_pmids)))

    for key, (pmid, expect) in NAMED.items():
        rec = recs.get(str(pmid))
        if not rec:
            raise SystemExit(f"{key}: PMID {pmid} not found; fix before formatting references")
        if expect.lower() not in rec.get("title", "").lower():
            raise SystemExit(f"{key}: PMID {pmid} title is '{rec.get('title')}', expected to contain "
                             f"'{expect}'. Refusing to cite the wrong paper.")

    order: list[str] = []

    def number(tag: str) -> int:
        if tag not in order:
            order.append(tag)
        return order.index(tag) + 1

    def repl(m: re.Match) -> str:
        """Number BOTH citation styles in a single pass, so numbering follows document order.

        Until 21 September 2026 this ran two passes, every [PMID ...] first and every [REF:...]
        second, which numbered a Discussion PMID ahead of an Introduction named reference and put
        reference 8 before reference 7 in the text.
        """
        if m.group("pmids"):
            nums = sorted(number(p) for p in re.findall(r"\d+", m.group("pmids")))
            return "[" + ",".join(str(n) for n in nums) + "]"
        key = m.group("named")
        pmid = str(NAMED[key][0]) if key in NAMED else key
        return f"[{number(pmid)}]"

    out = re.sub(r"\[PMID\s+(?P<pmids>[0-9;\s]+)\]|\[REF:(?P<named>[A-Z0-9_]+)\]", repl, text)

    lines = ["# References", "",
             "Generated by analysis/16_references.py from NCBI E-utilities records on "
             f"{time.strftime('%d %B %Y')}. Vancouver style, formatted from the returned fields.", ""]
    for i, tag in enumerate(order, 1):
        if tag in MANUAL:
            lines.append(f"{i}. {MANUAL[tag]['text']}")
        else:
            rec = recs.get(tag)
            if not rec:
                lines.append(f"{i}. [PMID {tag}: no record returned; complete by hand]")
            else:
                lines.append(f"{i}. {vancouver(rec)} PMID: {tag}.")
    (PUB / "references.md").write_text("\n".join(lines) + "\n")
    # Replace ONLY the reference block: from the "## References" heading to the next heading of the
    # same level. Splitting on a horizontal rule silently deleted every section after References
    # (declarations, figure legends) in manuscripts that do not use one.
    head, marker, tail = out.partition("## References")
    if not marker:
        raise SystemExit("no '## References' heading found")
    m = re.search(r"\n## ", tail)
    rest = tail[m.start():] if m else ""
    body = "\n".join(lines[4:])
    out = head + "## References\n\n" + body + "\n" + rest
    Path(str(ms).replace(".md", ".referenced.md")).write_text(out)
    print(f"{len(order)} references; wrote publication/references.md and "
          f"{Path(str(ms).replace('.md', '.referenced.md')).name}")
    for ln in lines[4:10]:
        print("  ", ln[:150])


if __name__ == "__main__":
    main()
