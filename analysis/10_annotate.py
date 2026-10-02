#!/usr/bin/env python3
"""
10_annotate.py

Two annotations the analysis plan requires for every truth-set variant.

1. CODING CONSEQUENCE CLASS (amendment A1), from the VariantValidator protein
   prediction already cached by 05_verify.py: noncoding (5'UTR / intronic, p.?
   or p.(=) at an intronic or UTR position), synonymous, missense, nonsense,
   other_coding.

2. CLINVAR STATUS (plan section 6, leakage). NCBI E-utilities search of ClinVar
   by GRCh38 position, then an exact match on the transcript-level c. notation
   in the record title. Records presence, accession, germline classification.

Output: data/curated/annotations.tsv

Author: Christopher Lawrence
"""
from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.parse
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
CUR = ROOT / "data" / "curated"
CACHE = ROOT / "data" / "raw" / "vv_cache"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
VV = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{}/all?content-type=application/json"


def vv_record(transcript: str, hgvs_c: str) -> dict:
    url = VV.format(urllib.parse.quote(f"{transcript}:{hgvs_c}", safe=""))
    key = CACHE / (hashlib.sha1(url.encode()).hexdigest() + ".json")
    data = json.loads(key.read_text())
    return next(v for k, v in data.items() if k not in ("flag", "metadata"))


def consequence(hgvs_c: str, rec: dict) -> tuple[str, str]:
    p = ((rec.get("hgvs_predicted_protein_consequence") or {}).get("tlr") or "")
    body = p.split(":", 1)[-1]
    c = hgvs_c
    # A range that starts in an exon and ends in an intron (or the reverse)
    # removes coding bases AND intronic bases, so it is not purely non-coding.
    if re.search(r"c\.\d+_\d+[+-]\d+", c) or re.search(r"c\.\d+[+-]\d+_\d+(?![-+\d])", c):
        return "other_coding", body
    intronic_or_utr = bool(re.search(r"c\.-|c\.\*|\d[+-]\d", c))
    if intronic_or_utr and not re.search(r"_\d+[^+-]", c.replace("c.", "")):
        return "noncoding", body
    if "Ter" in body or "*" in body.replace("fs*", "fs"):
        if "fs" not in body:
            return "nonsense", body
    if body.endswith("=)") or body.endswith("=") or "(=)" in body:
        return "synonymous", body
    if re.search(r"p\.\([A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)$", body):
        return "missense", body
    if intronic_or_utr:
        return "noncoding", body
    return "other_coding", body


def clinvar(chrom: str, pos: str, hgvs_c: str, transcript: str) -> dict:
    term = f"{chrom}[chr] AND {pos}[chrpos38]"
    r = requests.get(f"{EU}/esearch.fcgi", params={"db": "clinvar", "term": term,
                                                  "retmode": "json", "retmax": 50}, timeout=60)
    ids = r.json()["esearchresult"]["idlist"]
    time.sleep(0.4)
    out = {"clinvar_present": False, "clinvar_accession": "", "clinvar_classification": "",
           "clinvar_title": ""}
    if not ids:
        return out
    s = requests.get(f"{EU}/esummary.fcgi", params={"db": "clinvar", "id": ",".join(ids),
                                                   "retmode": "json"}, timeout=60).json()["result"]
    time.sleep(0.4)
    want = hgvs_c.replace(" ", "")
    for i in ids:
        d = s.get(i, {})
        title = d.get("title", "")
        m = re.search(r"\):(c\.[^\s(]+)", title)
        if m and m.group(1) == want and transcript.split(".")[0] in title:
            gc = d.get("germline_classification", {}) or {}
            out.update(clinvar_present=True, clinvar_accession=d.get("accession", ""),
                       clinvar_classification=gc.get("description", ""), clinvar_title=title)
            break
    return out


def main() -> None:
    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"].drop_duplicates(["chrom", "pos", "ref", "alt"])
    rows = []
    for r in ts.to_dict("records"):
        rec = vv_record(r["transcript"], r["hgvs_c"])
        cls, p = consequence(r["hgvs_c"], rec)
        cv = clinvar(r["chrom"], r["pos"], r["hgvs_c"], r["transcript"])
        rows.append(dict(entry_id=r["entry_id"], gene=r["gene"], hgvs_c=r["hgvs_c"], label=r["label"],
                         protein=p, consequence_class=cls, **cv))
        print(f"{r['entry_id']} {r['gene']} {r['hgvs_c']:32s} {cls:12s} {p:28s} "
              f"{'ClinVar ' + cv['clinvar_classification'] if cv['clinvar_present'] else '-'}")
    pd.DataFrame(rows).to_csv(CUR / "annotations.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
