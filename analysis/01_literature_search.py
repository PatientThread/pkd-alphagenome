#!/usr/bin/env python3
"""
01_literature_search.py

Recorded, reproducible PubMed search for functionally validated non-coding and
splice-acting variants in PKD1 and PKD2. The truth set must not depend on which
papers happened to surface in ad hoc web searching, so the query is fixed here,
run through NCBI E-utilities, and every hit is written to disk with its date.

Output: data/raw/pubmed_search_<date>.tsv, data/raw/pubmed_search_<date>.json

Author: Christopher Lawrence
"""
from __future__ import annotations

import datetime as dt
import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = {
    # Q1: gene-named papers on splicing / non-coding mechanisms
    "Q1_splice_noncoding": (
        '("PKD1"[tiab] OR "PKD2"[tiab] OR "polycystin"[tiab]) AND '
        '(intronic[tiab] OR "non-coding"[tiab] OR noncoding[tiab] OR "deep intronic"[tiab] '
        'OR "untranslated"[tiab] OR "UTR"[tiab] OR splicing[tiab] OR "splice"[tiab] '
        'OR minigene[tiab] OR pseudoexon[tiab] OR "pseudo-exon"[tiab] OR branchpoint[tiab] '
        'OR "branch point"[tiab] OR "cryptic"[tiab] OR "RNA analysis"[tiab] OR "RT-PCR"[tiab] '
        'OR promoter[tiab] OR "upstream open reading frame"[tiab] OR uORF[tiab])'),
    # Q2: disease-named papers with functional RNA work
    "Q2_adpkd_functional": (
        '("polycystic kidney disease"[tiab] OR ADPKD[tiab]) AND '
        '(minigene[tiab] OR pseudoexon[tiab] OR "deep intronic"[tiab] OR "intronic variant"[tiab] '
        'OR "splicing variant"[tiab] OR "splice variant"[tiab] OR "aberrant splicing"[tiab] '
        'OR "5\' UTR"[tiab] OR "5\'UTR"[tiab] OR "3\'UTR"[tiab] OR "3\' UTR"[tiab] '
        'OR "transcriptome sequencing"[tiab] OR "RNA sequencing"[tiab] OR "allele-specific expression"[tiab])'),
}


def esearch(term: str) -> list[str]:
    r = requests.get(f"{EU}/esearch.fcgi", params={
        "db": "pubmed", "term": term, "retmax": 5000, "retmode": "json"}, timeout=60)
    r.raise_for_status()
    return r.json()["esearchresult"]["idlist"]


def esummary(pmids: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i + 200]
        r = requests.post(f"{EU}/esummary.fcgi", data={
            "db": "pubmed", "id": ",".join(chunk), "retmode": "json"}, timeout=120)
        r.raise_for_status()
        res = r.json()["result"]
        for p in chunk:
            d = res.get(p, {})
            ids = {a["idtype"]: a["value"] for a in d.get("articleids", [])}
            out.append({
                "pmid": p, "pmcid": ids.get("pmc", ""), "doi": ids.get("doi", ""),
                "year": (d.get("pubdate") or "")[:4], "journal": d.get("source", ""),
                "title": d.get("title", ""),
                "first_author": (d.get("authors") or [{}])[0].get("name", ""),
            })
        time.sleep(0.4)
    return out


def main() -> None:
    today = dt.date.today().isoformat()
    hits: dict[str, set[str]] = {}
    for name, q in QUERIES.items():
        ids = esearch(q)
        print(f"{name}: {len(ids)} hits")
        for p in ids:
            hits.setdefault(p, set()).add(name)
        time.sleep(0.4)
    rows = esummary(sorted(hits, key=int))
    df = pd.DataFrame(rows)
    df["matched_queries"] = df["pmid"].map(lambda p: ";".join(sorted(hits[p])))
    df = df.sort_values("year", ascending=False)
    df.to_csv(RAW / f"pubmed_search_{today}.tsv", sep="\t", index=False)
    (RAW / f"pubmed_search_{today}.json").write_text(json.dumps(
        {"date": today, "queries": QUERIES, "n_unique": len(df)}, indent=2))
    print(f"unique PMIDs: {len(df)}; with PMC full text: {(df.pmcid != '').sum()}")


if __name__ == "__main__":
    main()
