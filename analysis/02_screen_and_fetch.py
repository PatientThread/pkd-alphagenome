#!/usr/bin/env python3
"""
02_screen_and_fetch.py

Title screen of the recorded PubMed search, then retrieval of every included
paper: PMC open-access full text as BioC JSON where it exists, the PubMed
abstract for everything. Every downstream curated variant must be anchored to a
verbatim quote in one of these files, so this is the evidence store.

INCLUSION: human PKD1 or PKD2 variants with RNA-level, reporter or regulatory
functional evidence, or papers defining PKD1/PKD2 regulatory elements.
EXCLUDED at title stage: animal-only models, cell signalling, therapeutics,
reviews, other genes (PKHD1, IFT140, etc.) and homonyms (protein kinase D,
Kluyveromyces pKD1 plasmid, PKD2L1).

Output: data/curated/screening_log.tsv, data/papers/<pmid>.{bioc.json,abstract.txt}

Author: Christopher Lawrence
"""
from __future__ import annotations

import glob
import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
PAPERS = ROOT / "data" / "papers"
PAPERS.mkdir(parents=True, exist_ok=True)
CUR = ROOT / "data" / "curated"
CUR.mkdir(parents=True, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/fullTextXML"
BIOC = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{}/unicode"

# pmid: (category, reason). Category A = variant-level functional evidence,
# B = regulatory element definition, C = context (checked, may yield variants).
INCLUDE = {
    "41006799": ("A", "PKD1 5'UTR variants, translation reporter"),
    "41242437": ("A", "de novo PKD1 splicing variant, molecular characterisation"),
    "41558825": ("B", "PKD1 3'UTR miR-17 six-nucleotide motif"),
    "41279527": ("B", "PKD1 3'UTR miR-17 six-nucleotide motif (second record)"),
    "42502691": ("A", "PKD1 intron 16 pseudoexon, complex allele, RNA-seq"),
    "42518289": ("B", "PKD1 uORFs and PC1 expression"),
    "42621998": ("A", "PKD1 c.7489+5G>A minigene"),
    "39757590": ("B", "PKD1 uORFs translational regulation"),
    "40069205": ("C", "targeted long-read PKD1, may report splice variants"),
    "41083429": ("A", "PKD1 splice variant, Chinese family"),
    "41141524": ("A", "PKD1 intronic variants, minigene"),
    "41659037": ("A", "PKD1 intronic micro-deletion, minigene"),
    "34749493": ("A", "PKD2 atypical splice mutation"),
    "38944240": ("A", "targeted RNA-seq of urinary cells, noncoding ADPKD variants"),
    "39169606": ("A", "PKD2 aberrant splicing, Korean family"),
    "39183159": ("C", "19 novel ADPKD variants, may include splicing assays"),
    "39200828": ("A", "PKD2 single base substitution dual exon skipping"),
    "37272738": ("A", "PKD1 variant causing alternative splicing"),
    "37419908": ("A", "PKD1 atypical splicing variants, RNA studies"),
    "37468838": ("A", "11 exonic PKD1/PKD2 variants alter splicing, minigene; negatives"),
    "34257392": ("A", "novel SNVs altering RNA splicing of PKD1 and PKD2"),
    "36029107": ("A", "common intronic SNV modifies PKD1 expression"),
    "34542828": ("A", "presumed synonymous PKD2 variant"),
    "30185468": ("C", "abnormal alternative splicing of wild-type PKD1"),
    "30424739": ("A", "novel PKD1 splicing mutation, case report"),
    "26692149": ("A", "three exonic PKD2 mutations alter splicing, minigene"),
    "27984604": ("A", "novel PKD1 splicing mutation, pedigree"),
    "25757501": ("A", "exonic PKD1 mutations causing splicing defects"),
    "24575920": ("A", "two PKD1 mutations affecting one splice site"),
    "24907393": ("A", "PKD1 missense and synonymous variants defective splicing"),
    "20950398": ("A", "PKD2 aberrant splicing from presumed missense"),
    "19158373": ("A", "atypical splice mutations in ADPKD, pathogenicity evidence"),
    "11800305": ("A", "PKD1 abnormal RNA processing, Thai family"),
    "10612835": ("A", "PKD1 IVS13-2A>T RNA processing defect"),
    "11058904": ("A", "PKD1 splicing mutations, expression of mutated genes"),
    "10541293": ("A", "PKD2 aberrant splicing"),
    "7633405":  ("A", "PKD1 intronic deletion causing splicing mutations"),
    "12070253": ("C", "PKD1 transcript screening by RT-PCR"),
    "36000373": ("C", "PKD1 pedigree analysis"),
    "31384335": ("C", "PKD1 rs201204878 co-segregation"),
    "35497784": ("C", "ADPKD genetic test case, rare genotype"),
    "29986647": ("B", "PKD2 long-range regulatory mechanisms"),
    "36283570": ("B", "enhancer landscape in PKD"),
    "15770226": ("B", "PKD1 and PKD2 promoter regulatory elements"),
    "17890878": ("B", "Pkd1 promoter Sp1 regulation"),
    "25173570": ("B", "PKD2 promoter characterisation"),
    "12832634": ("C", "PKD1 last intron conservation"),
}


def fetch_abstracts(pmids: list[str]) -> None:
    r = requests.post(f"{EU}/efetch.fcgi", data={
        "db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "text"},
        timeout=120)
    r.raise_for_status()
    # efetch text returns records separated by blank-line numbered blocks; also
    # fetch per PMID so each file is unambiguous.
    for p in pmids:
        rr = requests.get(f"{EU}/efetch.fcgi", params={
            "db": "pubmed", "id": p, "rettype": "abstract", "retmode": "text"}, timeout=60)
        rr.raise_for_status()
        (PAPERS / f"{p}.abstract.txt").write_text(rr.text)
        time.sleep(0.35)


def main() -> None:
    src = sorted(glob.glob(str(ROOT / "data/raw/pubmed_search_*.tsv")))[-1]
    df = pd.read_csv(src, sep="\t", dtype=str).fillna("")
    df["include"] = df["pmid"].isin(INCLUDE)
    df["category"] = df["pmid"].map(lambda p: INCLUDE.get(p, ("", ""))[0])
    df["reason"] = df["pmid"].map(lambda p: INCLUDE.get(p, ("", "title: out of scope"))[1])
    df.to_csv(CUR / "screening_log.tsv", sep="\t", index=False)
    missing = set(INCLUDE) - set(df["pmid"])
    assert not missing, f"included PMIDs not in search output: {missing}"
    inc = df[df.include]
    print(f"screened {len(df)}; included {len(inc)} "
          f"(A={sum(inc.category=='A')}, B={sum(inc.category=='B')}, C={sum(inc.category=='C')})")

    fetch_abstracts(list(inc.pmid))
    got = 0
    for r in inc.itertuples():
        if not r.pmcid:
            continue
        out = PAPERS / f"{r.pmid}.bioc.json"
        if (PAPERS / f"{r.pmid}.epmc.xml").exists():
            got += 1
            continue
        if out.exists():
            got += 1
            continue
        resp = requests.get(BIOC.format(r.pmcid), timeout=120)
        ok = False
        if resp.status_code == 200:
            # The BioC service returns HTTP 200 with an "[Error] : No result"
            # HTML page when it has no record. A prefix check on "[" passes
            # that page, so parse it: only valid JSON counts.
            try:
                json.loads(resp.text)
                ok = True
            except ValueError:
                ok = False
        if ok:
            out.write_text(resp.text)
            got += 1
        else:
            # Fallback: Europe PMC full-text XML for the same PMCID.
            x = requests.get(EPMC.format(r.pmcid), timeout=120)
            if x.status_code == 200 and x.text.lstrip().startswith("<"):
                (PAPERS / f"{r.pmid}.epmc.xml").write_text(x.text)
                got += 1
            else:
                print(f"  no OA full text for {r.pmid} {r.pmcid} "
                      f"(BioC {resp.status_code}, EuropePMC {x.status_code})")
        time.sleep(0.4)
    print(f"abstracts: {len(inc)}; full texts: {got}")


if __name__ == "__main__":
    main()
