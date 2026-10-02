#!/usr/bin/env python3
"""
03_extract_mentions.py

Mechanical extraction of every variant mention (HGVS c./g./r., legacy IVS
notation, rsIDs) from the evidence store, each with its verbatim sentence and
section. Curation (04) then works only from these quotes, and 05 re-checks that
every curated quote is present verbatim in the source text.

Also writes data/papers/<pmid>.txt: one plain-text rendering per paper with
section labels, which is the canonical text the quote check runs against.

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PAPERS = ROOT / "data" / "papers"
CUR = ROOT / "data" / "curated"

# c./n./g./r. HGVS, tolerant of the spacing and unicode dashes journals introduce
POS = r"[-−–*]?\d+(?:\s?[+\-−–]\s?\d+)?"
HGVS = re.compile(
    rf"\b[cgrn]\.\s?(?:\(?{POS}(?:\s?_\s?{POS})?\)?)"
    rf"\s?(?:[ACGTacgu]+\s?(?:>|&gt;)\s?[ACGTacgu]+|del[ACGT0-9]*(?:ins[ACGT0-9]+)?|dup[ACGT0-9]*|ins[ACGT0-9]+|inv|=)?")
IVS = re.compile(r"\bIVS\s?\d+\s?[+\-−–]\s?\d+\s?[ACGT]?\s?(?:>|&gt;)?\s?[ACGT]?")
RSID = re.compile(r"\brs\d{4,}\b")
SENT = re.compile(r"(?<=[.!?])\s+(?=[A-Z(])")


def bioc_text(path: Path) -> list[tuple[str, str]]:
    data = json.loads(path.read_text())
    docs = data if isinstance(data, list) else [data]
    out = []
    for coll in docs:
        for doc in coll.get("documents", []):
            for p in doc.get("passages", []):
                sec = p.get("infons", {}).get("section_type", "") or p.get("infons", {}).get("type", "")
                txt = p.get("text", "") or ""
                if txt.strip():
                    out.append((sec, txt))
    return out


def epmc_text(path: Path) -> list[tuple[str, str]]:
    root = ET.fromstring(path.read_text())
    out = []
    for sec in root.iter():
        if sec.tag in ("p", "title", "td", "th", "caption", "article-title"):
            txt = "".join(sec.itertext()).strip()
            if txt:
                out.append((sec.tag, txt))
    return out


def abstract_text(path: Path) -> list[tuple[str, str]]:
    return [("abstract", t) for t in path.read_text().split("\n\n") if t.strip()]


def main() -> None:
    rows = []
    pmids = sorted({p.name.split(".")[0] for p in PAPERS.iterdir() if p.suffix in (".json", ".xml", ".txt")
                    and not p.name.endswith(".fulltext.txt")})
    for pmid in pmids:
        passages = []
        if (PAPERS / f"{pmid}.bioc.json").exists():
            passages = bioc_text(PAPERS / f"{pmid}.bioc.json"); kind = "fulltext"
        elif (PAPERS / f"{pmid}.epmc.xml").exists():
            passages = epmc_text(PAPERS / f"{pmid}.epmc.xml"); kind = "fulltext"
        else:
            kind = "abstract"
        if (PAPERS / f"{pmid}.abstract.txt").exists() and kind == "abstract":
            passages = abstract_text(PAPERS / f"{pmid}.abstract.txt")
        canon = "\n\n".join(f"[{s}] {t}" for s, t in passages)
        (PAPERS / f"{pmid}.fulltext.txt").write_text(canon)
        for sec, txt in passages:
            for sent in SENT.split(txt):
                for rx, kindm in ((HGVS, "hgvs"), (IVS, "ivs"), (RSID, "rsid")):
                    for m in rx.finditer(sent):
                        rows.append({"pmid": pmid, "text_kind": kind, "section": sec,
                                     "mention_type": kindm, "mention": m.group(0).strip(),
                                     "sentence": sent.strip()[:700]})
    df = pd.DataFrame(rows).drop_duplicates()
    df.to_csv(CUR / "mentions.tsv", sep="\t", index=False)
    summ = df.groupby(["pmid", "text_kind"]).mention.nunique().reset_index(name="unique_mentions")
    print(summ.to_string(index=False))
    print(f"\ntotal mention rows {len(df)}, unique (pmid, mention) {df[['pmid','mention']].drop_duplicates().shape[0]}")


if __name__ == "__main__":
    main()
