#!/usr/bin/env python3
"""
05_verify.py

Two independent checks on every curated entry. An entry reaches the verified
truth set only if it passes both.

1. QUOTE CHECK. The entry's quote must occur verbatim in the stored source text
   (data/papers/<pmid>.fulltext.txt, falling back to the abstract). Only
   typography is normalised: whitespace runs, unicode minus and dashes, curly
   quotes. Wording is not.

2. VARIANTVALIDATOR. The HGVS string is validated against its RefSeq transcript
   on GRCh38. This confirms the stated reference base is really there and
   yields the hg38 VCF representation used for scoring. Any VariantValidator
   warning is recorded, and an entry with a reference mismatch fails.

rsIDs are resolved through Ensembl (hg38 position and alleles) and then passed
through VariantValidator like everything else.

Output: data/curated/truthset_verified.tsv (all entries, with pass/fail columns)
        data/raw/vv_cache/<hash>.json (raw responses, for audit)

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
PAPERS = ROOT / "data" / "papers"
CUR = ROOT / "data" / "curated"
CACHE = ROOT / "data" / "raw" / "vv_cache"
CACHE.mkdir(parents=True, exist_ok=True)
VV = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{}/all?content-type=application/json"
ENS = "https://rest.ensembl.org"


def norm(s: str) -> str:
    s = s.replace("−", "-").replace("–", "-").replace("‑", "-")
    s = s.replace("‘", "'").replace("’", "'").replace("′", "'")
    s = s.replace("“", '"').replace("”", '"').replace(" ", " ")
    s = s.replace("\ufb01", "fi").replace("\ufb02", "fl")  # PDF ligatures
    return re.sub(r"\s+", " ", s)


def source_text(pmid: str) -> str:
    parts = []
    for name in (f"{pmid}.fulltext.txt", f"{pmid}.abstract.txt", f"{pmid}.plain.txt", f"{pmid}.layout.txt", f"{pmid}.s2.txt"):
        p = PAPERS / name
        if p.exists():
            parts.append(p.read_text())
    return norm("\n".join(parts))


def cached_get(url: str, headers: dict | None = None) -> dict | list:
    key = CACHE / (hashlib.sha1(url.encode()).hexdigest() + ".json")
    if key.exists():
        return json.loads(key.read_text())
    for attempt in range(5):
        r = requests.get(url, headers=headers or {}, timeout=120)
        if r.status_code == 429:
            time.sleep(5 * (attempt + 1))
            continue
        r.raise_for_status()
        data = r.json()
        key.write_text(json.dumps(data))
        time.sleep(1.0)  # be polite to the public service
        return data
    raise RuntimeError(f"rate limited: {url}")


def resolve_rsid(rsid: str) -> str:
    """rsID -> NM_001009944.3 c. notation via Ensembl VEP (RefSeq transcripts)."""
    data = cached_get(f"{ENS}/vep/human/id/{rsid}?refseq=1&hgvs=1&content-type=application/json")
    hits = set()
    for rec in data:
        for tc in rec.get("transcript_consequences", []):
            h = tc.get("hgvsc", "")
            if h.startswith("NM_001009944.3:") or h.startswith("NM_000297.4:"):
                hits.add(h.split(":", 1)[1])
    if len(hits) != 1:
        raise ValueError(f"{rsid}: expected one RefSeq c. notation, got {sorted(hits)}")
    return hits.pop()


def validate(transcript: str, hgvs_c: str) -> dict:
    q = urllib.parse.quote(f"{transcript}:{hgvs_c}", safe="")
    data = cached_get(VV.format(q))
    out = {"vv_ok": False, "vv_hgvs": "", "vv_warnings": "", "chrom": "", "pos": "",
           "ref": "", "alt": "", "hgvs_g": ""}
    recs = {k: v for k, v in data.items() if k not in ("flag", "metadata")}
    warns = []
    for k, v in recs.items():
        warns += v.get("validation_warnings", []) or []
        loci = (v.get("primary_assembly_loci") or {}).get("grch38") or {}
        vcf = loci.get("vcf") or {}
        if vcf:
            out.update(vv_hgvs=v.get("hgvs_transcript_variant", ""),
                       chrom=vcf.get("chr", ""), pos=vcf.get("pos", ""),
                       ref=vcf.get("ref", ""), alt=vcf.get("alt", ""),
                       hgvs_g=loci.get("hgvs_genomic_description", ""))
    out["vv_warnings"] = " | ".join(dict.fromkeys(warns))
    bad = re.search(r"does not agree with reference|not match|invalid|Invalid|out of bounds|"
                    r"not in the reference", out["vv_warnings"])
    out["vv_ok"] = bool(out["pos"]) and data.get("flag") == "gene_variant" and not bad
    out["vv_flag"] = data.get("flag", "")
    return out


def main() -> None:
    df = pd.read_csv(CUR / "truthset_candidates.tsv", sep="\t", dtype=str).fillna("")
    rows = []
    for r in df.to_dict("records"):
        r["quote_ok"] = norm(r["quote"]) in source_text(r["pmid"])
        hg = r["hgvs_c"]
        if hg.startswith("rs"):
            r["rsid"] = hg
            hg = resolve_rsid(hg)
            r["hgvs_c"] = hg
        r.update(validate(r["transcript"], hg))
        r["verified"] = r["quote_ok"] and r["vv_ok"]
        rows.append(r)
        print(f"{r['entry_id']} {r['gene']} {hg:32s} quote={'OK ' if r['quote_ok'] else 'FAIL'} "
              f"vv={'OK ' if r['vv_ok'] else 'FAIL'} {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} "
              f"{('WARN: ' + r['vv_warnings'][:160]) if r['vv_warnings'] else ''}")
    out = pd.DataFrame(rows)
    out.to_csv(CUR / "truthset_verified.tsv", sep="\t", index=False)
    print(f"\nverified {out.verified.sum()} / {len(out)}; "
          f"quote failures {(~out.quote_ok).sum()}; VV failures {(~out.vv_ok).sum()}")


if __name__ == "__main__":
    main()
