#!/usr/bin/env python3
"""
09_score_cadd.py

CADD v1.7 (GRCh38) comparator.

  1. Locus background: every scored SNV across the same PKD1 and PKD2 loci used
     for the AVI background (06), read by remote tabix from the official CADD
     whole-genome SNV file, so CADD can be ranked exactly like AVI.
  2. Truth set: SNVs are taken from that same file. Indels and the
     dinucleotide allele are requested from the CADD web API; CADD only holds
     precomputed indels observed in population data, so private indels are
     expected to come back empty and are recorded as not scored rather than
     imputed.

Output: results/cadd_locus_<gene>.tsv.gz, results/cadd_truthset.tsv

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import pysam
import requests

ROOT = Path(__file__).resolve().parent.parent
CUR = ROOT / "data" / "curated"
RES = ROOT / "results"
URL = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
API = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{}:{}_{}_{}"
LOCI = {"PKD1": ("16", 2088708, 2135898 + 5000), "PKD2": ("4", 88007586 - 5000, 88077777)}


def main() -> None:
    tbx = pysam.TabixFile(URL)
    for gene, (chrom, start, end) in LOCI.items():
        out = RES / f"cadd_locus_{gene}.tsv.gz"
        if out.exists():
            continue
        rows = [line.split("\t") for line in tbx.fetch(chrom, start - 1, end)]
        df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "raw", "phred"])
        df.to_csv(out, sep="\t", index=False, compression="gzip")
        print(f"{gene}: {len(df):,} CADD SNVs (expected ~{(end - start + 1) * 3:,})")

    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"].drop_duplicates(["chrom", "pos", "ref", "alt"])
    extra = pd.DataFrame([
        dict(entry_id="X1", gene="PKD1", hgvs_c="c.7066-131G>A (alone)", label="NOT_IN_TRUTHSET",
             chrom="16", pos="2107079", ref="C", alt="T"),
        dict(entry_id="X2", gene="PKD1", hgvs_c="c.2853+5G>A", label="ALT_READING_OF_E011",
             chrom="16", pos="2114165", ref="C", alt="T")])
    todo = pd.concat([ts[["entry_id", "gene", "hgvs_c", "label", "chrom", "pos", "ref", "alt"]], extra])
    res = []
    for v in todo.to_dict("records"):
        r = dict(v, cadd_raw="", cadd_phred="", cadd_source="")
        if len(v["ref"]) == 1 and len(v["alt"]) == 1:
            p = int(v["pos"])
            for line in tbx.fetch(v["chrom"], p - 1, p):
                f = line.split("\t")
                if f[2] == v["ref"] and f[3] == v["alt"]:
                    r.update(cadd_raw=f[4], cadd_phred=f[5], cadd_source="whole_genome_SNVs")
        else:
            try:
                got = requests.get(API.format(v["chrom"], v["pos"], v["ref"], v["alt"]), timeout=60).json()
                if got:
                    r.update(cadd_raw=got[0].get("RawScore", ""), cadd_phred=got[0].get("PHRED", ""),
                             cadd_source="api_precomputed_indel")
                else:
                    r["cadd_source"] = "not_scored_private_indel"
            except (requests.RequestException, json.JSONDecodeError) as exc:
                r["cadd_source"] = f"api_error:{type(exc).__name__}"
            time.sleep(1)
        res.append(r)
    out = pd.DataFrame(res)
    out.to_csv(RES / "cadd_truthset.tsv", sep="\t", index=False)
    print(out[["entry_id", "gene", "hgvs_c", "cadd_phred", "cadd_source"]].to_string(index=False))


if __name__ == "__main__":
    main()
