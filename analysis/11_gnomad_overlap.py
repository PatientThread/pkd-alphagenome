#!/usr/bin/env python3
"""
11_gnomad_overlap.py

Amendment A2. AVI was trained on gnomAD v4.1 genome variants from chromosomes
1, 4, 7, 8, 10, 13, 15, and DeepMind's own evaluation dropped any benchmark
variant whose POSITION (or indel span) overlapped a training variant, allele
regardless. The same rule is applied here.

For every truth-set variant (and X1, X2) the gnomAD v4.1 GRCh38 region covering
the variant's reference span is queried. Recorded:
    gnomad_genome_overlap   any gnomAD v4.1 GENOME variant in the span
    gnomad_exact            the exact allele is in gnomAD v4.1 (genome or exome)
    gnomad_genome_af, gnomad_faf95_popmax for the exact allele where present
    avi_training_chrom      the chromosome is one of AVI's training chromosomes
    avi_leakage_flag        avi_training_chrom AND gnomad_genome_overlap

Only the GENOME dataset matters for the leakage rule, because AVI's training
pool was built from gnomAD v4.1 genomes. Exome data is recorded for context.

Output: data/curated/gnomad_overlap.tsv

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
CUR = ROOT / "data" / "curated"
API = "https://gnomad.broadinstitute.org/api"
AVI_TRAIN_CHROMS = {"1", "4", "7", "8", "10", "13", "15"}
AVI_VALID_CHROMS = {"2", "5", "11", "14", "17", "20", "22", "X"}

Q = """query ($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id pos
      genome { ac an af faf95 { popmax popmax_population } }
      exome { ac an af }
    }
  }
}"""


def query(chrom: str, start: int, stop: int) -> list[dict]:
    for attempt in range(8):
        r = requests.post(API, json={"query": Q, "variables": {"chrom": chrom, "start": start, "stop": stop}},
                          timeout=120)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(10 * (attempt + 1))
            continue
        r.raise_for_status()
        d = r.json()
        if "errors" in d:
            raise RuntimeError(json.dumps(d["errors"])[:400])
        return d["data"]["region"]["variants"]
    raise RuntimeError("gnomAD unavailable after retries")


def main() -> None:
    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"].drop_duplicates(["chrom", "pos", "ref", "alt"])
    todo = ts[["entry_id", "gene", "hgvs_c", "label", "chrom", "pos", "ref", "alt"]].to_dict("records")
    todo += [dict(entry_id="X1", gene="PKD1", hgvs_c="c.7066-131G>A (alone)", label="NOT_IN_TRUTHSET",
                  chrom="16", pos="2107079", ref="C", alt="T"),
             dict(entry_id="X2", gene="PKD1", hgvs_c="c.2853+5G>A", label="ALT_READING_OF_E011",
                  chrom="16", pos="2114165", ref="C", alt="T")]
    rows = []
    for v in todo:
        start = int(v["pos"])
        stop = start + len(v["ref"]) - 1
        hits = query(v["chrom"], start, stop)
        vid = f"{v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}"
        exact = next((h for h in hits if h["variant_id"] == vid), None)
        g_overlap = any(h.get("genome") for h in hits)
        g = (exact or {}).get("genome") or {}
        e = (exact or {}).get("exome") or {}
        rec = dict(v, gnomad_n_in_span=len(hits), gnomad_genome_overlap=g_overlap,
                   gnomad_exact=exact is not None,
                   gnomad_genome_af=g.get("af"),
                   gnomad_faf95_popmax=(g.get("faf95") or {}).get("popmax"),
                   gnomad_exome_af=e.get("af"),
                   avi_training_chrom=v["chrom"] in AVI_TRAIN_CHROMS,
                   avi_validation_chrom=v["chrom"] in AVI_VALID_CHROMS)
        rec["avi_leakage_flag"] = rec["avi_training_chrom"] and g_overlap
        rows.append(rec)
        print(f"{v['entry_id']:5s} {v['gene']} {v['hgvs_c']:32s} span_hits={len(hits)} "
              f"genome_overlap={g_overlap!s:5s} exact={rec['gnomad_exact']!s:5s} "
              f"af={rec['gnomad_genome_af']} leak={rec['avi_leakage_flag']}")
        time.sleep(0.6)
    pd.DataFrame(rows).to_csv(CUR / "gnomad_overlap.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
