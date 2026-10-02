#!/usr/bin/env python3
"""
06_score_atlas.py

Primary endpoint of docs/ANALYSIS_PLAN.md. Pulls precomputed AlphaGenome Atlas
scores:

  1. every possible SNV across each gene locus (gene span + 5 kb upstream,
     all three alternate alleles at every position): AVI_SCORE only. This is
     the background each truth-set variant is ranked against.
  2. every truth-set SNV: AVI_SCORE plus the 18 AVI feature attributions.

The Atlas holds SNVs, plus indels observed in gnomAD, UK Biobank and All of Us.
The truth set's indels and the dinucleotide allele are private variants, so
they are not expected here; 07_score_ondemand.py scores those.

Output: results/atlas_locus_<gene>.tsv.gz, results/atlas_truthset.tsv,
        results/atlas_log.json

AlphaGenome Output. Non-commercial use only; see LEGALLY_BINDING_TERMS_OF_USE.txt.

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import grpc
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CUR = ROOT / "data" / "curated"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

# GRCh38, Ensembl lookup 19 Sep 2026. PKD1 is on the minus strand, so its
# upstream 5 kb lies at HIGHER coordinates; PKD2 is on the plus strand.
CHUNK = 32          # server-side interval size used by the Atlas client
BASE_DELAY = 0.5    # seconds between requests; doubled on quota errors

LOCI = {
    "PKD1": ("chr16", 2088708, 2135898 + 5000),
    "PKD2": ("chr4", 88007586 - 5000, 88077777),
}


def ad_to_frame(ad) -> pd.DataFrame:
    x = np.asarray(ad.X)
    cols = list(ad.var["name"]) if "name" in ad.var else [str(c) for c in ad.var.index]
    df = pd.DataFrame(x, columns=cols)
    df.insert(0, "variant", [str(v) for v in ad.obs["variant"]])  # objects on some code paths
    return df


def main() -> None:
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--truthset-only", action="store_true",
                    help="skip the locus backgrounds; fetch only the truth-set SNVs")
    args = ap.parse_args()

    print("authenticating (key not logged)")
    client = atlas.create(get_key(), timeout=60)
    log = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "loci": {}}

    # 1. locus backgrounds. The Atlas serves 32 bp per request and enforces a
    # per-minute quota, and the client's own query_interval fans out ten
    # parallel workers, which trips the quota at once. So walk the locus
    # ourselves, one chunk at a time, checkpointing every chunk to disk so an
    # interrupted run resumes where it stopped, and back off on quota errors.
    for gene, (chrom, start, end) in ({} if args.truthset_only else LOCI).items():
        out = RES / f"atlas_locus_{gene}.tsv.gz"
        if out.exists():
            print(f"{gene}: cached")
            continue
        ckpt = RES / "atlas_chunks" / gene
        ckpt.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        starts = list(range(start - 1, end, CHUNK))
        done = 0
        for s0 in starts:
            f = ckpt / f"{s0}.tsv"
            if f.exists():
                done += 1
                continue
            iv = genome.Interval(chromosome=chrom, start=s0, end=min(s0 + CHUNK, end))
            delay = BASE_DELAY
            for attempt in range(12):
                try:
                    res = client.query_interval(iv, requested_scorers=["AVI_SCORE"],
                                                progress_bar=False, max_workers=1)
                    break
                except grpc.RpcError as e:  # quota and transient errors surface raw
                    if e.code() in (grpc.StatusCode.RESOURCE_EXHAUSTED, grpc.StatusCode.UNAVAILABLE):
                        time.sleep(min(60, delay)); delay *= 2
                        continue
                    raise
            else:
                raise RuntimeError(f"{gene} chunk {s0}: gave up after repeated quota errors")
            ad_to_frame(res["AVI_SCORE"]).rename(columns={"AVI_SCORE": "avi"}).to_csv(f, sep="\t", index=False)
            done += 1
            if done % 100 == 0:
                print(f"  {gene}: {done}/{len(starts)} chunks, {time.time()-t0:.0f}s", flush=True)
            time.sleep(BASE_DELAY)
        df = pd.concat([pd.read_csv(ckpt / f"{s0}.tsv", sep="\t") for s0 in starts], ignore_index=True)
        df = df.drop_duplicates("variant")
        df.to_csv(out, sep="\t", index=False, compression="gzip")
        expected = (end - start + 1) * 3
        log["loci"][gene] = {"interval": f"{chrom}:{start}-{end}", "n": len(df),
                             "expected_snvs": expected, "seconds": round(time.time() - t0, 1)}
        print(f"{gene}: {len(df):,} variants (expected ~{expected:,} SNVs) in {time.time()-t0:.0f}s")

    # 2. truth-set SNVs with attributions
    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"]
    snv = ts[(ts.ref.str.len() == 1) & (ts.alt.str.len() == 1)].drop_duplicates(["chrom", "pos", "ref", "alt"])
    variants = [genome.Variant(chromosome=f"chr{r.chrom}", position=int(r.pos),
                               reference_bases=r.ref, alternate_bases=r.alt)
                for r in snv.itertuples()]
    res = client.query_variants(variants, requested_scorers=["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"],
                                progress_bar=False, max_workers=1)
    avi = ad_to_frame(res["AVI_SCORE"]).rename(columns={"AVI_SCORE": "avi"})
    fi = ad_to_frame(res["AVI_SCORE_FEATURE_IMPORTANCE"])
    fi.columns = ["variant"] + [f"fi_{c}" for c in fi.columns[1:]]
    merged = avi.merge(fi, on="variant", how="outer")
    merged.to_csv(RES / "atlas_truthset.tsv", sep="\t", index=False)
    log["truthset_snvs_requested"] = len(variants)
    log["truthset_snvs_returned"] = len(merged)
    log["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (RES / ("atlas_truthset_log.json" if args.truthset_only else "atlas_log.json")).write_text(
        json.dumps(log, indent=2))
    print(f"truth-set SNVs: requested {len(variants)}, returned {len(merged)}")


if __name__ == "__main__":
    main()
