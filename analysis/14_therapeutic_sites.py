#!/usr/bin/env python3
"""
14_therapeutic_sites.py

Exploratory aim A4 (docs/ANALYSIS_PLAN.md). Two PKD1 regulatory sites that are
current drug targets, both shown experimentally to RAISE polycystin-1 when
disrupted:

  miR-17 seed match, 3'UTR c.*153_*158, chr16:2089569-2089574 (fwd AAAGTG)
      target of anti-miR-17 (farabursen) and steric-blocking oligonucleotides
  uORF1 / uORF2 start codons, 5'UTR c.-87 and c.-20 (chr16:2135776, 2135709)
      target of uORF-blocking oligonucleotides

Q1  gnomAD v4.1 variants at the six seed bases and the two uORF ATGs.
Q2  Atlas AVI for every SNV at those bases, ranked within its own UTR.
Q3  AlphaGenome on-demand PKD1 RNA-seq prediction for the full human 6-base
    seed edit (AAAGTG -> GGGACA, forward strand).
Q4  The same for the uORF start-loss edits c.-87A>T and c.-20A>T.

Output: results/therapeutic_sites.json

AlphaGenome Output is non-commercial; see LEGALLY_BINDING_TERMS_OF_USE.txt.

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "results"
sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

SEED = list(range(2089569, 2089575))
UORF_ATG = {"uORF1": [2135776, 2135775, 2135774], "uORF2": [2135709, 2135708, 2135707]}
# PKD1 is minus strand: 3'UTR = genomic 2088708 .. (CDS end), 5'UTR = (CDS start) .. 2135898
UTR3 = (2088708, 2089726)   # c.*1019 .. c.*1 (c.*153 = 2089574, confirmed by VariantValidator)
UTR5 = (2135690, 2135898)   # c.-1 .. c.-209
GQ = """query ($start: Int!, $stop: Int!) { region(chrom: "16", start: $start, stop: $stop,
  reference_genome: GRCh38) { variants(dataset: gnomad_r4) { variant_id pos
  genome { ac an af } exome { ac an af } } } }"""


def gnomad(start: int, stop: int) -> list[dict]:
    for attempt in range(8):
        r = requests.post("https://gnomad.broadinstitute.org/api",
                          json={"query": GQ, "variables": {"start": start, "stop": stop}}, timeout=120)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(10 * (attempt + 1)); continue
        r.raise_for_status()
        return r.json()["data"]["region"]["variants"]
    raise RuntimeError("gnomAD unavailable")


def main() -> None:
    out: dict = {}
    # Q1
    out["Q1_gnomad_seed"] = gnomad(SEED[0], SEED[-1])
    out["Q1_gnomad_uorf_atg"] = {k: gnomad(min(v), max(v)) for k, v in UORF_ATG.items()}

    # Q2: AVI from the complete locus background
    loc = pd.read_csv(RES / "atlas_locus_PKD1.tsv.gz", sep="\t")
    loc["pos"] = loc.variant.str.split(":").str[1].astype(int)
    def rank_in(region, positions):
        bg = loc[(loc.pos >= region[0]) & (loc.pos <= region[1])]
        sel = bg[bg.pos.isin(positions)].copy()
        sel["pct_within_utr"] = sel.avi.map(lambda x: round(100 * (bg.avi < x).mean(), 1))
        return {"utr_snvs": int(len(bg)), "utr_avi_median": round(float(bg.avi.median()), 3),
                "site_snvs": sel[["variant", "avi", "pct_within_utr"]].round(3).to_dict("records"),
                "site_median_pct": float(sel.pct_within_utr.median())}
    out["Q2_avi_seed_vs_3utr"] = rank_in(UTR3, SEED)
    out["Q2_avi_uorf_atg_vs_5utr"] = rank_in(UTR5, sum(UORF_ATG.values(), []))

    # Q3, Q4: on-demand PKD1 RNA-seq and splicing
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome
    client = dna_client.create(get_key())
    sc = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in ("RNA_SEQ", "SPLICE_SITE_USAGE")]
    edits = {"Q3_seed_6nt_edit": (2089569, "AAAGTG", "GGGACA"),
             "Q4_uORF1_c.-87A>T": (2135776, "T", "A"), "Q4_uORF2_c.-20A>T": (2135709, "T", "A")}
    kidney = {"cortex of kidney", "outer medulla of kidney", "kidney", "left kidney", "right kidney",
              "kidney epithelial cell", "kidney tubule cell", "epithelial cell of proximal tubule",
              "renal cortical epithelial cell", "glomerular visceral epithelial cell",
              "kidney capillary endothelial cell", "nephron progenitor cell"}
    for name, (pos, ref, alt) in edits.items():
        L = 1_048_576; st = pos - L // 2
        got = client.score_variant(interval=genome.Interval("chr16", st, st + L),
                                   variant=genome.Variant("chr16", pos, ref, alt), variant_scorers=sc)
        rec = {}
        for key, ad in zip(("rna_seq", "splice_site_usage"), got):
            rows = np.where(ad.obs["gene_id"].astype(str).str.split(".").str[0] == "ENSG00000008710")[0]
            X = np.asarray(ad.X)[rows, :]
            k = np.where(ad.var["biosample_name"].isin(kidney).values)[0]
            rec[key] = {"max_abs_all": round(float(np.abs(X).max()), 4),
                        "mean_all": round(float(X.mean()), 4),
                        "mean_kidney": round(float(X[:, k].mean()), 4) if len(k) else None,
                        "n_tracks": int(X.shape[1]), "n_kidney": int(len(k))}
        out[name] = rec
    (RES / "therapeutic_sites.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
