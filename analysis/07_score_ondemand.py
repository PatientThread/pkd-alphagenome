#!/usr/bin/env python3
"""
07_score_ondemand.py

AlphaGenome scored on demand, for every verified truth-set variant including the
indels and the dinucleotide allele that the Atlas cannot hold, plus two extra
variants the analysis plan needs:

    X1  c.7066-131G>A scored ALONE (chr16:2107079 C>T). Never observed without
        rs3874658, so it is not in the truth set, but H3 compares it with the
        dinucleotide allele.
    X2  c.2853+5G>A (chr16:2114165 C>T). The source paper is internally
        inconsistent between G>A and G>C; both alleles are scored.

SCORERS (AlphaGenome recommended defaults, 1,048,576 bp context centred on the
variant): SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, RNA_SEQ.

Every summary is restricted to the variant's own gene row (PKD1 or PKD2), and
reported twice: across ALL tracks, and across KIDNEY tracks only. The kidney
biosample list is the curated list from the companion kidney-representation
work, matched by exact biosample_name (never substring: "renal" would pull in
adrenal gland).

Result columns are positional and each result carries its own .var, so tracks
are always selected from that result's own biosample_name column.

Output: results/ondemand_scores.tsv, results/ondemand_log.json

AlphaGenome Output. Non-commercial use only; see LEGALLY_BINDING_TERMS_OF_USE.txt.

Author: Christopher Lawrence
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CUR = ROOT / "data" / "curated"
RES = ROOT / "results"
sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

SEQ_1MB = 1_048_576
GENE_ID = {"PKD1": "ENSG00000008710", "PKD2": "ENSG00000118762"}
KIDNEY_BIOSAMPLES = {
    "cortex of kidney", "outer medulla of kidney", "kidney", "left kidney",
    "right kidney", "kidney epithelial cell", "kidney tubule cell",
    "epithelial cell of proximal tubule", "renal cortical epithelial cell",
    "glomerular visceral epithelial cell",
    "kidney capillary endothelial cell", "nephron progenitor cell",
    "renal cortex interstitium", "left renal cortex interstitium",
    "right renal cortex interstitium",
}
SCORER_KEYS = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "RNA_SEQ")
EXTRAS = [
    dict(entry_id="X1", gene="PKD1", hgvs_c="c.7066-131G>A (alone)", label="NOT_IN_TRUTHSET",
         chrom="16", pos="2107079", ref="C", alt="T"),
    dict(entry_id="X2", gene="PKD1", hgvs_c="c.2853+5G>A", label="ALT_READING_OF_E011",
         chrom="16", pos="2114165", ref="C", alt="T"),
]


def gene_rows(ad, gene: str) -> np.ndarray:
    ids = ad.obs["gene_id"].astype(str).str.split(".").str[0].values
    return np.where(ids == GENE_ID[gene])[0]


def summarise(ad, gene: str, prefix: str, rec: dict) -> None:
    # The official recommendation aggregates "across all tracks and genes"; we record that as
    # <prefix>_maxabs_anygene and keep the gene-restricted value for the locus-specific analyses.
    Xall = np.asarray(ad.X)
    rec[f"{prefix}_maxabs_anygene"] = float(np.nanmax(np.abs(Xall))) if Xall.size else float("nan")
    rows = gene_rows(ad, gene)
    if len(rows) == 0:
        rec[f"{prefix}_missing_gene"] = True
        return
    X = np.asarray(ad.X)[rows, :]
    names = ad.var["name"].astype(str).values
    a = np.abs(X)
    j = np.unravel_index(np.nanargmax(a), a.shape)
    rec[f"{prefix}_maxabs_all"] = float(a[j])
    rec[f"{prefix}_signed_at_max_all"] = float(X[j])
    rec[f"{prefix}_track_at_max_all"] = names[j[1]]
    if "biosample_name" in ad.var.columns:
        k = np.where(ad.var["biosample_name"].isin(KIDNEY_BIOSAMPLES).values)[0]
        rec[f"{prefix}_n_kidney_tracks"] = int(len(k))
        if len(k):
            ak = a[:, k]
            jk = np.unravel_index(np.nanargmax(ak), ak.shape)
            rec[f"{prefix}_maxabs_kidney"] = float(ak[jk])
            rec[f"{prefix}_signed_at_max_kidney"] = float(X[:, k][jk])
            rec[f"{prefix}_track_at_max_kidney"] = names[k[jk[1]]]


def main() -> None:
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"].drop_duplicates(["chrom", "pos", "ref", "alt"])
    todo = ts[["entry_id", "gene", "hgvs_c", "label", "chrom", "pos", "ref", "alt"]].to_dict("records") + EXTRAS
    out_path = RES / "ondemand_scores.tsv"
    done = set()
    rows = []
    if out_path.exists():
        prev = pd.read_csv(out_path, sep="\t", dtype=str)
        rows = prev.to_dict("records")
        # Resume by VARIANT, never by entry_id: entry_ids are positional and shift whenever an entry
        # is inserted into 04_build_truthset.py, which silently skipped four variants on 20 Sep 2026.
        done = {(r["chrom"], str(r["pos"]), r["ref"], r["alt"]) for r in rows}

    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in SCORER_KEYS]
    t0, fails = time.time(), []
    for v in todo:
        if (v["chrom"], str(v["pos"]), v["ref"], v["alt"]) in done:
            continue
        pos = int(v["pos"])
        start = max(0, pos - SEQ_1MB // 2)
        iv = genome.Interval(chromosome=f"chr{v['chrom']}", start=start, end=start + SEQ_1MB)
        var = genome.Variant(chromosome=f"chr{v['chrom']}", position=pos,
                             reference_bases=v["ref"], alternate_bases=v["alt"])
        got = None
        for attempt in range(6):
            try:
                got = client.score_variant(interval=iv, variant=var, variant_scorers=scorers)
                break
            except Exception as exc:  # noqa: BLE001  quota / transient
                if attempt == 5:
                    fails.append((v["entry_id"], f"{type(exc).__name__}: {str(exc)[:120]}"))
                time.sleep(5 * (attempt + 1))
        if got is None:
            continue
        rec = dict(v)
        for key, ad in zip(SCORER_KEYS, got):
            summarise(ad, v["gene"], key.lower(), rec)
        rows.append(rec)
        pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)  # checkpoint every variant
        time.sleep(0.3)
    log = {"variants": len(todo), "scored": len(rows), "failures": fails,
           "seconds": round(time.time() - t0, 1),
           "finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (RES / "ondemand_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))
    if fails:
        raise SystemExit(f"{len(fails)} variants failed; rerun to retry (completed rows are kept)")


if __name__ == "__main__":
    main()
