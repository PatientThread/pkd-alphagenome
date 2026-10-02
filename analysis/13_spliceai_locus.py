#!/usr/bin/env python3
"""
13_spliceai_locus.py

SpliceAI for EVERY possible SNV across the PKD1 (then PKD2) locus, so SpliceAI
can be ranked on the same within-gene percentile scale as AVI and CADD
(docs/ANALYSIS_PLAN.md section 5: "if it is feasible, the percentile is
reported too").

The stock spliceai package scores one variant at a time (two sequences, five
models each), which would take days for 156,573 variants. This script runs the
SAME arithmetic as spliceai.utils.get_delta_scores for SNVs (distance 500,
unmasked, same gene-boundary N-padding, same strand handling), but predicts in
batches, and reuses one reference prediction for the three alternate alleles
at each position.

EQUIVALENCE IS TESTED, NOT ASSUMED. `--validate` scores the truth-set SNVs and
compares them with results/spliceai_scores.tsv from the stock code
(08_score_spliceai.py). The locus run refuses to start unless that validation
file exists and every difference is within rounding (0.005).

Variants outside the SpliceAI annotation of the gene (for example the 5 kb
upstream) get no score from the stock tool. They are recorded with
in_annotation = False and delta 0 for ranking, which is what a user of the
stock tool would see (no predicted effect). This convention is stated in the
output.

Output: results/spliceai_locus_<gene>.tsv.gz (chunked checkpoints in
        results/spliceai_chunks/<gene>/), results/spliceai_validation.tsv

Author: Christopher Lawrence
"""
from __future__ import annotations

import os

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")

import argparse  # noqa: E402
import logging  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import spliceai.utils as su  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FASTA = Path.home() / ".genomes" / "hg38" / "hg38_chr4_chr16.fa"
DIST = 500
COV = 2 * DIST + 1
WID = 10000 + COV
LOCI = {"PKD1": ("chr16", 2088708, 2135898 + 5000), "PKD2": ("chr4", 88007586 - 5000, 88077777)}
POS_PER_BATCH = 64          # positions per batch -> 256 sequences (ref + 3 alts)
CHUNK_POS = 1000            # positions per checkpoint file
_MAP = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.float32)


def onehot(seq: str) -> np.ndarray:
    s = (seq.upper().replace("A", "\x01").replace("C", "\x02")
         .replace("G", "\x03").replace("T", "\x04").replace("N", "\x00"))
    return _MAP[np.frombuffer(s.encode("latin-1"), np.int8) % 5]


def build(ann, gene: str, chrom: str, pos: int):
    """Return (ref_window, strand) exactly as the stock code builds it, or None."""
    genes, strands, idxs = ann.get_name_and_strand(chrom, pos)
    genes = list(genes)
    if gene not in genes:
        return None
    i = genes.index(gene)
    seq = ann.ref_fasta[chrom][pos - WID // 2 - 1:pos + WID // 2].seq
    if len(seq) != WID:
        return None
    d = ann.get_pos_data(idxs[i], pos)
    pad = [max(WID // 2 + d[0], 0), max(WID // 2 - d[1], 0)]
    x_ref = "N" * pad[0] + seq[pad[0]:WID - pad[1]] + "N" * pad[1]
    return x_ref, strands[i]


def predict(models, X: np.ndarray) -> np.ndarray:
    return np.mean([m.predict(X, batch_size=X.shape[0], verbose=0) for m in models], axis=0)


def score_positions(ann, gene: str, chrom: str, positions: list[int]) -> list[dict]:
    rows, batch, meta = [], [], []

    def flush():
        if not batch:
            return
        X = np.stack(batch)
        Y = predict(ann.models, X)
        for k, (pos, ref, alts, strand) in enumerate(meta):
            yr = Y[4 * k]
            for j, alt in enumerate(alts):
                ya = Y[4 * k + 1 + j]
                if strand == "-":
                    a, r = ya[::-1], yr[::-1]
                else:
                    a, r = ya, yr
                ag = float((a[:, 1] - r[:, 1]).max()); al = float((r[:, 1] - a[:, 1]).max())
                dg = float((a[:, 2] - r[:, 2]).max()); dl = float((r[:, 2] - a[:, 2]).max())
                rows.append(dict(chrom=chrom, pos=pos, ref=ref, alt=alt, in_annotation=True,
                                 ag=ag, al=al, dg=dg, dl=dl, delta_max=max(ag, al, dg, dl)))
        batch.clear(); meta.clear()

    for pos in positions:
        b = build(ann, gene, chrom, pos)
        ref = ann.ref_fasta[chrom][pos - 1:pos].seq.upper()
        alts = [x for x in "ACGT" if x != ref]
        if b is None or ref not in "ACGT":
            for alt in alts:
                rows.append(dict(chrom=chrom, pos=pos, ref=ref, alt=alt, in_annotation=False,
                                 ag=0.0, al=0.0, dg=0.0, dl=0.0, delta_max=0.0))
            continue
        x_ref, strand = b
        seqs = [x_ref] + [x_ref[:WID // 2] + alt + x_ref[WID // 2 + 1:] for alt in alts]
        enc = [onehot(s) for s in seqs]
        if strand == "-":
            enc = [e[::-1, ::-1] for e in enc]
        batch.extend(enc)
        meta.append((pos, ref, alts, strand))
        if len(meta) >= POS_PER_BATCH:
            flush()
    flush()
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--genes", default="PKD1,PKD2")
    ap.add_argument("--sample", type=int, default=0,
                    help="score an evenly spaced sample of this many positions per locus instead of "
                         "every position (background distribution only; truth-set variants are always "
                         "scored exactly, by 08_score_spliceai.py)")
    args = ap.parse_args()
    logging.getLogger().setLevel(logging.ERROR)
    su.one_hot_encode = onehot
    ann = su.Annotator(str(FASTA), "grch38")

    if args.validate:
        stock = pd.read_csv(RES / "spliceai_scores.tsv", sep="\t", dtype={"chrom": str})
        stock = stock[(stock.ref.str.len() == 1) & (stock.alt.str.len() == 1) & stock.spliceai_max.notna()]
        out = []
        for gene, g in stock.groupby("gene"):
            chrom = f"chr{g.chrom.iloc[0]}"
            rows = pd.DataFrame(score_positions(ann, gene, chrom, sorted(set(g.pos.astype(int)))))
            m = g.merge(rows, left_on=["pos", "ref", "alt"], right_on=["pos", "ref", "alt"], how="left")
            out.append(m)
        v = pd.concat(out)
        v["abs_diff"] = (v.spliceai_max - v.delta_max.round(2)).abs()
        v.to_csv(RES / "spliceai_validation.tsv", sep="\t", index=False)
        print(f"validated {len(v)} SNVs; max |stock - batched| = {v.abs_diff.max():.4f}; "
              f"all within 0.005: {bool((v.abs_diff <= 0.005).all())}")
        return

    val = RES / "spliceai_validation.tsv"
    if not val.exists():
        raise SystemExit("run --validate first")
    v = pd.read_csv(val, sep="\t")
    if not (v.abs_diff <= 0.005).all():
        raise SystemExit("validation failed: batched scores do not reproduce stock SpliceAI")

    for gene in args.genes.split(","):
        chrom, start, end = LOCI[gene]
        tag = f"sample{args.sample}" if args.sample else "locus"
        out = RES / (f"spliceai_locusbg_{gene}.tsv.gz" if args.sample else f"spliceai_locus_{gene}.tsv.gz")
        if out.exists():
            print(f"{gene}: cached"); continue
        ck = RES / f"spliceai_chunks_{tag}" / gene
        ck.mkdir(parents=True, exist_ok=True)
        if args.sample:
            step = max(1, (end - start + 1) // args.sample)
            allpos = list(range(start, end + 1, step))
        else:
            allpos = list(range(start, end + 1))
        starts = list(range(0, len(allpos), CHUNK_POS))
        t0 = time.time()
        for n, s0 in enumerate(starts, 1):
            f = ck / f"{s0}.tsv.gz"
            if f.exists():
                continue
            pos = allpos[s0:s0 + CHUNK_POS]
            pd.DataFrame(score_positions(ann, gene, chrom, pos)).to_csv(f, sep="\t", index=False, compression="gzip")
            el = time.time() - t0
            print(f"{gene} chunk {n}/{len(starts)}  {el/60:.1f} min elapsed", flush=True)
        df = pd.concat([pd.read_csv(ck / f"{s0}.tsv.gz", sep="\t") for s0 in starts], ignore_index=True)
        df.to_csv(out, sep="\t", index=False, compression="gzip")
        print(f"{gene}: {len(df):,} SNVs from {len(allpos):,} positions "
              f"(locus spans {end - start + 1:,} positions); in annotation {df.in_annotation.sum():,}")


if __name__ == "__main__":
    main()
