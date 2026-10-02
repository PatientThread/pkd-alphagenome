#!/usr/bin/env python3
"""
08_score_spliceai.py

SpliceAI comparator, run locally (the public SpliceAI-lookup service timed out
repeatedly on 19 Sep 2026). Settings fixed by docs/ANALYSIS_PLAN.md: distance
500, unmasked. Reference: UCSC hg38 chr4 + chr16. Environment:
~/.venvs/spliceai (spliceai 1.3.1, TensorFlow 2.21, tf-keras).

TWO DOCUMENTED DEVIATIONS FROM THE STOCK PACKAGE, both reproducing its logic:

1. NumPy 2 removed np.fromstring, which spliceai.utils.one_hot_encode uses. It
   is replaced with a byte-identical np.frombuffer version.
2. The stock package refuses any variant where ref and alt are both longer than
   one base, so it returns nothing for the dinucleotide allele
   c.7066-131_7066-130delinsAG. For an EQUAL-LENGTH substitution the stock
   code path needs no length adjustment, so the only change is to let equal-
   length multi-base substitutions through that guard. Unequal complex
   variants are still refused. Rows scored this way are flagged
   spliceai_mnv_extension = True.

VALIDATION. Before use, this installation reproduced the published SpliceAI
scores for c.7066-131G>A (acceptor gain 0.12, donor gain 0.06; PMID 42502691).

Output: results/spliceai_scores.tsv

Author: Christopher Lawrence
"""
from __future__ import annotations

import os

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")

import logging  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import spliceai.utils as su  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CUR = ROOT / "data" / "curated"
RES = ROOT / "results"
FASTA = Path.home() / ".genomes" / "hg38" / "hg38_chr4_chr16.fa"
DIST, MASK = 500, 0

_MAP = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def one_hot_encode(seq: str) -> np.ndarray:
    seq = (seq.upper().replace("A", "\x01").replace("C", "\x02")
           .replace("G", "\x03").replace("T", "\x04").replace("N", "\x00"))
    return _MAP[np.frombuffer(seq.encode("latin-1"), np.int8) % 5]


su.one_hot_encode = one_hot_encode


class Rec:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom, self.pos, self.ref, self.alts = chrom, pos, ref, [alt]


def delta(rec: Rec, ann) -> str:
    """Stock get_delta_scores, except equal-length multi-base substitutions are allowed."""
    if len(rec.ref) > 1 and len(rec.alts[0]) > 1 and len(rec.ref) == len(rec.alts[0]):
        # Present it to the stock function as if it were a 1-base-alt record is
        # NOT equivalent, so run the stock arithmetic directly on the substitution.
        cov, wid = 2 * DIST + 1, 10000 + 2 * DIST + 1
        genes, strands, idxs = ann.get_name_and_strand(rec.chrom, rec.pos)
        chrom = su.normalise_chrom(rec.chrom, list(ann.ref_fasta.keys())[0])
        seq = ann.ref_fasta[chrom][rec.pos - wid // 2 - 1:rec.pos + wid // 2].seq
        assert seq[wid // 2:wid // 2 + len(rec.ref)].upper() == rec.ref, "reference mismatch"
        out = []
        for i in range(len(idxs)):
            d = ann.get_pos_data(idxs[i], rec.pos)
            pad = [max(wid // 2 + d[0], 0), max(wid // 2 - d[1], 0)]
            x_ref = "N" * pad[0] + seq[pad[0]:wid - pad[1]] + "N" * pad[1]
            x_alt = x_ref[:wid // 2] + rec.alts[0] + x_ref[wid // 2 + len(rec.ref):]
            xr, xa = one_hot_encode(x_ref)[None, :], one_hot_encode(x_alt)[None, :]
            if strands[i] == "-":
                xr, xa = xr[:, ::-1, ::-1], xa[:, ::-1, ::-1]
            yr = np.mean([ann.models[m].predict(xr, verbose=0) for m in range(5)], axis=0)
            ya = np.mean([ann.models[m].predict(xa, verbose=0) for m in range(5)], axis=0)
            if strands[i] == "-":
                yr, ya = yr[:, ::-1], ya[:, ::-1]
            y = np.concatenate([yr, ya])
            ipa = (y[1, :, 1] - y[0, :, 1]).argmax(); ina = (y[0, :, 1] - y[1, :, 1]).argmax()
            ipd = (y[1, :, 2] - y[0, :, 2]).argmax(); ind = (y[0, :, 2] - y[1, :, 2]).argmax()
            out.append("{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}".format(
                rec.alts[0], genes[i],
                y[1, ipa, 1] - y[0, ipa, 1], y[0, ina, 1] - y[1, ina, 1],
                y[1, ipd, 2] - y[0, ipd, 2], y[0, ind, 2] - y[1, ind, 2],
                ipa - cov // 2, ina - cov // 2, ipd - cov // 2, ind - cov // 2))
        return ";".join(out)
    return ";".join(su.get_delta_scores(rec, ann, DIST, MASK))


def main() -> None:
    logging.getLogger().setLevel(logging.ERROR)
    ann = su.Annotator(str(FASTA), "grch38")
    for m in ann.models:  # silence per-call progress bars
        _p = m.predict
        m.predict = (lambda f: (lambda x, **k: f(x, verbose=0)))(_p)

    ts = pd.read_csv(CUR / "truthset_verified.tsv", sep="\t", dtype=str).fillna("")
    ts = ts[ts.verified == "True"].drop_duplicates(["chrom", "pos", "ref", "alt"])
    todo = ts[["entry_id", "gene", "hgvs_c", "label", "chrom", "pos", "ref", "alt"]].to_dict("records")
    todo += [dict(entry_id="X1", gene="PKD1", hgvs_c="c.7066-131G>A (alone)", label="NOT_IN_TRUTHSET",
                  chrom="16", pos="2107079", ref="C", alt="T"),
             dict(entry_id="X2", gene="PKD1", hgvs_c="c.2853+5G>A", label="ALT_READING_OF_E011",
                  chrom="16", pos="2114165", ref="C", alt="T")]
    rows = []
    for v in todo:
        rec = Rec(f"chr{v['chrom']}", int(v["pos"]), v["ref"], v["alt"])
        raw = delta(rec, ann)
        r = dict(v, spliceai_raw=raw,
                 spliceai_mnv_extension=len(v["ref"]) > 1 and len(v["alt"]) > 1 and len(v["ref"]) == len(v["alt"]))
        # keep the annotation for the variant's own gene
        best = None
        for part in raw.split(";"):
            f = part.split("|")
            if len(f) == 10 and f[1] == v["gene"] and f[2] != ".":
                best = f
        if best:
            ag, al, dg, dl = map(float, best[2:6])
            r.update(spliceai_ag=ag, spliceai_al=al, spliceai_dg=dg, spliceai_dl=dl,
                     spliceai_max=max(ag, al, dg, dl))
        rows.append(r)
        print(f"{v['entry_id']:5s} {v['gene']} {v['hgvs_c']:32s} {r.get('spliceai_max', 'NA')}", flush=True)
    pd.DataFrame(rows).to_csv(RES / "spliceai_scores.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
