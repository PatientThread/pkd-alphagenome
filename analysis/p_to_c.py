"""
p_to_c.py  -  derive a c. SNV from a missense protein change, mechanically.

Reads the coding sequence from the RefSeq GenBank record, checks the stated
reference amino acid is at that codon, enumerates all nine single-base changes
of the codon and returns the c. notation ONLY if exactly one produces the
stated alternate amino acid. Otherwise raises, so nothing ambiguous is guessed.

Author: Christopher Lawrence
"""
from __future__ import annotations

import re
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
TX = {"PKD1": "NM_001009944.3", "PKD2": "NM_000297.4"}
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter"}
BASES = "TCAG"
CODE = dict(zip([a + b + c for a in BASES for b in BASES for c in BASES],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def cds(gene: str) -> str:
    gb = (RAW / f"{TX[gene]}.gb").read_text()
    s, e = map(int, re.search(r"^\s{5}CDS\s+(\d+)\.\.(\d+)", gb, re.M).groups())
    seq = "".join(re.sub(r"[\s\d]", "", l) for l in gb.split("ORIGIN", 1)[1].splitlines() if l.strip() and l.strip() != "//").upper()
    return seq[s - 1:e]


def p_to_c(gene: str, p: str) -> str:
    m = re.fullmatch(r"p\.([A-Z])(\d+)([A-Z*])", p)
    if not m:
        raise ValueError(f"not a one-letter missense: {p}")
    ref, n, alt = m.group(1), int(m.group(2)), m.group(3)
    c = cds(gene)
    codon = c[3 * (n - 1):3 * n]
    if CODE[codon] != ref:
        raise ValueError(f"{gene} {p}: codon {n} is {codon} = {CODE[codon]}, not {ref}")
    hits = []
    for i in range(3):
        for b in "ACGT":
            if b == codon[i]:
                continue
            mut = codon[:i] + b + codon[i + 1:]
            if CODE[mut] == alt:
                hits.append(f"c.{3 * (n - 1) + i + 1}{codon[i]}>{b}")
    if len(hits) != 1:
        raise ValueError(f"{gene} {p}: {len(hits)} possible SNVs {hits}; refusing to guess")
    return hits[0]
