#!/usr/bin/env python3
"""
04_build_truthset.py

The curated truth set. Every entry was read from a primary source in
data/papers/ (full text where open access, otherwise the PubMed abstract) and
carries a short VERBATIM quote from that source. 05_verify.py rejects any row
whose quote is not literally present in the source text, so a curation error
cannot pass silently.

LABELS (what the experiment showed, never what a predictor said):
    POS        functional assay showed an effect: aberrant splicing, or a change in transcript level
               or translation. The direction is recorded in the mechanism field, because two
               translation-class variants INCREASE output rather than reduce it.
    NEG        tested in a functional assay and indistinguishable from wild type
    INDET      tested, but the assay could not be interpreted
Variants that were named but never functionally tested are NOT in this table.
They are listed in EXCLUDED at the bottom with the reason.

ORIGIN distinguishes a variant carried by a patient from one made in the lab
by site-directed mutagenesis. Both are valid molecular truth; only patient
variants can speak to disease.

Legacy "IVS" notation is converted mechanically from the RefSeq GenBank exon
table (data/raw/<transcript>.gb), never by hand.

Output: data/curated/truthset_candidates.tsv, data/curated/excluded.tsv

Author: Christopher Lawrence
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent))
from p_to_c import p_to_c  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
CUR = ROOT / "data" / "curated"

TX = {"PKD1": "NM_001009944.3", "PKD2": "NM_000297.4"}


def exon_table(acc: str) -> tuple[list[tuple[int, int]], int]:
    """Exon (start, end) in mRNA coordinates, and CDS start, from GenBank."""
    gb = (RAW / f"{acc}.gb").read_text()
    exons = [(int(a), int(b)) for a, b in re.findall(r"^\s{5}exon\s+(\d+)\.\.(\d+)", gb, re.M)]
    cds = int(re.search(r"^\s{5}CDS\s+(\d+)\.\.", gb, re.M).group(1))
    return exons, cds


def ivs_to_c(gene: str, intron: int, offset: int, ref: str, alt: str) -> str:
    """IVS<n>+k -> c.<last base of exon n>+k; IVS<n>-k -> c.<first base of exon n+1>-k."""
    exons, cds = exon_table(TX[gene])
    if offset > 0:
        anchor = exons[intron - 1][1] - cds + 1
        return f"c.{anchor}+{offset}{ref}>{alt}"
    anchor = exons[intron][0] - cds + 1
    return f"c.{anchor}{offset}{ref}>{alt}"


# (gene, hgvs_c as published or converted, label, mechanism, assay, origin,
#  pmid, quote, note)
E = []


def add(gene, hgvs, label, mechanism, assay, origin, pmid, quote, note=""):
    E.append(dict(gene=gene, transcript=TX[gene], hgvs_c=hgvs, label=label,
                  mechanism=mechanism, assay=assay, origin=origin, pmid=pmid,
                  quote=quote, note=note))


# ---------------------------------------------------------------- 5'UTR
add("PKD1", "c.-69dupG", "POS", "translation (uORF1 extension)", "5'UTR Gaussia luciferase reporter normalised to secreted alkaline phosphatase; no RNA study",
    "patient", "41006799",
    "reduced the translation efficiency of the main PKD1 open reading frame by ~87% compared to wildtype",
    "Paper gives hg38 NC_000016.10:g.2135757_2135758insC. SpliceAI <0.05: not a splice variant.")
add("PKD1", "c.-52C>T", "NEG", "none", "5'UTR Gaussia luciferase reporter normalised to secreted alkaline phosphatase; no RNA study", "population", "41006799",
    "variant selected as a negative control (mean difference 0.008",
    "Common in gnomAD with homozygotes; measured equal to wild type.")
add("PKD1", "c.-209G>A", "NEG", "none", "5'UTR Gaussia luciferase reporter normalised to secreted alkaline phosphatase; no RNA study", "patient", "41006799",
    "did not find a significant difference between the NHS1 variant and wildtype",
    "Patient variant (Genomics England NHS1) that proved functionally inert.")

# ---------------------------------------------------------------- npj Genom Med 2023
add("PKD1", "c.11017-25A>G", "POS", "branchpoint loss; exon 38 skipping + partial intron 37 retention",
    "patient blood RT-PCR", "patient", "37419908",
    "a novel variant was identified in intron 37 of PKD1 (c.11017-25 A > G)",
    "Patient RBW403. Detected on two earlier diagnostic tests and called benign in silico.")
add("PKD1", "c.11017-25A>C", "POS", "branchpoint loss; exon 38 skipping + partial intron 37 retention",
    "patient blood RT-PCR", "patient", "37419908",
    "RPA028 had a different nucleotide substitution (PKD1 c.11017-25 A > C)",
    "Patient RPA028. v0.1 wrongly merged this with RBW403 as an 'allele conflict'.")
add("PKD1", "c.11017-10C>A", "POS", "AG-exclusion-zone acceptor disruption; exon 38 skipping",
    "patient blood RT-PCR", "patient", "37419908",
    "RPA014 had a variant 10 base pairs from the start of exon 38 (PKD1 c.11017-10 C > A)",
    "Independently reported co-segregating in a Chinese pedigree (PMID 36000373, no RNA).")
add("PKD1", "c.1991C>T", "POS", "new cryptic donor; in-frame deletion p.Ala664_Ser699del",
    "patient blood RT-PCR", "patient", "37419908",
    "Introme predicted that the c.1991C > T variant would generate a new cryptic donor site",
    "Exonic (missense by annotation) but acts through splicing.")
add("PKD1", "c.10167+25_10167+43del", "POS", "critical intron shortening",
    "patient blood RT-PCR", "patient", "37419908",
    "c.10167+25_10167+43del",
    "19-bp intronic deletion; mechanism is intron length, not a splice motif.")
add("PKD1", "c.7489+5G>A", "POS", "donor +5; 93-bp intron 18 retention (blood) / exon 18 skipping (minigene)",
    "patient blood RT-PCR; minigene", "patient", "37419908",
    "resulted in the retention of 93 base pairs of intron 18",
    "DISCORDANT MECHANISM across sources: PMID 42621998 minigene shows exon 18 skipping (c.7210_7489del).")

# ---------------------------------------------------------------- KI Reports 2025 minigene
add("PKD1", "c.2097+5G>A", "POS", "donor +5; 2-bp intron 10 retention", "minigene (pSPL3, HEK293T)",
    "patient", "41141524", "a 2 bp GT was confirmed to remain in the 5' splice site region of intron 10")
add("PKD1", "c.2853+5G>C", "POS", "donor +5; intron 11 retention and c.2707_2852 deletion",
    "minigene (pSPL3, HEK293T)", "patient", "41141524",
    "Minigene assay results of NM_001009944.3(PKD1):c.2853+5G>C",
    "SOURCE INCONSISTENT: Results text says c.2853+5G>A; both figure legends and Discussion say G>C. "
    "PMID 42621998 cites it as G>A. Score both alleles; resolve with authors.")
add("PKD1", "c.7490-3C>G", "POS", "acceptor -3; 93-bp intron 18 retention", "minigene (pSPL3, HEK293T)",
    "engineered", "41141524",
    "The c.7490-3C>G variant was generated by site-directed mutagenesis, as no patient-derived sample was available")
add("PKD1", "c.7490-11C>G", "POS", "acceptor -11; 93-bp intron 18 retention", "minigene (pSPL3, HEK293T)",
    "patient", "41141524", "Both variants (c.7490-3C>G and c.7490-11C>G) were 850 bp longer than the normal control")
add("PKD1", "c.288-12C>A", "INDET", "unclear", "minigene (pSPL3, HEK293T)", "patient", "41141524",
    "the splicing outcome of the remaining variant (c.288-12C>A) could not be clearly interpreted",
    "All in silico tools except NNSPLICE predicted aberrant splicing; the paper's comparison section "
    "calls the minigene pattern normal. Hold out of scoring as a calibration case.")

# ---------------------------------------------------------------- pseudoexon 2026
add("PKD1", "c.7066-131_7066-130delinsAG", "POS", "new donor; 114-bp pseudoexon in intron 16",
    "patient fibroblast RNA-seq", "patient", "42502691",
    "the NM_001009944.3:c.7066-131_7066-130delinsAG dinucleotide variant creates a GA dinucleotide",
    "De novo SNV on the haplotype of common rs3874658. SpliceAI: each SNV alone predicted no impact.")
add("PKD1", "c.7066-130A>G", "NEG", "none (alone)", "fibroblast RNA-seq of carrier father", "population",
    "42502691", "confirmed the activation of pseudoexon in the proband only, which is absent in the father who carried only the common sequence variant",
    "rs3874658, NC_000016.10:g.2107078T>C per paper. Common variant; functional negative when alone.")

# ---------------------------------------------------------------- other PKD1 RNA / minigene
add("PKD1", "c.12445-34_12445-10del", "POS", "two aberrant transcripts in exon 46", "minigene",
    "patient", "41659037", "two abnormal splicing alterations with the c.12445-34_12445-10del variant")
add("PKD1", "c.7489+5G>A", "POS", "exon 18 skipping", "minigene", "patient", "42621998",
    "the c.7489 + 5G>A variant induced exon 18 skipping",
    "Second, independent source for the same variant. This paper cites NM_001009944.2.")
add("PKD1", "c.11538-2A>G", "POS", "canonical acceptor; complete exon 42 skipping", "minigene",
    "patient", "41083429", "the c.11538-2A>G variant could cause complete exon 42 skipping")
add("PKD1", "c.9202-16G>A", "POS", "acceptor region; 60-base excision", "minigene", "patient",
    "41242437", "confirmed this variant causes aberrant splicing with a 60-base excision")
add("PKD1", "c.12444G>A", "POS", "last exonic base of exon 45; weakened donor, less normal PC1",
    "cDNA cloning and sequencing", "patient", "24575920",
    "PKD1:c.12444G > A mainly weakened the site and decreased the expression of normal PC1")
add("PKD1", "c.12444+1G>A", "POS", "canonical donor destroyed; novel donor used",
    "cDNA cloning and sequencing", "patient", "24575920",
    "PKD1:c.12444 + 1G > A definitely destroyed the native splice site and created a novel donor site")
add("PKD1", "c.11156G>A", "POS", "exonic; abolishes intron 38 donor", "minigene (kidney cell lines)",
    "engineered", "24907393", "The substitution c.11156G>A, previously predicted as missense mutation p.R3719Q, abolished the donor splice site of intron 38",
    "Reported in literature/databases; tested by site-directed mutagenesis.")
add("PKD1", "c.327A>T", "POS", "synonymous; creates strong exonic donor, incomplete exon 3",
    "minigene (kidney cell lines)", "engineered", "24907393",
    "Two synonymous variants, c.327A>T (p.G109G) and c.11257C>A (p.R3753R), generated strong donor splice sites")
add("PKD1", "c.11257C>A", "POS", "synonymous; creates strong exonic donor, incomplete exon 39",
    "minigene (kidney cell lines)", "engineered", "24907393",
    "Two synonymous variants, c.327A>T (p.G109G) and c.11257C>A (p.R3753R), generated strong donor splice sites")
add("PKD1", "c.8791+1_8791+5del", "POS", "canonical donor deletion; 8 bases inserted at exon 23 end",
    "patient RT-PCR", "patient", "27984604",
    "Sequencing of RNA has confirmed that there were 8 bases inserted in the 3' end of exon 23",
    "Published as c.8791+1_8791+5delGTGCG.")
add("PKD1", "c.2854-3C>G", "POS", "acceptor -3; 29 bases inserted at exon 11 end", "patient RT-PCR",
    "patient", "30424739", "mRNA sequencing showed that 29 bases inserted into the 3")
add("PKD1", "c.1717_1722+11del", "POS", "exon-intron boundary deletion; two abnormal transcripts",
    "patient cDNA sequencing", "patient", "37272738",
    "Two abnormal transcription products were formed")
add("PKD1", ivs_to_c("PKD1", 13, -2, "A", "T"), "POS", "canonical acceptor; cryptic acceptor in exon 14, 74-nt deletion",
    "patient lymphocyte RNA", "patient", "10612835",
    "The IVS13-2A>T substitution resulted in an inactivation of this splice site",
    "Converted from legacy IVS13-2A>T via RefSeq exon table.")
add("PKD1", "rs3874648", "POS", "expression modifier; Tra2-beta site, cryptic acceptor, less full-length PKD1",
    "patient blood RT-PCR; allele-specific expression; binding assays", "population", "36029107",
    "rs3874648-A allele increased Tra2-β binding affinity and activated a cryptic acceptor splice-site",
    "COMMON variant. c. notation resolved in 05 via Ensembl. Full text not open access.")

# BMC Genomics 2023: exonic variants, minigene; positives are statistically significant
# shifts that are often modest (e.g. 40.2% vs 32.9% exon exclusion)
for v in ["c.7866C>A", "c.7960A>G", "c.7979A>T", "c.7987C>T", "c.11248C>G", "c.11251C>T",
          "c.11257C>G", "c.11257C>T", "c.11346C>T", "c.11393C>G"]:
    add("PKD1", v, "POS", "exonic; partial exon skipping (modest)", "minigene (pSPL3)", "engineered",
        "37468838", "were identified to result in exon skipping",
        "Reported ADPKD variants, tested by site-directed mutagenesis. Effect sizes often small.")
add("PKD1", "c.7866C>T", "NEG", "none", "minigene (pSPL3)", "engineered", "37468838",
    "it was concluded that the variant c.7866C > T did not alter mRNA splicing",
    "Borderline: 45.3% vs 33.0% exon 21 skipping, P = 0.0618. Same position as POS c.7866C>A.")

# ---------------------------------------------------------------- PKD2
add("PKD2", "c.1480G>T", "POS", "nonsense; new alternative splicing + partial exon 6 skipping",
    "minigene (pSPL3)", "engineered", "37468838", "were identified to result in exon skipping")
for v in ["c.741C>G", "c.796G>T", "c.1546G>T"]:
    add("PKD2", v, "NEG", "none", "minigene (pSPL3)", "engineered", "37468838",
        "were considered to have no influence on pre-mRNA splicing")
add("PKD2", "c.1717-2A>G", "POS", "canonical acceptor; skipping of exons 8 and 9 together",
    "patient RT-PCR", "patient", "39200828",
    "revealing the simultaneous skipping of both exons 8 and 9 of the PKD2 gene")
add("PKD2", ivs_to_c("PKD2", 6, 5, "G", "C"), "POS", "donor +5; abnormal exon 6 splicing", "patient RT-PCR",
    "patient", "34749493", "The RT-PCR revealed the abnormal splicing of exon 6",
    "Converted from legacy IVS6+5G>C via RefSeq exon table.")
add("PKD2", ivs_to_c("PKD2", 8, 5, "G", "A"), "POS", "donor +5; six cryptic-site transcripts",
    "patient leukocyte RT-PCR", "patient", "19158373",
    "The data provide strong evidence that IVS8 + 5G-->A is a pathogenic mutation for PKD2",
    "Converted from legacy IVS8+5G>A via RefSeq exon table.")
add("PKD2", "c.1320G>T", "POS", "first base of exon 6; abolishes intron 5 acceptor", "patient leukocyte RNA",
    "patient", "20950398", "This mutation abolishes a conserved acceptor splice site of intron 5",
    "The same aberrant transcript occurs at low level in controls (paper's own caveat).")
add("PKD2", "c.1532A>T", "POS", "exonic; exon 6 skipping and truncated exon", "minigene", "engineered",
    "26692149", "Mutation c.1532A>T resulted in skipping of exon 6")
add("PKD2", "c.1716G>A", "POS", "synonymous, last base of exon 7; exon 7 skipping", "minigene; patient RT-PCR",
    "patient", "26692149", "c.1716G>A led to skipping of exon 7",
    "Independently confirmed in patient RNA by PMID 34542828.")
add("PKD2", "c.1716G>A", "POS", "exon 7 skipping", "patient RT-PCR", "patient", "34542828",
    "this heterozygous synonymous mutation led to exon7 skipping in PKD2 gene")
add("PKD2", "c.2657A>G", "POS", "exonic; incomplete exon 14", "minigene; patient lymphoblast RNA",
    "patient", "26692149", "Mutation c.2657A>G resulted in incorporation of an incomplete exon 14")
add("PKD2", "c.2020-5A>G", "POS", "acceptor -5; intron retention with PTC", "patient RT-PCR",
    "patient", "39169606", "NM_000297.4:c.2020-5A>G")


# ---------------------------------------------------------------- amendment A3: hand-downloaded papers
# PMID 34257392 (Xie 2021, J Hum Genet), Table 2: eight minigene positives. Variants were drawn from
# the ADPKD database and ClinVar and tested by site-directed mutagenesis, so origin = engineered.
XIE = [
    ("c.1202-8G>A", "c.1202-8 G > A PKD1 AL Insertion of 118 nucleotides prior to exon 6 (fs)"),
    ("c.7210-9T>A", "c.7210-9 T > A PKD1 AL Insertion of 7 nucleotides prior to exon 18 (fs)"),
    ("c.10407T>A", "c.10407 T > A PKD1 AL Insertion of 77 nucleotides prior to exon 34 (fs)"),
    ("c.11270-10T>A", "c.11270-10 T > A PKD1 AL 28 bp deletion in exon 40 (fs)"),
    ("c.11412-3C>G", "c.11412-3 C > G PKD1 AL Insertion of 140 nucleotides prior to exon 41 (fs)"),
    ("c.2985+5G>A", "c.2985 + 5 G > A PKD1 DL 130 bp deletion in exon 12 (fs)"),
    ("c.7209+5G>A", "c.7209 + 5 G > A PKD1 DL 69 bp deletion in exon 17"),
    ("c.9397+3G>C", "c.9397 + 3 G > C PKD1 DL 41 bp deletion in exon 26 (fs)"),
]
for hg, q in XIE:
    add("PKD1", hg, "POS", "splice-region; aberrant transcript in minigene", "minigene", "engineered",
        "34257392", q, "A3. Negatives from this paper (29) are in Supplementary Dataset S2, not yet obtained.")

# PMID 26692149 (Gonzalez-Paredes 2016, Gene), Table 1: PKD2 missense variants with effect "None".
for prot, q in [("p.R306Q", "p.R306Q 4 251 73 (5')"), ("p.R322W", "p.R322W 4 251 121 (5')"),
                ("p.R322Q", "p.R322Q 4 251 122 (5')"), ("p.R325Q", "p.R325Q 4 251 121 (3')"),
                ("p.A356P", "p.A356P 4 251 29 (3')"), ("p.G390S", "p.G390S 5 225 73 (5')"),
                ("p.W414G", "p.W414G 5 225 80 (3')"), ("p.C632R", "p.C632R 8 182 5 (3')"),
                ("p.R807Q", "p.R807Q 13 164 62 (5')"), ("p.R848Q", "p.R848Q 14 148 21 (5')")]:
    add("PKD2", p_to_c("PKD2", prot), "NEG", "none (missense)", "minigene (293T and COS7)", "engineered",
        "26692149", q, f"A3. Published as {prot}; c. derived by analysis/p_to_c.py (unique SNV). "
        "Table 1 'Effect observed in minigene assay' = None.")

# PMID 24907393 (Gonzalez-Paredes 2014, Gene): negatives named in the main text.
Q14 = ("included in the study for being relatively closed to exon ends, but not selected by ESEfinder "
       "or RESCUE-ESE")
for prot in ["p.M1092T", "p.R2765C", "p.H2921P", "p.R3277C", "p.V3285I", "p.T3382M", "p.A3391V"]:
    add("PKD1", p_to_c("PKD1", prot), "NEG", "none (missense)", "minigene (kidney cell lines)", "engineered",
        "24907393", Q14, f"A3. Published as {prot}; c. derived by analysis/p_to_c.py (unique SNV). "
        "The paper states these had no impact on pre-mRNA splicing.")
add("PKD1", "c.11258G>A", "NEG", "none (missense p.R3753Q)", "minigene (kidney cell lines)", "engineered",
    "24907393", "did not activate this cryptic donor site, and did not alter pre-mRNA splicing in our minigene analysis",
    "A3.")
add("PKD1", "c.11257C>T", "NEG", "none (missense p.R3753W)", "minigene (kidney cell lines)", "engineered",
    "24907393", "one affecting the same nucleotide c.11257C N T (p.R3753W)",
    "A3. CONFLICTS with PMID 37468838 (POS). Resolved as DISCORDANT by the A3 rule in 12_analyse.py. "
    "'N' is how the PDF renders '>'.")

# PMID 24907393 Supplementary Table S3 (image supplied by Dr Lawrence, transcribed to
# data/papers/24907393_TableS3_transcribed.tsv). Rows whose variant is unique by codon and whose
# reference amino acid is confirmed; the paper states all variants other than the three positives
# "did not affect splicing". Ambiguous rows (p.Q987H, p.L3675L) are excluded; rows already added
# from the main text are not repeated.
Q14S3 = "The rest of the mutations analyzed did not affect splicing"
S3_NOTES = {"p.C508R": "Source exon length 218 vs 221 in NM_001009944.3 (older reference); position from 3' end matches.",
            "p.V3375M": "Source distance 43 from 3' end vs 45 computed; consistent with an older reference."}
# Completed 20 Sep 2026 from the full Table S3 image (all 33 tested variants). The table's three
# positives (p.G109G, p.R3719Q, p.R3753R) are already in the set from the main text, and the derived
# c. notation for p.R3719Q, p.R3753Q and p.R3753W matches the main text exactly, which is a further
# check on the conversion method.
for prot in ["p.S75F", "p.C210G", "p.C508R", "p.L845S", "p.H1093Y", "p.R2761P", "p.L2763V", "p.R2791Q",
             "p.R2985G", "p.Q3016R", "p.R3105W", "p.V3138M", "p.N3188S", "p.V3375M",
             "p.M3678T", "p.L3682P", "p.L3731Q", "p.Q3751R"]:
    add("PKD1", p_to_c("PKD1", prot), "NEG", "none (missense)", "minigene (kidney cell lines)", "engineered",
        "24907393", Q14S3, f"A3, Table S3. Published as {prot}; c. derived by analysis/p_to_c.py (unique SNV); "
        "exon and position cross-checked against Table S3. " + S3_NOTES.get(prot, ""))

# PMID 34257392 Dataset S2 (supplied by Dr Lawrence; transcribed to
# data/papers/34257392_DatasetS2_transcribed.tsv from a 220 dpi rendering). The 29 variants with
# minigene result "No impact". The paper's main text independently fixes the count ("eight of 30 PKD1
# splicing candidates ... and zero of seven PKD2 splicing candidates influenced RNA splicing") and its
# Table 2 names the eight positives, which match Dataset S2 exactly.
_S2 = pd.read_csv(ROOT / "data" / "papers" / "34257392_DatasetS2_transcribed.tsv", sep="\t", comment="#")
for r in _S2[_S2.minigene == "No impact"].itertuples():
    note = f"A3, Dataset S2. Published SpliceAI {r.spliceai_published}."
    if r.location == "c.11713-3C>G":
        note += (" Designation column says IVS42-3C>T; allele resolved as C>G because local SpliceAI "
                 "reproduces the published 0.95 for C>G and gives 0.03 for C>T.")
    add(r.gene, r.location, "NEG", "none", "minigene", "engineered", "34257392",
        f"{r.location} {r.designation}", note)

# PMID 39757590 (Chen 2025, RNA Biology): PKD1 5'UTR variants tested by luciferase reporter with
# RT-qPCR showing transcript level unchanged. Added after external review, which noted this study was
# cited but its tested variants were absent. Translation class: these are reported descriptively and
# are not used in any splicing statistic. Two of the three INCREASE translation.
Q_CHEN = ("Luciferase reporter assays and RT-qPCR results reveal that rs2092942382 and rs1596636969 "
          "increase, while rs2092942900 decreases main gene translation without affecting transcription")
add("PKD1", "rs1596636969", "POS", "translation INCREASED; disrupts a uORF start codon",
    "5'UTR Renilla/firefly dual-luciferase reporter; reporter RNA measured by RT-qPCR", "population", "39757590", Q_CHEN,
    "A7. Direction: increase. Transcript level unchanged by RT-qPCR.")
add("PKD1", "rs2092942382", "POS", "translation INCREASED about 42%; changes the -3 nucleotide of uORF2",
    "5'UTR Renilla/firefly dual-luciferase reporter; reporter RNA measured by RT-qPCR", "population", "39757590",
    "rs2092942382 variant increased the relative Rluc/Fluc ratio by approximately 42%",
    "A7. Direction: increase. Transcript level unchanged by RT-qPCR.")
add("PKD1", "rs2092942900", "POS", "translation DECREASED; removes the uORF1 stop codon",
    "5'UTR Renilla/firefly dual-luciferase reporter; reporter RNA measured by RT-qPCR", "population", "39757590", Q_CHEN,
    "A7. Direction: decrease. Transcript level unchanged by RT-qPCR.")


# ---------------------------------------------------------------- not in the truth set
EXCLUDED = [
    ("PKD1", "c.1202C>T; c.1248C>T; c.11025G>A; c.11119C>T", "37468838",
     "Selected for minigene testing but 'not verifiable for technical reasons'. UNTESTED, not negative."),
    ("PKD1", "c.7066-131G>A (alone)", "42502691",
     "De novo component never observed without rs3874658; single-variant function unknown."),
    ("PKD1", "c.-59C>A; c.-8C>A (in cis with c.-69dupG)", "41006799",
     "Engineered rescue constructs, not single variants."),
    ("PKD1", "uORF start-loss constructs", "41006799",
     "Engineered; exact edits in figure only. Candidate positive-direction test once extracted."),
    ("PKD1", "IVS45+56del25", "11058904",
     "Legacy notation; deleted bases not resolvable from the abstract. Needs full text."),
    ("PKD1", "20-bp deletion in intron 43", "11800305",
     "Position not stated in abstract. Needs full text."),
    ("PKD1", "18-bp and 20-bp deletions in a 75-bp intron", "7633405",
     "Positions not stated in abstract. Needs full text."),
    ("PKD1", "rs201204878", "31384335", "In silico prediction only, no functional assay."),
    ("PKD1", "3'UTR miR-17 six-nucleotide motif", "41558825",
     "Engineered in MOUSE Pkd1; human coordinates only in supplement. Moves to the therapeutic-site aim."),
    ("PKD1", "rs1016465177 (c.-67C>T)", "39757590",
     "Described as a uORF1 synonymous variant but not functionally tested in that paper."),
    ("PKD1", "p.Q987H; p.L3675L; p.L3754L", "24907393",
     "Table S3 negatives that more than one SNV could produce; not guessed."),
    ("PKD2", "~10 unnamed exonic variants tested negative", "26692149",
     "Negatives not named in abstract; full text not open access. MANUAL DOWNLOAD."),

    ("PKD1/PKD2", "urinary-cell targeted RNA-seq variants", "38944240",
     "Kidney Int research letter; no abstract text. MANUAL DOWNLOAD."),
    ("PKD2", "exon 14 nonsense and missense with cryptic splicing", "10541293",
     "Variants not named in abstract. MANUAL DOWNLOAD."),
]


def main() -> None:
    df = pd.DataFrame(E)
    df.insert(0, "entry_id", [f"E{i:03d}" for i in range(1, len(df) + 1)])
    df.to_csv(CUR / "truthset_candidates.tsv", sep="\t", index=False)
    pd.DataFrame(EXCLUDED, columns=["gene", "variant", "pmid", "reason"]).to_csv(
        CUR / "excluded.tsv", sep="\t", index=False)
    u = df.drop_duplicates(["gene", "hgvs_c"])
    print(f"entries {len(df)} (unique variants {len(u)})")
    print(u.groupby(["gene", "label"]).size().to_string())
    print("legacy conversions:", [e["hgvs_c"] for e in E if "legacy" in e["note"]])


if __name__ == "__main__":
    main()
