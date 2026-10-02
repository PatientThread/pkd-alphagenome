# Papers to download by hand

These could not be retrieved automatically: they are not in the PubMed Central open-access subset,
and PMC article pages sit behind a reCAPTCHA for scripts. Everything else in the truth set came from
full text or abstracts retrieved by `analysis/02_screen_and_fetch.py`.

**How to use this list.** Download each PDF through your browser (university or RCP library access,
or the publisher), and save it into `data/papers/` named by PMID, for example `34257392.pdf`. Tell me
when they are there and I will extract the variants, add them to `analysis/04_build_truthset.py` with
verbatim quotes, and re-run verification.

## Status, 19 September evening

**All three papers and both supplements have now been supplied and used** (amendment A3). Received and
filed in `data/papers/`: PMIDs 24907393, 34257392, 26692149, plus Table S3 (image), Dataset S2 (PDF),
the NAR miR-17 supplement and Figure 4A of the 5'UTR paper. The truth set now holds 121 variants with
68 negatives.

The only thing still missing:

| PMID | What is missing | Why it matters |
|---|---|---|
| 24907393 | About 9 rows of Table S3 below p.L3675L (the image supplied shows 24 of the 33 tested variants) | Would add roughly 7 more negatives, including the exon 38 and 39 variants |

Nothing else is outstanding. Two variants in Table S3 (p.Q987H, p.L3675L) can never be used: more than
one DNA change could produce them, so they are excluded rather than guessed.

Save them into `data/papers/` with the PMID at the start of the file name.

## Priority 1 (original list): these hold the missing negatives

The truth set has 44 positives and only 7 negatives. These three papers tested roughly 80 variants
between them and found most of them to be normal. Recovering them is the single biggest improvement
available, because it is what would allow a proper discrimination analysis with a confidence
interval.

| PMID | Paper | What it holds | Journal |
|---|---|---|---|
| 34257392 | Identification of novel single-nucleotide variants altering RNA splicing of PKD1 and PKD2 (2022) | 37 splice-region variants tested by minigene: 8 positive, **29 negative**. None named in the abstract. | not open access |
| 24907393 | Defective pre-mRNA splicing in PKD1 due to presumed missense and synonymous mutations (2014) | 33 PKD1 exonic variants by minigene: 3 positive (already in), **about 30 negative**. | not open access |
| 26692149 | Three exonic mutations in PKD2 alter splicing of its pre-mRNA in a minigene system (2016) | 13 PKD2 variants by minigene: 3 positive (already in), **10 negative**. | not open access |

## Priority 2: positives and key sources with detail missing

| PMID | Paper | Why |
|---|---|---|
| 38944240 | Targeted RNAseq from patients' urinary cells to validate pathogenic noncoding variants in ADPKD genes (Kidney Int 2024) | Research letter with no abstract text; likely several RNA-validated non-coding variants. |
| 36029107 | A common intronic single nucleotide variant modifies PKD1 expression level (Clin Genet 2022) | rs3874648 is in the truth set from the abstract; the full text has the effect sizes. |
| 19158373 | Evidence for pathogenicity of atypical splice mutations in ADPKD (CJASN 2009) | PKD2 IVS8+5G>A is in from the abstract; the paper also reviews 17 atypical splice mutations. |
| 10541293 | Aberrant splicing in the PKD2 gene as a cause of polycystic kidney disease (1999) | Two PKD2 exon 14 variants with RNA evidence, not named in the abstract. |
| 11058904 | Novel splicing and missense mutations in PKD1 (2000) | IVS45+56del25: legacy notation, deleted bases not resolvable without the full text. |
| 11800305 | Molecular defect of PKD1 gene resulting in abnormal RNA processing in a Thai family (2001) | 20-bp intron 43 deletion; position not given in the abstract. |
| 7633405 | Splicing mutations of PKD1 induced by intronic deletion (1995) | Two intronic deletions; positions not given in the abstract. |

## Also useful, lower priority

| PMID | Paper | Why |
|---|---|---|
| 25757501 | Splicing defects caused by exonic mutations in PKD1 (RNA Biol 2015) | Review of the 2014 work; may list the negatives in a table. |
| 30185468 | Human-specific abnormal alternative splicing of wild-type PKD1 (2018) | Context for PKD1 splicing background, not variants. |
| 12070253 | Mutation screening of the PKD1 transcript by RT-PCR (2002) | May contain RT-PCR-confirmed splice variants. |

## Supplementary tables already identified as needed

- **PMID 41558825 / 41279527** (miR-17 six-nucleotide motif, NAR 2026): Supplementary Table S4 gives
  the exact edited nucleotides. Needed for the therapeutic-site aim, not the truth set.
- **PMID 41006799** (PKD1 5'UTR variants, EJHG 2025): Figure 4A shows the exact uORF start-loss edits.
  These would add a rare *positive-direction* test (engineered variants that raised PC1 by up to 130%).
