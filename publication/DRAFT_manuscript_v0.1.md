# Where sequence-to-function models see, and do not see, proven PKD1 and PKD2 variants

**DRAFT v0.1, 19 September 2026.** Numbers marked PENDING await the gene-wide SpliceAI run.
Nothing in this draft should be submitted before the checks listed at the end are done.

Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP
Consultant Nephrologist, 9 Harley Street, London W1G 9QY, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net

---

## Abstract

**Background.** Sequence-to-function models are entering variant interpretation, and the AlphaGenome
Atlas now provides a precomputed impact score (AVI) for every possible single-nucleotide variant in
the genome. Their performance at *PKD1*, a locus that defeats conventional pipelines, is unknown.

**Methods.** A recorded PubMed search identified papers reporting functional tests of *PKD1* or *PKD2*
variants. A truth set of 121 unique variants was curated from 24 primary sources: 51 shown to disrupt
splicing or expression, 68 tested and indistinguishable from wild type, 1 inconclusive and 1
discordant between laboratories. Every entry carries a verbatim quote from its source, checked
mechanically against the stored text, and every variant was confirmed against GRCh38 by
VariantValidator. The analysis plan was fixed and hash-stamped before scoring. Variants were scored
with the AlphaGenome Atlas (AVI and its feature attributions), AlphaGenome on demand (splicing and
RNA-seq), SpliceAI run locally, and CADD v1.7.

**Results.** For splice-acting variants, AlphaGenome and SpliceAI were indistinguishable (AUC 0.84,
95% CI 0.77 to 0.92 for AlphaGenome splice-site usage; 0.84, 0.76 to 0.91 for SpliceAI). At the
conventional SpliceAI threshold of 0.2, 16 of 68 experimentally normal variants (24%) were flagged.
AVI ranked non-coding and synonymous positives at a median 94.2nd percentile of the *PKD1* locus, with
5 of 25 in the top 1%. Three classes of proven variant were missed by every tool: a 5'UTR variant that
reduces translation by 87%, a common intronic variant that lowers full-length *PKD1*, and engineered
disruption of the two elements that are current drug targets (the 3'UTR miR-17 seed and the 5'UTR
upstream open reading frames). Where AVI did rank the drug-target sites highly, evolutionary
conservation was the largest contributor for all 36 possible changes, not any functional prediction.

**Conclusions.** At *PKD1* and *PKD2*, a multimodal sequence model matches but does not beat a
dedicated splicing model, and both carry a substantial false-alarm burden. Post-transcriptional
mechanisms, including the two elements currently being drugged, are outside what these models
predict. PENDING: the gene-wide false-alarm comparison.

---

## Introduction

[To write. Points to cover, with sources already in the truth set:
 - ADPKD remains genetically unexplained in about 10% of patients after coding-region testing.
 - *PKD1* is hard: six pseudogenes over exons 1 to 33, GC-rich, poorly covered by capture methods.
 - ADPKD is a dosage disease with a functional threshold, so quantitative regulatory prediction has a
   clear biological job.
 - The AlphaGenome Atlas has just made a genome-wide precomputed score available, and laboratories
   will use it. No locus-specific evaluation exists for the PKD genes.
 - Two *PKD1* regulatory elements are drug targets in current trials, which makes their predictability
   a practical question, not an academic one.]

## Methods

### Literature search and screening
Two fixed PubMed queries (`analysis/01_literature_search.py`, run 19 September 2026) returned 391
unique records. Titles were screened for papers reporting human *PKD1* or *PKD2* variants with
RNA-level, reporter or regulatory functional evidence; 47 were included and retrieved as open-access
full text (26) or abstract. Three further papers and two supplementary files, not available through
automated retrieval, were obtained manually. Screening decisions are recorded in
`data/curated/screening_log.tsv`.

### Truth set
Entries were labelled by what the experiment showed, never by what a predictor said: POS (aberrant
splicing or reduced expression demonstrated), NEG (tested and indistinguishable from wild type),
INDET (assay uninterpretable). Variants named but never tested were excluded and listed with reasons
(`data/curated/excluded.tsv`). Where two sources reported opposite results for the same variant it
was labelled DISCORDANT and excluded from both classes.

Each entry carries a short verbatim quote from its source. `analysis/05_verify.py` rejects any entry
whose quote is not literally present in the stored source text, normalising only typography. Every
HGVS description was validated against the MANE Select transcripts NM_001009944.3 (*PKD1*) and
NM_000297.4 (*PKD2*) on GRCh38 by VariantValidator, which also supplied genomic coordinates; 124 of
124 entries passed both checks. Legacy IVS notation was converted mechanically from the RefSeq exon
table. Variants published only as protein changes were converted by reading the codon from the
reference coding sequence and accepting the result only when exactly one single-base substitution
produced the stated amino-acid change (`analysis/p_to_c.py`); all 17 such conversions were unique,
and the method reproduced exactly the five variants for which the same papers also printed c.
notation.

### Scoring
AlphaGenome Atlas AVI scores and the 18 AVI feature attributions were retrieved for every truth-set
SNV and for every possible SNV across each locus (gene span plus 5 kb upstream: 156,573 SNVs for
*PKD1*, 225,576 for *PKD2*). AlphaGenome on-demand predictions (splice sites, splice-site usage,
splice junctions, RNA-seq; 1,048,576 bp context) were obtained for all variants including indels,
summarised over all tracks and over 11 kidney tracks. SpliceAI 1.3.1 was run locally (distance 500,
unmasked) against GRCh38; the installation reproduced published SpliceAI scores exactly for a
reference variant, and agreed with a further 37 published scores in 29 of 37 cases to two decimal
places. CADD v1.7 PHRED scores were read from the official whole-genome file.

### Statistics and pre-specification
The analysis plan (`docs/ANALYSIS_PLAN.md`) was written and SHA-256 hash-stamped before scoring, with
four dated amendments, each recording whether the scores it could affect already existed. The primary
endpoint was the percentile rank of AVI within the variant's own locus. Discrimination is reported as
the Mann-Whitney AUC with a stratified bootstrap 95% interval (2,000 resamples). Composite scores
(AVI, CADD) are evaluated primarily on non-coding and synonymous variants, where a high score cannot
arise from a protein consequence.

## Results

### The truth set
[Table 1: composition by gene, label, mechanism class, assay type and source.]

### Discrimination
[Figure 1.] AlphaGenome splice-site usage AUC 0.84 (0.77 to 0.92); SpliceAI 0.84 (0.76 to 0.91);
AlphaGenome splice junctions 0.83 (0.75 to 0.91); splice sites 0.83 (0.74 to 0.90). Restricting
AlphaGenome to kidney tracks changed nothing (0.84). On non-coding and synonymous variants, AVI
reached 0.77 (0.64 to 0.89) and CADD 0.72 (0.58 to 0.84).

### False alarms
At SpliceAI 0.2, 41 of 51 positives were flagged and 16 of 68 negatives (24%); at 0.5, 35 of 51 and 7
of 68 (10%). PENDING: the same burden expressed across all 156,573 possible *PKD1* variants, and the
AVI threshold that matches it.

### Ranking within the locus
For 25 non-coding or synonymous *PKD1* positives, AVI's median locus percentile was 94.2 (5 of 25 in
the top 1%, none in the top 0.1%); CADD 88.2 (2 of 25). About 8,000 other possible *PKD1* changes
outrank a typical proven splice variant.

### What the models do not see
[Figure 2.] c.-69dupG (87% less translation), c.-87A>T and c.-20A>T (9% and 34% more translation),
the six-base miR-17 seed edit, and rs3874648. AlphaGenome predicted a change in *PKD1* RNA of 0.047
or less (log2) for all five.

### What drives AVI
[Figure 3.] For all 36 possible changes at the two drug-target sites, conservation was the largest
contributor. Among coding positives, 8 of 15 were driven by protein damage and only 4 by splicing.

### Adversarial cases
The intron 16 pseudoexon allele was detected by AlphaGenome and SpliceAI when presented as a single
two-base allele (splice-junction score 7.4 against 0.17 and 0.93 for its parts; SpliceAI 0.78 against
0.12 and 0.24), but CADD cannot score it, and conventional pipelines separate the two changes.

## Discussion

[To write. The argument:
 1. For the commonest testable class, the new model matches the old one. Laboratories gain little by
    switching, and AVI has no validated threshold.
 2. The false-alarm burden is the practical limit, not sensitivity.
 3. The blind spots are systematic, not incidental: anything acting after transcription is invisible,
    which includes both elements currently being drugged in ADPKD and the only proven common dosage
    modifier. For a disease whose severity varies with polycystin dose, that is the interesting class.
 4. AVI's high scores at conserved regulatory sites come from conservation, so the score can look
    informative where the model is not.
 5. The benchmark itself is assay-dependent: one variant has opposite results in two laboratories.
 6. Pipeline design matters as much as model choice, as the intron 16 case shows.]

## Limitations

Most positives were demonstrated in blood RNA or in minigenes in non-kidney cell lines. Published
variants are enriched for those a predictor once found interesting, which favours all predictors and
particularly the older tools used to choose what to test. Negatives are enriched for missense
variants, which are valid negatives for splicing but not for composite scores. *PKD2* lies on an AVI
training chromosome, so its AVI results carry a leakage caveat and are secondary; *PKD1* lies on
chromosome 16, used for neither training nor validation. About 21 further tested variants from one
supplementary table were not available.

## Data and code availability

Analysis code, the curated truth set with its quotes, the screening log and the pre-specified plan
with its hash record are in the project repository. AlphaGenome predictions are redistributed under
the AlphaGenome Output Terms of Use (non-commercial), as set out in the repository's
LEGALLY_BINDING_TERMS_OF_USE.txt. Publisher full texts are not redistributed.

## Declarations

No funding. No institutional or commercial affiliation. The AlphaGenome API was accessed under a
personal, individual, non-commercial account. Competing interests: [to complete].

## References

Generated mechanically from the stored PubMed records in `publication/references_auto.md`; two
entries need their author list completed by hand, and every entry should be checked against the
journal before submission.

---

## Before submission

1. Add the gene-wide SpliceAI results and Figure 4.
2. Decide *PKD1* only, or *PKD1* and *PKD2*.
3. Recover the remaining supplementary negatives if possible.
4. Write Introduction and Discussion.
5. Choose the journal and reformat. Candidates discussed: Clinical Kidney Journal, Kidney
   International Reports, European Journal of Human Genetics, a letter to AJKD.
6. Check every reference against the journal record.
7. Consider sending the discordant result (c.11257C>T) to both original groups before publishing it.
