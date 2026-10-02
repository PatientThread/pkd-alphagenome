# Where sequence-to-function models see, and do not see, functionally proven variants in PKD1 and PKD2

**DRAFT v1.3, 20 September 2026.** Complete draft, all results in for both genes, truth set complete.
Submission checklist at the end.

Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP
Consultant Nephrologist, 9 Harley Street, London W1G 9QY, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net

Running title: Sequence models at the PKD loci
Key words: autosomal dominant polycystic kidney disease, PKD1, PKD2, RNA splicing, variant
interpretation, machine learning

---

## Abstract

**Background.** Sequence-to-function models are entering clinical variant interpretation. The
AlphaGenome Atlas now provides a precomputed impact score for every possible single-nucleotide
variant in the human genome, including a combined AlphaGenome Variant Impact (AVI) score. Their
behaviour at the polycystic kidney disease loci, where about one patient in ten remains genetically
unexplained after coding-region testing, has not been examined.

**Methods.** A recorded PubMed search (391 records) identified papers reporting functional tests of
*PKD1* or *PKD2* variants. We curated 125 unique variants from 24 primary sources: 51 shown
experimentally to disrupt splicing or expression, 72 tested and indistinguishable from wild type, one
uninterpretable and one discordant between laboratories. Every entry carries a verbatim quote from
its source, machine-checked against the stored text, and every variant was confirmed on GRCh38 by
VariantValidator. The analysis plan was hash-stamped before scoring, with five dated amendments.
Variants were scored with the AlphaGenome Atlas (AVI and its 18 feature attributions), AlphaGenome on
demand (splice sites, splice-site usage, splice junctions and RNA-seq), SpliceAI [7] run locally, and CADD v1.7 [8]. Discrimination is reported as the Mann-Whitney area under the curve (AUC) with a
stratified bootstrap 95% confidence interval.

**Results.** For splice-acting variants, AlphaGenome and SpliceAI were statistically indistinguishable
(AlphaGenome splice-site usage AUC 0.85, 95% CI 0.77 to 0.92; SpliceAI 0.84, 0.76 to 0.91). Applying
SpliceAI's conventional threshold of 0.2 flagged 41 of 51 positives but also 17 of 72 experimentally
normal variants (24%). For non-coding and synonymous variants, AVI achieved AUC 0.77 (0.64 to 0.89)
and CADD 0.72 (0.58 to 0.84). Ranked against all 156,573 possible single-nucleotide variants in
*PKD1*, proven non-coding and synonymous positives sat at a median 94.2nd percentile, with 5 of 25 in
the top 1%. Three classes of proven variant were missed by every tool: a 5'UTR variant reducing
translation by 87%, a common intronic variant that lowers full-length *PKD1*, and disruption of the
two elements currently being targeted therapeutically. Where AVI ranked those drug-target sites
highly, evolutionary conservation was the largest contributor for all 36 possible changes.

**Conclusions.** At *PKD1* and *PKD2* a multimodal sequence model matches, but does not exceed, a
dedicated splicing model, and both carry a substantial false-alarm burden at conventional thresholds.
Mechanisms acting after transcription, including both elements under therapeutic development and the
only proven common dosage modifier, lie outside what these models predict. Laboratories should treat
a high score at a conserved regulatory element as a prompt for RNA studies, not as evidence of
mechanism.

---

## Introduction

Autosomal dominant polycystic kidney disease (ADPKD) affects approximately 1 in 400 to 1,000 people
and nearly 12 million worldwide, and is caused principally by variants in *PKD1* and *PKD2* [1,2]. The disease is a disorder of gene dosage: heterozygous loss reduces functional
polycystin, and cyst formation follows when the amount of functional protein in a tubular cell falls
below a threshold [2]. This is why therapeutic attention has turned to raising polycystin
production from the intact allele, and two regulatory elements in *PKD1* are now targets. A
microRNA-17 binding site in the 3' untranslated region represses *PKD1*, and disrupting it, either by
base editing or transiently with a steric-blocking oligonucleotide, stabilises the transcript, raises
polycystin-1 and attenuates cyst growth [2]. Two upstream open reading frames in the 5'
untranslated region likewise restrain translation, and blocking them raises polycystin-1 expression
[3,4].

Genetic diagnosis remains incomplete. Approximately 10% of individuals with an ADPKD phenotype have no
causal variant identified after coding-region testing [4,5]. *PKD1* is unusually
hostile to sequencing and to interpretation: exons 1 to 32 are repeated six times on chromosome 16 as
pseudogenes, so mapping and variant calling are error-prone and gene conversion confuses
interpretation [6]. Non-coding variants in regulatory regions are not an established
cause of disease in ADPKD, yet the unexplained fraction is precisely where they would be found [4].

Against this background, sequence-to-function models have become directly available to diagnostic
laboratories. AlphaGenome predicts molecular readouts from DNA sequence [9], and the AlphaGenome Atlas
provides precomputed predictions for all possible single-nucleotide variants together with a combined
AlphaGenome Variant Impact score, which in the developers' own rare-disease evaluation placed the
causal variant among the top 50 candidates 29.5% of the time against 12.5% for CADD [AlphaGenome
Atlas preprint]. Its developers are explicit that these are research tools that "can only act as part
of the evidence chain leading to clinical diagnoses, and are not sufficient evidence on their own"
[10]. No locus-specific evaluation exists for the PKD genes, and general
benchmarks cannot answer a locus-specific question: performance depends on the local density of
exons, on conservation, and on which mechanisms are common in that gene.

We therefore asked three questions. Does AlphaGenome identify *PKD1* and *PKD2* variants that have
been functionally tested in the laboratory? Does it do so better than the tools laboratories already
use? And what does it miss?

## Methods

### Search and screening
Two fixed PubMed queries were run on 19 September 2026 and returned 391 unique records. Titles were
screened for reports of human *PKD1* or *PKD2* variants with RNA-level, reporter or regulatory
functional evidence, and 47 papers were included and retrieved, 26 as open-access full text and the
remainder as abstracts. Three further papers and three supplementary files were obtained manually.
Queries, hits and screening decisions are recorded in the repository, and the 24 primary sources with what each contributed are listed in Supplementary Table S1.

### Truth set construction
Variants were labelled by experimental outcome, never by predictor output: POS (aberrant splicing or
reduced expression demonstrated), NEG (tested by the same class of assay and indistinguishable from
wild type), INDET (assay uninterpretable). Variants named but not tested were excluded and listed with
reasons. Where two sources reported opposite outcomes for the same variant it was labelled DISCORDANT
and excluded from both classes.

Two mechanical checks protect the set. First, each entry carries a short verbatim quote from its
source; the verification script rejects any entry whose quote is not literally present in the stored
source text, normalising only typography. Second, every HGVS description was validated against the
MANE Select transcripts NM_001009944.3 (*PKD1*) and NM_000297.4 (*PKD2*) on GRCh38 using
VariantValidator [11], which also supplied genomic coordinates. All 124 entries passed both checks. Legacy
IVS notation was converted from the RefSeq exon table. Variants published only as protein changes were
converted by reading the codon from the reference coding sequence and accepting a result only where
exactly one single-base substitution produced the stated amino-acid change; all 17 such conversions
were unique, and as a control the method reproduced exactly the five variants for which the same
papers also printed cDNA notation.

### Scoring
AVI scores and the 18 AVI feature attributions were retrieved from the AlphaGenome Atlas for every
truth-set single-nucleotide variant, and for every possible single-nucleotide variant across each
locus (gene span plus 5 kb upstream; 156,573 for *PKD1* and 225,576 for *PKD2*) to provide the
background distribution. AlphaGenome on-demand predictions (splice sites, splice-site usage, splice
junctions, RNA-seq; 1,048,576 bp context) were obtained for all variants including insertions,
deletions and one dinucleotide allele, summarised as the maximum absolute value across all tracks and
across 11 kidney tracks. SpliceAI 1.3.1 [7] was run locally (distance 500, unmasked) on GRCh38. CADD v1.7
[8] scores were read from the official whole-genome file, which also provided the CADD background.

### Validation of the comparator
The local SpliceAI installation reproduced the published scores for a reference variant exactly
(acceptor gain 0.12, donor gain 0.06) and, against 37 SpliceAI scores published in a source paper,
agreed exactly to two decimal places in 29 and within 0.02 in 31, the largest difference being 0.12.
For the locus-wide background the same arithmetic was implemented in batches; the batched
implementation reproduced the stock package exactly on 48 test variants (maximum difference 0.0000).

### Pre-specification and statistics
The analysis plan was written and SHA-256 hash-stamped before any truth-set variant was scored, with
one disclosed exception (a single variant used to probe the interface). Five dated amendments each
record whether the results they could affect already existed. Two amendments matter for
interpretation. First, AVI and CADD are composite pathogenicity scores, so a protein-truncating
variant that does not affect splicing is correctly scored high; their primary evaluation therefore
uses non-coding and synonymous variants only. Second, AVI was trained on gnomAD [12] variants from
chromosomes 1, 4, 7, 8, 10, 13 and 15, so *PKD1* (chromosome 16) is held out whereas *PKD2*
(chromosome 4) is not; AVI results are reported primarily for *PKD1*, and *PKD2* variants overlapping
gnomAD positions are flagged.

The primary endpoint was the percentile rank of AVI within the variant's own locus. Discrimination is
the Mann-Whitney AUC with a stratified bootstrap 95% interval (2,000 resamples). No AUC was computed
before the truth set contained sufficient negatives.

## Results

### The truth set
The final set contains 125 unique variants from 24 primary sources: 51 positive, 72 negative, one
uninterpretable (*PKD1* c.288-12C>A) and one discordant (*PKD1* c.11257C>T, positive in one minigene
study and negative in another). Positives comprise 30 non-coding, 10 missense, 5 synonymous, 5
nonsense and one exon-boundary deletion; negatives comprise 43 missense, 24 non-coding, 3 synonymous
and 2 nonsense. Twenty-seven positives were carried by patients and the remainder engineered from
reported variants (Supplementary Table S1); 32 of 51 positives were demonstrated by minigene assay and the others by patient
RNA. (Table 1.)

**Table 1. Composition of the truth set.**

| Gene | Variant class | Positive | Negative | Patient-derived | Demonstrated by minigene |
|---|---|---|---|---|---|
| PKD1 | noncoding | 26 | 18 | 18 | 30 |
| PKD1 | synonymous | 4 | 2 | 1 | 5 |
| PKD1 | missense | 7 | 32 | 1 | 38 |
| PKD1 | nonsense | 4 | 0 | 0 | 4 |
| PKD1 | other coding | 1 | 0 | 1 | 0 |
| PKD2 | noncoding | 4 | 6 | 4 | 6 |
| PKD2 | synonymous | 1 | 1 | 1 | 2 |
| PKD2 | missense | 3 | 11 | 2 | 13 |
| PKD2 | nonsense | 1 | 2 | 0 | 3 |
| Total |  | 51 | 72 | 28 | 101 |

Counts are entries per class; a variant demonstrated both in patient RNA and in a minigene
contributes to both assay columns.

### AlphaGenome and SpliceAI perform equivalently on splice-acting variants
Across all labelled variants (51 positive, 72 negative), AlphaGenome splice-site usage achieved an AUC
of 0.85 (95% CI 0.77 to 0.92), splice junctions 0.83 (0.75 to 0.91) and splice sites 0.82 (0.74 to
0.89), against SpliceAI 0.84 (0.76 to 0.91). Restricting AlphaGenome to its 11 kidney tracks changed
nothing material (0.84, 0.77 to 0.91). Restricted to *PKD1*, the held-out gene, the ordering was
unchanged (AlphaGenome splice junctions 0.82, 0.73 to 0.91; SpliceAI 0.81, 0.71 to 0.91). The confidence intervals overlap almost completely, and the result is
unchanged whether the discordant variant is excluded or counted as positive. (Figure 1, Table 2.)

**Table 2. Discrimination between variants shown to disrupt splicing or expression and variants tested and found normal.**

| Score | Variants evaluated | Positives | Negatives | AUC | 95% CI |
|---|---|---|---|---|---|
| SpliceAI (delta score) | all labelled | 51 | 72 | 0.84 | 0.76 to 0.91 |
| AlphaGenome splice sites | all labelled | 51 | 72 | 0.82 | 0.74 to 0.89 |
| AlphaGenome splice-site usage | all labelled | 51 | 72 | 0.85 | 0.77 to 0.92 |
| AlphaGenome splice junctions | all labelled | 51 | 72 | 0.83 | 0.75 to 0.91 |
| AlphaGenome splice-site usage, kidney tracks only | all labelled | 51 | 72 | 0.84 | 0.77 to 0.91 |
| AlphaGenome splice junctions, kidney tracks only | all labelled | 51 | 72 | 0.83 | 0.75 to 0.91 |
| AVI (Atlas composite) | non-coding and synonymous | 30 | 27 | 0.77 | 0.64 to 0.89 |
| CADD v1.7 | non-coding and synonymous | 31 | 27 | 0.72 | 0.58 to 0.84 |

For composite scores on non-coding and synonymous variants, AVI reached 0.77 (0.64 to 0.89) and CADD
0.72 (0.58 to 0.84).

### False alarms at conventional thresholds
At the SpliceAI threshold of 0.2 commonly used for screening, 41 of 51 positives (80%) were flagged,
together with 17 of 72 experimentally normal variants (24%). At 0.5, 35 of 51 (69%) and 7 of 72 (10%)
were flagged. Individual cases are striking: *PKD1* c.11713-3C>G scores 0.95 and c.9201+5G>A 0.96,
and both were normal in minigene assays. This is consistent with the specificity of 0.55 reported by
one source group on their own 37 variants.

Extending this across the whole gene shows what the burden means in practice. SpliceAI at 0.5 flags
1.04% of all possible *PKD1* single-nucleotide variants. Setting AVI to the threshold that flags the
same number of variants in the gene, SpliceAI recovers 20 of the 25 non-coding and synonymous
positives and AVI recovers 1. Ranked against the whole gene, the median proven positive sits at the
99.4th percentile by SpliceAI and the 94.2nd by AVI (Figure 4).

That comparison is unfavourable to AVI by construction, because AVI is a whole-genome pathogenicity
score and correctly ranks protein-truncating variants at the top of a gene, whereas the question here
concerns non-coding variants. Restricting the background to non-coding positions (a post hoc
sensitivity analysis, amendment A6) narrows the gap considerably: SpliceAI at 0.5 flags 0.56% of
non-coding positions and recovers 17 of 21 non-coding positives, against 8 of 21 for AVI at the same
burden, with both flagging 3 of 18 non-coding negatives. On this restricted background the median
positive ranks at the 99.9th percentile by SpliceAI and the 99.8th by AVI, that is, almost
identically. Both analyses are reported because they answer different questions: the first asks how
a score behaves when pointed at a whole gene, the second how it behaves once coding variants have
been set aside.

*PKD2* shows the same pattern with much smaller numbers. Against the whole gene, SpliceAI at 0.5
flags 0.28% of possible variants and recovers 5 of 5 non-coding and synonymous positives, against 0
of 5 for AVI at matched burden; against a non-coding background the burden falls to 0.12% and the
counts are 4 of 4 and 3 of 4, with median ranks of 99.99th percentile for both tools. With four
positives these figures are descriptive only, and *PKD2* lies on an AVI training chromosome.

### Ranking within the locus
(Figure 4.) For the 25 non-coding or synonymous *PKD1* positives, AVI's median rank was the 94.2nd percentile of
all 156,573 possible variants in the gene, with 5 in the top 1% and none in the top 0.1%; CADD's
median was the 88.2nd percentile with 2 in the top 1%. In practical terms, roughly 8,000 other
possible changes in *PKD1* outrank a typical proven splice variant. Attribution was nonetheless
appropriate: splicing was the largest contributor for 18 of 21 splice-acting *PKD1* positives.

In *PKD2* all five non-coding or synonymous positives ranked in the top 1% (median 99.8th percentile),
but this is not evidence of better performance. The *PKD2* positives are canonical or near-canonical
splice-site variants, whereas the *PKD1* set includes branchpoint, AG-exclusion-zone and deeper
intronic variants; and the background differs markedly, with 4.9% of all possible *PKD1* variants
scoring above an AVI of 1.0 against 1.3% in *PKD2*, because *PKD1* packs 46 exons into 47 kb. Percentile
ranks are therefore comparable within a gene but not between genes.

### Mechanisms that no tool detects
Five variants with demonstrated function act after transcription, and every tool missed all five
(Figure 2). *PKD1* c.-69dupG extends an upstream open reading frame and reduces translation of the
main reading frame by 87%; AlphaGenome predicted a change in *PKD1* RNA of at most 0.047 (log2) and
SpliceAI scored it 0.01. The engineered start-loss variants c.-87A>T and c.-20A>T raise translation by
9% and 34%; predicted changes were 0.029 and 0.012, the former in the wrong direction. The six-base
edit of the miR-17 seed, which stabilises the transcript and raises polycystin-1, produced a predicted
change of 0.003. The common intronic variant rs3874648 (c.10051-239G>A), shown to activate a cryptic
acceptor site and lower full-length *PKD1*, scored at the 10.9th percentile of the locus by AVI, 0.05
by SpliceAI and PHRED 0.3 by CADD.

### High scores at conserved sites come from conservation
AVI ranked the two drug-target sites highly within their own untranslated regions (median 92nd
percentile for the miR-17 seed, 95th for the uORF start codons). For all 36 possible changes across
both sites, evolutionary conservation was the largest contributor to the score, with AlphaGenome's own
functional predictions contributing 0.01 to 0.05 (Figure 3). The same caution applies to coding
positives: of 15 with AVI scores, 8 were driven by protein damage and only 4 by splicing, so a high
composite score at such a variant is not evidence that the splicing effect was detected.

### An adversarial case that models pass and pipelines fail
A reported *PKD1* pseudoexon is activated only by two adjacent changes acting together. Presented as a
single two-base allele, both models detected it clearly (AlphaGenome splice-junction score 7.4 against
0.17 and 0.93 for the constituent changes alone; SpliceAI 0.78 against 0.12 and 0.24). CADD cannot
score the combined allele at all. The diagnostic failure in the original case was therefore a pipeline
failure, since conventional pipelines separate the two changes and filter the common one as benign.

### Leakage and population checks
*PKD1* positives already recorded in ClinVar did not score higher by AVI than those absent from
ClinVar (median 94.5th against 97.4th percentile), consistent with AVI having been trained on
population frequency rather than clinical assertions. At the therapeutic target sites, gnomAD v4.1
records rare changes in the miR-17 seed and start-loss changes at both upstream reading frames,
including the exact uORF2 change used in the laboratory (two alleles) and c.-87A>G at uORF1 (twelve
alleles); whether such individuals have higher polycystin-1 is untested.

## Discussion

Three findings bear on practice.

First, for the class of variant that laboratories most often need to adjudicate in ADPKD, a splicing
effect in or near *PKD1*, the new multimodal model performs like the tool already in use. There is no
accuracy argument for changing practice, and AVI carries the additional difficulty that no validated
threshold exists for it, so a laboratory cannot state in advance how many variants per gene will clear
a given cut-off. The percentile framing used here is one workable substitute. Used across a whole gene it shows a
clear limit: at an equal alarm budget SpliceAI recovered 20 of 25 proven positives and AVI 1.
Restricted to non-coding positions, where the comparison is fairer, AVI ranks proven positives at the
99.8th percentile and the difference falls to roughly twofold in recovery. The practical reading is
that AVI is informative once the search has been narrowed to non-coding candidates, and is not a
substitute for a splicing-specific score when the question is a splicing one.

Second, the practical limit is specificity rather than sensitivity. Nearly a quarter of
experimentally normal variants exceeded SpliceAI's usual screening threshold, and individual variants
with scores of 0.95 and 0.96 were shown to be normal in minigene assays. In a gene where most variants
are private and family segregation is often uninformative, a false alarm has a cost: it invites RNA
studies that may not be available, or worse, an unjustified reclassification. Any prediction at these
loci should be treated as a hypothesis for an RNA study rather than as evidence in itself, which is
also the position taken by the Atlas developers.

Third, the blind spots are systematic rather than incidental. AlphaGenome models what is made from
DNA, not what happens to the message afterwards. Every mechanism we tested that acts after
transcription was invisible: translational repression by upstream reading frames, microRNA-mediated
repression, and a common variant that lowers full-length transcript. This matters more in ADPKD than
in most diseases, for two reasons. The two *PKD1* elements under therapeutic development are exactly
of this kind, so the model cannot be used to find or interpret natural variation at the sites being
drugged. And because ADPKD severity tracks polycystin dose, the variants most likely to explain why
relatives with the same pathogenic allele differ in severity are quantitative modifiers, of which the
one proven common example, rs3874648, was missed by every tool. A search for dosage modifiers using
these scores would not find variants of that kind.

A fourth observation is methodological. Where AVI did rank regulatory sites highly, the score came
almost entirely from evolutionary conservation. Conservation is genuine evidence, but it is not
mechanism, and a score that looks informative for the wrong reason is a trap in variant
interpretation, particularly when feature attributions are not routinely inspected. The AVI
attribution breakdown should be read whenever the score is used.

Finally, the benchmark itself deserves scepticism. One variant in this set, c.11257C>T, altered
splicing in one laboratory's minigene and not in another's. Minigene assays lack the native genomic
context, use non-kidney cell lines, and can disagree. Benchmarks built from published functional data
inherit that variability, and reported accuracies should be read with it in mind.

## Limitations

The truth set is literature-derived, so it is enriched for variants that some earlier predictor found
interesting enough to test. This favours all predictors, and particularly the older tools that were
used to select candidates. Most positives were demonstrated in blood RNA or minigene assays in
non-kidney cell lines rather than in kidney tissue. Negatives are enriched for missense variants,
which are valid negatives for splicing predictions but not for composite pathogenicity scores, which
is why the two classes of score were evaluated on different subsets. *PKD2* lies on an AVI training
chromosome, so its AVI results carry a leakage caveat and are secondary, and its sample is small.
The gene-wide SpliceAI background is an evenly spaced sample of each locus rather than every position, which estimates the flagged proportion to within a few tenths of a percentage point. Finally,
the locus percentile depends on the composition of the gene, so it should not be compared across
genes without stating the background.

## Conclusion

At *PKD1* and *PKD2*, AlphaGenome matches but does not surpass SpliceAI for splice-acting variants,
and both flag a substantial minority of experimentally normal variants at conventional thresholds.
The model does not see mechanisms acting after transcription, which include both *PKD1* elements
currently being targeted therapeutically and the only proven common dosage modifier. These scores are
useful for prioritising RNA studies at the PKD loci, and are not sufficient evidence of mechanism.

## Declarations

**Funding.** No funding was received in support of this study.
**Competing interests.** The author is a director and shareholder of Lawrence Medical Ltd, Patient
Thread Ltd and Scira Diagnostics Ltd, companies developing clinical software and genomic diagnostics,
and is named on United Kingdom patent applications relating to genomic diagnostics, one of which
concerns the use of sequence-to-function model outputs for expression-weighted assessment in
transplantation. None of these interests relates to *PKD1* or *PKD2* variant interpretation, to
polycystic kidney disease, or to any of the tools evaluated here, and none funded this work. The
AlphaGenome application programming interface was accessed under a personal, individual,
non-commercial account and not on behalf of any company.
**Ethics.** No patients, samples or identifiable data were involved; all data are from published
literature and public databases.
**Author contributions.** Sole author: conception, curation, analysis, drafting.
**Acknowledgements.** None.
**Data and code availability.** Analysis code, the curated truth set with its source quotes, the
screening log, and the pre-specified analysis plan with its hash record are available in the project
repository. AlphaGenome predictions are redistributed under the AlphaGenome Output Terms of Use
(non-commercial use only), as set out in the repository. Publisher full texts are not redistributed.
**Use of artificial intelligence.** The author used an AI coding assistant (Claude, Anthropic) to
write the analysis code, to assemble the curated data set from published sources, and to draft the
manuscript. Every variant assignment carries a verbatim quotation from its source that is checked
mechanically against the stored text, and every variant description was validated independently
against the reference genome, as described in Methods. The author reviewed and verified all content,
takes full responsibility for it, and confirms that the AI tool is not an author and is not cited as
one.

## References

1. Xu P, Huang S, Li J, Zou Y, Gao M, Kang R, et al. A novel splicing mutation in the PKD1 gene causes autosomal dominant polycystic kidney disease in a Chinese family: a case report. BMC Med Genet. 2018;19(1):198. doi:10.1186/s12881-018-0706-6. PMID: 30424739.
2. Lakhia R, Song C, Biggers L, Zumwalt M, Alvarez J, Somasundaram A, et al. Disruption of a six-nucleotide miRNA motif improves PKD1 dosage and ameliorates polycystic kidney disease. bioRxiv. 2025. doi:10.1101/2025.10.22.683929. PMID: 41279527.
3. Chen L, Gao X, Liu X, Zhu Y, Wang D. Translational regulation of PKD1 by evolutionarily conserved upstream open reading frames. RNA Biol. 2025;22(1):1-12. doi:10.1080/15476286.2024.2448387. PMID: 39757590.
4. Wedd L, Hort Y, Patel C, Sayer JA, Rius R, Mallett AJ, et al. PKD1 5'UTR variants are a rare cause of disease in ADPKD and suggest a new focus for therapeutic development. Eur J Hum Genet. 2026;34(1):61-69. doi:10.1038/s41431-025-01949-z. PMID: 41006799.
5. Zhang Z, Blumenfeld J, Ramnauth A, Barash I, Zhou P, Levine D, et al. A common intronic single nucleotide variant modifies PKD1 expression level. Clin Genet. 2022;102(6):483-493. doi:10.1111/cge.14214. PMID: 36029107.
6. Pan Q, Liu Y, Sun X, Lu S, Li L, Shen J. Case Report: Functional validation of a PKD1 c.7489 + 5G>A variant in an ADPKD family. Front Genet. 2026;17:1850997. doi:10.3389/fgene.2026.1850997. PMID: 42621998.
7. Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, Darbandi SF, Knowles D, Li YI, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015. PMID: 30661751.
8. Schubach M, Maass T, Nazaretyan L, Röner S, Kircher M. CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Res. 2024;52(D1):D1143-D1154. doi:10.1093/nar/gkad989. PMID: 38183205.
9. Avsec Ž, Latysheva N, Cheng J, Novati G, Taylor KR, Ward T, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature. 2026;649(8099):1206-1218. doi:10.1038/s41586-025-10014-0. PMID: 41606153.
10. AlphaGenome Atlas team; Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. Google DeepMind; 2026. Available from: https://deepmind.google/science/alphagenome/atlas
11. Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R. VariantValidator: Accurate validation, mapping, and formatting of sequence variation descriptions. Hum Mutat. 2018;39(1):61-68. doi:10.1002/humu.23348. PMID: 28967166.
12. Chen S, Francioli LC, Goodrich JK, Collins RL, Kanai M, Wang Q, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625(7993):92-100. doi:10.1038/s41586-023-06045-0. PMID: 38057664.

---

## Submission checklist

1. Write nothing further into Results without re-running `analysis/12_analyse.py`.
3. Complete competing interests, acknowledgements and the AI-use statement to the journal's wording.
4. Format references and verify each against the journal record.
5. Decide whether to report *PKD2* or restrict to *PKD1*.
6. Consider contacting both groups about the discordant variant before publication.
7. Confirm the AlphaGenome terms permit the final form of data release (publication is expressly
   carved out; the repository must carry the terms file).
