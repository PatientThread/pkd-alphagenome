# Functional benchmarking of AlphaGenome and SpliceAI at PKD1 and PKD2

Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP

Consultant Nephrologist, 9 Harley Street, London W1G 9QY, United Kingdom

ORCID 0000-0002-8159-0879

Correspondence: Christopher.lawrence3@nhs.net

Running title: Benchmarking sequence models at the PKD loci

Keywords: autosomal dominant polycystic kidney disease; PKD1; PKD2; RNA splicing; variant interpretation; machine learning

---

## Abstract

Sequence-to-function models are entering variant interpretation, and the AlphaGenome Atlas now
provides precomputed predictions and a combined impact score (AVI) for every possible
single-nucleotide variant in the genome. Their behaviour at the polycystic kidney disease loci is
unknown, although about one patient in ten remains genetically unexplained after coding-region
testing. From a recorded search of 391 records we curated 128 unique *PKD1* and *PKD2* variants from
25 primary sources that had been tested experimentally, each carrying a verbatim quotation from its
source that is checked mechanically and a description validated against GRCh38. Variants were
assigned to the class of experiment that produced their evidence, because a translation reporter is
not evidence about splicing. The splicing class contained 49 variants shown to alter splicing and 70
tested with no effect detected. AlphaGenome's documented splicing composite and SpliceAI performed
similarly (AUC 0.86, 95% CI 0.79 to 0.93 and 0.85, 0.77 to 0.93); the paired difference was 0.009
(-0.021 to 0.042), so these data do not resolve which is better. Resampling whole source studies
widened both intervals to approximately 0.72 to 1.00, showing that the controls are concentrated in
few experiments. The composite impact score AVI was weaker for this question (0.73, 0.64 to 0.82) and
CADD weaker still (0.58, 0.46 to 0.68). Where AVI ranked two *PKD1* regulatory elements highly,
conservation contributed most of the score. Predictions at these loci are best treated as grounds for
an RNA study rather than as evidence of mechanism.

---

## Introduction

Autosomal dominant polycystic kidney disease (ADPKD) affects approximately 1 in 400 to 1,000 people
and is caused principally by variants in *PKD1* and *PKD2* [1]. Disease follows a
reduction in functional polycystin below a threshold, so the amount of protein produced by the intact
allele matters [2]. Two regulatory elements in *PKD1* are consequently under therapeutic
investigation: a microRNA-17 binding site in the 3' untranslated region, disruption of which
stabilises the transcript and raises polycystin-1 in model systems [2], and two upstream
open reading frames in the 5' untranslated region that restrain translation [3,4].

Genetic diagnosis remains incomplete. About 10% of individuals with an ADPKD phenotype have no causal
variant identified after coding-region testing, for reasons that include assay scope, structural
variants and mosaicism as well as regulatory variation [4,5]. *PKD1* is difficult
to sequence and to interpret because exons 1 to 32 are repeated six times on chromosome 16 as
pseudogenes [6].

Sequence-to-function models are now directly available to diagnostic laboratories. AlphaGenome
predicts molecular readouts from DNA sequence [8]; the AlphaGenome Atlas stores
precomputed predictions for all possible single-nucleotide variants and adds AVI, a supervised
composite that integrates those predictions with protein and conservation features [9]. These
are distinct objects, and a result about AVI is not necessarily a result about the underlying
sequence model. The developers state that they "can only act as part of the evidence chain leading to
clinical diagnoses, and are not sufficient evidence on their own" [9].

General benchmarks cannot answer a locus-specific question, because performance depends on the local
density of exons, on conservation, and on which mechanisms are common in that gene. We therefore
assembled a locus-specific benchmark containing both variants shown experimentally to disrupt
splicing or expression and variants tested with no effect detected, and asked three questions: how
well these tools separate the two classes at *PKD1* and *PKD2*; whether AlphaGenome performs
differently from the splicing predictor already in use; and what the scores represent when they are
high.

## Materials and Methods

### Search, screening and curation
Two fixed PubMed queries run on 19 September 2026 returned 391 unique records. Titles were screened
for reports of human *PKD1* or *PKD2* variants with functional evidence; 47 papers were included, of
which 26 were available as open-access full text. Three further papers and three supplementary files
were obtained directly. Screening decisions are in the repository.

Variants were labelled by experimental outcome, never by predictor output. POS denotes a demonstrated
effect; NEG denotes "no effect detected in the reported assay", which is not a claim of clinical
benignity; INDET denotes an uninterpretable assay. Variants named but not tested were excluded, with
reasons. Where two sources reported opposite outcomes for the same variant it was labelled DISCORDANT
and excluded from both classes.

Two mechanical checks protect the set. Each entry carries a short verbatim quotation from its source,
and the verification script rejects any entry whose quotation is not literally present in the stored
source text. Every description was validated against the MANE Select transcripts NM_001009944.3
(*PKD1*) and NM_000297.4 (*PKD2*) on GRCh38 with VariantValidator [10], which also
supplied genomic coordinates. All 131 source-level entries passed both checks; they describe 128
unique variants, 3 of which were reported independently by two papers. Legacy IVS notation was
converted from the RefSeq exon table. Variants published only as protein changes were converted by
reading the codon from the reference coding sequence and accepting a result only where exactly one
single-base substitution produced the stated amino-acid change; all 17 such conversions were unique,
and the method reproduced exactly the eight variants for which the same papers also printed cDNA
notation.

### Assignment to an evidence class
Each variant was assigned to the class of experiment that produced its evidence: splicing (patient
RNA or minigene read out as splicing), abundance (transcript-level measurement), or translation (5'
untranslated region reporter, measuring protein output with transcript level shown to be unchanged).
The primary analysis uses the splicing class alone. The abundance class contains one variant and the
translation class six; these are described individually and no summary statistic is computed for
them. Because two translation-class variants increase output, a demonstrated effect is recorded with
its direction rather than being assumed to be a reduction.

### Scoring
SpliceAI 1.3.1 [11] was run locally with the stock implementation (distance 500, unmasked)
on GRCh38. For AlphaGenome we used the documented splicing composite, max(splice sites) + max(splice
site usage) + max(splice junctions)/5, each component the maximum absolute score across all tracks and
genes [12]. AVI scores and the 18 AVI feature attributions were taken from the Atlas, as
were background distributions for every possible single-nucleotide variant across each locus (gene
span plus 5 kb upstream; 156,573 variants for *PKD1*, 225,576 for *PKD2*). CADD v1.7 [13]
PHRED scores came from the official whole-genome file. A gene-wide SpliceAI background was computed on
an evenly spaced sample of 12,000 positions per locus; the batched implementation used for that
background reproduces the stock package exactly on 48 test variants. Coverage is reported per tool,
and a variant a tool cannot score is counted as unscored, never as a low score.

### Statistics
Discrimination is the Mann-Whitney area under the curve (AUC) with a stratified bootstrap 95%
interval (2,000 resamples). The comparison between AlphaGenome and SpliceAI is a paired difference on
identical variants with its own bootstrap interval; no equivalence margin is claimed and overlapping
marginal intervals are not treated as evidence of equivalence. Because the controls are concentrated
in a few papers, we also report leave-one-study-out AUCs and a bootstrap that resamples whole source
studies. Composite scores (AVI, CADD) are additionally reported on non-coding and synonymous variants,
where a high score cannot arise from a protein consequence.

An analysis plan was written and hash-stamped before scoring, with seven dated amendments recording
whether the results each could affect already existed. Following external review, the analysis was
rebuilt to use one reference standard per assay class and the documented splicing composite; that
rebuild is post hoc and is identified as such in the plan.

## Results

### Composition
The benchmark contains 128 unique variants from 25 sources: 54 with a demonstrated effect, 72 tested
with no effect detected, one uninterpretable (*PKD1* c.288-12C>A) and one discordant (*PKD1*
c.11257C>T, positive in one minigene study and negative in another). Assigned by evidence class, 119
variants are splicing, one abundance and six translation (Table 1). Thirty-one variants with a
demonstrated effect were carried by patients or observed in population data; the remainder were
engineered from reported variants.

### Discrimination within the splicing class
Among 49 positives and 70 negatives, the AlphaGenome splicing composite achieved an AUC of 0.861
(95% CI 0.789 to 0.927) and SpliceAI 0.852 (0.766 to 0.926) (Table 2, Figure 1). The paired difference
was 0.009 (-0.021 to 0.042): these data do not resolve which score is better, and the interval
excludes differences larger than about 0.04 in either direction. AVI reached 0.733 (0.643 to 0.821)
and CADD 0.577 (0.463 to 0.683) on the same variants, with five variants unscored by both because
they are indels absent from the precomputed set. Restricted to non-coding and synonymous variants,
AVI achieved 0.785 (0.644 to 0.903) and CADD 0.749 (0.604 to 0.870).

### Robustness
Leave-one-study-out AUCs ranged from 0.812 to 0.924 for SpliceAI and 0.813 to 0.933 for AlphaGenome,
so no single source determines the result. Resampling whole studies rather than variants widened the
intervals to 0.724 to 1.000 and 0.719 to 0.999 respectively, which is the more honest expression of
uncertainty given that 22 of the negatives come from one paper and 29 from another. Performance was
lower for variants 3 to 20 nucleotides from a splice site (SpliceAI 0.791, AlphaGenome 0.789) than
for exonic and untranslated-region variants (0.819 and 0.836), and lower for engineered than for
patient-derived variants (Table 3).

### Ranking within the locus
Ranked against all 156,573 possible *PKD1* variants, proven non-coding and synonymous positives sat
at a median 99.4th percentile by SpliceAI and 94.2nd by AVI (Figure 3). At a matched burden, defined
as the threshold at which each score flags 1.04% of the locus, SpliceAI recovered 20 of 25 such
positives and AVI 5. Restricting the background to non-coding positions, which is the fairer
comparison for a whole-genome pathogenicity score asked a non-coding question, the burden falls to
0.56% and the counts become 17 of 21 for SpliceAI and 18 of 21 for AVI, with both flagging 3 of 18
negatives. The median rank is then 99.9th percentile for SpliceAI and 99.8th for AVI.

### Operating characteristics
At the conventional SpliceAI threshold of 0.2, 41 of 51 positives were flagged together with 17 of 72
variants in which no effect was detected; at 0.5 the counts were 35 and 7. Several individual
examples are informative: *PKD1* c.11713-3C>G scores 0.95 and c.9201+5G>A 0.96, and neither altered
splicing in a minigene. AVI has no validated threshold, so no equivalent operating point can be
quoted for it.

### What a high AVI score represents
At the two *PKD1* regulatory elements under therapeutic investigation, AVI ranked changes highly
within their own untranslated regions (median 92nd percentile for the microRNA-17 site and 95th for
the upstream open reading frame start codons), but for all 36 possible substitutions the largest
contributor was evolutionary conservation, with AlphaGenome's own functional predictions contributing
0.01 to 0.05 (Figure 2). Among coding positives with AVI scores, 8 of 15 were driven by protein
features and 4 by splicing, so a high composite score at such a variant is not evidence that a
splicing effect was detected.

### Variants whose evidence is not a splicing assay
Seven variants are reported individually (Table 4). rs3874648, which activates a cryptic acceptor and
lowers full-length *PKD1*, received low scores from both splicing predictors; this is a miss on a task
both models address. The six 5' untranslated region variants, from two independent laboratories, were
assayed by luciferase reporters measuring protein output with transcript level shown to be unchanged,
an endpoint neither model predicts. Their measured effects span roughly 42% higher to 87% lower
translation, yet the AlphaGenome splicing composite varied only between 0.030 and 0.046 across all
six, and SpliceAI scored five of six at zero. The two variants in which no effect was detected did
not score lower than those with demonstrated effects: c.-209G>A, which was indistinguishable from
wild type, scored 0.044, above the variant that reduced translation by 87% (0.033). Two of these
variants sit at the uORF2 start codon region tested independently by both laboratories, a natural
variant at c.-20 and an engineered change at the same position.

## Discussion

For the class of variant that dominates unresolved ADPKD, a splicing effect in or near *PKD1*, these
data do not distinguish AlphaGenome's documented splicing composite from SpliceAI. The paired interval
excludes differences beyond about 0.04, which is a more useful statement than either an equivalence
claim or a bare comparison of point estimates. Laboratories using SpliceAI have no accuracy reason to
change, and those adopting AlphaGenome should expect comparable behaviour rather than improvement.

The composite impact score behaves differently from the splicing-specific outputs, and the difference
matters in practice. AVI was clearly weaker for this question, and when it did rank regulatory sites
highly the score came mainly from conservation. Conservation is legitimate evidence, but it is not
mechanism, and it is not specific to the effect a laboratory is trying to establish. Feature
attributions should be inspected whenever a composite score contributes to an interpretation. The
converse also holds: a protein-truncating variant correctly ranked at the top of a gene by AVI says
nothing about splicing, which is why the two classes of score were evaluated on different subsets
here.

Specificity, not sensitivity, is the practical limit. Almost a quarter of variants in which no effect
was detected exceeded SpliceAI's usual screening threshold, and two variants scoring above 0.9 were
normal in minigene assays. In a gene where most variants are private and segregation is often
uninformative, each such flag invites an RNA study that may not be available. This supports treating
a prediction as grounds for a functional test, in the manner of the ClinGen splicing recommendations
[14], rather than as evidence in itself.

Three caveats limit how far these results generalise. The controls are concentrated in a small number
of experiments, and the study-cluster intervals are correspondingly wide. Most evidence comes from
minigene constructs in non-kidney cell lines, which may not reproduce native splicing. And the
benchmark inherits the variability of its sources: one variant altered splicing in one laboratory and
not in another, which places a floor under the accuracy any benchmark of this kind can measure.

Two observations concern pipelines rather than models. A reported *PKD1* pseudoexon is activated only
by two adjacent changes acting together [7]; presented as a single two-base allele both
models detected it clearly, but conventional pipelines separate the changes and filter the common one.
Five indels in this set have no precomputed Atlas score, so any workflow relying on precomputed
values must handle them separately rather than treating absence as a low score.

The translation class is small but internally consistent. Six variants from two laboratories, with
effects from about 42% higher to 87% lower protein output and transcript level shown to be unchanged,
produced predicted values within a range of 0.016 on the AlphaGenome composite, and the two inert
variants were not separated from the four active ones. These models predict RNA, and the experiment
measures protein from an unchanged transcript, so this is a statement about endpoints rather than an
error by either model. It nonetheless bears on practice, because the two *PKD1* elements under
therapeutic investigation act through exactly this mechanism.

Finally, the variants that would matter most for explaining why relatives with the same pathogenic
allele differ in severity are quantitative modifiers. Only one such variant with functional evidence
exists at these loci, and both splicing predictors scored it low. That is a single observation, not a
general claim, but it indicates where a benchmark of this kind is currently thinnest.

## Data Availability Statement

The curated benchmark, its source quotations, the screening log, all score tables and the analysis
plan with its hash record are available in the project repository [URL and archived release
identifier to be inserted on acceptance]. AlphaGenome predictions are redistributed under the
AlphaGenome Output Terms of Use, which permit non-commercial use and scientific publication.
Publisher full texts are not redistributed.

## Code Availability

All analysis code is in the same repository, under an MIT licence that covers the code only, with
the commands needed to regenerate every result, table and figure.

## References

1. Xu P, Huang S, Li J, Zou Y, Gao M, Kang R, et al. A novel splicing mutation in the PKD1 gene causes autosomal dominant polycystic kidney disease in a Chinese family: a case report. BMC Med Genet. 2018;19(1):198. doi:10.1186/s12881-018-0706-6. PMID: 30424739.
2. Lakhia R, Song C, Biggers L, Zumwalt M, Alvarez J, Somasundaram A, et al. Disruption of a six-nucleotide miRNA motif improves PKD1 dosage and ameliorates polycystic kidney disease. bioRxiv. 2025. doi:10.1101/2025.10.22.683929. PMID: 41279527.
3. Chen L, Gao X, Liu X, Zhu Y, Wang D. Translational regulation of PKD1 by evolutionarily conserved upstream open reading frames. RNA Biol. 2025;22(1):1-12. doi:10.1080/15476286.2024.2448387. PMID: 39757590.
4. Wedd L, Hort Y, Patel C, Sayer JA, Rius R, Mallett AJ, et al. PKD1 5'UTR variants are a rare cause of disease in ADPKD and suggest a new focus for therapeutic development. Eur J Hum Genet. 2026;34(1):61-69. doi:10.1038/s41431-025-01949-z. PMID: 41006799.
5. Zhang Z, Blumenfeld J, Ramnauth A, Barash I, Zhou P, Levine D, et al. A common intronic single nucleotide variant modifies PKD1 expression level. Clin Genet. 2022;102(6):483-493. doi:10.1111/cge.14214. PMID: 36029107.
6. Pan Q, Liu Y, Sun X, Lu S, Li L, Shen J. Case Report: Functional validation of a PKD1 c.7489 + 5G>A variant in an ADPKD family. Front Genet. 2026;17:1850997. doi:10.3389/fgene.2026.1850997. PMID: 42621998.
7. Han ST, Rehman AU, Srinivasa ST, Lam A, Ficurilli M, Felice V, et al. Genome and transcriptome sequencing reveal pathogenic activation of a pseudoexon in PKD1 via a de novo-common variant complex allele. Genet Med Open. 2026;4:104407. doi:10.1016/j.gimo.2026.104407. PMID: 42502691.
8. Avsec Ž, Latysheva N, Cheng J, Novati G, Taylor KR, Ward T, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature. 2026;649(8099):1206-1218. doi:10.1038/s41586-025-10014-0. PMID: 41606153.
9. AlphaGenome Atlas team; Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. Google DeepMind; 2026. Available from: https://deepmind.google/science/alphagenome/atlas
10. Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R. VariantValidator: Accurate validation, mapping, and formatting of sequence variation descriptions. Hum Mutat. 2018;39(1):61-68. doi:10.1002/humu.23348. PMID: 28967166.
11. Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, Darbandi SF, Knowles D, Li YI, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015. PMID: 30661751.
12. Google DeepMind. AlphaGenome documentation: recommended splicing score and variant scoring definitions. Available from: https://www.alphagenomedocs.com/faqs.html (accessed 20 September 2026).
13. Schubach M, Maass T, Nazaretyan L, Röner S, Kircher M. CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Res. 2024;52(D1):D1143-D1154. doi:10.1093/nar/gkad989. PMID: 38183205.
14. Walker LC, Hoya M, Wiggins GAR, Lindy A, Vincent LM, Parsons MT, et al. Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing: Recommendations from the ClinGen SVI Splicing Subgroup. Am J Hum Genet. 2023;110(7):1046-1067. doi:10.1016/j.ajhg.2023.06.002. PMID: 37352859.

## Acknowledgments

None.

## Author Contribution Statement

CL conceived the study, curated the data, wrote the code, performed the analysis and wrote the
manuscript.

## Ethical Approval

Ethical approval was not required. The study used published summary data and public databases only.
No human participants were recruited, no individual-level or identifiable data were accessed, and no
new samples were obtained.

## Competing Interests

The author is a director and shareholder of Lawrence Medical Ltd, Patient Thread Ltd and Scira
Diagnostics Ltd, companies developing clinical software and genomic diagnostics, and is named on
United Kingdom patent applications relating to genomic diagnostics, one of which concerns the use of
sequence-to-function model outputs for expression-weighted assessment in transplantation. None of
these interests relates to *PKD1* or *PKD2* variant interpretation, to polycystic kidney disease or to
any tool evaluated here, and none funded this work. The AlphaGenome interface was accessed under a
personal, individual, non-commercial account. No funding was received for this study.

The author used an AI coding assistant (Claude, Anthropic) to write the analysis code, to assemble
the curated data set from published sources and to draft the manuscript. Every variant assignment
carries a verbatim quotation from its source that is checked mechanically against the stored text,
and every variant description was validated independently against the reference genome, as described
in Materials and Methods. The author reviewed and verified all content and takes full responsibility
for it. The tool is not an author.

## Figure Legends

**Figure 1.** Discrimination within the splicing class. Each point is one variant; the black bar is
the median. Positive denotes a demonstrated splicing effect and negative denotes no effect detected in
the reported assay. The AVI panel is restricted to non-coding and synonymous variants, where a high
score cannot arise from an amino-acid change.

**Figure 2.** Mean contribution of each feature group to the AVI score, by variant class. Conservation
accounts for most of the score at the two regulatory elements under therapeutic investigation, whereas
protein features account for most of it among coding positives.

**Figure 3.** Rank within the *PKD1* locus. For each variant, the percentage of all possible *PKD1*
single-nucleotide variants scoring higher, for SpliceAI and for AVI, on logarithmic axes. Points below
the diagonal are ranked higher by SpliceAI. The SpliceAI background is an evenly spaced sample of the
locus; values are clipped at 0.001%.
