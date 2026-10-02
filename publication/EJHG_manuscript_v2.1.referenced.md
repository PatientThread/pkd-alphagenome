# Functional benchmarking of AlphaGenome and SpliceAI at PKD1 and PKD2

Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP

Consultant Nephrologist, 9 Harley Street, London W1G 9QY, United Kingdom

ORCID 0000-0002-8159-0879

Correspondence: Christopher.lawrence3@nhs.net

Running title: Benchmarking sequence models at the PKD loci

Keywords: autosomal dominant polycystic kidney disease; PKD1; PKD2; RNA splicing; variant interpretation; machine learning

---

## Abstract

Sequence-to-function models are entering variant interpretation, and the AlphaGenome Atlas provides
precomputed predictions and a combined impact score (AVI) for every single-nucleotide variant in the
genome. Their behaviour at the polycystic kidney disease loci is unknown, though about
one patient in ten remains unexplained after coding-region testing. From a recorded search of 391
records we curated 128 experimentally tested *PKD1* and *PKD2* variants from 25 sources, each
carrying a verbatim quotation checked against its source and a description validated against GRCh38.
Variants were assigned to the class of experiment producing their evidence, since a reporter of
protein output is not evidence about splicing. The splicing class contained 50 variants shown to
alter splicing and 70 tested with no effect detected. AlphaGenome's documented splicing composite and
SpliceAI had the same observed AUC (0.85 each, 95% CI 0.77 to 0.92 and 0.76 to 0.92) and a paired
difference of 0.000 (-0.033 to 0.037), which does not establish that their errors coincide.
Resampling whole source studies widened both intervals to about 0.71 to 1.00, the controls being
concentrated in few experiments. The composite scores were weaker on the same variants, with AUCs of
0.72 (AVI) and 0.57 (CADD). Where AVI ranked two *PKD1* regulatory elements highly, conservation
contributed most of the score. Ranking comparisons depended on the candidate background. All but one
splicing-class variant with no detected effect was engineered, so these estimates describe assay
behaviour, not clinical specificity. Predictions here are best treated as grounds for an RNA study.

---

## Introduction

Autosomal dominant polycystic kidney disease (ADPKD) affects about 1 in 400 to 1,000 people and is
caused principally by variants in *PKD1* and *PKD2* [1]. Disease follows a reduction in
functional polycystin below a threshold, so the protein produced by the intact allele matters
[2]. Two *PKD1* regulatory elements are consequently under therapeutic investigation: a
microRNA-17 binding site in the 3' untranslated region, whose disruption stabilises the transcript
and raises polycystin-1 in model systems [2], and two upstream open reading frames in the
5' untranslated region that restrain translation [3,4].

Genetic diagnosis remains incomplete: about 10% of individuals with an ADPKD phenotype have no causal
variant identified after coding-region testing, for reasons including assay scope, structural
variants and mosaicism as well as regulatory variation [4,5]. *PKD1* is hard to
sequence and interpret because exons 1 to 32 are repeated six times on chromosome 16 as pseudogenes
[6].

Sequence-to-function models are now directly available to diagnostic laboratories. AlphaGenome
predicts molecular readouts from DNA sequence [7]; the AlphaGenome Atlas stores
precomputed predictions for all possible single-nucleotide variants and adds AVI, a supervised
composite integrating those predictions with protein and conservation features [8]. These are
distinct objects, and a result about AVI is not necessarily a result about the underlying sequence
model. The developers state that they "can only act as part of the evidence chain leading to clinical
diagnoses, and are not sufficient evidence on their own" [8].

General benchmarks cannot answer a locus-specific question, since performance depends on local exon
density, conservation and which mechanisms are common in that gene. We therefore
assembled a locus-specific benchmark of variants shown experimentally to disrupt splicing or
expression together with variants tested with no effect detected, and asked how well these tools
separate the two classes at *PKD1* and *PKD2*, whether AlphaGenome performs differently from the
splicing predictor already in use, and what the scores represent when they are high.

## Materials and Methods

### Search, screening and curation
Two fixed PubMed queries run on 19 September 2026 returned 391 unique records. Titles were screened
for reports of human *PKD1* or *PKD2* variants with functional evidence, and 47 papers were included,
26 as open-access full text, with three further papers and three supplementary files obtained
directly. Of the 47, 25 contributed a variant passing both verification checks, and the 120 labelled
splicing-class variants come from 21 of those. Screening decisions and exclusion reasons are in the
repository.

Variants were labelled by experimental outcome, never by predictor output: a demonstrated effect, "no
effect detected in the reported assay", which is not a claim of clinical benignity, or an
uninterpretable assay. Variants named but not tested, and one given opposite outcomes by two sources,
were excluded with reasons.

Two mechanical checks protect the set. Each entry carries a verbatim quotation, and the verification
script rejects any whose quotation is absent from the stored source text. Every description was validated against the MANE Select transcripts NM_001009944.3
(*PKD1*) and NM_000297.4 (*PKD2*) on GRCh38 with VariantValidator [9], which also
supplied genomic coordinates. All 131 source-level entries passed both checks, describing 128 unique
variants, three reported independently by two papers. Legacy IVS notation was converted from the RefSeq exon
table. Variants published only as protein changes were converted by reading the codon from the
reference coding sequence, accepting a result only where exactly one single-base substitution
produced the stated amino-acid change; all 17 conversions were unique, and the method exactly
reproduced the eight variants whose papers also printed cDNA notation.

### Assignment to an evidence class
Each variant was assigned to the class of experiment that produced its evidence: splicing (patient
RNA, cDNA sequencing, RNA-seq or minigene read out as splicing) or reporter output (5' untranslated
region luciferase reporter measuring protein output). The second is named for its endpoint, not a
mechanism, because for half its members the RNA-level mechanism is unresolved. Where a source
reported both a splicing experiment and a downstream consequence, the splicing experiment determines
the class; rs3874648 is assigned to the splicing class on this rule, reverse-transcription PCR of
patient leukocytes having identified an atypical *PKD1* splice form from a cryptic acceptor
[5]. The primary analysis uses the splicing class alone. The reporter-output class holds
six variants, described individually with no summary statistic. Its membership rule deliberately
differs: it admits only naturally observed human alleles that were tested, excluding engineered
substitutions and alleles described but not assayed, whereas the splicing class applies no such restriction and
contains 92 engineered variants. Because two increase output, a demonstrated effect is recorded with
its direction, not assumed to be a reduction. Transcript abundance was
measured only for the three variants of Chen et al. [3], who reported reporter RNA
unchanged by RT-qPCR; the other three authors state that without RNA studies they cannot exclude
effects on pre-mRNA splicing or transcript stability [4]. Transcript abundance was measured only in the three variants studied by Chen et
al. [3], who reported reporter RNA unchanged by RT-qPCR; the authors of the other three
state explicitly that without RNA studies they cannot exclude effects on pre-mRNA splicing or
transcript stability [4], so no transcript-level claim is made for those variants.

### Scoring
SpliceAI 1.3.1 [10] ran locally with the stock implementation (distance 500, unmasked) on
GRCh38. For AlphaGenome we used the documented splicing composite, max(splice sites) + max(splice
site usage) + max(splice junctions)/5, each the maximum absolute score across all tracks and genes
[11]. AVI scores and its 18 feature attributions came from the Atlas, as did background
distributions for every single-nucleotide variant across each locus (gene span plus 5 kb upstream;
156,573 variants for *PKD1*, 225,576 for *PKD2*). CADD v1.7 [12]
PHRED scores came from the official whole-genome file. The AVI training chromosome list is 1, 4, 7,
8, 10, 13 and 15 and the validation list 2, 5, 11, 14, 17, 20, 22 and X. *PKD1* lies on chromosome
16, which appears on neither list, so *PKD1* results are independent of that split. *PKD2* lies on
chromosome 4, which is a training chromosome; we could not check whether these particular positions
or alleles entered training, so overlap is not excluded rather than demonstrated. AVI results are
therefore reported primarily for *PKD1*, with the pooled values given alongside. A gene-wide SpliceAI background was computed on an evenly
spaced sample of positions across each locus. The step is the locus span divided by a target of
12,000, rounded down to whole bases, so the realised sample slightly exceeds the target: at *PKD1*
the span is 52,189 bp and the step 4 bp, giving 13,048 positions and, with all three substitutions at
each, 39,144 alleles, 25.00% of the 52,191 positions the Atlas covers; at *PKD2* the step is 6 bp,
giving 12,532 positions and 37,596 alleles. The batched implementation reproduces the stock package
exactly on 48 test variants. Coverage is reported per tool, and a variant a tool cannot
score counts as unscored, never as a low score.

### Statistics
Discrimination is the Mann-Whitney area under the curve (AUC) with a stratified bootstrap 95%
interval (2,000 resamples). AlphaGenome and SpliceAI are compared by a paired difference on identical
variants with its own bootstrap interval; no equivalence margin is claimed. Because the controls are
concentrated in a few papers, we also report leave-one-study-out AUCs, a bootstrap resampling whole
source studies, and a paired difference under those same study draws. Because five variants have no
value in the precomputed AVI and CADD resources used, all four tools are compared a second time on
the subset every tool scores, and the variants lost are listed. Composite scores (AVI, CADD) are
additionally reported on non-coding and synonymous variants, a subset excluding variants that
directly change the amino-acid sequence; it does not exclude protein consequences altogether, since a
splicing change can itself alter the protein. Sensitivity and the proportion of no-effect variants
flagged are given with full counts and Wilson intervals. Seeds, tie handling, the treatment of
resamples lacking an outcome class and of the two variants reported by two sources are recorded in
the reviewer package. Variant descriptions follow Human Genome Variation Society nomenclature on the
MANE Select transcripts named above.

An analysis plan was written and hash-stamped before scoring, with eight dated amendments recording
whether the results each could affect already existed. During methodological revision, prompted by an
AI-assisted critique of the draft rather than by journal peer review, the analysis was rebuilt to use
one reference standard per assay class and the documented splicing composite, and one variant was
reassigned between classes. Those changes are post hoc and are identified as such in the plan, which
records the numbers before and after.

The author used an AI coding assistant (Claude, Anthropic) to write the analysis code, to assemble
the curated data set from published sources, to critique successive drafts and to draft the
manuscript. Every variant assignment carries a verbatim quotation from its source that is checked
mechanically against the stored text, and every variant description was validated independently
against the reference genome. The author made and verified all final analytical and
interpretive decisions, reviewed all content and takes full responsibility for it. The tool is not an
author.

## Results

### Composition
The benchmark contains 128 unique variants from 25 sources: 54 with a demonstrated effect, 72 with no
effect detected, one uninterpretable (*PKD1* c.288-12C>A) and one discordant (*PKD1* c.11257C>T,
positive in one minigene study and negative in another). By evidence class, 122 are splicing and six
reporter output (Table 1); the uninterpretable and discordant variants are splicing-class and
excluded, leaving 120 labelled splicing-class variants.

Provenance is strongly associated with outcome, which limits what the comparison can support. Of the
50 splicing-class positives, 26 were carried by patients, one came from population data and 23 were
engineered; of the 70 with no detected effect, 69 were engineered and one was a population variant.
There are no patient-derived negatives. Selective publication is the obvious
candidate explanation, a laboratory finding nothing having little reason to report it, but we did not
test that. Every
estimate below therefore describes how these tools behave against tested sequence, not clinical
specificity in patients.

### Discrimination within the splicing class
Among 50 positives and 70 negatives, the AlphaGenome splicing composite and SpliceAI achieved the
same observed AUC of 0.846 (95% CI 0.770 to 0.916 and 0.761 to 0.919 respectively; Table 2,
Figure 1). The paired difference on identical variants was 0.000 (-0.033 to 0.037), and 0.000 (-0.042
to 0.017) when the bootstrap resampled whole source studies. These data do not resolve which score is
better; the intervals exclude a difference beyond 0.037 favouring AlphaGenome or 0.042 favouring
SpliceAI. Equal observed AUCs do not imply that the two tools err on the same variants. The paired
interval is narrow while the marginal intervals are wide because the two tools' scores are correlated
across variants, a property of the comparison rather than a contradiction.

AVI reached 0.718 (0.620 to 0.808) and CADD 0.565 (0.454 to 0.668), each on the 45 positives and 70
negatives they score. Restricted to non-coding and synonymous splicing-class variants, which excludes
direct amino-acid-changing variants, AVI achieved 0.763 (0.624 to 0.891) on 30 positives and 25
negatives and CADD 0.727 (0.591 to 0.852) on the same variants.

AVI and CADD are reported primarily for *PKD1*, whose chromosome is on neither AVI list. Restricted
to *PKD1*, AVI gives 0.722 (0.609 to 0.827) on 36 positives and 50 negatives against a pooled 0.718;
Table 2 holds every restricted analysis, and the two agree.

Five variants have no AVI or CADD value here: they are insertions or deletions absent from the
precomputed resources used, the Atlas single-nucleotide release and the CADD whole-genome SNV file,
which is a property of those files rather than of either tool. All five have a demonstrated effect,
so the absent values are not missing at random. On the 115 variants scored by all four tools the ordering is unchanged: SpliceAI 0.832 (0.746
to 0.913), AlphaGenome 0.833 (0.748 to 0.911), AVI 0.718 and CADD 0.565. Absent scores are counted as
unscored throughout and never as low scores.

### Robustness
Leave-one-study-out AUCs ranged from 0.805 to 0.920 for SpliceAI and 0.797 to 0.916 for
AlphaGenome, so no single source determines it. Resampling whole studies rather than variants widened the
intervals to 0.715 to 0.999 and 0.707 to 0.993, the more honest expression of uncertainty, since
22 negatives come from one paper and 29 from another. Performance was lower for variants 3 to 20 nucleotides from a splice
site (0.791 and 0.789) than for exonic and untranslated-region variants (0.819 and 0.836), and lower
at *PKD1* than *PKD2*, though that subgroup holds nine positives (Table 3). No AUC is defined within the patient-derived
variants, which contain no negatives; the engineered subgroup, carrying almost all controls, gives
0.734 and 0.756.

### Operating characteristics
Within the splicing class, at the conventional SpliceAI threshold of 0.2, 41 of 50 positives were
flagged and 9 missed, and 17 of 70 no-effect variants were flagged and 53 not:
sensitivity 82.0% (95% CI 69.2 to 90.2) with 24.3% (15.8 to 35.5) of no-effect variants flagged. At
0.5 the counts are 35 flagged and 15 missed among positives, and 7 flagged and 63 not among
no-effect variants: sensitivity 70.0% (56.2 to 80.9) with 10.0% (4.9 to 19.2) flagged. Because all but one no-effect variant was
engineered, the second figure is the rate at which the tool flags constructs in which no effect was
detected, not a false-positive rate in practice. Two examples are informative:
*PKD1* c.11713-3C>G scores 0.95 and c.9201+5G>A 0.96, and neither altered splicing in a minigene.
AVI has no validated threshold, so no equivalent operating point can be quoted for it.

### What a high AVI score represents
At the two *PKD1* regulatory elements under therapeutic investigation, AVI ranked changes highly
within their own untranslated regions (median 92nd percentile for the microRNA-17 site and 95th for
the upstream open reading frame start codons), but for all 36 possible substitutions the largest
contributor was evolutionary conservation, with AlphaGenome's own functional predictions contributing
0.01 to 0.05 (Figure 2). Among coding positives with AVI scores, 8 of 15 were driven by protein
features and 4 by splicing, so a high composite score at such a variant is not evidence that a
splicing effect was detected.

### Ranking within the locus
These analyses are restricted to the splicing class. The unit is the alternate allele: both
backgrounds carry all three substitutions at each position, and both tools are matched on one rule,
the proportion of their own background flagged. Counts, thresholds and achieved proportions are in
the reviewer package.

Against the whole *PKD1* locus, SpliceAI at 0.5 flags 1.04% of 39,144 sampled alleles, and the AVI
threshold flagging 1.04% of all 156,573 alleles is 1.550. Non-coding and synonymous positives sat at
a median 99.4th percentile by SpliceAI and 94.2nd by AVI; SpliceAI recovered 20 of 25 such positives
and AVI 5, flagging 3 and 2 of 18 negatives respectively (Figure 3a).

Against a background of non-coding positions only, the same proportion rule gives 0.56% of 28,401
sampled alleles for SpliceAI and an AVI threshold of 0.665. SpliceAI then recovered 17 of 21
non-coding positives and AVI 18 of 21, both flagging 3 of 16 negatives, with median ranks of 99.9th and
99.8th percentile (Figure 3b). The two negative denominators differ because the whole-locus
comparison admits synonymous variants and the non-coding comparison does not.

These are not a wrong value and a corrected one; they answer different questions, and the difference
is itself informative. AVI is a whole-genome pathogenicity score and ranks amino-acid-changing
variants highly, so its top slice of the whole locus is filled by coding positions: of the 1,632
alleles above the whole-locus threshold only 2.1% are non-coding, although 72.6% of the locus is.
Ranking a non-coding candidate against the whole gene makes it compete with coding variants it would
never be confused with; restricting the background removes that competition, and the two tools then
rank these variants alike. Separately, and recorded in the amendment log, both comparisons matched
raw counts rather than proportions before 20 September 2026, setting a fourfold stricter AVI
threshold; every figure above is post-correction.

### The reporter-output class
Six 5' untranslated region variants, from two laboratories, were assayed by luciferase reporters
measuring protein output (Table 4). This endpoint neither model predicts, so the comparison tests the
boundary of what these tools claim, not their accuracy. Measured effects span roughly 42%
higher to 87% lower output, yet the AlphaGenome splicing composite varied only between 0.030 and
0.046 across all six, and SpliceAI scored five of six at zero. The two variants with no detected
effect did not score lower than those with one: c.-209G>A scored 0.044, above c.-69dupG, which
reduced output (0.033). Chen et al. accompanied their three with RT-qPCR showing reporter RNA
unchanged, so for those the effect is attributable to translation; Wedd et al. measured no RNA and
state that splicing and transcript-stability effects cannot be excluded, so their three act on
protein output by an undetermined route. Both laboratories examined the uORF2 initiation region, Chen
et al. reporting the natural allele c.-20A>G, included here, and Wedd et al. constructing engineered
substitutions there, excluded under the membership rule in Materials and Methods.

rs3874648, which activates a cryptic acceptor and lowers full-length *PKD1*, belongs to the splicing
class and is in the analyses above. Both splicing predictors scored it low, and it is one of the nine
positives missed at a SpliceAI threshold of 0.2. It is the only variant here whose source reports it
as an expression modifier rather than a loss-of-function allele.

## Discussion

For the class of variant represented in this benchmark, a splicing effect in or near *PKD1*, these
data do not distinguish AlphaGenome's splicing composite from SpliceAI: the paired difference is
0.000 and the interval excludes differences beyond about 0.04. How common such variants are among
unexplained cases is not established here, since the variants were selected for having been tested.
This is a benchmark on curated experimental evidence, not a clinical validation: it licenses neither
tool for diagnostic use, and a laboratory changing tools would weigh curation, thresholds and
reporting as well as rank-order accuracy.

The composite score behaves differently from the splicing-specific outputs, and the difference
matters in practice. AVI was weaker here, and when it did rank regulatory sites highly the score came
mainly from conservation, for the individual variant and not only on average: conservation was the
largest of the four feature groups for all 36 substitutions at the two elements, exceeding the next
by at least 0.31 in every case, while no other group reached 0.08 in any. Conservation is legitimate evidence, but it is not mechanism and is not specific to the
effect a laboratory seeks to establish, so feature attributions should be inspected whenever a
composite score informs an interpretation. The converse holds too: a protein-truncating variant
ranked top of a gene by AVI says nothing about splicing.

Both error types are substantial. At 0.2, SpliceAI missed 9 of 50 demonstrated effects, and about a
quarter of the constructs with no detected effect exceeded that threshold, two scoring above 0.9
while normal in minigene assays. Since the no-effect variants are almost all engineered, that second
figure describes behaviour on tested sequence, not the rate at which a laboratory would be misled;
establishing the latter needs patient-derived variants tested and found negative, seldom published. In a gene where most variants are private and segregation is often uninformative,
each flag nonetheless invites an RNA study that may not be available. The ClinGen recommendations
provide for exactly this: a prediction may contribute PP3 or BP4 at calibrated strength and may
justify a functional study, which our results support rather than treating a score as a statement of
mechanism [13].

Four caveats limit how far these results generalise. The controls are concentrated in a few
experiments, so the study-cluster intervals are wide. All but one no-effect variant was engineered,
so clinical specificity is not estimated. Most evidence comes from minigene constructs in non-kidney
cell lines, which may not reproduce native splicing. And the benchmark inherits its sources'
variability: one variant altered splicing in one laboratory and not another, placing a floor under
the accuracy any such benchmark can measure. For AVI, *PKD2* lies on a training chromosome, so only
the *PKD1* results are independent of that split.

Two observations concern pipelines rather than models. A reported *PKD1* pseudoexon is activated
only by the de novo c.7066-131G>A together with the adjacent common rs3874658 in cis; the original
report found the de novo change alone scored weakly by SpliceAI, 0.12 acceptor gain and 0.06 donor
gain, while the dinucleotide allele gave a donor gain of 0.78 [14]. Our local SpliceAI
reproduces both exactly, and AlphaGenome behaves the same way: composite 3.21 for the two-base allele
against 0.12 for the de novo change alone. A workflow annotating the two substitutions separately, or
not reporting deep intronic positions at adequate coverage, never puts the combined allele in front
of either model. Separately, five indels here have no precomputed Atlas score, so any workflow
relying on precomputed values must handle them rather than treat absence as a low score.

The reporter-output class is small but consistent. Six variants from two laboratories, with effects
from 42% higher to 87% lower protein output, produced predictions within a range of 0.016 on the
AlphaGenome composite, and the two with no detected effect were not separated from the four with one.
These models predict RNA and the experiment measures protein, so this concerns endpoints rather than
model error, and it is clean only for the three variants whose reporter RNA was shown unchanged. It nonetheless bears on practice, because the two *PKD1*
elements under therapeutic investigation act through different mechanisms, and only one of them is
outside what these models address. The upstream open reading frames restrain translation of a message
whose abundance is unchanged, an endpoint no sequence-to-function model evaluated here predicts.
Disrupting the microRNA-17 motif, by contrast, stabilises *Pkd1* messenger RNA and raises
polycystin-1 [2]; that is an RNA-abundance effect, and RNA abundance is within
AlphaGenome's scope, so it is worth asking what the model predicts there. This is exploratory rather than a
replication: the published genomic edit was made in mouse *Pkd1* while the model was queried on human
sequence, and that study's human experiments used oligonucleotide targeting rather than an edit. We
replaced the six bases of the human seed match, GRCh38 chr16:2,089,569 AAAGTG>GGGACA in the *PKD1* 3'
untranslated region, in the recommended 1,048,576 bp context. The predicted change is small: the largest absolute change across 371 RNA-seq
tracks was 0.017 on a log2 scale, about 1%, with a mean of 0.003 across all tracks and across the
eight tracks matching our curated kidney biosample list. The engineered uORF start-codon changes gave
0.093 and 0.069. These are the model's own outputs rather than feature attributions, because a small
contribution to a composite score is not the same statement as a small molecular prediction. Given the differences in species and perturbation, this is a prompt for an
experiment, not a negative result about the model.

Finally, one variant here, rs3874648, is an expression modifier: it lowers full-length *PKD1*
without abolishing it [5]. Both splicing predictors scored it low, and it is the only
graded expression effect in the set. Its source also examined clinical severity and found no association
between genotype and total kidney volume, Mayo imaging class, age at kidney failure or hypertension,
but held eight homozygotes, so the null does not exclude an effect.
Three of those eight carried no detectable *PKD1* or *PKD2* mutation, and homozygotes were
over-represented among mutation-negative patients (p = 0.009). That is why variants of this kind matter for
unexplained disease, and why it is unwelcome that both tools missed this one. A single observation,
but it marks where a benchmark of this kind is thinnest.

## Data Availability Statement

The curated benchmark, its source quotations, the screening log, all score tables, every figure and
the analysis plan with its eight dated amendments and hash record are available in the project
repository [URL and archived release identifier to be inserted on acceptance]. A reviewer package
(PKD_AlphaGenome_reviewer_package.zip) is supplied with this submission and contains the labelled
data set, per-variant analysis-set membership flags, all stored scores, the ranking manifest with
every background, threshold and recovery count, the amendment record and the commands needed to
regenerate every number, table and figure in this manuscript without repeating any model query. AlphaGenome predictions are redistributed under the
AlphaGenome Output Terms of Use, which permit non-commercial use and scientific publication.
Publisher full texts are not redistributed.

## Code Availability

All analysis code is in the same repository, under an MIT licence that covers the code only, with
the commands needed to regenerate every result, table and figure.

## References

1. Xu P, Huang S, Li J, Zou Y, Gao M, Kang R, et al. A novel splicing mutation in the PKD1 gene causes autosomal dominant polycystic kidney disease in a Chinese family: a case report. BMC Med Genet. 2018;19(1):198. doi:10.1186/s12881-018-0706-6. PMID: 30424739.
2. Lakhia R, Song C, Biggers L, Zumwalt M, Alvarez J, Somasundaram A, et al. Disruption of a six-nucleotide miRNA motif improves PKD1 dosage and ameliorates polycystic kidney disease. Nucleic Acids Res. 2026;54(2):gkaf1538. doi:10.1093/nar/gkaf1538. PMID: 41558825.
3. Chen L, Gao X, Liu X, Zhu Y, Wang D. Translational regulation of PKD1 by evolutionarily conserved upstream open reading frames. RNA Biol. 2025;22(1):1-12. doi:10.1080/15476286.2024.2448387. PMID: 39757590.
4. Wedd L, Hort Y, Patel C, Sayer JA, Rius R, Mallett AJ, et al. PKD1 5'UTR variants are a rare cause of disease in ADPKD and suggest a new focus for therapeutic development. Eur J Hum Genet. 2026;34(1):61-69. doi:10.1038/s41431-025-01949-z. PMID: 41006799.
5. Zhang Z, Blumenfeld J, Ramnauth A, Barash I, Zhou P, Levine D, et al. A common intronic single nucleotide variant modifies PKD1 expression level. Clin Genet. 2022;102(6):483-493. doi:10.1111/cge.14214. PMID: 36029107.
6. Pan Q, Liu Y, Sun X, Lu S, Li L, Shen J. Case Report: Functional validation of a PKD1 c.7489 + 5G>A variant in an ADPKD family. Front Genet. 2026;17:1850997. doi:10.3389/fgene.2026.1850997. PMID: 42621998.
7. Avsec Ž, Latysheva N, Cheng J, Novati G, Taylor KR, Ward T, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature. 2026;649(8099):1206-1218. doi:10.1038/s41586-025-10014-0. PMID: 41606153.
8. AlphaGenome Atlas team; Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. Google DeepMind; 2026. Available from: https://deepmind.google/science/alphagenome/atlas
9. Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R. VariantValidator: Accurate validation, mapping, and formatting of sequence variation descriptions. Hum Mutat. 2018;39(1):61-68. doi:10.1002/humu.23348. PMID: 28967166.
10. Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, Darbandi SF, Knowles D, Li YI, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015. PMID: 30661751.
11. Google DeepMind. AlphaGenome documentation: recommended splicing score and variant scoring definitions. Available from: https://www.alphagenomedocs.com/faqs.html (accessed 20 September 2026).
12. Schubach M, Maass T, Nazaretyan L, Röner S, Kircher M. CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Res. 2024;52(D1):D1143-D1154. doi:10.1093/nar/gkad989. PMID: 38183205.
13. Walker LC, Hoya M, Wiggins GAR, Lindy A, Vincent LM, Parsons MT, et al. Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing: Recommendations from the ClinGen SVI Splicing Subgroup. Am J Hum Genet. 2023;110(7):1046-1067. doi:10.1016/j.ajhg.2023.06.002. PMID: 37352859.
14. Han ST, Rehman AU, Srinivasa ST, Lam A, Ficurilli M, Felice V, et al. Genome and transcriptome sequencing reveal pathogenic activation of a pseudoexon in PKD1 via a de novo-common variant complex allele. Genet Med Open. 2026;4:104407. doi:10.1016/j.gimo.2026.104407. PMID: 42502691.

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
personal, individual, non-commercial account.

## Funding

This study received no funding from any public, commercial or not-for-profit body. No article
processing charge was paid from a grant.

## Use of Artificial Intelligence

The author's use of an AI assistant, and the mechanical checks that constrain it, are described in
Materials and Methods. The tool does not meet authorship criteria and is not an author.

---

## Figure Legends

**Figure 1.** Discrimination within the splicing class. Each point is one variant; the black bar is
the median. Positive denotes a demonstrated splicing effect and negative denotes no effect detected in
the reported assay. The SpliceAI and AlphaGenome panels show 50 positives and 70 negatives. The AVI
panel is restricted to non-coding and synonymous variants, which excludes direct amino-acid-changing
variants, and shows 30 positives and 25 negatives; it is a different subset and is labelled as such,
so the panels should not be read as a direct three-way comparison on identical variants. The
four-tool comparison on identical variants is in Table 2.

**Figure 2.** Contribution of each feature group to the AVI score, by variant class, as the mean of
the grouped feature attributions returned by the Atlas. Sixteen of the 18 attributions are grouped by summation
into splicing (merged splicing), expression and chromatin (nine features), protein (four) and
conservation (two); the two remaining indicators flag insertions and deletions and are not applicable
to substitutions. Values are the signed attributions as returned, not absolute values. Conservation accounts for most of the score at the two regulatory
elements under therapeutic investigation, and is the largest group for every one of the 36
substitutions individually, not only on average. Protein features account for most of the score among
coding positives.

**Figure 3.** Rank within the *PKD1* locus, for splicing-class variants only. **(a)** Against the whole
locus: for each variant, the percentage of all possible *PKD1* single-nucleotide variants scoring
higher, for SpliceAI and for AVI, on logarithmic axes. Points below the diagonal are ranked higher by
SpliceAI. **(b)** Against a non-coding background, matched on the proportion of positions flagged
rather than the raw count, which is the comparison we regard as fair to a whole-genome pathogenicity
score asked a non-coding question. The SpliceAI background is an evenly spaced sample of the locus;
values below the resolution of the sampled
background, one allele in 39,144 or 0.0026%, are censored at that floor and plotted there.
