<!-- medRxiv preprint version. Free to deposit; no peer review certification.
     Requirements checked 1 Oct 2026 at medrxiv.org/submit-a-manuscript: title,
     author names and affiliations and abstract on the title page; a competing
     interest declaration covering the past 36 months; a funding statement.
     Key learning points removed (an NDT-specific box). -->

# Predicting splicing variants in ADPKD: a functional benchmark of AlphaGenome and SpliceAI at PKD1 and PKD2

Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP

Consultant Nephrologist, 9 Harley Street, London W1G 9QY, United Kingdom

ORCID 0000-0002-8159-0879

Correspondence: Christopher.lawrence3@nhs.net

Running head: Predicting splicing variants in ADPKD

Keywords: autosomal dominant polycystic kidney disease; PKD1; RNA splicing; variant interpretation; machine learning

---

## Abstract

**Background and hypothesis.** Genetic testing in autosomal dominant polycystic kidney disease
(ADPKD) informs prognosis, trial eligibility and living-donor assessment, yet about one patient in
ten has no causal variant found after coding-region testing. Laboratories rely on computational
predictors to decide which non-coding variant deserves an RNA study. SpliceAI is in routine use;
AlphaGenome and its precomputed variant impact score (AVI) are newly available, and their behaviour
at *PKD1* and *PKD2* is unknown.

**Methods.** From a recorded search of 391 records we curated 128 *PKD1* and *PKD2* variants tested
experimentally in 25 published studies, each carrying a verbatim quotation from its source, checked
mechanically, and a description validated against GRCh38. Variants were assigned to the class of
experiment producing their evidence, because a reporter of protein output is not evidence about
splicing. The primary analysis compared AlphaGenome's documented splicing composite with
SpliceAI on identical variants, as a paired difference with bootstrap intervals.

**Results.** The splicing class contained 50 variants shown to alter splicing and 70 tested with no
effect detected. AlphaGenome and SpliceAI had the same observed area under the curve, 0.85 each (95%
confidence interval 0.77 to 0.92 and 0.76 to 0.92), with a paired difference of 0.000 (-0.033 to
0.037). AVI was weaker (0.72) and CADD weaker still (0.57). At the usual SpliceAI threshold of 0.2,
41 of 50 demonstrated effects were flagged, together with 17 of 70 variants in which no effect was
detected. Neither tool separated six 5' untranslated region variants whose effects on protein output
ranged from 42% higher to 87% lower. Which tool ranked candidates better depended on the comparison
set used.

**Conclusions.** AlphaGenome offers no accuracy gain over SpliceAI at these loci. Both miss about one
demonstrated effect in five. Predictions should prompt an RNA study, not stand as evidence of
mechanism.

---

## Introduction

Autosomal dominant polycystic kidney disease (ADPKD) affects about 1 in 400 to 1,000 people and is
caused principally by variants in *PKD1* and *PKD2*.[PMID 30424739] Genetic diagnosis now carries
practical weight in nephrology: it refines prognosis alongside imaging, determines eligibility for
trials and for tolvaptan in some health systems, and resolves ambiguous living-donor assessment in
affected families.

Diagnosis is nonetheless incomplete. About 10% of individuals with an ADPKD phenotype have no causal
variant identified after coding-region testing, for reasons including assay scope, structural
variants and mosaicism as well as regulatory variation.[PMID 36029107; 41006799] *PKD1* is difficult to sequence
and to interpret because exons 1 to 32 are repeated six times on chromosome 16 as pseudogenes.[PMID 42621998]
Disease follows a reduction in functional polycystin below a threshold, so the amount of protein made
by the intact allele matters,[PMID 41558825] and two *PKD1* regulatory elements are consequently under
therapeutic investigation: a microRNA-17 binding site in the 3' untranslated region[PMID 41558825] and
two upstream open reading frames that restrain translation.[PMID 39757590; 41006799]

The practical question for a laboratory is narrower. Faced with a non-coding variant of uncertain
significance in a patient with cystic kidney disease, which candidates justify the cost and delay of
an RNA study? SpliceAI is in routine use for this purpose. AlphaGenome, which predicts molecular
readouts directly from sequence,[REF:ALPHAGENOME] and its companion Atlas, which stores precomputed
predictions for every possible single-nucleotide variant and adds a supervised combined impact score
(AVI),[REF:ATLAS] are newly available. These are distinct objects, and a result about AVI is not
necessarily a result about the underlying model; the developers state that such predictions "can only
act as part of the evidence chain leading to clinical diagnoses, and are not sufficient evidence on
their own".[REF:ATLAS]

General benchmarks cannot answer a locus-specific question, because performance depends on local exon
density, conservation and which mechanisms are common in that gene. We therefore assembled a
benchmark of *PKD1* and *PKD2* variants that had actually been tested in a laboratory, and asked how
well these tools separate variants with a demonstrated effect from those in which none was detected.

## Materials and Methods

### Search, screening and curation

Two fixed PubMed queries run on 19 September 2026 returned 391 unique records. Titles were screened
for reports of human *PKD1* or *PKD2* variants with functional evidence, and 47 papers were included,
26 as open-access full text, with three further papers and three supplementary files obtained
directly. Of the 47, 25 contributed a variant passing both verification checks.

Variants were labelled by experimental outcome, never by predictor output: a demonstrated effect, "no
effect detected in the reported assay", which is not a claim of clinical benignity, or an
uninterpretable assay. Variants named but not tested, and one given opposite outcomes by two sources,
were excluded with reasons recorded.

Two mechanical checks protect the set. Each entry carries a verbatim quotation from its source, and a
verification script rejects any entry whose quotation is absent from the stored source text. Every
description was validated against the MANE Select transcripts NM_001009944.3 (*PKD1*) and NM_000297.4
(*PKD2*) on GRCh38 with VariantValidator,[REF:VARIANTVALIDATOR] which also supplied genomic coordinates. All
131 source-level entries passed both checks, describing 128 unique variants.

### Assignment to an evidence class

Each variant was assigned to the class of experiment that produced its evidence: splicing (patient
RNA, cDNA sequencing, RNA sequencing or minigene read out as splicing) or reporter output (5'
untranslated region luciferase reporter measuring protein output). Where a source reported both a
splicing experiment and a downstream consequence, the splicing experiment determines the class. The
primary analysis uses the splicing class alone; the six reporter-output variants are described
individually, with no summary statistic. The reporter-output class admits only naturally observed
human alleles that were tested, whereas the splicing class applies no such restriction and contains
92 engineered variants; it is therefore not a random sample of tested 5' untranslated region alleles.
Full class definitions and membership rules are in Supplementary Methods.

### Scoring and statistics

SpliceAI 1.3.1[REF:SPLICEAI] ran locally with the stock implementation (distance 500, unmasked) on
GRCh38. For AlphaGenome we used the documented splicing composite, max(splice sites) + max(splice
site usage) + max(splice junctions)/5, each the maximum absolute score across all tracks and
genes.[REF:ATLASDOCS] AVI scores, feature attributions and locus backgrounds came from the Atlas; CADD
v1.7[REF:CADD17] PHRED scores from the official whole-genome file. AVI was trained on chromosomes
that include chromosome 4 but not chromosome 16, so AVI results are reported primarily for *PKD1*.

Discrimination is the Mann-Whitney area under the curve (AUC) with a stratified bootstrap 95%
interval (2,000 resamples). AlphaGenome and SpliceAI are compared by a paired difference on identical
variants. Because the controls are concentrated in a few papers, we also report leave-one-study-out
AUCs and a bootstrap resampling whole source studies. Coverage is reported per tool, and a variant a
tool cannot score counts as unscored, never as a low score. Sampling, resampling seeds and full
specifications are in Supplementary Methods.

The author used an AI coding assistant (Claude, Anthropic) to write the analysis code, assemble the
curated data set from published sources, critique successive drafts and draft the manuscript. Every
variant assignment carries a verbatim quotation checked mechanically against the stored text, and
every description was validated independently against the reference genome. The author made and
verified all final analytical and interpretive decisions and takes full responsibility for the work.

## Results

### Composition

The benchmark contains 128 unique variants from 25 sources: 54 with a demonstrated effect, 72 with no
effect detected, one uninterpretable and one discordant. By evidence class, 122 are splicing and six
reporter output (Table 1), leaving 120 labelled splicing-class variants for analysis.

Provenance is strongly associated with outcome, which limits what the comparison can support. Of the
50 splicing-class positives, 26 were carried by patients, one came from population data and 23 were
engineered;
of the 70 with no detected effect, 69 were engineered and one was a population variant. There are no
patient-derived negatives. Selective publication is the obvious candidate explanation, a laboratory
finding nothing having little reason to report it, but we did not test that. Every estimate below
therefore describes how these tools behave against experimentally tested sequence, not clinical
specificity in patients.

### Discrimination

Among 50 positives and 70 negatives, the AlphaGenome splicing composite and SpliceAI had the same
observed AUC of 0.846 (95% CI 0.770 to 0.916 and 0.761 to 0.919; Table 2, Figure 1). The paired
difference on identical variants was 0.000 (-0.033 to 0.037), and 0.000 (-0.042 to 0.017) when the
bootstrap resampled whole source studies. Equal observed AUCs do not imply that the two tools err on
the same variants.

AVI reached 0.718 (0.620 to 0.808) and CADD 0.565 (0.454 to 0.668), each on the 45 positives and 70
negatives they score. Restricted to *PKD1*, AVI gives 0.722 (0.609 to 0.827). Five variants have no
AVI or CADD value because they are insertions or deletions absent from the precomputed resources
used; all five have a demonstrated effect, so the absent values are not missing at random. On the 115 variants scored by all four tools the
ordering is unchanged: SpliceAI 0.832, AlphaGenome 0.833, AVI 0.718
and CADD 0.565.

Leave-one-study-out AUCs ranged from 0.805 to 0.920 for SpliceAI and 0.797 to 0.916 for AlphaGenome.
Resampling whole studies widened the intervals to 0.715 to 0.999 and 0.707 to 0.993, the more honest
expression of uncertainty given that 22 negatives come from one paper and 29 from another (Table 3).

### Operating characteristics

At the conventional SpliceAI threshold of 0.2, 41 of 50 positives were flagged and 9 missed, and 17
of 70 no-effect variants were flagged and 53 not: sensitivity 82.0% (95% CI 69.2 to 90.2) with 24.3%
(15.8 to 35.5) of no-effect variants flagged. At 0.5 the counts are 35 and 15 among positives, and 7
and 63 among no-effect variants: sensitivity 70.0% (56.2 to 80.9) with 10.0% (4.9 to 19.2) flagged.
Because all but one no-effect variant was engineered, the second figure is the rate at which the tool
flags constructs in which no effect was detected, not a false-positive rate in practice. Two examples
are informative: *PKD1* c.11713-3C>G scores 0.95 and c.9201+5G>A 0.96, and neither altered splicing
in a minigene.

### What a high composite score represents

At the two *PKD1* regulatory elements under therapeutic investigation, AVI ranked changes highly
within their own untranslated regions, but for all 36 possible substitutions the largest contributor
was evolutionary conservation, exceeding the next feature group by at least 0.31 in every case, while
AlphaGenome's own functional predictions contributed 0.01 to 0.05 (Figure 2). Among coding positives
with AVI scores, 8 of 15 were driven by protein features and 4 by splicing, so a high composite score
at such a variant is not evidence that a splicing effect was predicted.

### Ranking depends on the comparison set

Ranked against the whole *PKD1* locus, SpliceAI at 0.5 flags 1.04% of 39,144 sampled alleles;
matching that proportion, SpliceAI recovered 20 of 25 non-coding and synonymous positives and AVI 5
(Figure 3a). Against a background of non-coding positions only, the same proportion rule gives 0.56%,
and SpliceAI recovered 17 of 21 and AVI 18 of 21 (Figure 3b).

These answer different questions, and the difference is informative. AVI ranks amino-acid-changing
variants highly, so its top slice of the whole locus is filled by coding positions: of the 1,632
alleles above the whole-locus threshold only 2.1% are non-coding, although 72.6% of the locus is.
Ranking a non-coding candidate against the whole gene makes it compete with coding variants it would
never be confused with.

### Variants acting after transcription

Six 5' untranslated region variants, from two laboratories, were assayed by luciferase reporters
measuring protein output (Table 4). This endpoint neither model predicts. Measured effects span
roughly 42% higher to 87% lower output, yet the AlphaGenome splicing composite varied only between
0.030 and 0.046 across all six, and SpliceAI scored five of six at zero. The two variants with no
detected effect did not score lower than those with one. Reporter RNA was measured in only one of the
two studies, which reported it unchanged; the other performed no RNA study, so for its three variants
the route to altered protein output is undetermined.

## Discussion

For the class of variant represented in this benchmark, a splicing effect in or near *PKD1*, these
data do not distinguish AlphaGenome's splicing composite from SpliceAI: the paired difference is
0.000 and the interval excludes differences beyond about 0.04. Neither tool was more accurate here, so accuracy alone does not
argue for changing tool, and a laboratory adopting AlphaGenome should expect equivalent rather than
improved discrimination. This is a benchmark on curated experimental evidence, not a clinical
validation, and it licenses neither tool for diagnostic use.

Both error types are substantial and relevant to practice. At 0.2, SpliceAI missed 9 of 50
demonstrated effects, and about a quarter of the constructs with no detected effect exceeded that
threshold, two scoring above 0.9 while normal in minigene assays. In a gene where most variants are
private and segregation is often uninformative, each flag invites an RNA study that may not be
available, and each miss leaves a family without a genetic diagnosis. The ClinGen recommendations
provide for exactly this: a prediction may contribute evidence at calibrated strength and may justify
a functional study,[REF:CLINGEN] which our results support, rather than treating a score as a
statement of mechanism.

The composite score behaves differently from the splicing-specific outputs, and the difference
matters at the bench. AVI was weaker here, and when it did rank regulatory sites highly the score
came mainly from conservation, for every individual variant and not only on average. Conservation is
legitimate evidence but it is not mechanism, so feature attributions should be inspected whenever a
composite score informs an interpretation. Equally, apparent superiority in ranking depended on the
comparison set rather than on the tool, which matters for any workflow that prioritises a candidate
list.

Two observations concern pipelines rather than models. A reported *PKD1* pseudoexon is activated only
by a de novo change together with an adjacent common variant in cis; the original report found the de
novo change alone scored 0.12 by SpliceAI while the dinucleotide allele gave 0.78.[PMID 42502691] Our
local SpliceAI reproduces both exactly, and AlphaGenome behaves the same way. A workflow annotating
the two substitutions separately, or not reporting deep intronic positions at adequate coverage,
never puts the combined allele in front of either model. Separately, five insertions or deletions
here have no precomputed score, so any workflow relying on precomputed values must handle them rather
than treat absence as a low score.

The reporter-output class defines a blind spot with direct therapeutic relevance. Six variants with
effects from 42% higher to 87% lower protein output produced predictions within a range of 0.016, and
the two with no detected effect were not separated from the four with one. These models predict RNA
and the experiment measures protein, so this concerns endpoints rather than model error. It bears on
practice because the two *PKD1* elements under drug development act through post-transcriptional
mechanisms, and these are not the same as one another: the upstream open reading frames restrain
translation of a message whose abundance is unchanged, whereas disrupting the microRNA-17 motif
stabilises *Pkd1* messenger RNA.[PMID 41558825] The second is an RNA-abundance effect and therefore
within AlphaGenome's scope; an exploratory query replacing the six bases of the human seed match
predicted a change of 0.017 on a log2 scale, about 1%. That query is on human sequence whereas the
published edit was made in mouse, so it is an analogous construct rather than a replication.

Four caveats limit generalisation. The controls are concentrated in few experiments, so the
study-cluster intervals are wide. All but one no-effect variant was engineered, so clinical
specificity is not estimated. Most evidence comes from minigene constructs in non-kidney cell lines,
which may not reproduce native splicing. And the benchmark inherits its sources' variability: one
variant altered splicing in one laboratory and not another, placing a floor under the accuracy any
such benchmark can measure.

Finally, one variant here, rs3874648, is an expression modifier: it lowers full-length *PKD1* without
abolishing it.[PMID 36029107] Both splicing predictors scored it low, and it is the only graded
expression effect in the set. Its source examined clinical severity and found no association with
total kidney volume, Mayo imaging class, age at kidney failure or hypertension, but held eight
homozygotes, so the null does not exclude an effect. Three of those eight carried no detectable
*PKD1* or *PKD2* mutation, and homozygotes were over-represented among mutation-negative
patients.[PMID 36029107] That is why variants of this kind matter for unexplained cystic kidney disease,
and why it is unwelcome that both tools missed this one.

## Acknowledgements

None.

## Declaration of generative AI use

During the preparation of this work the author used Claude (Opus 5, Anthropic) to curate and
cross-check the literature search, to write and revise the analysis code, and to draft and revise the
text. Every truth-set entry carries a verbatim quotation from its source and was checked
mechanically against that source and against GRCh38 by VariantValidator; every reported number is
regenerated from the deposited score tables by a script, independently of how that code was written.
After using this tool the author reviewed and edited the content and takes full responsibility for
the content of the publication. The tool is not listed as an author and does not meet authorship
criteria.

## Conflict of Interest Statement

In the past 36 months the author has been a director and shareholder of Lawrence Medical Ltd, Patient Thread Ltd and Scira
Diagnostics Ltd, companies developing clinical software and genomic diagnostics, and is named on
United Kingdom patent applications relating to genomic diagnostics, one of which concerns the use of
sequence-to-function model outputs for expression-weighted assessment in transplantation. None of
these interests relates to *PKD1* or *PKD2* variant interpretation, to polycystic kidney disease or
to any tool evaluated here, and none funded this work. The AlphaGenome interface was accessed under a
personal, individual, non-commercial account. The results presented in this paper have not been
published previously in whole or part.

## Authors' Contributions

CL is the sole author. CL conceived the study, curated and verified the data, directed the analysis,
verified every reported result against the deposited score tables, and prepared and approved the
final manuscript. Generative AI assistance is declared above; CL takes full responsibility for the
content.

## Funding

This study received no funding from any public, commercial or not-for-profit body.

## Data Availability Statement

The curated benchmark, its source quotations, the screening log, all score tables, every figure and
the analysis plan with its dated amendments are archived at Zenodo, doi:10.5281/zenodo.23103370,
and in the project repository at https://github.com/PatientThread/pkd-alphagenome (release
v1.0-jmg). The archive includes the labelled data set, per-variant analysis-set membership flags,
all stored scores, the ranking manifest and a script that regenerates every number in this
manuscript without repeating any model query. AlphaGenome predictions are redistributed under the
AlphaGenome Output Terms of Use, which permit non-commercial use and scientific publication;
SpliceAI-derived and CADD-derived files carry their own upstream terms, set out in the archive.
Publisher full texts are not redistributed: each truth-set entry carries a verbatim quotation and a
PMID so a reader can verify it against the source.

## References

[Generated by analysis/16_references.py]

## Figure Legends

**Figure 1.** Discrimination within the splicing class. Each point is one variant; the black bar is
the median. Positive denotes a demonstrated splicing effect and negative denotes no effect detected
in the reported assay. The SpliceAI and AlphaGenome panels show 50 positives and 70 negatives. The
AVI panel is restricted to non-coding and synonymous variants, which excludes variants directly
changing the amino-acid sequence, and shows 30 positives and 25 negatives; it is a different subset,
so the panels should not be read as a direct three-way comparison.

**Figure 2.** Contribution of each feature group to the AVI score, by variant class, as the mean of
the grouped feature attributions returned by the Atlas. Conservation accounts for most of the score
at the two regulatory elements under therapeutic investigation, and is the largest group for every
one of the 36 substitutions individually. Protein features account for most of the score among coding
positives.

**Figure 3.** Rank within the *PKD1* locus, splicing-class variants only. **(a)** Against the whole
locus and **(b)** against a non-coding background, in each case matching the two tools on the
proportion of their own background flagged. Axes give the percentage of alleles scoring higher, on
logarithmic scales; points below the diagonal are ranked higher by SpliceAI. Each panel is censored
at the resolution of its own background.
