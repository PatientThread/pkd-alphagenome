# Supplementary Methods

Companion to "Predicting splicing variants in ADPKD: a functional benchmark of AlphaGenome and
SpliceAI at *PKD1* and *PKD2*". Everything here is also generated or asserted by scripts in the
reviewer package, which reproduces every number in the main text without repeating any model query.

---

## S1. Search flow, reconciled

| Stage | Count |
|---|---|
| Records returned by two fixed PubMed queries, 19 September 2026 | 391 |
| Papers included after title and abstract screening | 47 |
| Of those, retrieved as open-access full text | 26 |
| Further papers and supplementary files obtained directly | 3 + 3 |
| Papers contributing a variant that passed both verification checks | 25 |
| Papers contributing the 120 labelled splicing-class variants | 21 |

The gap between 47 and 25 is papers included on screening that contributed no variant passing both
checks, most often because a variant was named but not experimentally tested. Every exclusion carries
a recorded reason.

## S2. Evidence classes and membership rules

**Splicing class.** Evidence is a splicing experiment: patient RNA, cDNA sequencing, RNA sequencing
or a minigene read out as splicing. Where a source reports both a splicing experiment and a
downstream consequence, the splicing experiment determines the class. rs3874648 is assigned here on
that rule: reverse-transcription PCR of leukocytes from a homozygous patient identified an atypical
*PKD1* splice form arising from a cryptic acceptor. No restriction is placed on provenance, and the
class contains 92 engineered variants.

**Reporter-output class.** Evidence is a 5' untranslated region luciferase reporter measuring protein
output. The class is named for its endpoint rather than a mechanism, because for three of its six
members the RNA-level mechanism is unresolved. Its membership rule deliberately differs from the
splicing class: it admits only naturally observed human alleles that were tested, excluding the
engineered substitutions both laboratories also constructed at these elements and alleles described
but not assayed. It is therefore not a random sample of tested 5' untranslated region alleles, and no
summary statistic is computed from it.

**Transcript abundance.** Measured only for the three variants of Chen et al., who reported reporter
RNA unchanged by RT-qPCR. The authors of the other three state that without RNA studies they cannot
exclude effects on pre-mRNA splicing or transcript stability, so no transcript-level claim is made
for those variants.

**Labelling.** By experimental outcome only, never by predictor output. "No effect detected in the
reported assay" is not a claim of clinical benignity. A variant given opposite outcomes by two
sources is excluded from both classes.

## S3. Scoring

| Component | Version |
|---|---|
| `alphagenome` Python client | 0.9.0 |
| SpliceAI | 1.3.1, stock implementation, run locally |
| TensorFlow / Keras (SpliceAI environment) | 2.21.0 / 3.15.1 |
| CADD | v1.7, `whole_genome_SNVs.tsv.gz` plus the precomputed-indel lookup |
| Reference genome | GRCh38 / hg38 |

AlphaGenome scorers `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS` and `RNA_SEQ` were taken
from `RECOMMENDED_VARIANT_SCORERS` with a 1,048,576 bp context centred on the variant. Gene-restricted
values were selected on the Ensembl gene identifier (*PKD1* ENSG00000008710, *PKD2* ENSG00000118762),
never by substring matching on a name. Kidney tracks were selected by exact `biosample_name` against
a curated list, never by substring, since "renal" alone would admit adrenal gland; eight RNA
sequencing tracks matched.

SpliceAI 1.3.1 required a one-line patch for NumPy 2 and an extension for equal-length
multi-nucleotide substitutions. The patched build reproduces the stock package exactly on 48 test
variants, and reproduces both published values for the complex allele discussed in the main text
(0.12 and 0.78).

AVI was trained on chromosomes 1, 4, 7, 8, 10, 13 and 15, with validation on 2, 5, 11, 14, 17, 20, 22
and X. *PKD1* lies on chromosome 16, which appears on neither list, so *PKD1* results are independent
of that split. *PKD2* lies on chromosome 4, a training chromosome; we could not check whether these
particular positions entered training, so overlap is not excluded rather than demonstrated.

## S4. Background sampling

The step is the locus span divided by a target of 12,000 and rounded **down** to whole bases, so the
realised sample slightly exceeds the target.

| Locus | Span (bp) | Step | Positions | Alleles | Atlas positions | Sampling fraction |
|---|---|---|---|---|---|---|
| *PKD1* | 52,189 | 4 | 13,048 | 39,144 | 52,191 | 25.00% |
| *PKD2* | 75,187 | 6 | 12,532 | 37,596 | 75,192 | 16.67% |

Every position contributes all three substitutions, so the unit is the alternate allele in both the
SpliceAI and the Atlas backgrounds, and matching the two tools on the proportion of their own
background flagged compares like with like. Resolution is one allele in 39,144 (0.0026%) for the
whole-locus background and one in 28,401 (0.0035%) for the non-coding background; Figure 3 censors
each panel at its own floor.

## S5. Statistics

Discrimination is the Mann-Whitney AUC, ties counted as one half. Intervals are percentile bootstrap,
2,000 resamples, 2.5th and 97.5th centiles. Seeds: 20260920 for the variant-level AUC, study-cluster
AUC and paired difference; 20260921 for the study-cluster paired difference.

The variant-level bootstrap resamples positives and negatives independently with replacement,
preserving observed class sizes. The cluster bootstrap resamples whole source studies with
replacement; a resample containing only one outcome class yields no AUC and is discarded, and 1,994
of 2,000 were retained for the paired study-cluster difference. The paired comparison evaluates both
tools on the same resampled indices, which is why the paired interval is narrow while the marginal
intervals are wide.

Three variants appear in more than one source. One is the excluded discordant variant, leaving two in
the analysed set; each is assigned to its first listed PMID for clustering, so every variant belongs
to exactly one of the 21 clusters and none is counted twice.

No equivalence margin is claimed, and overlapping marginal intervals are not treated as evidence of
equivalence. Coverage is reported per tool, and a variant a tool cannot score counts as unscored,
never as a low score.

## S6. Subgroup and restricted analyses

| Analysis | SpliceAI | AlphaGenome | n (pos/neg) |
|---|---|---|---|
| Primary, all splicing-class variants | 0.846 | 0.846 | 50/70 |
| Study-cluster bootstrap (95% interval) | 0.715 to 0.999 | 0.707 to 0.993 | 50/70 |
| Leave-one-study-out (range) | 0.805 to 0.920 | 0.797 to 0.916 | 21 studies |
| *PKD1* | 0.824 | 0.835 | 41/50 |
| *PKD2* | 0.928 | 0.917 | 9/20 |
| Exonic or untranslated region | 0.819 | 0.836 | 20/48 |
| Near-splice, 3 to 20 nucleotides | 0.791 | 0.789 | 19/21 |
| Engineered variants | 0.734 | 0.756 | 23/69 |
| Patient-derived variants | undefined, no controls | undefined, no controls | 26/0 |

An AUC is undefined where a subgroup contains no variants of one outcome class; this applies to the
canonical splice-site and patient-derived subgroups. Where both classes are present but controls
number one or two, an AUC is defined but too unstable to report and is withheld rather than shown.

Composite scores restricted to *PKD1*: AVI 0.722 (0.609 to 0.827) on 36 positives and 50 negatives,
CADD 0.557 (0.431 to 0.684). Restricted further to *PKD1* non-coding and synonymous variants, which
excludes variants directly changing the amino-acid sequence: AVI 0.793 (0.627 to 0.936) and CADD
0.731 (0.555 to 0.876).

## S7. Ranking manifest

| Background | Tool | Alleles scored | Threshold | Proportion flagged | Positives recovered | Negatives flagged |
|---|---|---|---|---|---|---|
| Whole locus | SpliceAI | 39,144 | 0.5 | 1.0423% | 20/25 | 3/18 |
| Whole locus | AVI | 156,573 | 1.550 | 1.0423% | 5/25 | 2/18 |
| Non-coding only | SpliceAI | 28,401 | 0.5 | 0.5598% | 17/21 | 3/16 |
| Non-coding only | AVI | 113,607 | 0.665 | 0.5598% | 18/21 | 3/16 |

Both comparisons match on proportion, not raw count. The negative denominators differ because the
whole-locus comparison admits synonymous variants and the non-coding comparison does not. The two
comparisons answer different questions: of the 1,632 alleles above the whole-locus AVI threshold only
2.1% are non-coding, although 72.6% of the locus is.

## S8. Analysis plan and amendments

An analysis plan was written and hash-stamped before scoring. Eleven dated amendments record, for
each change, whether the results it could affect already existed. Substantive amendments were: use of
one reference standard per assay class and the documented splicing composite; correction of a
matched-burden comparison that had matched raw counts rather than proportions; reassignment of
rs3874648 from an abundance class to the splicing class on the evidence above; withdrawal of a
transcript-abundance claim for three variants whose authors performed no RNA study; and restriction
of the locus-rank analyses to the splicing class. The plan, with every amendment, is in the reviewer
package.

## S9. Nomenclature

Gene symbols follow HGNC: *PKD1* (HGNC:9008) and *PKD2* (HGNC:9009). Variant descriptions follow
Human Genome Variation Society recommendations on the MANE Select transcripts NM_001009944.3
(NP_001009944.3) and NM_000297.4 (NP_000288.1), on GRCh38 (NC_000016.10, NC_000004.12). Every
description was validated with VariantValidator, which also supplied genomic coordinates; all 131
source-level entries passed. Legacy IVS notation was converted from the RefSeq exon table. Variants
published only as protein changes were converted by reading the codon from the reference coding
sequence and accepting a result only where exactly one single-base substitution produced the stated
amino-acid change; all 17 conversions were unique, and the method exactly reproduced the eight
variants whose papers also printed coding-sequence notation.
