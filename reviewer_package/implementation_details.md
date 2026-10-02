# Implementation details

Everything a reader needs to know about how the scores were produced. None of it belongs in the
main text, and all of it is needed to repeat the work exactly.

## Software versions

| Component | Version |
|---|---|
| `alphagenome` Python client | 0.9.0 |
| SpliceAI | 1.3.1 (stock implementation, run locally) |
| TensorFlow / Keras (SpliceAI environment) | 2.21.0 / 3.15.1 |
| NumPy | 2.5.3 |
| pandas | 3.0.6 |
| CADD | v1.7, `whole_genome_SNVs.tsv.gz` plus the precomputed-indel lookup |
| Reference genome | GRCh38 / hg38 |

SpliceAI 1.3.1 required a one-line patch for NumPy 2 (`one_hot` indexing) and an extension for
equal-length multi-nucleotide substitutions. The patched build reproduces the stock package exactly
on 48 test variants, and reproduces the two published values in PMID 42502691 (0.12 and 0.78).

## AlphaGenome query settings

- Scorers: `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`, `RNA_SEQ`, taken from
  `RECOMMENDED_VARIANT_SCORERS`.
- Sequence context: 1,048,576 bp, centred on the variant.
- Splicing composite: `max(splice_sites) + max(splice_site_usage) + max(splice_junctions)/5`, each
  the maximum absolute score across all tracks and genes, as documented by the developers.
- Gene-restricted values are also stored, selected on the Ensembl gene identifier
  (PKD1 ENSG00000008710, PKD2 ENSG00000118762), never by substring matching on a name.
- Kidney track selection is by exact `biosample_name` against a curated list, never by substring:
  "renal" alone would pull in adrenal gland. Eight RNA-seq tracks matched.
- On-demand scoring run finished 2026-09-20T16:11:54Z; 130 variants scored, 0 failures.
- AVI scores, the 18 AVI feature attributions and the locus backgrounds are the Atlas precomputed
  single-nucleotide release, not recomputed here.

## Background sampling

The step is the locus span divided by a target of 12,000, rounded **down** to whole bases, so the
realised sample exceeds the target:

| Locus | Span (bp) | Step | Positions | Alleles | Atlas positions | Sampling fraction |
|---|---|---|---|---|---|---|
| *PKD1* | 52,189 | 4 | 13,048 | 39,144 | 52,191 | 25.00% |
| *PKD2* | 75,187 | 6 | 12,532 | 37,596 | 75,192 | 16.67% |

Every position contributes all three substitutions, so the unit is the alternate allele in both the
SpliceAI and the Atlas backgrounds and proportion matching compares like with like. The resolution of
the *PKD1* whole-locus background is one allele in 39,144 (0.0026%) and of the non-coding background
one in 28,401 (0.0035%); Figure 3 censors each panel at its own floor.

## Resampling

- 2,000 resamples throughout; percentile intervals at the 2.5th and 97.5th centiles.
- Seeds: 20260920 for the variant-level AUC, study-cluster AUC and paired difference; 20260921 for
  the study-cluster paired difference.
- Variant-level bootstrap resamples positives and negatives independently with replacement,
  preserving the observed class sizes.
- Cluster bootstrap resamples whole source studies with replacement. A resample containing only one
  outcome class yields no AUC and is discarded; the number retained is reported beside the interval
  (1,994 of 2,000 for the paired study-cluster difference).
- The paired comparison evaluates both tools on the **same** resampled indices, which is why the
  paired interval is narrow while the marginal intervals are wide.
- Ties count as one half in the AUC, so identical scores neither help nor harm a tool.
- Three variants appear in more than one source. One is the excluded discordant variant
  (*PKD1* c.11257C>T), leaving two in the analysed set (*PKD1* c.7489+5G>A, *PKD2* c.1716G>A). Each
  is assigned to its **first listed PMID** for clustering, so every variant belongs to exactly one
  cluster and none is counted twice. The 21 study clusters are counted on that basis.

## The exploratory miR-17 query

The published genomic edit was made in **mouse** *Pkd1*; that study's human experiments used
oligonucleotide targeting rather than an edit. Our query is therefore an analogous construct on
human sequence, not a replication: GRCh38 chr16:2,089,569, reference `AAAGTG`, alternate `GGGACA`,
forward strand, in the *PKD1* 3' untranslated region, 1,048,576 bp context. Results are in
`therapeutic_sites.json` under `Q3_seed_6nt_edit`, with the uORF start-loss edits under `Q4_*`.
