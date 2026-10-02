# Reviewer package

Submitted with "Predicting splicing variants in ADPKD: a functional benchmark of AlphaGenome and
SpliceAI at PKD1 and PKD2" (Nephrology Dialysis Transplantation). `NDT_supplementary_methods.md`
in this folder is the manuscript's Supplementary Methods and is the best single description of how
the numbers were produced.

Everything needed to check every number in the manuscript without an AlphaGenome account and
without repeating any model query. All scores are stored; nothing here calls a paid service.

## Search flow, reconciled

| Stage | Count |
|---|---|
| Records returned by the two fixed PubMed queries (19 September 2026) | 391 |
| Papers included after title and abstract screening | 47 |
| Of those, retrieved as open-access full text | 26 |
| Further papers and supplementary files obtained directly | 3 + 3 |
| Papers that ultimately contributed a variant passing both verification checks | 25 |
| Papers contributing the 120 labelled splicing-class variants | 21 |

The gap between 47 and 25 is papers that were included on screening but contributed no variant that
passed both the quotation check and the description check, most often because the variant was named
but not experimentally tested. Every exclusion carries a reason in the screening log.

## Files

| File | What it is |
|---|---|
| `truthset_verified.tsv` | the curated benchmark: one row per source-level entry with the verbatim quotation, the outcome label and the VariantValidator coordinates. All 131 entries pass both checks. |
| `per_variant.tsv` | one row per unique variant with every score: SpliceAI, the AlphaGenome splicing composite and its components, AVI, CADD, locus percentiles. |
| `analysis_set_membership.tsv` | **the file to start from.** Every variant with its class, provenance, source, raw scores, and a boolean flag for each analysis population used in the manuscript, so any reported denominator can be rebuilt by filtering one column. |
| `excluded_indels.tsv` | the five variants absent from the precomputed files, all of which have a demonstrated effect. |
| `stratified.json` | the primary analysis: class counts, AUCs with intervals, paired differences (variant-level and study-cluster), leave-one-study-out, subgroups, the common four-tool subset, threshold confusion matrices, provenance by outcome, PKD1-only composites, and `resampling_details` (seeds, resample count, interval method, tie handling, single-class resamples, and the two variants with two source papers). |
| `ranking_manifest.json` | for each background and tool: alleles scored, positions sampled, threshold, achieved flag proportion, the eligible benchmark variants by identifier, and recovery counts. Also the quantitative reason the two backgrounds disagree. |
| `summary.json` | locus-wide background statistics and the matched-burden comparisons. |
| `coverage_matrix.csv` | which tool scored which variants, so an unscored variant is never read as a low score. |
| `noncoding_background_ranks.tsv` | per-variant ranks against the non-coding background: the source of Figure 3b. |
| `therapeutic_sites.json`, `therapeutic_sites_attribution.tsv` | raw AlphaGenome RNA-seq and splice outputs at the miR-17 and uORF sites, and the per-variant AVI feature attributions behind Figure 2. These are the model's own outputs, distinct from the composite attributions. |
| `ANALYSIS_PLAN.md`, `PLAN_FREEZE.txt` | the pre-scoring plan, its hash stamp, and nine dated amendments, each recording whether the results it could affect already existed. |
| `tables_ejhg.md`, `references.md`, `supp_table_S1_sources.md`, `nomenclature_statement.md` | manuscript tables, reference list, source table, nomenclature statement. |
| `figure*.png` | the figures. |

## Start here

    python3 verify_from_package.py

This rebuilds the manuscript's headline numbers from this folder alone, using only the stored
scores. It imports nothing from the analysis code, touches no network and makes no model query; it
needs only pandas and numpy, and exits non-zero on any mismatch. It also confirms that both ranking
comparisons flag the same PROPORTION of their own background, which is the point on which the two
comparisons turn. `implementation_details.md` holds the software versions, query settings, sampling
arithmetic and resampling specification.

## Reproducing the numbers from source

    python3 analysis/12_analyse.py             # locus background and burden statistics
    python3 analysis/18_stratified_analysis.py # the primary analysis
    python3 analysis/23_ranking_manifest.py    # ranking backgrounds, thresholds, recovery
    python3 analysis/21_build_tables.py        # regenerates the tables from the results
    python3 analysis/22_check_manuscript.py    # asserts the manuscript matches the results

`22_check_manuscript.py` is the check to run first: it exits non-zero if any number quoted in the
manuscript disagrees with the result files, and it also fails on a list of phrases withdrawn in
earlier revisions.

## Two points a reader should know before interpreting the results

**There are no patient-derived negatives.** Of 50 splicing-class positives, 26 came from patients; of
70 variants with no detected effect, 69 were engineered. Nothing here estimates clinical specificity.

**The two ranking comparisons answer different questions.** Both match the two tools on the
proportion of their own background flagged, and both are post-correction. They differ because the
candidate background differs: AVI ranks amino-acid-changing variants highly, so its top slice of the
whole locus is 97.9% coding although the locus is 72.6% non-coding.

## Licence

Analysis code: MIT. AlphaGenome predictions are redistributed under the AlphaGenome Output Terms of
Use (non-commercial; scientific publication permitted): see `LEGALLY_BINDING_TERMS_OF_USE.txt`.
Publisher full texts are not redistributed.
