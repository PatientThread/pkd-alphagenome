# Analysis plan: AlphaGenome at PKD1 and PKD2

**Status:** pre-specified. Written and committed on 19 September 2026, after the truth set was
verified and **before any variant in it was scored.** Changes after this point are recorded in
the amendments section at the end, with the date and the reason.

**Author:** Christopher Lawrence

---

## 1. Question

Does AlphaGenome, used through the AlphaGenome Atlas (AVI score and its feature attributions)
and through on-demand variant scoring, correctly identify experimentally validated non-coding
and splice-acting variants in *PKD1* and *PKD2*, and does it do better than SpliceAI and CADD at
these two loci?

*PKD1* is chosen because it is a hard locus for general pipelines (six pseudogenes over exons
1 to 33, high GC content) and because ADPKD is a dosage disease with a known functional
threshold, so a quantitative regulatory predictor has a clear biological job to do there.

## 2. Truth set

`data/curated/truthset_verified.tsv`, built by `analysis/04_build_truthset.py` and verified by
`analysis/05_verify.py`. At the time of writing: 54 entries, 52 unique variants.

| Label | Meaning | n (unique) |
|---|---|---|
| POS | functional assay showed aberrant splicing or reduced expression | 44 |
| NEG | functionally tested, indistinguishable from wild type | 7 |
| INDET | tested, assay could not be interpreted | 1 |

Every entry passed a verbatim quote check against its source and VariantValidator reference
confirmation on GRCh38. Untested variants are excluded, not treated as negative.

**The negative set is small (7).** That is the single largest limitation and it fixes what
statistics are appropriate (section 5).

## 3. Scores

1. **AVI** (Atlas, SNVs only): the raw AVI value and the 18 feature attributions.
2. **AlphaGenome on demand** (all variants, including indels and the dinucleotide allele):
   recommended scorers for splice sites, splice-site usage, splice junctions and RNA-seq, with
   1,048,576 bp context. Summarised two ways: across all tracks, and restricted to kidney
   tracks (the same curated kidney biosample list as the companion kidney-representation paper).
3. **SpliceAI** delta score (maximum of the four), run locally, distance 500, unmasked.
4. **CADD** v1.7 PHRED.

## 4. Pre-specified hypotheses

- **H1.** Splice-acting positives (branchpoint, AG-exclusion zone, near-splice-site, intronic
  deletions, pseudoexon) will rank in the top 1% of all possible SNVs in their locus by AVI,
  with MERGED_SPLICING as the dominant attribution.
- **H2 (expected failure).** c.-69dupG acts by lengthening a uORF and reducing *translation*.
  AlphaGenome predicts RNA, not protein, so it has no modality that can see this mechanism.
  **Prediction: AlphaGenome will not flag c.-69dupG.** Recorded in advance so that a miss is a
  finding about the model's scope, not a post hoc excuse.
- **H3 (haplotype).** For the intron 16 pseudoexon, the dinucleotide allele
  c.7066-131_7066-130delinsAG will score far above either constituent SNV scored alone. The
  common constituent (c.7066-130A>G, rs3874658) is a functional negative and should score low.
- **H4 (allele resolution).** At c.7866, C>A is positive and C>T is negative in the same assay.
  A model with allele-level resolution should separate them. At c.11257, C>A, C>G and C>T are
  all positive.
- **H5 (modest effects).** The exonic positives from the BMC Genomics minigene study are mostly
  small shifts in exon exclusion. Expect weaker, less consistent scores there than for the
  near-splice-site positives, for every tool.
- **H6 (common modifier).** rs3874648 (c.10051-239G>A) lowers full-length *PKD1* through a cryptic
  acceptor. Being a common variant with a modest quantitative effect, it is expected to score
  well below the pathogenic positives. This is recorded as a description, not a test.

## 5. Endpoints and statistics

**Primary endpoint.** For each POS SNV, its percentile rank by AVI among **all possible SNVs in
its gene locus** (gene span from the RefSeq transcript plus 5 kb upstream, all three alternate
alleles at every position), obtained from the Atlas by interval query. This is scale-free, so it
does not depend on how raw AVI maps to DeepMind's published percentile scale, and it mirrors the
diagnostic task of picking the causal variant from everything a patient could carry at that gene.

Reported as: the per-variant percentile, the median across POS variants, and the proportion in
the top 1% and top 0.1% of the locus.

**Secondary endpoints.**
- The same locus percentile for CADD, which is available genome-wide. SpliceAI is compared on
  its standard delta-score thresholds (0.2 and 0.5), because a locus-wide SpliceAI run may not
  be computationally feasible; if it is, the percentile is reported too.
- POS versus NEG, per variant, for every tool.
- Dominant AVI attribution versus the experimentally shown mechanism.
- Kidney-track versus all-track on-demand scores (link to the kidney-representation paper).

**Not done, deliberately.** No ROC AUC is reported as a headline number with 7 negatives: the
confidence interval would be too wide to mean anything and would invite over-reading. If
additional negatives are recovered from the manual-download papers before submission, an AUC
with a bootstrap interval becomes a secondary endpoint and this plan is amended.

**Thresholds.** SpliceAI 0.2 and 0.5 (conventional). CADD PHRED 20. AVI has no validated
threshold, so none is invented; the locus percentile is used instead.

## 6. Leakage and circularity checks

- **ClinVar.** If AVI or its components were trained or calibrated on ClinVar pathogenic
  variants, any truth-set variant already in ClinVar may be scored well for that reason alone.
  Every variant's ClinVar status and submission date is recorded, and results are reported
  stratified by in-ClinVar versus not in ClinVar.
- **AlphaMissense.** AVI includes AlphaMissense. Several positives are missense by annotation but
  act through splicing. Their AVI may be driven by protein damage, not by the mechanism that
  was shown. The attribution breakdown is reported for these specifically.
- **Conservation.** AVI includes conservation features. High conservation at a splice site can
  flag a variant without the model "understanding" the splicing mechanism. Reported alongside.

## 7. Sources of bias acknowledged in advance

- Published positives are enriched for variants that some earlier predictor found interesting
  enough to test. This favours all predictors, and favours older tools (MaxEntScan, SpliceAI)
  that were used to choose what to test.
- Minigene results in HEK293T are not kidney RNA.
- Engineered variants are molecular truth but not disease truth. Results are reported for
  patient-carried and engineered variants separately.

## 8. What would make this a negative paper, and that is acceptable

If AlphaGenome does not outperform SpliceAI on the splice-acting variants, and misses
c.-69dupG as predicted, the conclusion is that at *PKD1* the multimodal model adds nothing over a
dedicated splicing model for the variants that are currently testable. That is reportable and
useful to diagnostic laboratories.

---

## Amendments

### A1. Coding consequence stratum (19 September 2026)

**Timing, stated plainly.** Written after SpliceAI and CADD had been run on the truth set, and
**before any AVI score for a truth-set variant existed** (the Atlas pull was still in progress; the
single pre-freeze probe of c.11017-25A>G is disclosed in `PLAN_FREEZE.txt`). It was prompted by the
CADD output, not by AVI.

**Problem.** The truth-set labels describe a molecular effect on splicing or expression. They are
not pathogenicity labels. Some variants carry a protein consequence that is independent of that
molecular effect:

- Two NEG variants are nonsense (PKD2 c.741C>G p.Tyr247\*, c.796G>T p.Glu266\*). They are negative
  for splicing, but as variants they are truncating and almost certainly damaging. A composite
  pathogenicity score (CADD, AVI) that ranks them high is behaving correctly, so they are not valid
  negatives for those tools.
- Several POS variants from the BMC Genomics minigene study are nonsense or missense. AVI and CADD
  can rank them high through protein features (AVI's PROTEIN_TERMINATION and ALPHAMISSENSE
  attributions) without detecting the splicing effect at all.

**Change.** Every variant is assigned a coding-consequence class from its VariantValidator protein
prediction: non-coding (5'UTR, intronic), synonymous, missense, nonsense, or other coding.

1. **Splicing-specific outputs** (SpliceAI; AlphaGenome splice-site, splice-site-usage and
   splice-junction scores) are evaluated on the full truth set as planned, because for them the
   molecular label is exactly the right target.
2. **Composite scores** (AVI, CADD) are evaluated primarily on the **non-coding and synonymous**
   variants, where a high score cannot come from a protein consequence. Results on the full set
   are still reported, labelled secondary, and the AVI attribution for every nonsense and missense
   POS variant is shown so a reader can see which feature drove it.
3. For AVI, a positive counts as detected "for the right reason" only when MERGED_SPLICING (or an
   expression attribution, for the regulatory variants) is its largest attribution. Detection
   driven by PROTEIN_TERMINATION or ALPHAMISSENSE is reported separately.

The primary endpoint (locus percentile of AVI for POS SNVs) is otherwise unchanged; it is now
computed on the non-coding and synonymous POS SNVs, with the all-variant figure as secondary.

### A2. Training-chromosome leakage for AVI (19 September 2026)

**Timing.** Written after reading the AlphaGenome Atlas preprint's methods, and **before any AVI
score for a truth-set variant existed** (confirmed: `results/atlas_truthset.tsv` absent).

**Fact, from the preprint's Data Split section.** AVI was trained on gnomAD v4.1 variants from
chromosomes 1, 4, 7, 8, 10, 13 and 15, with labels from allele frequency (FAF95_GRPMAX below 0.001
= proxy impactful, at or above = proxy neutral). Validation used chromosomes 2, 5, 11, 14, 17, 20,
22 and X. DeepMind's own benchmark evaluation excluded any variant whose genomic position overlapped
a training or validation variant, allele regardless.

**Consequences.**
1. *PKD1* is on chromosome 16, which was in neither the training nor the validation set. The *PKD1*
   truth set is therefore a genuinely held-out test for AVI. **AVI results are reported primarily on
   *PKD1*.**
2. *PKD2* is on chromosome 4, a training chromosome. The same exclusion rule DeepMind applied is
   applied here: any *PKD2* truth-set variant whose position (or, for an indel, spanned interval)
   overlaps a variant in gnomAD v4.1 genomes is flagged as potentially seen in training. AVI results
   for *PKD2* are reported separately, split by that flag, and are secondary.
3. This affects AVI only. The AlphaGenome base model (used for on-demand scores) was trained on
   functional genomics tracks, not on variant labels, so the variant-level leakage argument does not
   apply to it in the same way. SpliceAI and CADD have their own, different training data and are
   not adjusted.

### A3. New negatives from hand-downloaded papers; discordant results (19 September 2026)

**Timing, stated plainly.** Written after all scores for the original 52 variants had been seen
(results note of 19 September). None of the variants added here had been scored when this was
written. The changes below follow rules fixed in advance or are mechanical; none was chosen by
looking at how a tool performed.

**New data.** Dr Lawrence obtained three papers that were not open access (PMIDs 24907393, 34257392,
26692149). From them:
- 8 new positives with full HGVS (PMID 34257392, Table 2).
- New negatives named in the main text: 10 *PKD2* missense (PMID 26692149, Table 1) and 8 *PKD1*
  (PMID 24907393: seven missense named in protein notation, plus c.11258G>A).
- The remaining negatives of PMIDs 24907393 (about 21) and 34257392 (29) are only in online
  supplementary files and are **not** added until those files are obtained.

**Protein-to-DNA derivation.** Variants given only in protein notation are converted by
`analysis/p_to_c.py`: the codon is read from the RefSeq coding sequence, the stated reference amino
acid is checked, and a DNA change is accepted only if exactly one single-base substitution produces
the stated amino-acid change. All 17 such conversions were unique. For *PKD2*, every derived position
matches the paper's own "distance from exon end" column. As a control, the same method reproduces
exactly the five variants for which these two papers also print c. notation (three in *PKD2*, two in
*PKD1*).

**Discordant results (rule).** Where two sources report opposite functional results for the same
variant, the variant is labelled DISCORDANT and excluded from both POS and NEG. This applies now to
*PKD1* c.11257C>T: positive in PMID 37468838 (minigene, 33.8% vs 27.7% exon exclusion, P < 0.05),
negative in PMID 24907393 (minigene, no effect). Because this changes one original label after its
scores were seen, every affected endpoint is reported twice: with the rule applied and with the
original label.

**AUC.** Section 5 anticipated this: with the added negatives, an AUC with a bootstrap 95% interval
becomes a secondary endpoint for each tool. It is computed on the variant classes each tool can
legitimately be judged on (amendment A1): all variants for the splicing-specific scores; for AVI and
CADD, non-coding and synonymous variants only, and separately on all variants, labelled.

**Caveat for the new negatives.** Almost all are missense. They are valid negatives for splicing
scores, but not for AVI or CADD, which will rightly score damaging missense variants high (A1).
Several are listed as pathogenic in the ADPKD mutation database by the paper's own account; one,
*PKD1* p.R3277C, is a well-known hypomorphic allele.

### A4. Exploratory aim: the miR-17 drug-target site in the PKD1 3'UTR (19 September 2026)

**Status and timing.** A new, exploratory, descriptive aim. Written after obtaining Supplementary
Table 4 of Lakhia et al., NAR 2026 (PMID 41558825), and **before any score or population data for
the site was looked at**. It is not part of the benchmark and cannot change any earlier endpoint.

**The site.** The mouse edit (allele-specific primers in Table 4) changes the six seed-matching bases
CACTTT to TGTCCC. The human orthologue was located by sequence: CACTTT occurs exactly once in the
1,019-nt *PKD1* 3'UTR (NM_001009944.3 c.\*153_\*158), with 31 of 35 flanking bases identical to the
mouse context. GRCh38 chr16:2,089,569-2,089,574 (forward strand AAAGTG), reference bases confirmed by
VariantValidator. Experimental truth, from the paper: disrupting this site (mouse, by base editing;
human cells, by steric-blocking oligonucleotide) stabilises *PKD1* mRNA and raises polycystin-1.

**Questions and expectations, fixed now.**
- **Q1, natural variation.** Which variants does gnomAD v4.1 record at the six bases, at what
  frequency and in which populations? A variant disrupting the seed would in principle mimic the
  drug. Recorded as description and hypothesis only; nothing here tests protection.
- **Q2, does AlphaGenome see the site?** AVI for all 18 possible SNVs at the six bases, ranked against
  every other SNV in the 3'UTR. **Expectation: no distinguishable signal**, because AlphaGenome has no
  microRNA modality and miR-17 repression is post-transcriptional and cell-type dependent.
- **Q3, the edit itself.** AlphaGenome on-demand RNA-seq prediction for *PKD1* for the full human
  6-base edit (AAAGTG to GGGACA on the forward strand). Experimental direction: *more* PKD1 mRNA.
  **Expectation: a predicted change near zero.** This is the only positive-direction test available.
- **Q4, uORF start-loss (added the same evening, before scoring).** Figure 4A of PMID 41006799,
  supplied by Dr Lawrence, gives the engineered edits c.-87A>T (uORF1 ATG to TTG) and c.-20A>T (uORF2
  ATG to TTG); both start codons confirmed in the RefSeq sequence, and the flanking sequence matches
  the figure. Measured effect (paper text): +9% translation (uORF1), +34% (uORF2), +130% (both).
  These act on translation, like c.-69dupG. **Expectation: AlphaGenome predicts near-zero change.**
  The double edit is two non-adjacent changes and is not scored as a single allele.

### A5. Gene-wide SpliceAI background by evenly spaced sample (19 September 2026)

**Timing.** Written before any gene-wide SpliceAI result was analysed (3 of 53 chunks of the
exhaustive *PKD1* run had been computed and none had been read).

**Reason.** Exhaustive scoring of every position takes about 22 minutes per 1,000 positions on this
machine, which is roughly 19 hours for *PKD1* and a further 28 for *PKD2*. The background is used only
to convert a score into a percentile and to count what fraction of a locus a threshold flags. An
evenly spaced sample of positions estimates both to within a few tenths of a percentage point.

**Change.** The SpliceAI locus background is computed on an evenly spaced sample of 12,000 positions
per locus (about one position in four for *PKD1*, one in six for *PKD2*), all three alternate alleles
at each. This affects the background distribution only: every truth-set variant is still scored
exactly, by `08_score_spliceai.py`, and AVI and CADD backgrounds remain exhaustive. Sample size and
spacing are recorded in the output. The exhaustive run can be completed later; the partial results
already computed are kept.

### A6. Non-coding background as a sensitivity analysis (20 September 2026)

**Timing, stated plainly.** Devised and run **after** seeing the first burden result, which compared
AVI and SpliceAI at an equal number of flagged variants against a background of every possible
variant in *PKD1*. It is therefore exploratory and is reported as a sensitivity analysis, not as a
pre-specified endpoint. Both results are reported together and neither is presented alone.

**Why.** Ranking a non-coding candidate against a background that includes every coding variant
penalises AVI unfairly: AVI is a whole-genome pathogenicity score and is correct to place
protein-truncating variants at the top of a gene. A laboratory searching for a non-coding cause,
having already assessed the coding sequence, would compare candidates against other non-coding
variants.

**Method.** The background is restricted to positions outside the exons of the canonical Ensembl
transcript, padded by 2 bp so that canonical splice-site positions are also excluded, using exon
coordinates fetched from Ensembl (`data/raw/exons_grch38.json`). The comparison is otherwise
identical: the AVI threshold is set so that AVI flags the same proportion of the restricted
background as SpliceAI flags at 0.5, and recovery of proven positives is compared at that matched
burden. Only variants classified as non-coding are counted.

### A7. Rebuild after external review (20 September 2026)

**Timing.** After the full analysis was complete and after an external editorial review of draft v1.3.
Everything below is therefore post hoc. The earlier analyses are retained in the repository and the
change is described, not hidden.

**Errors found by the review and corrected.**
1. *Burden matching used counts, not proportions.* The AVI threshold was set to flag the same NUMBER
   of variants as SpliceAI at 0.5, but the SpliceAI background is an evenly spaced sample of the
   locus (A5) while the AVI background is complete, so the same count was a four times stricter
   threshold. Corrected to match proportions. Whole locus: AVI recovers 5 of 25, not 1 of 25.
   Non-coding background: AVI 18 of 21, not 8 of 21. The corrected result reverses the direction of
   the earlier claim.
2. *The documented AlphaGenome splicing composite was not tested.* The published recommendation is
   max(splice_sites) + max(splice_site_usage) + max(splice_junctions)/5, each aggregated across all
   tracks and genes. Only the separate components, gene-restricted, had been used. The composite is
   now the primary AlphaGenome score, and the all-gene aggregation is recorded alongside the
   gene-restricted one.
3. *Equivalence was inferred from overlapping intervals.* Replaced by a paired difference in AUC with
   a bootstrap interval. No equivalence margin is claimed.
4. *Counts were inconsistent* between the abstract, Results and Supplementary Table S1, because
   entries and unique variants were mixed. Reconciled explicitly: 128 source-level entries describe
   125 unique variants.

**Change of design: one reference standard per assay class.** Variants are assigned to the class of
experiment that produced their evidence, and each class is analysed separately:
  splicing (patient RNA or minigene read out as splicing), abundance (transcript-level measurement),
  translation (5'UTR reporter, protein output with RNA unchanged).
The primary comparison uses the splicing class only. The abundance class contains one variant and
the translation class three; both are reported descriptively and no AUC is computed for them.

**Added robustness analyses**, because the controls are concentrated in a few papers: leave-one-study
-out AUC, a bootstrap that resamples whole source studies rather than variants, and subgroups by
gene, by distance from the splice site, and by patient RNA versus minigene evidence.

**Coverage.** A variant a tool cannot score is reported as unscored, never as a low score, in a
per-tool coverage matrix.

**Claims withdrawn.** The general statement that these models are blind to mechanisms acting after
transcription is withdrawn. rs3874648 acts through a cryptic splice acceptor, so a near-zero score
there is a splicing miss, not evidence about scope. Statements about translation-acting variants are
confined to the endpoint actually measured.

---

## Amendment A8 — 21 September 2026 (after second methodological critique)

Made before the revised numbers were read into the manuscript; the analysis code was changed
first and every downstream figure and table regenerated from one snapshot.

**A8.1 rs3874648 moves from the abundance class to the splicing class.** The classifier keyed on
the mechanism string, which begins "expression modifier", and so missed that Zhang et al.
(PMID 36029107) demonstrated the mechanism directly: "Reverse-transcription PCR (RT-PCR) of
peripheral blood leukocytes (PBL) from an ADPKD patient homozygous for rs3874648-A identified an
atypical PKD1 splice form... the rs3874648-A allele increased Tra2-beta binding affinity and
activated a cryptic acceptor splice-site." That is splicing evidence. `assay_class()` now tests
for splicing evidence before falling through to abundance. Primary set: 49 -> 50 positives.

**A8.2 The claim that transcript abundance was unchanged is withdrawn for the three Wedd variants.**
Wedd et al. (PMID 41006799) state: "without performing RNA studies, we cannot rule out additional
effects the variant may have on pre-mRNA splicing and transcript stability." Only Chen et al.
measured reporter RNA. The translation class is therefore described as a protein-output class in
which RNA was controlled for in one study of two.

**A8.3 Assay labels corrected.** Wedd used Gaussia luciferase normalised to secreted alkaline
phosphatase; Chen used Renilla/firefly dual luciferase. Table 4 now names each.

**A8.4 The miR-17 site is separated from the uORF mechanism.** Lakhia et al. report that base
substitution of the motif "stabilizes Pkd1 messenger RNA and increases polycystin-1 (PC1) protein
levels", i.e. an RNA-abundance mechanism, not protein output from an unchanged transcript.

**A8.5 Pre-specified additions**, all reported whatever they showed: a common analysis set of
variants scored by all four tools; a study-cluster paired difference; full confusion matrices at
the two conventional SpliceAI thresholds on the splicing class; and provenance cross-tabulated
against outcome.

---

## Amendment A9 — 21 September 2026 (after third methodological critique)

**A9.1 The ranking comparison is reported as two questions, not a correction.** The previous text
attributed the difference between the whole-locus and non-coding comparisons to the count-versus-
proportion bug fixed in A7. That was wrong: both comparisons have used proportions since A7, and the
difference between them is a change of candidate background, which changes the question. The reason
they differ is now reported and quantified: AVI's top 1.042% slice of the whole PKD1 locus contains
1,632 alleles, of which only 2.1% are non-coding although 72.6% of the locus is. `analysis/23_ranking
_manifest.py` writes the allele counts, thresholds, achieved proportions, eligible variant
identifiers and recovery counts for every background and tool.

**A9.2 The unit is stated: the alternate allele.** Both backgrounds carry all three substitutions at
each position (SpliceAI 39,144 scored alleles at 13,048 sampled positions; Atlas 156,573 alleles at
52,191 positions, a 25.0% sampling fraction). Proportion matching is therefore between like and like.

**A9.3 AVI and CADD are now reported for PKD1 alone as well as pooled.** The Methods said AVI was
reported primarily for PKD1, but every AVI figure was pooled across both genes, which the scored
counts exposed (45 positives against 41 PKD1 positives in total). PKD1-only AVI is 0.722 (0.609 to
0.827) on 36 positives and 50 negatives, against a pooled 0.718; the two agree.

**A9.4 Training exposure is stated at the strength the evidence supports.** Chromosome 4 appears on
the AVI training list, which makes overlap possible for PKD2; it does not demonstrate that these
positions entered training, and no position-level overlap check was available. The claim that PKD2
"was seen in training" is withdrawn.

**A9.5 The class formerly called "translation" is renamed "reporter output".** Naming it for a
mechanism asserted something unresolved for the three variants whose authors performed no RNA study.
Its membership rule, which admits only naturally observed human alleles, is now stated, because it
differs from the splicing class, which contains 92 engineered variants.

**A9.6 rs3874648 is described as an expression modifier only.** The source reports an effect on
PKD1 expression. It reports nothing this benchmark can use about clinical severity, and the stored
text available to us is the abstract, so no severity claim is made in either direction.

**A9.7 Indel coverage is described accurately.** CADD scores short indels on request and the Atlas
includes indel predictions; the five variants are absent from the precomputed files used, which is a
property of those files.

**A9.8 Subgroups distinguish two reasons for withholding an AUC**: no variants of one outcome class
(arithmetically undefined) versus both classes present but controls too sparse to report.

**A9.9 Resampling details recorded** in `results/stratified.json` under `resampling_details`: seeds,
2,000 resamples, percentile intervals, handling of single-class resamples, tie handling as one half,
and the two variants with two source papers assigned to their first listed PMID for clustering.

**A9.10 Reference numbering fixed.** `16_references.py` numbered all PMID citations in one pass and
all named citations in a second, so a Discussion PMID preceded an Introduction reference. It is now a
single document-order pass. Journals numbering articles rather than pages now get the article
identifier from the DOI.

---

## Amendment A10 — 21 September 2026 (fourth critique; reporting only, no result changed)

No analysis was rerun and no reported value changed. These are reporting corrections.

**A10.1 The sampling target and the realised sample were conflated.** Methods said 12,000 sampled
positions per locus and 39,144 scored *PKD1* alleles, which cannot both be true at three alleles per
position. The step is the locus span divided by 12,000 and rounded DOWN to whole bases, so the
realised sample exceeds the target: *PKD1* span 52,189 bp, step 4 bp, 13,048 positions, 39,144
alleles, 25.00% of the 52,191 positions the Atlas covers; *PKD2* step 6 bp, 12,532 positions, 37,596
alleles, 16.67%. The stored backgrounds already support every reported rank, so this is a Methods
correction, made by reporting the realised sample rather than by altering a number to fit.

**A10.2 Two Statistics sentences still contradicted the corrected Results.** "Because two tools
cannot score insertions and deletions" and "where a high score cannot arise from a protein
consequence" both survived amendment A9 because a line break fell inside the phrase and the
banned-phrase check was not whitespace-insensitive. The check now normalises whitespace, and also
reads the tables file, since a number stated only in a table is still reported.

**A10.3 Figure 3 is censored per panel.** One allele in 39,144 is 0.0026% for the whole-locus
background but one in 28,401 is 0.0035% for the non-coding background. Each panel now carries its
own floor and states it, rather than sharing one that overstates the resolution of the second.

**A10.4 The miR-17 query is labelled exploratory.** Lakhia et al. made the genomic edit in mouse
*Pkd1* and used oligonucleotide targeting in their human experiments, whereas we queried the model
on human sequence. The query is now fully specified in the text (GRCh38 chr16:2,089,569
AAAGTG>GGGACA, forward strand, 1,048,576 bp context) and described as an analogous construct rather
than a replication. The phrase "a fair test" is withdrawn.

**A10.5 Table 4 retitled** to the reporter-output class, matching Methods and Results, and stating
that reporter RNA was measured for three of the six.

**A10.6 rs3874648 and clinical severity.** A reviewer stated that Zhang et al. examined clinical
severity associations and found none significant. Only the abstract is held locally and the article
is in PMC but not open access, so this could not be verified. No claim is made in either direction;
the variant is described solely as an expression modifier. FOR THE AUTHOR TO CHECK with journal
access, since mentioning a null severity analysis would strengthen the paragraph if correct.

**A10.7 Reproducibility package made self-verifying.** `verify_from_package.py` rebuilds every
headline number from the package alone, importing none of the analysis code and making no network
or model call. `implementation_details.md` records software versions, AlphaGenome query settings,
the sampling arithmetic and the full resampling specification, including that the two analysed
variants with two source papers are clustered on their first listed PMID and that 1,994 of 2,000
cluster resamples were retained. The archive has been extracted to a clean directory and the
verification run there.

---

## Amendment A11 — 21 September 2026 (Zhang severity analysis resolved)

A10.6 recorded that the clinical-severity analysis in PMID 36029107 could not be verified because
the article is in PMC but not open access and only the abstract was held locally. The author
supplied the relevant passage from the published text on 21 September 2026. It is stored verbatim,
with provenance, in `data/papers/36029107.author_supplied.txt`, kept separate from
machine-retrieved text because it did not come through the retrieval pipeline.

The passage supports three statements, all now in the Discussion:

1. Zhang et al. did test clinical severity. Multivariate analysis against height-adjusted total
   kidney volume, GFR and hypertension found no association with rs3874648 genotype for total
   kidney volume, Mayo Clinic Classification, age at onset of end-stage kidney disease, or
   hypertension.
2. The cohort contained **eight** homozygotes in total. A null from a group that size does not
   exclude an effect, and the manuscript says so rather than reporting the null as if it settled
   the question.
3. Of those eight, three carried no pathogenic *PKD1* or *PKD2* mutation detectable by LR-PCR-NGS
   and MLPA, and homozygotes were over-represented among mutation-negative patients (p = 0.009).

Point 3 is the addition of substance. It links the one graded expression modifier in the benchmark
to the unexplained fraction of ADPKD that the Introduction opens with, and therefore gives a
concrete reason why a tool missing this class of variant matters. It is reported as the source
reports it, with the small numbers visible, because three subjects is a fragile basis for a
prevalence claim.

No analysis was rerun. rs3874648 remains in the splicing class on the evidence recorded in
amendment A8, and no result changed.
