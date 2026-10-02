# AlphaGenome at PKD1 and PKD2: first results

**Date:** 19 September 2026. **Status:** working research note, not a manuscript. *PKD1* results are
complete. *PKD2* results were added the same evening (see Addendum).

**Author:** Christopher Lawrence

---

## What was done

1. **A recorded literature search** (PubMed, two fixed queries, 391 papers), title-screened to 47,
   read in full text where open access (26) and as abstracts otherwise.
2. **A verified truth set** of 52 variants from 23 papers: 44 where a lab experiment showed abnormal
   splicing or reduced expression (positives), 7 where the same kind of experiment showed nothing
   (negatives), 1 inconclusive. Every entry carries a verbatim quote from its paper, checked by
   machine against the stored text, and every variant was confirmed against the reference genome by
   VariantValidator.
3. **An analysis plan written and hash-stamped before scoring**, with two amendments made and
   timestamped before the AlphaGenome scores they affect existed.
4. **Scoring with four tools**: AlphaGenome Atlas (the AVI score and its breakdown by feature),
   AlphaGenome on demand (splicing and RNA predictions, including kidney-only tracks), SpliceAI
   (the standard splicing predictor, run locally and validated against a published value), and CADD
   (the standard general-purpose score).

**Why *PKD1* is the fair test.** DeepMind trained AVI on chromosomes 1, 4, 7, 8, 10, 13 and 15 and
validated on 2, 5, 11, 14, 17, 20, 22 and X. *PKD1* is on chromosome 16, which was used for neither,
so AVI has never seen it. *PKD2* is on chromosome 4, a training chromosome, so its AVI results carry
a leakage caveat and are secondary.

---

## Main result

**For the 18 non-coding or synonymous positives in *PKD1*:**

| | AVI | CADD | SpliceAI |
|---|---|---|---|
| Median position among all 156,573 possible SNVs in *PKD1* | 94.8th percentile | 88.7th percentile | not yet computed |
| In the top 1% of the gene | 5 of 18 | 2 of 18 | |
| In the top 0.1% | 0 of 18 | 0 of 18 | |
| Flagged at the usual clinical cut-off | no validated cut-off | 8 of 18 (PHRED 20) | 16 of 18 (0.2 and 0.5) |

**In plain terms.** AVI points the right way, and for the right reason: 17 of these 18 variants have
splicing as the main driver of their score, which matches what the experiments showed. But it is not
sharp. A typical proven splice variant sits around the 95th percentile of the gene, which means
roughly 8,000 other possible changes in *PKD1* score higher. As a tool for picking the one causal
variant out of a patient's genome, that is not enough on its own at this gene.

**The comparison with SpliceAI is not yet fair, and should not be read as "SpliceAI wins".** SpliceAI
flags 16 of 18 at its standard threshold, but we do not yet know how many of the other 156,000
possible *PKD1* changes it would also flag. A threshold count and a percentile are different
yardsticks. Running SpliceAI across the whole gene would settle it (see next steps).

---

## The pre-specified questions

**H1. Do splice-acting variants reach the top 1% by AVI? Mostly no.** 3 of 14 did. 13 of 14 were
correctly attributed to splicing. Right mechanism, modest magnitude.

**H2. c.-69dupG (the 5'UTR variant that lowers translation by 87%). Predicted to be missed. It was.**
It is not in the Atlas at all (it is a rare private insertion, and the Atlas only holds insertions
already seen in population data). Scored on demand, AlphaGenome predicts almost no change in RNA
(largest effect 0.047 on a log scale) and no splicing change. SpliceAI 0.01. This is a limit of scope,
not a bug: AlphaGenome models RNA, and this variant acts on translation of an unchanged RNA. **Any
pipeline that relies on these tools will not find this class of *PKD1* variant.**

**H3. The intron 16 pseudoexon (two adjacent changes that only work together). Supported, with a
catch.** Given the two changes as one allele, AlphaGenome's splice-junction score is 7.4, against
0.17 for the de novo change alone and 0.93 for the common change alone; SpliceAI gives 0.78 against
0.12 and 0.24. Both tools see the combined allele clearly. **The catch is upstream of the model:** a
standard pipeline splits them into two separate SNVs, filters out the common one as benign, and
never asks the question. The failure in the published case was a pipeline failure, and a model
cannot rescue that. CADD cannot score the combined allele at all.

**H4. Can the tools tell alleles apart at the same position? Yes, but look at why.** At c.7866,
C>A (positive) scores 97.8th percentile by AVI and C>T (negative) 1.6th. But C>A is also a stop
codon, and AVI's score for it is driven by protein truncation, not splicing. SpliceAI separates the
pair on splicing grounds (0.74 against 0.01). At c.11257, all three positives score high by AVI, but
two of them through conservation rather than splicing; SpliceAI catches only C>A (0.53).

**H5. Weaker signal for the exonic minigene positives. Supported.** SpliceAI median 0.17 against
0.84 for near-splice-site variants; AlphaGenome splice-site usage 0.04 against 0.93. AVI and CADD
still score these highly, but mostly for protein reasons: of 15 coding positives AVI scored, only 4
were driven by splicing, 8 by protein damage and 3 by conservation.

**H6. rs3874648, the common variant shown to lower full-length *PKD1*. Missed by every tool.** AVI
10.9th percentile, CADD PHRED 0.3, SpliceAI 0.05. This matters for the modifier question you started
with: **a real, functionally demonstrated, common dosage modifier of *PKD1* is invisible to all
current predictors.** A modifier study cannot use these scores to find candidates of this kind.

---

## Negatives, and the one false alarm

The four clean *PKD1* negatives (non-coding or synonymous, experimentally normal) scored 1.6, 69.9,
76.3 and 84.9th percentile by AVI. So several true positives (for example c.7490-11C>G at 80.3, and
rs3874648 at 10.9) score *below* a true negative (c.-209G>A at 84.9). AVI's separation at *PKD1* is
partial.

SpliceAI flagged 1 of 7 negatives above 0.2: rs3874658, the common variant in the pseudoexon case,
at 0.24. AlphaGenome also gives it a modest splice-site score (0.23). It is functionally inert on its
own, proven by the unaffected father's RNA.

---

## Checks that came back clean

- **ClinVar leakage.** *PKD1* positives already in ClinVar score no differently by AVI from those not
  in ClinVar (median 98.0 against 97.8th percentile). AVI was trained on population frequency, not
  ClinVar, and this is consistent with that.
- **Kidney tracks.** Restricting AlphaGenome to its 11 kidney tracks makes splicing scores slightly
  lower, not higher (splice-site usage median 0.80 against 0.86 across all tracks). Nothing here
  suggests kidney-specific tracks rescue performance; it fits your Paper 1 finding that kidney is
  thinly represented.
- **SpliceAI installation.** Reproduced the published scores for c.7066-131G>A exactly (0.12 and
  0.06).

---

## Addendum: *PKD2* (Atlas download completed the same evening)

All 225,576 possible *PKD2* SNVs were retrieved, complete. The five non-coding or synonymous *PKD2*
positives all rank in the **top 1% of the gene** by AVI (median 99.8th percentile), against 5 of 18
for *PKD1*. All five are splicing-attributed.

**This does not mean AVI works better on *PKD2*, and must not be written up that way.** Three things
differ between the genes:

1. **The background, which moves every percentile.** *PKD1* packs 46 exons into 47 kb, so its locus is
   crowded with high-scoring coding changes: 4.9% of all possible *PKD1* SNVs have AVI above 1.0,
   against 1.3% in *PKD2*. The same raw score therefore ranks lower in *PKD1*.
2. **The variants.** The *PKD2* positives are canonical or near-canonical (-2, +5, -5). The *PKD1* set
   includes harder cases: branchpoint, AG-exclusion zone, -10, -11, -16.
3. **Training exposure.** *PKD2* is on an AVI training chromosome; one of the five
   (c.1898+5G>A) overlaps a gnomAD variant under DeepMind's own exclusion rule.

Like-for-like, the raw scores still lean the same way: donor +5 variants score 1.40 to 1.71 in *PKD1*
and 1.85 to 1.90 in *PKD2*. With five variants this is a description, not a finding.

**Lesson for the paper.** A within-gene percentile depends on what else is in the gene. It is the
right yardstick for "how hard is it to find the causal variant here", but it cannot be compared
across genes without saying so. Report raw AVI alongside it.

The two *PKD2* nonsense negatives rank at the very top of the gene (99.9th and 100th percentile),
driven by protein truncation, exactly as amendment A1 anticipated.

---

## Update 2: three hand-downloaded papers (evening, 19 September)

Dr Lawrence supplied the three priority papers. What they yielded (amendment A3):

- **26 new variants**, all quote-checked and VariantValidator-confirmed: 8 positives (Xie 2021) and
  18 negatives (10 *PKD2* and 8 *PKD1*). The truth set is now **51 positive, 25 negative, 1
  inconclusive, 1 discordant**.
- Seventeen negatives were published only as protein changes (for example p.R322W). Each was
  converted to its DNA change by reading the codon from the reference sequence, and accepted only
  when exactly one single-base change could produce it. All 17 were unique, and the method
  reproduces exactly the five variants where the authors also printed DNA notation.
- **Two laboratories disagree about one variant.** *PKD1* c.11257C>T changed splicing in one
  minigene study (2023) and did not in another (2014). It is now labelled DISCORDANT and left out;
  the results are shown both ways and it makes no difference. It is a useful point for the paper:
  the "truth" in these benchmarks depends on the assay.
- About 50 further negatives from these papers are only in their online supplements (see
  `docs/MANUAL_DOWNLOADS.md`).

**Discrimination, now that there are enough negatives (AUC with 95% bootstrap interval).**

| Score | All labelled variants (51 vs 25) | *PKD1* only (42 vs 12) |
|---|---|---|
| SpliceAI | 0.90 (0.83 to 0.97) | 0.88 (0.78 to 0.96) |
| AlphaGenome splice-site usage | 0.92 (0.85 to 0.97) | |
| AlphaGenome splice junctions | 0.91 (0.84 to 0.97) | 0.91 (0.83 to 0.97) |
| AlphaGenome splice sites | 0.90 (0.82 to 0.96) | 0.85 (0.74 to 0.94) |

**In plain terms: AlphaGenome's splicing outputs and SpliceAI perform the same at these two genes.**
The intervals overlap almost completely. Restricting AlphaGenome to kidney tracks changes nothing
that matters (0.92 against 0.92).

**AVI and CADD cannot be judged on the full set**, and their full-set figures (0.70 and 0.49) should
not be quoted as performance. The new negatives are damaging missense changes that are normal only
for splicing; a pathogenicity score is right to rate them highly. On the fair subset (non-coding and
synonymous variants) AVI reaches 0.92, but with only 4 negatives the interval (0.80 to 1.00) is too
wide to support any claim.

**What this changes in the story.** The earlier "SpliceAI flags more" impression is gone once both
are judged the same way. The overnight gene-wide SpliceAI run will add the last piece: how many of
the ~156,000 possible *PKD1* changes each tool would flag, which is the false-alarm burden a
laboratory would actually face.

---

## Update 3: the two drug-target sites (exploratory, amendment A4)

Two *PKD1* control switches are current drug targets, and in both cases disrupting them *raises*
polycystin-1 in the laboratory:

- **The miR-17 "brake" in the 3' end of the message** (target of farabursen and of blocking
  oligonucleotides). Supplementary Table 4 of the NAR paper gave the exact six-letter edit in mouse.
  The human equivalent was found by sequence: the six letters occur exactly once in the human 3' UTR,
  surrounded by near-identical sequence (31 of 35 letters), at chr16:2,089,569-2,089,574.
- **The two "uORF" decoy start signals in the 5' end** (target of uORF-blocking oligonucleotides).
  Figure 4A of the 5'UTR paper, supplied by Dr Lawrence, gave the exact edits (c.-87A>T and
  c.-20A>T); both were confirmed against the reference sequence. Measured effect: +9%, +34%, and
  +130% for both together.

**Does AlphaGenome see what these edits do? No, as predicted before scoring.**

| Edit | What the experiments showed | AlphaGenome's predicted change in *PKD1* RNA |
|---|---|---|
| miR-17 seed, six-letter edit | more *PKD1* message and protein | essentially none (at most 0.017, log scale) |
| uORF1 start loss (c.-87A>T) | +9% translation | essentially none (-0.029 average, wrong direction) |
| uORF2 start loss (c.-20A>T) | +34% translation | essentially none (+0.012 average) |

**AVI does rank these sites highly, but for the wrong reason.** Within their own UTRs, changes at the
seed rank around the 92nd percentile and changes at the uORF start codons around the 95th. For every
one of the 36 possible changes, the largest contributor is *evolutionary conservation*; AlphaGenome's
own functional predictions contribute almost nothing. AVI flags these bases because evolution has
kept them unchanged, not because it models what they do.

**Natural carriers exist.** gnomAD records six rare changes inside the miR-17 seed (1 to 9 carriers
each, all in exome data), and people carrying start-loss changes at both uORFs, including the exact
uORF2 change the laboratory used (c.-20A>T, 2 carriers) and c.-87A>G at uORF1 (12 carriers). In
principle these people carry a natural version of the drug. That is a hypothesis, not a finding:
nothing here tests it, and exome calls in this GC-rich stretch of *PKD1* need care.

**Why this matters for the paper.** These are the two best-characterised *PKD1* regulatory elements,
both are drug targets, and the leading sequence model cannot see either mechanism; its headline
score picks them up only through conservation. That sharpens the "blind spots" message: the model
covers splicing well and post-transcriptional control not at all.

---

## Update 4: the full truth set, rescored (late evening)

With both supplements in, the truth set is **121 unique variants: 51 positive, 68 negative**, plus 1
inconclusive and 1 discordant. Everything was rescored with every tool. The negative set has gone
from 7 this afternoon to 68, so the comparison below is the first one that carries weight.

### Discrimination (AUC, 95% bootstrap interval)

| Score | All labelled variants (51 vs 68) | *PKD1* only (42 vs 48) |
|---|---|---|
| AlphaGenome splice-site usage | 0.84 (0.77 to 0.92) | |
| SpliceAI | 0.84 (0.76 to 0.91) | 0.82 (0.71 to 0.90) |
| AlphaGenome splice junctions | 0.83 (0.75 to 0.91) | 0.83 (0.73 to 0.91) |
| AlphaGenome splice sites | 0.83 (0.74 to 0.90) | 0.80 (0.70 to 0.89) |
| AVI, non-coding and synonymous only | 0.77 (0.64 to 0.89), 30 vs 27 | |
| CADD, non-coding and synonymous only | 0.72 (0.58 to 0.84), 31 vs 27 | |

The earlier figures (around 0.90) came from a set with only 25 negatives, most of them easy. With the
harder, more realistic negatives from the two supplements, every tool drops by about 0.06, and
**AlphaGenome's splicing outputs and SpliceAI remain indistinguishable.**

### What this costs a laboratory in false alarms

Across the 119 labelled variants:

| SpliceAI cut-off | Positives caught | Negatives wrongly flagged |
|---|---|---|
| 0.2 (the usual screening cut-off) | 41 of 51 (80%) | **16 of 68 (24%)** |
| 0.5 | 35 of 51 (69%) | 7 of 68 (10%) |

Nearly a quarter of experimentally normal variants clear the usual threshold. This is consistent
with the original authors' own assessment (they reported SpliceAI specificity 0.55 on their 37
variants) and with individual cases in the new data: c.11713-3C>G scores 0.95 and c.9201+5G>A 0.96,
and both are normal in a minigene.

### The primary endpoint, unchanged in character

For the 25 non-coding or synonymous *PKD1* positives, AVI's median rank is the 94.2nd percentile of
the gene (5 of 25 in the top 1%), against CADD's 88.2nd (2 of 25). AVI still points the right way
without being sharp enough to pick the culprit on its own.

### An independent check on the comparator

The 2021 paper printed its own SpliceAI score for all 37 variants. My local SpliceAI reproduces **29
of 37 exactly to two decimal places and 31 of 37 within 0.02**; the largest difference is 0.12
(c.12003+5G>T, 0.01 published against 0.13 here). The residual differences are most likely the
SpliceAI version and gene annotation used, and they do not affect any conclusion. Recorded in
`results/summary.json` under `validation_spliceai_vs_published`.

### Figures

Three figures are drafted in `publication/figures/` (PNG and PDF, each with a CSV of its values):

1. **fig1_scores_by_label**: what each tool scores tested-positive and tested-negative variants.
2. **fig2_blind_spots**: the five variants that act after the message is made. Only three of the five
   have a percentage in their source, so the other two are shown as text rather than an invented bar.
3. **fig3_avi_attribution**: what drives AVI by variant class, showing conservation carrying the
   drug-target sites and protein damage carrying the coding positives.

A fourth figure (AVI against SpliceAI on the same within-gene scale) is waiting on the gene-wide run.

---

## Update 5: the gene-wide comparison (overnight, 20 September)

The *PKD1* background finished at 02:38. Both tools can now be judged on the same scale, and the
answer depends on which question is asked.

**Across the whole gene.** SpliceAI at 0.5 flags 1.04% of all possible *PKD1* changes. Setting AVI to
the threshold that flags the same number, SpliceAI recovers **20 of 25** proven non-coding and
synonymous positives and AVI recovers **1**. Median rank of a proven positive: 99.4th percentile by
SpliceAI, 94.2nd by AVI.

**That comparison is unfair to AVI, and the fair version matters.** AVI scores the whole genome for
pathogenicity, so it correctly fills the top of a gene with protein-truncating variants. A laboratory
hunting a non-coding cause has already dealt with the coding sequence. Restricting the background to
non-coding positions: SpliceAI at 0.5 flags 0.56% of non-coding positions and recovers **17 of 21**,
AVI recovers **8 of 21** at the same burden, and both flag 3 of 18 non-coding negatives. Median rank
is then 99.9th percentile for SpliceAI and **99.8th for AVI**, essentially identical.

So the twenty-fold gap was an artefact of the background, and the real difference is about twofold in
recovery at a matched alarm budget. I have recorded the restricted analysis as a post hoc sensitivity
analysis (amendment A6) because I devised it after seeing the first result, and both versions are
reported in the paper.

**What this means in practice.** Use a splicing-specific score for a splicing question. AVI is
informative once candidates have been narrowed to non-coding variants, and is not a replacement for
SpliceAI at these loci.

Figure 4 (`fig4_locus_rank`) shows this per variant: the share of the gene outranking each variant,
on a log scale, for both tools. Most positives sit below the diagonal, meaning SpliceAI places them
higher in the gene-wide list.

**PKD2 (finished 06:33)** shows the same pattern with much smaller numbers. Whole gene: SpliceAI at
0.5 flags 0.28% of possible variants and recovers 5 of 5, AVI 0 of 5 at matched burden. Non-coding
background: burden 0.12%, SpliceAI 4 of 4, AVI 3 of 4, median rank 99.99th percentile for both. With
four positives this is descriptive only, and PKD2 sits on an AVI training chromosome.

**All scoring is now complete for both genes.**

---

## Update 6: complete truth set, and two bugs found and fixed (20 September, morning)

Dr Lawrence supplied the full Table S3, closing the last literature gap. It added four negatives
(p.M3678T, p.L3682P, p.L3731Q, p.Q3751R). **The truth set is now complete at 125 unique variants: 51
positive, 72 negative**, plus one inconclusive and one discordant.

The table also validated the method. Three of its nine new rows are variants already in the set from
the paper's main text, and the codon-based conversion independently reproduced the authors' own DNA
notation for all three (p.R3719Q, p.R3753Q, p.R3753W). Two rows (p.R3753R, p.L3754L) could each arise
from four different single-base changes and were excluded rather than guessed.

**Two bugs, found by checking counts rather than trusting them.**

1. *Four variants were silently skipped.* The on-demand scorer resumed by entry number, and inserting
   rows mid-file shifted every later number, so the four new variants were treated as already done.
   It now resumes by variant identity. Every truth-set variant now has a score from every tool
   (verified: zero missing).
2. *The false-alarm counts were computed on the wrong subset.* A loop added on 19 September reused
   the variable names holding the positive and negative sets, so the threshold counts printed after
   it came from a seven-variant subset. Renamed. The corrected figures are almost identical in
   proportion (24% of negatives flagged at 0.2, as before), so no conclusion changes, but the
   numbers in the paper have been updated to match the data exactly.

**Final numbers (51 positive, 72 negative).**

| Score | AUC (95% CI) |
|---|---|
| AlphaGenome splice-site usage | 0.85 (0.77 to 0.92) |
| SpliceAI | 0.84 (0.76 to 0.91) |
| AlphaGenome splice junctions | 0.83 (0.75 to 0.91) |
| AlphaGenome splice sites | 0.82 (0.74 to 0.89) |
| AVI, non-coding and synonymous only | 0.77 (0.64 to 0.89) |

SpliceAI at 0.2: 41 of 51 positives and 17 of 72 negatives (24%). At 0.5: 35 of 51 and 7 of 72 (10%).

Everything downstream was regenerated: analysis, all four figures, the manuscript (v1.3) and its
Word copies.

---

## What this does NOT show yet

- **Whether AlphaGenome or SpliceAI is better at *PKD1*.** Needs SpliceAI run across all 156,573
  possible *PKD1* SNVs so both can be judged on the same percentile scale. Estimated 10 to 17 hours
  of local computing; it can run overnight.
- **Kidney tissue.** Most underlying experiments were in blood RNA or minigenes in non-kidney cells.
- **Anything about kidney tissue specifically.** Most positives were shown in blood RNA or in
  minigenes in HEK293T cells, not kidney.
- **Clinical validity.** These are research tools and AlphaGenome's terms forbid using the outputs
  for decisions about individual patients.

---

## What I think the paper is

A short report, not a benchmark race: **"Where sequence-to-function models see, and do not see,
proven *PKD1* variants."** The finding is a map of blind spots with clinical consequences: splice
variants are seen (weakly by AVI, strongly by SpliceAI); translation-level 5'UTR variants and common
dosage modifiers are not seen by anything; and the one haplotype case is a pipeline problem, not a
model problem. That framing survives the small negative set, and it survives being scooped on a
generic benchmark, because it rests on a curated, verified, mechanism-labelled truth set that nobody
else has built.

## Next steps, in order

1. Download the priority-1 papers in `docs/MANUAL_DOWNLOADS.md` (the negatives).
2. Run SpliceAI across the whole *PKD1* locus overnight, for a like-for-like comparison.
3. Report raw AVI alongside locus percentiles in all figures.
4. Figures, then a draft.

---

## Update 7 — 21 September 2026: second methodological critique, amendment A8

A second round of simulated editorial review raised three biological criticisms. All three were
checked against the stored primary sources before any code was changed, and all three were correct.

**1. rs3874648 was in the wrong class.** Zhang et al. (36029107) did not merely observe an
expression difference; they did RT-PCR on leukocytes from a homozygous ADPKD patient, found an
atypical *PKD1* splice form, and showed the variant allele activates a cryptic acceptor. That is
splicing evidence. My classifier keyed on the mechanism string, which begins "expression modifier",
and so routed it to the abundance class. Corrected: splicing evidence now takes precedence.

**Effect on the results.** The splicing class went from 49 to 50 positives. Every estimate was
recomputed from one snapshot:

| | before A8 | after A8 |
|---|---|---|
| AlphaGenome splicing composite | 0.861 (0.789-0.927) | 0.846 (0.770-0.916) |
| SpliceAI | 0.852 (0.766-0.926) | 0.846 (0.761-0.919) |
| paired difference | +0.009 (-0.021 to 0.042) | 0.000 (-0.033 to 0.037) |
| AVI | 0.733 | 0.718 |
| CADD | 0.577 | 0.565 |

The two splicing tools now perform identically, which strengthens rather than weakens the paper's
central claim. The abundance class is now empty and has been removed from the manuscript.

**2. The transcript-abundance claim was overstated.** The draft said all six 5'UTR reporter
variants had transcript level "shown to be unchanged". Only Chen et al. (39757590) measured
reporter RNA. Wedd et al. (41006799) state in terms: "without performing RNA studies, we cannot
rule out additional effects the variant may have on pre-mRNA splicing and transcript stability."
Withdrawn for those three variants.

**3. The assay labels were wrong.** Wedd used Gaussia luciferase normalised to secreted alkaline
phosphatase; Chen used Renilla/firefly dual luciferase. Both now named per row in Table 4.

**4. Two therapeutic mechanisms were conflated.** Lakhia et al. report that disrupting the miR-17
motif "stabilizes *Pkd1* messenger RNA and increases polycystin-1 protein levels". That is RNA
abundance, which is within AlphaGenome's scope; the uORF mechanism is translation of an unchanged
message, which is not. The Discussion now separates them.

### Why the AVI recovery numbers changed twice

The reviewers asked for this explicitly, and it is worth stating plainly because I reported the
earlier version to Dr Lawrence twice before catching it.

The matched-burden comparison sets an AVI threshold that flags "as many" locus variants as
SpliceAI at 0.5. The AVI background is every position in the locus; the SpliceAI background is an
evenly spaced 1-in-4 sample. Matching raw *counts* therefore made the AVI threshold roughly four
times too strict. Matching *proportions* is correct. Whole-locus recovery went from 1 of 25 to 5
of 25, and non-coding recovery from 8 of 21 to 18 of 21 — which reverses the conclusion from "AVI
performs poorly at ranking" to "the two perform alike on the fair comparison". Both numbers are now
generated from `results/summary.json` and asserted by `analysis/22_check_manuscript.py`.

### Other changes in this round

- Locus-rank and matched-burden analyses restricted to the splicing class (they were splicing-only
  by accident until the Chen variants were added on 20 September).
- Provenance reported against outcome: 26 patient-derived positives, 0 patient-derived negatives.
  All clinical-specificity language withdrawn; the flag rate is described as a rate on engineered
  constructs.
- Full confusion matrices with Wilson intervals at SpliceAI 0.2 and 0.5, on the splicing class.
- Common four-tool subset (115 variants); the 5 variants lost are all indels and all positives, so
  the missingness is not at random and is stated as such.
- Study-cluster paired difference added alongside the variant-level one.
- AVI training-chromosome qualification restored (*PKD1* held out, *PKD2* in training).
- ClinGen use corrected: predictions can contribute PP3/BP4, so "rather than as evidence in itself"
  was wrong and is gone.
- Reference 2 updated from the bioRxiv record to the published Nucleic Acids Res paper (41558825).
- "Following external review" changed to name what it was: an AI-assisted critique during
  methodological revision, not journal peer review.
- Results subsections reordered so figures are cited 1, 2, 3.
- Tables are now generated by `analysis/21_build_tables.py` instead of being maintained by hand,
  which is how they drifted.
- Figure 3 gained panel (b), the non-coding-background comparison.
- Reviewer package built at `publication/PKD_AlphaGenome_reviewer_package.zip`.

### Still outstanding, needs Dr Lawrence

- The repository URL and archived release identifier for the Data Availability statement.
- Confirmation that the competing-interests paragraph is complete.

---

## Update 8 — 21 September 2026: third critique, amendment A9

The previous round's corrections were accepted. This round found one error of exposition that
changed what a principal result meant, one internal inconsistency, and a list of precision fixes.

### The exposition error: I blamed a bug for a change of question

The v2.1 text said the whole-locus ranking comparison was "unfair to AVI" and that the direction
reversed "when it is made fair", attributing the difference to the count-versus-proportion bug fixed
in A7. That was wrong, and I checked the code before accepting it: `12_analyse.py` has used
proportions in **both** comparisons since A7. So the two results differ for a different reason
entirely, and the real reason is more interesting than the bug.

| Background | SpliceAI | AVI | negatives flagged |
|---|---|---|---|
| whole locus, 1.042% matched | 20/25 | 5/25 | 3/18 and 2/18 |
| non-coding only, 0.560% matched | 17/21 | 18/21 | 3/16 both |

**Why they differ:** AVI is a whole-genome pathogenicity score and ranks amino-acid-changing variants
highly. Its top 1.042% slice of *PKD1* is 1,632 alleles, of which only **2.1% are non-coding**,
although **72.6%** of the locus is. Ranking a non-coding candidate against the whole gene makes it
compete with coding variants it would never be confused with. That is now the reported explanation,
quantified, with the bug relegated to the amendment log where it belongs.

Also settled: the unit is the alternate allele in both backgrounds (SpliceAI 39,144 alleles at 13,048
sampled positions; Atlas 156,573 at 52,191, a 25.0% sampling fraction), so proportion matching
compares like with like. `analysis/23_ranking_manifest.py` writes all of it.

### The inconsistency: AVI was pooled while the Methods claimed PKD1

Caught by arithmetic: the AVI row showed 45 scored positives, but *PKD1* has only 41 positives in
total, so it could not have been a *PKD1*-only analysis. Now reported both ways. *PKD1*-only AVI is
0.722 (0.609 to 0.827) on 36 positives and 50 negatives against a pooled 0.718. They agree, so the
conclusion is unchanged, but the reader can now see which analysis is which.

### Claims checked before acceptance, and what the sources actually said

- **Conservation dominance at all 36 substitutions** — SUPPORTED, and now reported at variant level:
  conservation is the largest of the four groups for every one of the 36, exceeding the next group by
  at least 0.307, while no other group reaches 0.076 in any.
- **miR-17 predicted RNA change is small** — SUPPORTED, and now given as the model's raw output
  rather than an attribution: largest absolute predicted change 0.017 log2 across 371 RNA-seq tracks,
  about 1%, mean 0.003 overall and across the eight kidney tracks.
- **Complex allele behaviour** — SUPPORTED and strengthened. Han reports SpliceAI 0.12 for the de
  novo change alone and 0.78 for the dinucleotide allele; our local SpliceAI reproduces both exactly.
  AlphaGenome composite 3.21 versus 0.12, a twenty-six-fold difference.
- **Zhang and clinical severity** — NOT VERIFIABLE. Only the abstract is stored, so the reviewer's
  statement that Zhang found no severity association could not be confirmed either way. Resolved the
  safe way: rs3874648 is described as an expression modifier and every severity implication is gone.
- **CADD and the Atlas can score indels** — CORRECT. I used `whole_genome_SNVs.tsv.gz` plus the
  precomputed-indel lookup, so the five variants are absent from the files used, not unscoreable.
  Wording fixed.

### Other changes

- Reporter class renamed from "translation" to "reporter output"; its membership rule (naturally
  observed human alleles only) is stated, against 92 engineered variants in the splicing class.
- Training exposure softened: chromosome 4 eligibility is not demonstrated exposure for *PKD2*.
- Subgroups now distinguish "undefined, no controls" from "withheld, sparse controls".
- Resampling details recorded: seeds, 2,000 resamples, percentile intervals, ties as one half,
  single-class resamples discarded, and the 2 variants with two source papers (of 3 total, one being
  the excluded discordant variant) assigned to their first listed PMID.
- Figure 3 floor moved from 0.001% to 0.0026%, the actual resolution of one allele in 39,144.
- Search flow reconciled in the package: 391 records, 47 included, 26 open-access, 25 contributing,
  21 behind the splicing class.
- Nomenclature statement written as a separate file; HGNC IDs, MANE transcripts and protein
  accessions all verified live against genenames.org and NCBI.
- **Reference numbering bug found and fixed**: `16_references.py` numbered all `[PMID ...]` tags in
  one pass and all `[REF:...]` tags in a second, so a Discussion PMID was numbered before an
  Introduction reference. Now a single document-order pass. Citation order is 1-14 and figure order
  1-3, both asserted mechanically by `22_check_manuscript.py`, which now runs 29 numeric checks plus
  order checks and a banned-phrase list.
- Tables merged so Table 2 carries the restricted analyses as rows: four tables, within the limit.
- Abstract 244 words, body 3,997, four tables, three figures, 14 references.

### Still outstanding, needs Dr Lawrence

- Repository URL and archived release identifier for the Data Availability statement.
- Confirmation that the competing-interests paragraph is complete.

---

## Update 9 — 21 September 2026: fourth critique, amendment A10

Editorial recommendation moved from major revision to short final revision before submission. No
result changed this round. Two things were genuinely wrong in the reporting, and one claim I had to
decline to make.

### The arithmetic that did not add up

Methods said 12,000 sampled positions per locus and 39,144 scored *PKD1* alleles. At three alleles
per position that is 36,000, not 39,144. I checked the code rather than adjusting a number: the step
is `span // 12000`, floor division, so the step is 4 bp at *PKD1* and the realised sample is 13,048
positions, hence 39,144 alleles, which is exactly 25.00% of the Atlas positions. 12,000 was a target,
not the realised count. The stored backgrounds support every rank already reported, so this was a
Methods correction. *PKD2*: step 6 bp, 12,532 positions, 37,596 alleles, 16.67%.

### A gap in my own checker

Two sentences withdrawn in the previous round were still in Statistics: "two tools cannot score
insertions and deletions" and "cannot arise from a protein consequence". My banned-phrase check
missed both because a line break fell inside the phrase. The check now normalises whitespace before
matching, and reads `tables_ejhg.md` as well, since a number stated only in a table is still
reported. That second change also fixed two checks that failed once I moved figures out of prose
and into Table 2.

### The claim I did not make

A reviewer said Zhang et al. examined clinical severity and found nothing significant, and suggested
saying so. I only hold the abstract; the paper is in PMC but not open access, so I could not verify
it. I have left the manuscript making no severity claim in either direction. **Worth Dr Lawrence
checking with journal access**, because if correct it would strengthen that paragraph.

### Other corrections

- **miR-17 query relabelled exploratory.** Lakhia edited MOUSE *Pkd1*; their human work used
  oligonucleotide targeting. I queried human sequence. The query is now fully specified (GRCh38
  chr16:2,089,569 AAAGTG>GGGACA, forward strand, 1,048,576 bp context) and called an analogous
  construct, not a replication. "A fair test" withdrawn.
- **Figure 3 censored per panel**: 0.0026% for the whole-locus background, 0.0035% for the
  non-coding one. A shared floor overstated the resolution of the second.
- Table 4 retitled to the reporter-output class.
- Abstract reworded so "0.72 (AVI) and 0.57 (CADD)" cannot be read as interval bounds or as results
  of the clustering analysis.
- Word counts: abstract 247, body 3,993, both with margin. Four tables, three figures, 14 references.

### The package now proves itself

`verify_from_package.py` rebuilds every headline number from the package alone: it imports none of
the analysis code, makes no network or model call, needs only pandas and numpy, and exits non-zero
on any mismatch. I extracted the archive to a clean directory and ran it there: all 14 checks pass,
and it confirms both ranking comparisons flag the same proportion (1.0423% and 0.5598%) of their own
backgrounds, which is the point the two comparisons turn on. `implementation_details.md` adds
software versions, AlphaGenome query settings, the sampling arithmetic and the resampling
specification.

### Still outstanding, needs Dr Lawrence

- Repository URL and archived release identifier for the Data Availability statement.
- Confirmation that the competing-interests paragraph is complete.
- Optional: confirm the Zhang severity analysis if he has journal access.

---

## Update 10 — 21 September 2026: the Zhang severity question, answered

Dr Lawrence supplied the passage from the published text. Stored verbatim with provenance in
`data/papers/36029107.author_supplied.txt`, separate from machine-retrieved text.

The reviewer was right that Zhang tested severity and found nothing. But the passage carries two
things the reviewer did not mention, and Dr Lawrence spotted both immediately: **there were only
eight homozygotes in the entire cohort, and three of those eight had no pathogenic PKD1 or PKD2
mutation at all.**

That changes how the finding should be used. Reporting "they tested severity and found nothing"
as a clean null would have been misleading from a group of eight. The Discussion now says they
tested it, names the endpoints, and states the cohort size so the reader can weigh the null
themselves.

**The genuinely useful part is the third sentence of that passage**, which the reviewer passed over:
homozygotes were significantly over-represented among patients with no detectable *PKD1* or *PKD2*
mutation (p = 0.009), three of the eight. That connects the one graded expression modifier in the
benchmark directly to the unexplained 10% the Introduction opens with. It gives the closing
paragraph a concrete reason why a tool that misses this class of variant matters, instead of a vague
appeal to phenotypic variability. It is reported with the small numbers on show, because three
subjects will not carry a strong prevalence claim.

Net effect on the paper: the closing paragraph went from a hedge to an argument. No analysis rerun,
no result changed. Body 3,998 words after trimming elsewhere to make room; abstract 247.

---

## Update 11 — 21 September 2026: retargeted to NDT

EJHG declined the companion paper (Paper 1, kidney representation in AlphaGenome, manuscript
1401-26-EJHG) on 21 September, not on quality but on fit: "better placed in a renal journal". This
paper is more kidney-specific still, so it was retargeted to Nephrology Dialysis Transplantation
before submission rather than after collecting the same decision.

**No result changed.** Both manuscript versions are driven by the same result files and both pass
all 29 numeric checks plus the citation and figure order checks. The truth set still verifies
131/131.

### What NDT requires that EJHG did not

Requirements read live from the NDT Author Guidelines on 21 September 2026. The binding constraint
is the word limit, which is 3,500 words **including** the abstract, where EJHG allowed 4,000
excluding it. That is roughly 750 words tighter.

- Structured abstract, 300 words, with the four headings Background and hypothesis / Methods /
  Results / Conclusions.
- Mandatory Key Learning Points: What was known, This study adds, Potential impact, at most three
  bullets and 50 words each.
- At most 5 keywords (EJHG version had 6), running head at most 50 characters.
- Superscript Vancouver citations, section order ending References, Tables, Figure Legends, Figures.

### What I changed, beyond reformatting

The paper is reframed for nephrologists rather than human geneticists. The Introduction now opens on
why genetic diagnosis matters in ADPKD clinically (prognosis, trial eligibility, living-donor
assessment) and on the laboratory's real question, which non-coding variant justifies an RNA study,
instead of on sequence-to-function models as a class of tool.

Methodological specification that no longer fits went to **Supplementary Methods S1-S9**: search-flow
reconciliation, class membership rules, versions and query settings, sampling arithmetic, the
resampling specification, subgroup and restricted analyses, the ranking manifest, the amendment
summary and nomenclature. Nothing was dropped. Every item a previous reviewer demanded is either in
the main text or in S1-S9, and `22_check_manuscript.py` now reads the supplementary file so those
numbers are still asserted.

Current lengths: abstract 298/300, abstract plus main text 2,745/3,500, so there is real headroom if
reviewers ask for additions. I used some of it to keep the common four-tool comparison, the
membership rule and the Wedd RNA caveat in the main text, since reviewers had specifically asked for
each.

### Two things worth recording

**Citations are now mechanical in this version too.** I first drafted the NDT version with typed
superscripts, then converted all 18 to `[PMID ...]` / `[REF:...]` tags so `16_references.py` numbers
them in document order. Hand numbering is exactly what produced the reference-order error fixed in
amendment A10, and there was no reason to reintroduce it.

**The checker caught two banned phrases I had reinstated by accident.** Writing fresh prose for the
new framing, I wrote "performed identically" and "no accuracy reason to change", both withdrawn in
earlier rounds for good reason. The banned-phrase list caught both immediately. That list is now
earning its keep across manuscript versions, not just across drafts of one.

### Files

`NDT_submission_manuscript.docx`, `NDT_review_copy.docx`, `NDT_supplementary_methods.docx`,
`NDT_cover_letter.docx` (carrying the AI disclosure NDT requires there as well as in Methods), and
the reviewer package. `analysis/24_build_ndt_docx.py` refuses to build if any stated NDT limit is
breached, so the limits are enforced rather than eyeballed. The EJHG version is retained unchanged.

### Still outstanding, needs Dr Lawrence

- Repository URL and archived release identifier for the Data Availability statement.
- Confirmation that the competing-interests paragraph is complete.
- **Paper 1 also needs a new home** after the EJHG decision.
- Optional: a graphical abstract, which NDT encourages but does not require.
