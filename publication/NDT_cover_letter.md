# Cover letter

The Editors
*Nephrology Dialysis Transplantation*

21 September 2026

Dear Editors,

I am pleased to submit "Predicting splicing variants in ADPKD: a functional benchmark of AlphaGenome
and SpliceAI at *PKD1* and *PKD2*" for consideration as an Original Article.

About one patient in ten with an ADPKD phenotype has no causal variant identified after
coding-region testing, and laboratories increasingly lean on computational predictors to decide
which non-coding variant of uncertain significance justifies an RNA study. SpliceAI is already in
routine use for this. AlphaGenome, with precomputed predictions for every possible single-nucleotide
variant in the genome, has just become available to diagnostic laboratories. Nobody has asked how
either performs at the polycystic kidney disease loci specifically, and general benchmarks cannot
answer that, because performance depends on local exon density, conservation and which mechanisms
are common in a given gene.

I assembled a benchmark of 128 *PKD1* and *PKD2* variants that had actually been tested in a
laboratory, and compared the tools on it. The findings a nephrologist can act on are these.
AlphaGenome offers no accuracy gain over SpliceAI at these loci: the paired difference in
discrimination was exactly zero. Both tools miss about one demonstrated effect in five at the usual
threshold. Both are blind to variants acting on protein output rather than splicing, which includes
the two *PKD1* regulatory elements currently under drug development. And apparent superiority in
ranking candidates depended on the comparison set used rather than on the tool, which matters for
any laboratory prioritising a candidate list.

I have tried to make the limitations as prominent as the results. Almost every variant in which no
effect was detected was engineered rather than patient-derived, because laboratories rarely publish
a negative result on a patient's variant. The paper therefore states plainly that it does not and
cannot estimate clinical specificity.

Two features may be of particular interest to your reviewers. Every variant assignment carries a
verbatim quotation from its source that is checked mechanically against the stored text, and every
variant description was validated independently against the reference genome; all 131 source-level
entries pass both checks. Second, a reviewer package accompanies this submission containing the
labelled data set, per-variant analysis-set membership flags, all stored scores and a script that
regenerates every number in the manuscript from those scores alone, with no network access and no
model query. Running `verify_from_package.py` reproduces the headline results in a few seconds.

**Use of artificial intelligence.** As disclosed in Materials and Methods, I used an AI coding
assistant (Claude, Anthropic) to write the analysis code, to assemble the curated data set from
published sources, to critique successive drafts and to draft the manuscript. The mechanical
quotation and coordinate checks described above exist precisely to constrain that use. I made and
verified all final analytical and interpretive decisions and I take full responsibility for the
work. The tool does not meet authorship criteria and is not an author.

This work has not been published previously in whole or in part and is not under consideration
elsewhere. There are no co-authors. My competing interests are declared in full in the manuscript: I
hold directorships and shareholdings in companies developing clinical software and genomic
diagnostics and am named on patent applications relating to genomic diagnostics, none of which
relates to *PKD1* or *PKD2* variant interpretation, to polycystic kidney disease, or to any tool
evaluated here. No funding was received for this study.

I would be grateful for your consideration.

Yours faithfully,

Dr Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP
Consultant Nephrologist
ORCID 0000-0002-8159-0879
Christopher.lawrence3@nhs.net
