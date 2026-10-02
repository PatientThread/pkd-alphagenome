Christopher Lawrence, BSc, MBBS, LLM, MD(Res), FRCP
Consultant Nephrologist
9 Harley Street, London W1G 9QY, United Kingdom
Christopher.lawrence3@nhs.net
ORCID 0000-0002-8159-0879

2 October 2026

The Editors
Journal of Medical Genetics

Dear Editors,

I am submitting "Predicting splicing variants in ADPKD: a functional benchmark of AlphaGenome and
SpliceAI at PKD1 and PKD2" for consideration as Original research.

Diagnostic laboratories now use computational predictors to decide which non-coding variant in a
patient with autosomal dominant polycystic kidney disease deserves an RNA study. About one patient
in ten has no causal variant found after coding-region testing, so that decision is made often.
SpliceAI is the incumbent. AlphaGenome arrived this year with a precomputed variant impact score,
and laboratories are already asking whether to adopt it. Nobody had tested either tool at PKD1 or
PKD2 against variants with published experimental evidence.

This paper does that. From a recorded search of 391 records I curated 128 PKD1 and PKD2 variants
tested experimentally in 25 published studies. Every entry carries a verbatim quotation from its
source, checked mechanically against the stored text, and a description validated against GRCh38 by
VariantValidator. Variants are assigned to the class of experiment that produced their evidence,
because a reporter of protein output is not evidence about splicing, and conflating the two inflates
any benchmark built from them.

The answer is that AlphaGenome offers no accuracy gain over SpliceAI at these loci. The two had the
same area under the curve, 0.85 each, with a paired difference of 0.000 (-0.033 to 0.037) on
identical variants. Both miss about one demonstrated effect in five, and both miss every variant
that acts after transcription, including a 5' untranslated region variant that lowers protein output
by 87%. Which tool ranks candidates better depends on the comparison set used, which is itself worth
knowing before a laboratory quotes a rank to a clinician.

I believe this belongs in your journal rather than a nephrology one. The question is a variant
interpretation question: how much weight can a laboratory place on a computational prediction when
no RNA study has been done? The answer bears on practice now, and the curated set is a resource
others can reuse and extend for these two genes.

The work reports no new patient data. It is a secondary analysis of published experimental results
together with model predictions obtained under AlphaGenome's non-commercial research terms. The
truth set and analysis plan were frozen, with a recorded hash, before scoring began; two amendments
made after seeing results are labelled exploratory and dated. The full archive, including the
labelled data set, every stored score and a script that regenerates every number without repeating a
model query, is openly available at doi:10.5281/zenodo.23103370.

In line with your policy on generative AI, I disclose that I used Claude (Opus 5, Anthropic) to
curate and cross-check the literature search, to write and revise the analysis code, and to draft
and revise the text. Every truth-set entry is verified against its source and every number is
regenerated from the deposited score tables, independently of how that code was written. A
declaration appears in the manuscript. The tool is not listed as an author, and I take full
responsibility for the content.

The manuscript is original, is not under consideration elsewhere, and has not been submitted to
another journal. I am the sole author. My competing interests are declared in full in the
manuscript: I am a director and shareholder of companies developing clinical software and genomic
diagnostics and am named on patent applications in that field, none of which relates to PKD1, PKD2,
polycystic kidney disease or any tool evaluated here, and none of which funded this work. No funding
was received.

Yours faithfully,

Christopher Lawrence
