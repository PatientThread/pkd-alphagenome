# AlphaGenome at PKD1 and PKD2

A locus-specific evaluation of AlphaGenome (Atlas AVI scores and on-demand predictions) against a
curated, literature-derived truth set of functionally tested non-coding and splice-acting variants in
*PKD1* and *PKD2*, with SpliceAI and CADD as comparators.

**Author:** Christopher Lawrence, Consultant Nephrologist, 9 Harley Street, London. ORCID
0000-0002-8159-0879. This is individual, non-commercial research. The AlphaGenome API key used is a
personal key; no company account or resource was used.

**Status (19 September 2026, late):** truth set built and verified (121 variants, 51 positive and 68
negative); analysis plan frozen with five dated amendments; Atlas, on-demand, SpliceAI and CADD
scoring complete for both genes; three figures drafted; manuscript draft v0.1 started. Outstanding:
the gene-wide SpliceAI background (running), Introduction and Discussion.

---

## Read these first

| File | What it is |
|---|---|
| `docs/ANALYSIS_PLAN.md` | Pre-specified question, hypotheses, endpoints, amendments A1 and A2 |
| `docs/PLAN_FREEZE.txt` | Hashes and timestamps fixing the plan before scoring, and one disclosed exception |
| `docs/MANUAL_DOWNLOADS.md` | Papers that could not be retrieved automatically, and why they matter |
| `data/curated/truthset_verified.tsv` | The truth set: every row quote-checked and VariantValidator-checked |
| `data/curated/excluded.tsv` | Variants deliberately left out, with the reason |

## Pipeline

Python 3.12. `~/.venvs/alphagenome` (alphagenome 0.9.0, pandas, numpy, scipy, pysam) for everything
except SpliceAI, which runs in `~/.venvs/spliceai` (spliceai 1.3.1, TensorFlow 2.21, tf-keras).

| Step | Script | Does |
|---|---|---|
| 01 | `01_literature_search.py` | Recorded PubMed search (two fixed queries, 391 unique papers) |
| 02 | `02_screen_and_fetch.py` | Title screen (47 included) and full-text / abstract retrieval |
| 03 | `03_extract_mentions.py` | Mechanical extraction of every variant mention with its sentence |
| 04 | `04_build_truthset.py` | The curated truth set, each entry with a verbatim quote |
| 05 | `05_verify.py` | Quote check against source text; VariantValidator on GRCh38 |
| 06 | `06_score_atlas.py` | Atlas AVI: locus backgrounds and truth-set attributions |
| 07 | `07_score_ondemand.py` | AlphaGenome on demand: splicing and RNA-seq, all and kidney tracks |
| 08 | `08_score_spliceai.py` | SpliceAI, local, validated against a published score |
| 09 | `09_score_cadd.py` | CADD v1.7: locus backgrounds (remote tabix) and truth set |
| 10 | `10_annotate.py` | Coding-consequence class and ClinVar status |
| 11 | `11_gnomad_overlap.py` | gnomAD v4.1 overlap, for the AVI training-leakage rule |
| 12 | `12_analyse.py` | The pre-specified analysis |
| 13 | `13_spliceai_locus.py` | SpliceAI across the locus, batched and validated against the stock tool |
| 14 | `14_therapeutic_sites.py` | The two drug-target sites (exploratory aim A4) |
| 15 | `15_figures.py` | Manuscript figures |
| - | `p_to_c.py` | Protein change to DNA change, only where the answer is unique |

## Truth set in one paragraph

121 unique variants from 24 primary sources: 51 positive (functional assay showed aberrant splicing
or reduced expression), 68 negative (functionally tested, indistinguishable from wild type), 1
indeterminate and 1 discordant between laboratories. Every entry carries a verbatim quote from its source that is machine-checked against
the stored text, and every HGVS description was validated by VariantValidator against the MANE
Select transcripts NM_001009944.3 (*PKD1*) and NM_000297.4 (*PKD2*). For two variants the authors
printed hg38 coordinates themselves; VariantValidator's coordinates match both exactly.

## Licensing, in brief

Analysis code: MIT (see `LICENSE`, whose scope note matters). AlphaGenome and Atlas predictions:
non-commercial only, governed by `LEGALLY_BINDING_TERMS_OF_USE.txt`, which lists every affected file.
CADD scores: CADD's non-commercial terms. Publisher full texts are not redistributed (`data/papers/`
is excluded from any repository).

The API key lives at `~/.alphagenome/api_key` (mode 600), outside this tree, and is loaded by
`analysis/ag_auth.py`. Never place it inside this directory.
