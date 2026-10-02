# Gene and variant nomenclature statement

Submitted separately, as requested by the journal, in support of the manuscript
"Predicting splicing variants in ADPKD: a functional benchmark of AlphaGenome and SpliceAI at *PKD1*
and *PKD2*".

## Genes

Gene symbols follow HGNC: *PKD1* (HGNC:9008) and *PKD2* (HGNC:9009). Gene symbols are italicised
throughout; protein products are given as polycystin-1 and polycystin-2 and are not italicised.

## Reference sequences

All variant descriptions are given on the MANE Select transcripts on GRCh38:

| Gene | MANE Select transcript | Protein | Genomic reference |
|---|---|---|---|
| *PKD1* | NM_001009944.3 | NP_001009944.3 | NC_000016.10 |
| *PKD2* | NM_000297.4 | NP_000288.1 | NC_000004.12 |

Coding-sequence numbering begins at the A of the initiator ATG, so 5' untranslated region positions
carry a minus sign (for example c.-20A>G) and 3' untranslated region positions an asterisk (for
example c.*153).

## Variant descriptions

Descriptions follow the Human Genome Variation Society recommendations. Every description in this
manuscript, its tables and its supplementary material was validated mechanically with
VariantValidator against the transcripts above, which also supplied the genomic coordinates used for
scoring. All 131 source-level entries passed. Those 131 entries collapse to the 128 unique variants reported
in the manuscript, because three variants were described independently in more than one source
publication and each description was verified separately. No description is carried over unchecked
from its source publication.

Two conversions were necessary and are documented in the repository:

1. **Legacy IVS notation.** Descriptions such as IVS30-25A>G were converted to coding-sequence
   notation using the RefSeq exon table for the transcript above.
2. **Protein-only descriptions.** Where a source gave only an amino-acid change, the codon was read
   from the reference coding sequence and a result accepted only where exactly one single-base
   substitution produced the stated change. All 17 such conversions were unique, and the method
   exactly reproduced the eight variants whose papers also printed coding-sequence notation.

Where a source publication's own description was internally inconsistent, this is recorded in the
curated data set with the verbatim quotation from that source, and the variant is either resolved
against an independent line of evidence or excluded. No inconsistency was resolved silently.

## Reference single-nucleotide variant identifiers

dbSNP identifiers are given only where the source publication gives them, and always alongside the
coding-sequence description: rs3874648 is *PKD1* c.10051-239G>A on NM_001009944.3.
