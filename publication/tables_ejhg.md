**Table 1. Composition of the benchmark by type of experimental evidence.**

| Evidence type | Variant class | Effect shown | No effect detected | From patients or populations | Genes |
| --- | --- | --- | --- | --- | --- |
| splicing | noncoding | 29 | 22 | 22 | PKD1, PKD2 |
| splicing | synonymous | 5 | 3 | 2 | PKD1, PKD2 |
| splicing | nonsense | 5 | 2 | 0 | PKD1, PKD2 |
| splicing | missense | 10 | 43 | 3 | PKD1, PKD2 |
| splicing | other coding | 1 | 0 | 1 | PKD1 |
| reporter output | noncoding | 4 | 2 | 6 | PKD1 |

Provenance is reported in the fifth column because it is strongly associated with outcome: in the splicing class there are no patient-derived variants in which no effect was detected.

**Table 2. Discrimination within the splicing class.**

| Score | Analysis population | Positives | Negatives | AUC | 95% CI |
| --- | --- | --- | --- | --- | --- |
| SpliceAI | all splicing-class variants | 50 | 70 | 0.846 | 0.76 to 0.92 |
| AlphaGenome splicing composite | all splicing-class variants | 50 | 70 | 0.846 | 0.77 to 0.92 |
| AVI | all splicing-class variants (5 unscored) | 45 | 70 | 0.718 | 0.62 to 0.81 |
| CADD | all splicing-class variants (5 unscored) | 45 | 70 | 0.565 | 0.45 to 0.67 |
| SpliceAI | PKD1 only | 41 | 50 | 0.824 | 0.73 to 0.91 |
| AlphaGenome splicing composite | PKD1 only | 41 | 50 | 0.835 | 0.74 to 0.91 |
| AVI | PKD1 only | 36 | 50 | 0.722 | 0.61 to 0.83 |
| CADD | PKD1 only | 36 | 50 | 0.557 | 0.43 to 0.68 |
| AVI | PKD1, non-coding and synonymous | 25 | 18 | 0.793 | 0.63 to 0.94 |
| CADD | PKD1, non-coding and synonymous | 25 | 18 | 0.731 | 0.56 to 0.88 |
| SpliceAI | scored by all four tools | 45 | 70 | 0.832 | 0.75 to 0.91 |
| AlphaGenome splicing composite | scored by all four tools | 45 | 70 | 0.833 | 0.75 to 0.91 |
| AVI | scored by all four tools | 45 | 70 | 0.718 | 0.62 to 0.81 |
| CADD | scored by all four tools | 45 | 70 | 0.565 | 0.45 to 0.67 |

Paired difference, AlphaGenome composite minus SpliceAI, on the 50 positives and 70 negatives both tools score: +0.000 (95% CI -0.033 to +0.037); resampling whole source studies, +0.000 (-0.042 to +0.017). Equal observed AUCs do not imply equal per-variant behaviour.

AVI and CADD are reported primarily for PKD1, which lies on a chromosome absent from AVI's training and validation lists; PKD2 lies on a training chromosome, so overlap is not excluded for it, though neither is it demonstrated. The pooled and PKD1-only analyses agree. The 5 variants absent from the four-tool comparison are insertions or deletions missing from the precomputed files used, and all have a demonstrated effect, so the absent values are not missing at random.

**Table 3. Robustness of the primary comparison.**

| Analysis | SpliceAI | AlphaGenome composite | n (pos/neg) |
| --- | --- | --- | --- |
| Primary (all splicing-class variants) | 0.846 | 0.846 | 50/70 |
| Study-cluster bootstrap (95% interval) | 0.71 to 1.00 | 0.71 to 0.99 | 50/70 |
| Leave-one-study-out (range) | 0.805 to 0.920 | 0.797 to 0.916 | 21 studies |
| Subgroup: PKD1 | 0.824 | 0.835 | 41/50 |
| Subgroup: PKD2 | 0.928 | 0.917 | 9/20 |
| Subgroup: canonical | undefined, no controls | undefined, no controls | 5/0 |
| Subgroup: deep intronic (>20 nt) | withheld, sparse controls | withheld, sparse controls | 6/1 |
| Subgroup: exonic or UTR | 0.819 | 0.836 | 20/48 |
| Subgroup: near-splice (3-20 nt) | 0.791 | 0.789 | 19/21 |
| Evidence from engineered variants | 0.734 | 0.756 | 23/69 |
| Evidence from patient variants | undefined, no controls | undefined, no controls | 26/0 |
| Evidence from population variants | withheld, sparse controls | withheld, sparse controls | 1/1 |

Two different reasons are distinguished. Where a subgroup contains no variants of one outcome class an AUC is undefined: this applies to the canonical splice-site and the patient-derived subgroups, the latter being the important one, since it means no discrimination estimate on patient-derived variants is available from these data. Where both classes are present but the controls number one or two, an AUC is defined but far too unstable to report, and is withheld rather than shown.

**Table 4. The reporter-output class: 5' untranslated region variants assayed by protein output. Reporter RNA was measured for three of the six.**

| Variant | Assay | Transcript measured | Measured effect | AlphaGenome splicing composite | AlphaGenome RNA (log2) | SpliceAI | Source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| c.-209G>A | Gaussia luciferase normalised to secreted alkaline phosphatase, protein output | no RNA study performed | no difference from wild type | 0.044 | 0.037 | 0.00 | 41006799 |
| c.-20A>G | Renilla/firefly dual-luciferase reporter, protein output | yes, RT-qPCR: unchanged | higher output (uORF2 start codon disrupted) | 0.041 | 0.079 | 0.00 | 39757590 |
| c.-23G>T | Renilla/firefly dual-luciferase reporter, protein output | yes, RT-qPCR: unchanged | about 42% higher output (uORF2 -3 nucleotide) | 0.046 | 0.075 | 0.00 | 39757590 |
| c.-52C>T | Gaussia luciferase normalised to secreted alkaline phosphatase, protein output | no RNA study performed | no difference from wild type | 0.030 | 0.019 | 0.00 | 41006799 |
| c.-66T>C | Renilla/firefly dual-luciferase reporter, protein output | yes, RT-qPCR: unchanged | lower output (uORF1 stop codon removed) | 0.044 | 0.070 | 0.00 | 39757590 |
| c.-69dupG | Gaussia luciferase normalised to secreted alkaline phosphatase, protein output | no RNA study performed | 87% lower output (uORF1 extended) | 0.033 | 0.047 | 0.01 | 41006799 |

The AlphaGenome RNA column is the RNA_SEQ variant score: the maximum ABSOLUTE log2 predicted change across all RNA-seq tracks, restricted to the PKD1 gene row, at the recommended 1,048,576 bp context. It is an unsigned magnitude, so it states how large a change is predicted, not its direction. Direction of the measured effect is stated because two variants increase output. Chen et al. [39757590] measured reporter RNA by RT-qPCR and found it unchanged, so for those three variants the effect is attributable to translation. Wedd et al. [41006799] performed no RNA study and state that effects on pre-mRNA splicing and transcript stability cannot be excluded, so for their three variants the mechanism of the change in protein output is undetermined. Predicted values are shown for the same variants; neither predictor addresses protein output.
