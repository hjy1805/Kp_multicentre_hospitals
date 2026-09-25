# Pangenome GWAS Output: Column Definitions

This document describes the columns contained in the **pangenome genome-wide association study (GWAS)** output files.

The files combine gene presence/absence association statistics generated using **Scoary** with additional quality-control, sequence-similarity, and Panaroo presence information.

## 1. Scoary pangenome GWAS columns

| Column | Description |
|---|---|
| `Gene` | Gene or pangenome feature identifier used in the association analysis. |
| `Non-unique Gene name` | Alternative or non-unique gene name assigned to the pangenome feature, where available. |
| `Annotation` | Functional annotation associated with the gene or pangenome feature. |
| `Number_pos_present_in` | Number of phenotype-positive isolates in which the gene is present. |
| `Number_neg_present_in` | Number of phenotype-negative isolates in which the gene is present. |
| `Number_pos_not_present_in` | Number of phenotype-positive isolates in which the gene is absent. |
| `Number_neg_not_present_in` | Number of phenotype-negative isolates in which the gene is absent. |
| `Sensitivity` | Sensitivity obtained when gene presence is used to predict phenotype positivity. |
| `Specificity` | Specificity obtained when gene absence is used to predict phenotype negativity. |
| `Odds_ratio` | Odds ratio for the association between gene presence and the phenotype. Values >1 indicate enrichment of the gene among phenotype-positive isolates, whereas values <1 indicate enrichment among phenotype-negative isolates. |
| `Naive_p` | Unadjusted p-value for the gene–phenotype association. |
| `Bonferroni_p` | P-value corrected for multiple testing using the Bonferroni method. |
| `Benjamini_H_p` | P-value corrected using the Benjamini–Hochberg false discovery rate (FDR) procedure. |
| `Max_Pairwise_comparisons` | Maximum number of non-intersecting isolate pairs contrasting in both gene presence/absence and phenotype status. |
| `Max_supporting_pairs` | Maximum number of contrasting pairs supporting the observed direction of association. |
| `Max_opposing_pairs` | Maximum number of contrasting pairs opposing the observed direction of association. |
| `Best_pairwise_comp_p` | Lowest possible p-value obtained from the maximum set of contrasting phylogenetic pairs. |
| `Worst_pairwise_comp_p` | Highest possible p-value obtained from the maximum set of contrasting phylogenetic pairs. |

The Scoary pairwise-comparison statistics provide information on whether a gene–phenotype association remains supported after considering the phylogenetic relationships among isolates.

## 2. Quality-control and annotation columns

| Column | Description |
|---|---|
| `QC_Correlation` | Correlation-based quality-control metric used to evaluate agreement or consistency of the gene presence/absence pattern with the corresponding validation data. |
| `QC_report` | Quality-control result or status assigned to the pangenome feature based on the QC procedure. |
| `blast` | Result of the BLAST-based sequence similarity search used to support or validate the identity/annotation of the pangenome feature. |
| `panaroo_present` | Indicator of whether the corresponding gene/pangenome feature is present according to the Panaroo pangenome analysis. |

## 3. Interpretation of association statistics

For each pangenome feature, Scoary compares gene presence/absence with the phenotype of interest.

An `Odds_ratio > 1` indicates that gene presence is associated with higher odds of the positive phenotype, whereas an `Odds_ratio < 1` indicates that the gene is more frequent among phenotype-negative isolates.

`Naive_p` represents the uncorrected association p-value. Because many pangenome features are tested simultaneously, the multiple-testing-adjusted values (`Bonferroni_p` and `Benjamini_H_p`) should be considered when identifying statistically supported associations.

The pairwise-comparison columns provide an additional phylogenetically informed assessment of the association.

## 4. Pangenome and validation information

The final four columns provide additional information for evaluating GWAS hits beyond the initial Scoary association.

`QC_Correlation` and `QC_report` provide quality-control information for the feature.

`blast` provides sequence-similarity information from BLAST and can be used to support gene identification or functional annotation.

`panaroo_present` provides the corresponding gene-presence information obtained from the Panaroo pangenome analysis.

Together, these fields allow statistically associated pangenome features to be cross-checked against sequence-level evidence and the original Panaroo gene presence/absence calls.

## 5. Missing values

`NA` indicates that a value, annotation, QC result, or sequence-similarity result was unavailable or not applicable to the corresponding pangenome feature.
