# SNP GWAS Output: Column Definitions

This document describes the columns contained in the SNP genome-wide association study (GWAS) output files for three clinical outcomes:

- `SNPs_Mortality_GWAS.csv` — mortality
- `SNPs_ICU_GWAS.csv` — ICU admission
- `SNPs_LOS_GWAS.csv` — length of stay (LOS)

All three files contain an initial set of SNP–phenotype association statistics generated using **Scoary**, followed by outcome-specific regression analyses and genomic annotation.

The regression analyses differ according to the clinical outcome:

| File | Outcome | Statistical analyses |
|---|---|---|
| `SNPs_Mortality_GWAS.csv` | Mortality | Cox proportional hazards, crude logistic, adjusted logistic |
| `SNPs_ICU_GWAS.csv` | ICU admission | Crude logistic, adjusted logistic |
| `SNPs_LOS_GWAS.csv` | Length of stay | Crude linear regression, adjusted linear regression |

## 1. Common Scoary GWAS columns

| Column | Description |
|---|---|
| `SNP` | SNP identifier. |
| `Number_pos_present_in` | Number of phenotype-positive isolates carrying the SNP. |
| `Number_neg_present_in` | Number of phenotype-negative isolates carrying the SNP. |
| `Number_pos_not_present_in` | Number of phenotype-positive isolates not carrying the SNP. |
| `Number_neg_not_present_in` | Number of phenotype-negative isolates not carrying the SNP. |
| `Sensitivity` | Sensitivity obtained when SNP presence is used to predict phenotype positivity. |
| `Specificity` | Specificity obtained when SNP absence is used to predict phenotype negativity. |
| `Odds_ratio` | Odds ratio for the association between SNP presence and the phenotype. |
| `Naive_p` | Unadjusted p-value for the SNP–phenotype association. |
| `Bonferroni_p` | P-value corrected for multiple testing using the Bonferroni method. |
| `Benjamini_H_p` | P-value corrected using the Benjamini–Hochberg false discovery rate (FDR) procedure. |
| `Max_Pairwise_comparisons` | Maximum number of non-intersecting isolate pairs contrasting in both SNP and phenotype status. |
| `Max_supporting_pairs` | Maximum number of contrasting pairs supporting the observed association. |
| `Max_opposing_pairs` | Maximum number of contrasting pairs opposing the observed association. |
| `Best_pairwise_comp_p` | Lowest possible p-value obtained from the maximum set of contrasting phylogenetic pairs. |
| `Worst_pairwise_comp_p` | Highest possible p-value obtained from the maximum set of contrasting phylogenetic pairs. |

The Scoary pairwise-comparison statistics provide information on whether SNP–phenotype associations remain supported when phylogenetic relationships among isolates are considered.

## 2. Mortality GWAS

**File:** `SNPs_Mortality_GWAS.csv`

Mortality was evaluated using Cox proportional hazards analysis and crude and confounder-adjusted logistic regression.

### Cox proportional hazards analysis

| Column | Description |
|---|---|
| `cox_coef` | Cox regression coefficient (log hazard ratio) associated with SNP presence. |
| `cox_HR` | Hazard ratio (HR), calculated as `exp(cox_coef)`. |
| `cox_SE` | Standard error of the Cox regression coefficient. |
| `cox_z` | Wald z statistic for the Cox coefficient. |
| `cox_p` | P-value for the SNP effect in the Cox proportional hazards model. |

`cox_HR > 1` indicates an increased hazard of the defined event, whereas `cox_HR < 1` indicates a decreased hazard, subject to the coding of the event variable.

### Crude logistic analysis

| Column | Description |
|---|---|
| `crude_OR` | Unadjusted odds ratio for the association between SNP presence and mortality. |
| `crude_lower` | Lower confidence limit of the crude odds ratio. |
| `crude_upper` | Upper confidence limit of the crude odds ratio. |
| `crude_midp_p` | Mid-P exact p-value for the crude SNP–mortality association. |
| `crude_fisher_p` | Fisher's exact test p-value for the crude association. |
| `crude_chisq_p` | Chi-square test p-value for the crude association. |

### Adjusted logistic regression

| Column | Description |
|---|---|
| `adjusted_logit_coef` | SNP coefficient (log odds ratio) after adjustment for the included confounder(s). |
| `adjusted_logit_SE` | Standard error of the adjusted SNP coefficient. |
| `adjusted_logit_z` | Wald z statistic for the adjusted SNP coefficient. |
| `adjusted_logit_p` | P-value for the SNP after adjustment for the included confounder(s). |

The adjusted odds ratio can be calculated as:

`adjusted OR = exp(adjusted_logit_coef)`

## 3. ICU admission GWAS

**File:** `SNPs_ICU_GWAS.csv`

ICU admission was analysed as a binary outcome using crude and confounder-adjusted logistic regression.

### Crude logistic regression

| Column | Description |
|---|---|
| `coef_crude` | SNP coefficient (log odds ratio) from the crude logistic regression model. |
| `exp(coef)_crude` | Exponentiated crude regression coefficient, corresponding to the crude odds ratio. |
| `se(coef)_crude` | Standard error of the crude logistic regression coefficient. |
| `z_crude` | Wald z statistic for the crude SNP coefficient. |
| `Pr(>|z|)_crude` | P-value for the SNP coefficient in the crude logistic regression model. |
| `estimate_crude` | Crude odds-ratio estimate from the corresponding association analysis. |
| `lower_crude` | Lower confidence limit for the crude odds ratio. |
| `upper_crude` | Upper confidence limit for the crude odds ratio. |
| `midp.exact_crude` | Mid-P exact p-value for the crude SNP–ICU association. |
| `fisher.exact_crude` | Fisher's exact test p-value for the crude association. |
| `chi.square_crude` | Chi-square test p-value for the crude association. |

### Adjusted logistic regression

| Column | Description |
|---|---|
| `Estimate_adjusted` | SNP coefficient (log odds ratio) after adjustment for the included confounder(s). |
| `Std. Error_adjusted` | Standard error of the adjusted SNP coefficient. |
| `z value_adjusted` | Wald z statistic for the adjusted SNP coefficient. |
| `Pr(>|z|)_adjusted` | P-value for the SNP coefficient in the adjusted logistic regression model. |

The adjusted odds ratio can be calculated as:

`adjusted OR = exp(Estimate_adjusted)`

## 4. Length-of-stay GWAS

**File:** `SNPs_LOS_GWAS.csv`

Length of stay (LOS) was analysed as a continuous outcome using crude and confounder-adjusted linear regression.

### Crude model

| Column | Description |
|---|---|
| `Df` | Degrees of freedom associated with the model term. |
| `Sum Sq_crude` | Sum of squares attributable to the SNP in the crude model. |
| `Mean Sq_crude` | Mean square for the SNP term in the crude model. |
| `F value_crude` | F statistic for the SNP effect in the crude model. |
| `Pr(>F)_crude` | P-value corresponding to the crude-model F statistic. |
| `Estimate_crude` | Estimated crude regression coefficient for SNP presence. |
| `Std. Error_crude` | Standard error of the crude regression coefficient. |
| `t value_crude` | t statistic for the crude SNP coefficient. |
| `Pr(>|t|)_crude` | P-value for the SNP coefficient in the crude linear regression model. |

A positive `Estimate_crude` indicates an increase in LOS associated with SNP presence, whereas a negative estimate indicates a decrease, according to the coding and scale used for LOS.

### Adjusted model

| Column | Description |
|---|---|
| `Estimate_adjusted` | Regression coefficient for SNP presence after adjustment for the included confounder(s). |
| `Std. Error_adjusted` | Standard error of the adjusted SNP coefficient. |
| `t value_adjusted` | t statistic for the adjusted SNP coefficient. |
| `Pr(>|t|)_adjusted` | P-value for the SNP coefficient after adjustment for the included confounder(s). |

`Estimate_adjusted` represents the estimated difference in LOS associated with SNP presence after controlling for the confounder(s) included in the model.

## 5. SNP annotation columns

| Column | Description |
|---|---|
| `CHROM` | Reference chromosome, contig, or sequence containing the SNP. |
| `POS` | Genomic position of the SNP. |
| `TYPE` | Type of sequence variant. |
| `REF` | Reference allele. |
| `ALT` | Alternative allele. |
| `EVIDENCE` | Supporting information associated with the variant annotation. |
| `FTYPE` | Type of genomic feature affected by the SNP. |
| `STRAND` | Strand of the affected genomic feature. |
| `NT_POS` | Nucleotide position of the SNP within the affected feature. |
| `AA_POS` | Amino-acid position affected by the SNP, where applicable. |
| `EFFECT` | Predicted functional consequence of the SNP. |
| `LOCUS_TAG` | Locus tag of the affected gene. |
| `GENE` | Gene associated with the SNP. |
| `PRODUCT` | Annotated product or function of the affected gene. |

## 6. Nucleotide-frequency columns

| Column | Description |
|---|---|
| `missing` | Proportion of isolates with missing nucleotide information at the corresponding position. |
| `A_ratio` | Proportion of nucleotide calls that are adenine (A). |
| `C_ratio` | Proportion of nucleotide calls that are cytosine (C). |
| `T_ratio` | Proportion of nucleotide calls that are thymine (T). |
| `G_ratio` | Proportion of nucleotide calls that are guanine (G). |

## Statistical interpretation

For **mortality**, Cox regression evaluates the association between SNP presence and the time-to-event outcome, while logistic regression evaluates the association between SNP presence and mortality as a binary outcome.

For **ICU admission**, logistic regression evaluates the association between SNP presence and the probability of ICU admission.

For **length of stay**, linear regression evaluates the association between SNP presence and LOS as a continuous outcome.

For each outcome, the **crude model** represents the unadjusted SNP–outcome association, whereas the **adjusted model** estimates the SNP association after accounting for the confounder(s) included in the corresponding regression model.

Scoary statistics represent the initial bacterial GWAS association analysis and should not be interpreted as equivalent to the outcome-specific regression estimates.

## Missing values

`NA` indicates that a statistic or genomic annotation was unavailable, could not be calculated, or was not applicable to the corresponding SNP.
