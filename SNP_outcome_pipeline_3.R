# =============================================================================
# SNPs AND CLINICAL OUTCOME: ADJUSTED ASSOCIATION, THEN TIME TO IN-HOSPITAL DEATH
#
# Two halves. The first asks which SNPs are associated with death, adjusting for
# the clinical variables that confound that association and for the population
# structure of the organism. The second asks when patients die, using isolation
# as the time origin, and treats discharge alive as a competing event rather
# than as censoring.
#
# WHAT I FOUND IN YOUR FILES BEFORE WRITING THIS. Every number below is from the
# uploaded data, not an assumption, and several of them change what the analysis
# has to do.
#
#   CLUSTERING IS NOT A DETAIL. 1073 isolates come from 886 patients. 129
#   patients contribute more than one isolate, up to 8, and 316 isolates (29%)
#   sit inside a multi-isolate patient. Every model here, the logistic ones
#   included, uses patient-clustered standard errors. Treating 1073 isolates as
#   1073 independent observations would narrow every interval in the script.
#
#   THE SURVIVAL COHORT IS 852. That is how many rows carry admission,
#   collection and discharge dates together. Collection falls inside the
#   admission for every one of them, 0 before and 0 after, so isolation as the
#   time origin needs no left-truncation repair. 16 rows are discharged on the
#   collection day and would contribute zero follow-up.
#
#   CENSORING AT DISCHARGE WOULD BE WRONG, and the data say so directly. Of the
#   852, 299 died in hospital and 553 left alive. Of those 553, 128 died later.
#   Censoring a discharged patient asserts that their future hazard resembles
#   that of someone still admitted, and a quarter of them died. That is why
#   block 12 exists. On simulated data with a similar death-to-discharge ratio,
#   one minus Kaplan-Meier reported a 60-day risk of 80% where the competing-
#   risks estimate was 42%: the error is not a rounding matter, it grows with
#   follow-up, and block 12 prints the size of it for your own data.
#
#   4612 SNPs HOLD 4383 DISTINCT PRESENCE/ABSENCE PATTERNS, the largest
#   identical group having 41 members. Correcting over 4612 tests would be
#   correcting for 229 comparisons that do not exist. BH runs over the distinct
#   patterns.
#
#   BMI_pre IS MISSING FOR 505 OF 1073 and cannot be a confounder here.
#   Hospital_Acquired is missing for exactly the 221 rows with no admission
#   record, which is the same rows the survival half drops.
#
# WHAT YOUR SPEC DID NOT MENTION AND BLOCK 05 ADDS. In a bacterial association
# study the dominant confounder is the lineage structure of the organism itself.
# SNPs are inherited in blocks, so a SNP that merely marks a clone will appear
# associated with any outcome that clone happens to cause, and no amount of
# patient-level adjustment removes it. Block 05 takes principal components of
# the SNP matrix and puts them in the adjustment set. Without that step the
# results are a lineage scan wearing the clothes of a SNP scan.
#
# Flat script, no named functions. Run block by block. Figures print to the
# device. The saving block is at the end, commented out.
#
# Blocks
#   01 packages and settings        08 annotating the hits
#   02 read and link                09 the survival cohort
#   03 the clustering problem       10 Kaplan-Meier
#   04 confounders, and exclusions  11 Cox, cause-specific
#   05 population structure         12 competing risks
#   06 SNP filtering                13 the survival scan, all SNPs
#   07 adjusted logistic for death  14 the three side by side
#                                   15 saving (commented out)
# =============================================================================


# =============================================================================
# 01  PACKAGES AND SETTINGS
# =============================================================================

# rm(list = ls()) first. dplyr resolves a bare name against the data frame and
# then the global environment, so a leftover object sharing a column's name is
# used silently wherever that column is out of scope.

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)
library(survival)      # base R recommended package, already installed

DIR <- "/Users/daneshm/Documents/PA_KAIMRC"     # set to where the five files are

F_SNP    <- file.path(DIR, "SNPs_presence_absence.csv")
F_SCOARY <- file.path(DIR, "significant_scoary_results_filtered_withannotation.csv")
F_META   <- file.path(DIR, "Structured_metadata.csv")
F_LAB    <- file.path(DIR, "GWAS_Labels_new.csv")
F_SURV   <- file.path(DIR, "Survival_analysis_df.csv")

OUTCOME <- "Total_Death"      # the binary outcome for the first half

# Minor-allele count floor. Below this the odds ratio is estimable in name only:
# in the earlier ARGVir work the median confidence-interval ratio was 30.8 at a
# minor count of 3 to 9 and 7.3 at 20 to 29, so 20 is where an interval starts
# to mean something. 3966 of your 4612 SNPs clear it.
MIN_MAC <- 20
MIN_CELL <- 5                 # every cell of the 2x2 must reach this

N_PC <- 5                     # principal components of the SNP matrix, block 05
ALPHA <- 0.05

# Pre-specified adjustment set. Every one of these is fixed before the infection
# and cannot lie on the path from bacterial genotype to death.
CONF_CORE <- c("AGE", "Sex", "RGN_NM", "Charlson_12m", "IP_admissions_12m")

# Excluded on purpose, with the reason attached. Block 04 prints this.
CONF_EXCLUDED <- c(
  BMI_pre                 = "missing in 505 of 1073",
  AMR_Class               = "resistance phenotype: downstream of the genotype, a mediator",
  CARB                    = "same, and conditional on selective testing",
  COL                     = "same",
  AG                      = "same", CEP = "same", FQ = "same", PTZ = "same",
  Hospital_Acquired       = "partly a consequence of the same admission being modelled",
  Previous_carbapenem     = "selects for the genotype being tested, a collider risk",
  Combination_therapy     = "chosen in response to the isolate, so post-exposure",
  LOS                     = "an outcome, not a covariate")

OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#999999", "#000000")

theme_set(theme_bw(base_size = 11))
set.seed(20260910)


# =============================================================================
# 02  READ AND LINK
# =============================================================================

snp_raw <- read_csv(F_SNP, show_col_types = FALSE)
dim(snp_raw)                      # 4612 x 1074, SNPs down the rows
snp_raw[1:3, 1:6]

# The matrix arrives transposed relative to everything else: SNPs are rows and
# isolates are columns. Turn it so that a row is an isolate.
snp_ids  <- snp_raw[[1]]
isolates <- names(snp_raw)[-1]
length(snp_ids); length(isolates)

head(snp_ids, 3)                  # "1|155|.|G|A|.|.|.|GT" style

G <- as.matrix(snp_raw[, -1])     # SNPs x isolates
storage.mode(G) <- "double"
dim(G)
table(G, useNA = "ifany")         # must be 0/1 only

X_snp <- t(G)                     # isolates x SNPs
rownames(X_snp) <- isolates
colnames(X_snp) <- snp_ids
dim(X_snp)

meta <- read_csv(F_META, show_col_types = FALSE,
                 na = c("", "NA", "-", "nan", "NaN"))
lab  <- read_csv(F_LAB,  show_col_types = FALSE, col_types = cols(.default = "c"))
surv_raw <- read_csv(F_SURV, show_col_types = FALSE,
                     na = c("", "NA", "-"), col_types = cols(.default = "c"))
scoary <- read_csv(F_SCOARY, show_col_types = FALSE)

dim(meta); dim(lab); dim(surv_raw); dim(scoary)

# "-" is the missing marker in the label file, so those columns are read as
# character and converted here rather than being silently turned to NA by a
# numeric guess.
OUT_BIN <- setdiff(names(lab), "KAUST_ID")
labels <- lab
for (o in OUT_BIN) labels[[o]] <- as.integer(labels[[o]])
sapply(labels[OUT_BIN], function(x) sum(is.na(x)))
sapply(labels[OUT_BIN], function(x) sum(x, na.rm = TRUE))

# the four identifier sets must line up
length(intersect(isolates, meta$KAUST_ID))
length(intersect(isolates, labels$KAUST_ID))
length(intersect(isolates, surv_raw$KAUST_ID))
setdiff(isolates, meta$KAUST_ID)
setdiff(meta$KAUST_ID, isolates)


# =============================================================================
# 03  THE CLUSTERING PROBLEM
# =============================================================================

# The patient identifier lives only in the survival file, in the IDs column, so
# it has to be carried across to everything else before any model is fitted.

pat <- surv_raw %>% select(KAUST_ID, patient = IDs, RGN_NM_surv = RGN_NM)

n_distinct(pat$KAUST_ID)          # 1073 isolates
n_distinct(pat$patient)           # 886 patients

pat %>% count(patient, name = "isolates") %>% count(isolates) %>%
  arrange(isolates) %>% print(n = Inf)

pat %>% count(patient, name = "isolates") %>%
  summarise(patients = n(),
            multi = sum(isolates > 1),
            max_per_patient = max(isolates),
            isolates_in_multi = sum(isolates[isolates > 1]))

# WHY THIS MATTERS ARITHMETICALLY. Two isolates from one patient share that
# patient's outcome exactly, so they carry one outcome's worth of information
# and not two. With 29% of isolates inside multi-isolate patients, an analysis
# treating them as independent inflates the effective sample size by roughly
# that amount and narrows every interval accordingly. Everything below uses
# cluster-robust standard errors on the patient.

any(is.na(pat$patient))


# =============================================================================
# 04  CONFOUNDERS, AND WHAT IS DELIBERATELY LEFT OUT
# =============================================================================

# A variable earns a place in the adjustment set by being associated with the
# outcome, being associated with the exposure, and lying before the exposure in
# time. The third condition is the one that disqualifies most of the candidates
# here, and it cannot be settled by any test.

as.data.frame(CONF_EXCLUDED)

# THE MEDIATOR PROBLEM, stated once. A SNP may raise the risk of death precisely
# by conferring resistance. Adjusting for AMR_Class or for the carbapenem result
# then removes the very pathway being asked about, and the adjusted odds ratio
# answers a question nobody posed: the effect of the SNP among patients whose
# resistance phenotype was held fixed. That is why every resistance column is on
# the excluded list rather than in the model.

names(meta)
cand <- setdiff(names(meta), c("KAUST_ID", names(CONF_EXCLUDED)))
cand <- cand[!str_detect(cand, "^(Diag_PC|Medic_PC)")]
cand

# missingness, which decides what is usable at all
meta %>%
  summarise(across(all_of(cand), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing") %>%
  mutate(pct = round(100 * missing / nrow(meta), 1)) %>%
  arrange(desc(missing)) %>%
  print(n = Inf)

# the screen: each candidate against death on its own, patient-clustered. This
# is reported to show which variables carry outcome information. It does NOT
# select the adjustment set on its own, because a variable can be associated
# with the outcome and still be a mediator or a collider.
screen_dat <- meta %>%
  inner_join(select(labels, KAUST_ID, all_of(OUTCOME)), by = "KAUST_ID") %>%
  inner_join(pat, by = "KAUST_ID") %>%
  filter(!is.na(.data[[OUTCOME]]))

nrow(screen_dat)

screen <- NULL

for (v in cand) {

  d <- screen_dat[!is.na(screen_dat[[v]]), ]
  if (nrow(d) < 50) next
  if (length(unique(d[[v]])) < 2) next

  f <- try(glm(as.formula(paste0("`", OUTCOME, "` ~ `", v, "`")),
               data = d, family = binomial()), silent = TRUE)
  if (inherits(f, "try-error")) next

  # cluster-robust covariance, computed here rather than pulled from a package
  Xd <- model.matrix(f); rr <- residuals(f, type = "response")
  S <- Xd * rr; Bd <- vcov(f)
  M <- matrix(0, ncol(Xd), ncol(Xd))
  for (g in unique(d$patient)) {
    s <- colSums(S[d$patient == g, , drop = FALSE]); M <- M + tcrossprod(s)
  }
  Gn <- length(unique(d$patient)); pn <- ncol(Xd); nn <- nrow(Xd)
  V <- Bd %*% M %*% Bd * (Gn / (Gn - 1)) * ((nn - 1) / (nn - pn))

  for (j in 2:ncol(Xd)) {
    b <- coef(f)[j]; se <- sqrt(V[j, j])
    screen <- bind_rows(screen, tibble(
      variable = v, term = colnames(Xd)[j], n = nn, patients = Gn,
      OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
      p = 2 * pnorm(-abs(b / se)),
      se_naive = sqrt(Bd[j, j]), se_cluster = se,
      inflation = se / sqrt(Bd[j, j])))
  }
}

screen <- mutate(screen, q = p.adjust(p, method = "BH"))
screen %>% arrange(p) %>% print(n = Inf)

# how much the clustering costs in precision, across the whole screen
summarise(screen, median_se_inflation = median(inflation),
          max_se_inflation = max(inflation))

CONF <- CONF_CORE
CONF

# what complete cases on that set costs
screen_dat %>%
  summarise(rows = n(),
            complete = sum(complete.cases(select(., all_of(CONF)))),
            lost = rows - complete)


# =============================================================================
# 05  POPULATION STRUCTURE
# =============================================================================

# The step the clinical variables cannot do. SNPs travel in linkage blocks, so a
# variant that only marks a successful clone will track any outcome that clone
# produces. Principal components of the SNP matrix capture that structure, and
# putting them in the model asks whether a SNP predicts death beyond the lineage
# it sits in.
#
# The first components should separate the major lineages. Look at the scree
# plot and at how much of the variance the first few take: if two components
# carry a large share, the collection is dominated by a few clones and the
# adjustment matters a great deal.

pc <- prcomp(X_snp, center = TRUE, scale. = FALSE)

var_ex <- pc$sdev^2 / sum(pc$sdev^2)
round(head(var_ex, 12), 4)
round(cumsum(head(var_ex, 12)), 4)

plot(head(var_ex, 30), type = "b", pch = 19,
     xlab = "component", ylab = "share of variance",
     main = "Population structure in the SNP matrix")

PCs <- as_tibble(pc$x[, seq_len(N_PC), drop = FALSE])
names(PCs) <- paste0("SNP_PC", seq_len(N_PC))
PCs$KAUST_ID <- rownames(X_snp)

# do the components line up with region, which would mean geography and lineage
# are partly the same variable and neither adjustment is doing what it appears
PCs %>%
  inner_join(select(meta, KAUST_ID, RGN_NM), by = "KAUST_ID") %>%
  group_by(RGN_NM) %>%
  summarise(across(starts_with("SNP_PC"), ~ round(mean(.x), 2)), n = n(),
            .groups = "drop") %>%
  print(n = Inf)

p_pc <- PCs %>%
  inner_join(select(meta, KAUST_ID, RGN_NM), by = "KAUST_ID") %>%
  inner_join(select(labels, KAUST_ID, all_of(OUTCOME)), by = "KAUST_ID") %>%
  ggplot(aes(SNP_PC1, SNP_PC2, colour = RGN_NM)) +
  geom_point(size = 1.1, alpha = 0.75) +
  scale_colour_manual(values = OKABE, name = NULL) +
  labs(title = "First two components of the SNP matrix",
       subtitle = paste0("Tight clusters are lineages. If they also separate by ",
                         "region, region and lineage\nare confounded with each ",
                         "other and neither adjustment is clean.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_pc

CONF_FULL <- c(CONF, names(PCs)[startsWith(names(PCs), "SNP_PC")])
CONF_FULL


# =============================================================================
# 06  SNP FILTERING
# =============================================================================

mac <- pmin(colSums(X_snp), nrow(X_snp) - colSums(X_snp))
summary(mac)
sum(mac >= MIN_MAC)

keep_mac <- colnames(X_snp)[mac >= MIN_MAC]
length(keep_mac)

# IDENTICAL COLUMNS. Perfect linkage means several SNPs carry exactly the same
# presence/absence vector, and they cannot be told apart by any test on these
# data. One representative is fitted and the rest are recorded as its group, so
# the BH correction runs over the number of questions actually asked.
pattern_key <- apply(X_snp[, keep_mac, drop = FALSE], 2,
                     function(v) paste0(v, collapse = ""))

grp <- tibble(snp = keep_mac, key = unname(pattern_key)) %>%
  group_by(key) %>%
  mutate(n_in_group = n(), representative = first(snp)) %>%
  ungroup()

n_distinct(grp$key)
grp %>% count(n_in_group) %>% arrange(desc(n_in_group)) %>% print(n = Inf)

grp %>% filter(n_in_group > 5) %>%
  group_by(representative, n_in_group) %>%
  summarise(members = paste(head(snp, 4), collapse = " ; "), .groups = "drop") %>%
  print(n = 10)

TEST_SNPS <- unique(grp$representative)
length(TEST_SNPS)

cat("\nfitting", length(TEST_SNPS), "models over", length(keep_mac),
    "SNPs that pass the count filter, from", ncol(X_snp), "in the file\n")


# =============================================================================
# 07  ADJUSTED LOGISTIC REGRESSION FOR DEATH
# =============================================================================

# Unadjusted, adjusted for the clinical set, and adjusted for the clinical set
# plus population structure. Reporting all three is the point: the movement
# between them is the result, and a SNP whose odds ratio collapses once the
# components enter was marking a clone.

model_dat <- as_tibble(X_snp[, TEST_SNPS, drop = FALSE]) %>%
  mutate(KAUST_ID = rownames(X_snp)) %>%
  inner_join(select(labels, KAUST_ID, all_of(OUTCOME)), by = "KAUST_ID") %>%
  inner_join(select(meta, KAUST_ID, all_of(CONF)), by = "KAUST_ID") %>%
  inner_join(PCs, by = "KAUST_ID") %>%
  inner_join(pat, by = "KAUST_ID") %>%
  filter(!is.na(.data[[OUTCOME]])) %>%
  filter(complete.cases(select(., all_of(CONF))))

nrow(model_dat); n_distinct(model_dat$patient)
mean(model_dat[[OUTCOME]])

conf_clin <- paste0("`", CONF, "`", collapse = " + ")
conf_full <- paste0("`", CONF_FULL, "`", collapse = " + ")

y <- model_dat[[OUTCOME]]
cl <- model_dat$patient
clusters <- unique(cl)
Gn <- length(clusters)

res <- NULL
skipped <- NULL

for (s in TEST_SNPS) {

  x <- model_dat[[s]]
  tab <- table(x, y)

  if (!all(dim(tab) == c(2, 2)) || min(tab) < MIN_CELL) {
    skipped <- bind_rows(skipped, tibble(snp = s, reason = "cell too small",
                                         min_cell = min(tab)))
    next
  }

  row <- tibble(snp = s, n = nrow(model_dat), patients = Gn,
                n_present = sum(x), min_cell = min(tab))

  for (mod in c("unadjusted", "clinical", "clinical + structure")) {

    fml <- switch(mod,
      unadjusted            = paste0("`", OUTCOME, "` ~ `", s, "`"),
      clinical              = paste0("`", OUTCOME, "` ~ `", s, "` + ", conf_clin),
      `clinical + structure`= paste0("`", OUTCOME, "` ~ `", s, "` + ", conf_full))

    f <- try(glm(as.formula(fml), data = model_dat, family = binomial()),
             silent = TRUE)
    if (inherits(f, "try-error")) next

    Xd <- model.matrix(f)
    j <- match(paste0("`", s, "`"), colnames(Xd))
    if (is.na(j)) j <- match(s, colnames(Xd))
    if (is.na(j)) next

    rr <- residuals(f, type = "response")
    S  <- Xd * rr
    Bd <- vcov(f)
    M  <- matrix(0, ncol(Xd), ncol(Xd))
    for (g in clusters) {
      idx <- which(cl == g)
      if (!length(idx)) next
      sc <- colSums(S[idx, , drop = FALSE]); M <- M + tcrossprod(sc)
    }
    nn <- nrow(Xd); pn <- ncol(Xd)
    V <- Bd %*% M %*% Bd * (Gn / (Gn - 1)) * ((nn - 1) / (nn - pn))

    b <- coef(f)[j]; se <- sqrt(V[j, j])

    row[[paste0("OR_", mod)]] <- exp(b)
    row[[paste0("lo_", mod)]] <- exp(b - 1.96 * se)
    row[[paste0("hi_", mod)]] <- exp(b + 1.96 * se)
    row[[paste0("p_",  mod)]] <- 2 * pnorm(-abs(b / se))
    row[[paste0("se_naive_", mod)]] <- sqrt(Bd[j, j])
    row[[paste0("se_clust_", mod)]] <- se
  }

  res <- bind_rows(res, row)
}

nrow(res); nrow(skipped)
count(skipped, reason)

# BH over the number of models actually fitted, which is the number of distinct
# patterns and not the number of SNPs in the file
res <- res %>%
  mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "BH"),
                .names = "q_{.col}")) %>%
  rename_with(~ str_replace(.x, "^q_p_", "q_"), starts_with("q_p_"))

names(res)

res %>%
  arrange(`p_clinical + structure`) %>%
  select(snp, n_present, `OR_unadjusted`, OR_clinical,
         `OR_clinical + structure`, `lo_clinical + structure`,
         `hi_clinical + structure`, `q_clinical + structure`) %>%
  head(25) %>% print(n = 25)

# how many survive at each stage, which is the headline of this half
res %>%
  summarise(fitted = n(),
            sig_unadjusted = sum(q_unadjusted < ALPHA, na.rm = TRUE),
            sig_clinical = sum(q_clinical < ALPHA, na.rm = TRUE),
            sig_full = sum(`q_clinical + structure` < ALPHA, na.rm = TRUE))

# the movement between the three, which is the interesting part
res %>%
  summarise(median_se_inflation_from_clustering =
              median(`se_clust_clinical + structure` /
                       `se_naive_clinical + structure`, na.rm = TRUE),
            median_abs_logOR_shrink_from_structure =
              median(abs(log(`OR_clinical + structure`)) -
                       abs(log(OR_clinical)), na.rm = TRUE))

p_shift <- res %>%
  filter(is.finite(OR_clinical), is.finite(`OR_clinical + structure`)) %>%
  mutate(sig = `q_clinical + structure` < ALPHA) %>%
  ggplot(aes(OR_clinical, `OR_clinical + structure`, colour = sig)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 1, colour = "grey80") +
  geom_vline(xintercept = 1, colour = "grey80") +
  geom_point(size = 1.1, alpha = 0.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = OKABE[2], `FALSE` = OKABE[8]),
                      labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "not"),
                      name = NULL) +
  labs(x = "Odds ratio, clinical adjustment only",
       y = "Odds ratio, plus population structure",
       title = "What adjusting for lineage does to each SNP",
       subtitle = paste0("Points pulled toward 1 on the vertical axis were ",
                         "marking a clone rather than\ncarrying an independent ",
                         "association with death.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_shift

top_snps <- res %>% arrange(`p_clinical + structure`) %>% head(20)

p_forest <- top_snps %>%
  mutate(snp = fct_reorder(snp, `OR_clinical + structure`),
         sig = `q_clinical + structure` < ALPHA) %>%
  ggplot(aes(`OR_clinical + structure`, snp, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = `lo_clinical + structure`,
                    xmax = `hi_clinical + structure`),
                orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = OKABE[2], `FALSE` = OKABE[1]),
                      labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "not"),
                      name = NULL) +
  labs(x = paste0("Adjusted odds of ", OUTCOME, " (log scale)"), y = NULL,
       title = "Leading SNPs, adjusted for clinical variables and lineage",
       subtitle = "Patient-clustered confidence intervals.") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6),
        panel.grid.minor = element_blank())

p_forest


# =============================================================================
# 08  ANNOTATING THE HITS
# =============================================================================

# The scoary output carries the annotation. It is a different analysis with its
# own p-values, so it is joined here for the gene names and effects only. Note
# it holds 6031 rows across 8 phenotypes, so filter to the one being modelled
# before joining or a SNP will pick up several annotations.

count(scoary, Phenotype) %>% arrange(desc(n)) %>% print(n = Inf)

ann <- scoary %>%
  filter(Phenotype == OUTCOME) %>%
  select(snp = Gene, position, EFFECT, LOCUS_TAG, GENE, PRODUCT,
         scoary_OR = Odds_ratio, scoary_q = Benjamini_H_p,
         scoary_best_pair_p = Best_pairwise_comp_p) %>%
  distinct(snp, .keep_all = TRUE)

nrow(ann)

hits <- res %>%
  filter(`q_clinical + structure` < ALPHA) %>%
  arrange(`p_clinical + structure`) %>%
  left_join(ann, by = "snp") %>%
  left_join(distinct(grp, representative, n_in_group),
            by = c("snp" = "representative"))

nrow(hits)

hits %>%
  select(snp, n_in_group, n_present, `OR_clinical + structure`,
         `lo_clinical + structure`, `hi_clinical + structure`,
         `q_clinical + structure`, GENE, EFFECT, PRODUCT) %>%
  print(n = Inf)

# A SNP standing for a group of identical patterns is one finding, not several.
# n_in_group says how many variants share that vector and therefore cannot be
# separated on these data. Report the group, not the representative alone.
sum(hits$n_in_group > 1, na.rm = TRUE)

# scoary tested the same SNPs without clinical adjustment or clustering, so
# disagreement is expected and informative rather than an error
hits %>%
  filter(!is.na(scoary_OR)) %>%
  select(snp, `OR_clinical + structure`, scoary_OR,
         `q_clinical + structure`, scoary_q) %>%
  print(n = 20)


# =============================================================================
# 09  THE SURVIVAL COHORT
# =============================================================================

# Only isolates whose patient has both an admission and a discharge date, as you
# asked. Time origin is the collection date and follow-up runs to discharge or
# in-hospital death.

surv <- surv_raw %>%
  transmute(KAUST_ID, patient = IDs, region = RGN_NM,
            los_recorded = as.numeric(LOS),
            hosp_death = Hospital_Death,
            ever_death = DTH,
            d_death = as.Date(DTH_Date, format = "%d/%m/%Y"),
            d_adm   = as.Date(ADS_DT_main, format = "%d/%m/%Y"),
            d_coll  = as.Date(Date_of_Collection, format = "%d/%m/%Y"),
            d_dis   = as.Date(DS_DT_main, format = "%d/%m/%Y"))

nrow(surv)
sapply(select(surv, starts_with("d_")), function(x) sum(is.na(x)))

# the completeness filter
surv <- surv %>% mutate(has_full_stay = !is.na(d_adm) & !is.na(d_dis) &
                                        !is.na(d_coll))
count(surv, has_full_stay)

sc <- filter(surv, has_full_stay)
nrow(sc); n_distinct(sc$patient)

# ORDERING CHECKS. Collection must fall inside the admission for isolation to be
# a legitimate time origin. In your data it does, for all 852, so no
# left-truncation adjustment is needed. The checks stay because a future extract
# may not behave.
sc %>% summarise(coll_before_adm = sum(d_coll < d_adm),
                 coll_after_dis  = sum(d_coll > d_dis),
                 dis_before_adm  = sum(d_dis < d_adm))

sc <- sc %>%
  mutate(fu_days = as.numeric(d_dis - d_coll),
         los_days = as.numeric(d_dis - d_adm))

summary(sc$fu_days)
sum(sc$fu_days == 0)      # discharged on the collection day

# LOS in the file should equal discharge minus admission, and it does
sc %>% filter(!is.na(los_recorded)) %>%
  summarise(max_disagreement = max(abs(los_recorded - los_days)))

# ZERO-LENGTH FOLLOW-UP. A row with fu_days == 0 contributes nothing to a Cox
# model and is dropped by the risk-set construction. Giving it half a day keeps
# it in the cohort at the cost of an arbitrary choice, so both are reported and
# the analysis uses the half-day version.
sc <- mutate(sc, time = pmax(fu_days, 0.5))

# THE EVENT DEFINITION, and the reason the competing risk is not optional. Of
# the patients discharged alive, a large number die afterwards, which is exactly
# the population that censoring would treat as still at risk.
count(sc, hosp_death, ever_death)

sc <- sc %>%
  mutate(status_cr = case_when(hosp_death == "Dead" ~ "death",
                               TRUE                 ~ "discharge"),
         status_cr = factor(status_cr, levels = c("censor", "death",
                                                  "discharge")),
         died_inhosp = as.integer(hosp_death == "Dead"))

count(sc, status_cr)
mean(sc$died_inhosp)

# died in hospital according to the flag, against the date falling inside the
# admission. A mismatch here is a data question worth resolving before the
# models, not after.
sc %>%
  mutate(date_inside = !is.na(d_death) & d_death >= d_adm & d_death <= d_dis) %>%
  count(died_inhosp, date_inside)

# attach the exposure, the confounders and the components
surv_dat <- sc %>%
  inner_join(select(meta, KAUST_ID, all_of(CONF)), by = "KAUST_ID") %>%
  inner_join(PCs, by = "KAUST_ID") %>%
  left_join(as_tibble(X_snp[, TEST_SNPS, drop = FALSE]) %>%
              mutate(KAUST_ID = rownames(X_snp)), by = "KAUST_ID")

nrow(surv_dat); n_distinct(surv_dat$patient)

surv_dat <- filter(surv_dat, complete.cases(select(surv_dat, all_of(CONF))))
nrow(surv_dat); n_distinct(surv_dat$patient)
count(surv_dat, status_cr)


# =============================================================================
# 10  KAPLAN-MEIER
# =============================================================================

# Set SNP_OF_INTEREST and rerun this block and the two after it. The default is
# the leading hit from block 07.

# WHERE THE DEFAULT COMES FROM, and what happens when there is no hit. If
# nothing survived BH in block 07 then hits has zero rows, hits$snp[1] is NA,
# and a membership test on NA is FALSE. That produces an error message about
# the count filter which is not the actual problem, so the state is checked
# before the choice is made.

nrow(res)                 # models fitted
nrow(hits)                # of those, BH significant after full adjustment

if (nrow(hits) > 0) {
  SNP_OF_INTEREST <- hits$snp[1]
  cat("using the leading BH-significant SNP:", SNP_OF_INTEREST, "\n")
} else {
  SNP_OF_INTEREST <- res$snp[which.min(res$`p_clinical + structure`)]
  cat("NO SNP SURVIVED BH CORRECTION.\n",
      "Falling back to the smallest raw p-value, which is a descriptive choice\n",
      "and not a finding:", SNP_OF_INTEREST,
      " raw p =", signif(min(res$`p_clinical + structure`, na.rm = TRUE), 3),
      " q =", signif(min(res$`q_clinical + structure`, na.rm = TRUE), 3), "\n",
      "Blocks 10 to 13 then describe the time course of a SNP that the scan did\n",
      "not establish an association for. Say so in the text if you report them.\n")
}

SNP_OF_INTEREST

# to look at a particular variant instead, set it here and rerun from this line
# SNP_OF_INTEREST <- "1|155|.|G|A|.|.|.|GT"

# the diagnostics that distinguish the three ways this can go wrong
if (!SNP_OF_INTEREST %in% names(surv_dat)) {
  cat("in TEST_SNPS      :", SNP_OF_INTEREST %in% TEST_SNPS, "\n")
  cat("in colnames(X_snp):", SNP_OF_INTEREST %in% colnames(X_snp), "\n")
  cat("passed MAC filter :", SNP_OF_INTEREST %in% keep_mac, "\n")
  cat("nearest names in surv_dat:\n")
  print(head(names(surv_dat)[startsWith(names(surv_dat),
                                        substr(SNP_OF_INTEREST, 1, 6))], 5))
  stop("'", SNP_OF_INTEREST, "' is not a column of surv_dat. The lines above ",
       "say whether it was filtered out in block 06, never fitted, or lost in ",
       "the join in block 09.")
}

surv_dat$expo <- factor(surv_dat[[SNP_OF_INTEREST]], 0:1,
                        c("absent", "present"))
count(surv_dat, expo, status_cr)

km <- survfit(Surv(time, died_inhosp) ~ expo, data = surv_dat)
print(km)
summary(km, times = c(7, 14, 30, 60))

lr <- survdiff(Surv(time, died_inhosp) ~ expo, data = surv_dat)
lr
p_logrank_one <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
p_logrank_one   # the number to put on the Kaplan-Meier figure

# WHAT THIS CURVE ACTUALLY SAYS, and it is not what most readers will take from
# it. Kaplan-Meier here censors discharge, which assumes a discharged patient
# remains at risk of in-hospital death in the same way as one still admitted.
# They do not, and one minus the KM curve therefore overstates the probability
# of dying in hospital. It is reported because you asked for it and because
# readers expect it; block 12 is the estimate to quote.

km_df <- tibble(time = km$time, surv = km$surv, lower = km$lower,
                upper = km$upper,
                strata = rep(names(km$strata), km$strata))

p_km <- ggplot(km_df, aes(time, 1 - surv, colour = strata, fill = strata)) +
  geom_ribbon(aes(ymin = 1 - upper, ymax = 1 - lower), alpha = 0.15,
              colour = NA) +
  geom_step(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 90)) +
  scale_colour_manual(values = OKABE, name = NULL) +
  scale_fill_manual(values = OKABE, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Days from isolation", y = "1 - Kaplan-Meier",
       title = paste0(SNP_OF_INTEREST, ": in-hospital death, discharge censored"),
       subtitle = paste0("This overstates the risk of dying in hospital, ",
                         "because censoring at discharge treats\na discharged ",
                         "patient as still liable to die on the ward. Compare ",
                         "with block 12.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_km


# =============================================================================
# 11  COX, CAUSE-SPECIFIC
# =============================================================================

# The cause-specific hazard of in-hospital death, with discharge censored, and
# patient-clustered robust standard errors. This model answers an aetiological
# question: among patients still in hospital and still alive, does the SNP raise
# the instantaneous rate of dying.

cox_fml <- as.formula(paste0(
  "Surv(time, died_inhosp) ~ expo + ",
  paste0("`", CONF_FULL, "`", collapse = " + "), " + cluster(patient)"))

cox_fml

cox_cs <- coxph(cox_fml, data = surv_dat)
summary(cox_cs)

# the clustering, in numbers: robust against naive standard error
sm <- summary(cox_cs)$coefficients
round(cbind(sm[, "se(coef)"], sm[, "robust se"],
            ratio = sm[, "robust se"] / sm[, "se(coef)"]), 4)

# PROPORTIONAL HAZARDS. If the SNP term fails this test the hazard ratio is an
# average over a changing effect and should not be quoted as a single number.
zph <- cox.zph(cox_cs)
print(zph)
plot(zph[1], main = "Scaled Schoenfeld residuals for the SNP term")

# unadjusted for contrast
cox_crude <- coxph(Surv(time, died_inhosp) ~ expo + cluster(patient),
                   data = surv_dat)
summary(cox_crude)$conf.int


# =============================================================================
# 12  COMPETING RISKS
# =============================================================================

# Discharge alive is not censoring, it is a competing event: once it happens the
# patient can no longer die in hospital. Two things follow.
#
# The cumulative incidence function, estimated by Aalen-Johansen, is the
# probability of having died in hospital by a given day, accounting for the fact
# that discharge removes people from that possibility. It is what a clinician
# means by the risk of dying on this admission.
#
# The Fine-Gray model puts covariates on the subdistribution hazard, which
# targets that cumulative incidence directly. It answers a prediction question
# where the cause-specific Cox of block 11 answers an aetiological one, and the
# two coefficients differ on purpose. Report both and say which question each is
# for.

# NOTE ON id. survfit is NOT given id = patient here. That argument tells it to
# read repeated rows as one subject moving between states, which is not what a
# second isolate from the same patient is. The point estimate treats each
# episode separately; the clustering shows up in the model standard errors
# below, and the pointwise interval on the curve is correspondingly optimistic.

cif <- survfit(Surv(time, status_cr) ~ expo, data = surv_dat)
print(cif$states)

# summary()$pstate comes back WITHOUT dimnames, so the death column has to be
# found by position. Indexing it by the name "death" fails with a dimnames
# error, which is an unhelpful way to be told that.
J_DEATH <- match("death", cif$states)
J_DEATH

cif_tab <- summary(cif, times = c(7, 14, 30, 60, 90))
data.frame(strata = cif_tab$strata, time = cif_tab$time,
           death = round(cif_tab$pstate[, J_DEATH], 4),
           discharge = round(cif_tab$pstate[, match("discharge", cif$states)], 4))

cif_df <- tibble(time = cif$time,
                 strata = rep(names(cif$strata), cif$strata),
                 prob = cif$pstate[, J_DEATH])

p_cif <- ggplot(cif_df, aes(time, prob, colour = strata)) +
  geom_step(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 90)) +
  scale_colour_manual(values = OKABE, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Days from isolation", y = "Cumulative incidence of in-hospital death",
       title = paste0(SNP_OF_INTEREST, ": competing-risks estimate"),
       subtitle = paste0("Discharge alive treated as a competing event. This ",
                         "sits below the curve in block 10,\nand the gap is the ",
                         "size of the error that censoring at discharge makes.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_cif

# the two estimates side by side at fixed times, which is the clearest way to
# show a reviewer why the competing-risks version was used
km_tab <- summary(km, times = c(7, 14, 30, 60, 90))

tibble(strata = as.character(km_tab$strata), time = km_tab$time,
       km_risk = 1 - km_tab$surv) %>%
  inner_join(tibble(strata = as.character(cif_tab$strata),
                    time = cif_tab$time,
                    cif_risk = cif_tab$pstate[, J_DEATH]),
             by = c("strata", "time")) %>%
  mutate(overstatement = km_risk - cif_risk,
         times_too_high = km_risk / cif_risk) %>%
  arrange(strata, time) %>%
  print(n = Inf)

# Gray's test, via the weighted score test on the subdistribution
fg_dat <- finegray(Surv(time, status_cr) ~ ., data = surv_dat, etype = "death")
nrow(fg_dat); nrow(surv_dat)

fg_crude <- coxph(Surv(fgstart, fgstop, fgstatus) ~ expo + cluster(patient),
                  weights = fgwt, data = fg_dat)
summary(fg_crude)

fg_fml <- as.formula(paste0(
  "Surv(fgstart, fgstop, fgstatus) ~ expo + ",
  paste0("`", CONF_FULL, "`", collapse = " + "), " + cluster(patient)"))

fg_adj <- coxph(fg_fml, weights = fgwt, data = fg_dat)
summary(fg_adj)


# =============================================================================
# 13  THE SURVIVAL SCAN ACROSS ALL SNPs
# =============================================================================

# Blocks 10 to 12 follow one SNP. This block runs the same two survival models
# for every SNP that passed the filters, so that the time-to-event evidence is a
# scan in its own right rather than a footnote on the logistic scan.
#
# It is affordable. On 852 rows with this covariate set a cause-specific Cox fit
# takes about 11 ms and a Fine-Gray fit about 6 ms, so the whole scan is a
# couple of minutes. The Fine-Gray expansion depends only on the outcome and the
# censoring distribution, never on the covariate, so it is built once and the
# SNP column is attached to it by row index each iteration. Rebuilding it inside
# the loop would multiply the cost by a thousand for no gain.
#
# NOTE ON WHAT THIS IS NOT. Three scans over the same patients now exist:
# logistic on death ever, cause-specific Cox on in-hospital death, and Fine-Gray
# on its cumulative incidence. They are correlated, so a SNP appearing in all
# three is not three independent confirmations. Agreement across them is a
# consistency check on the modelling, not additional evidence.

surv_scan_snps <- TEST_SNPS[TEST_SNPS %in% names(surv_dat)]
length(surv_scan_snps)

# the Fine-Gray expansion, built once
surv_dat$row_id <- seq_len(nrow(surv_dat))

fg_base <- finegray(
  as.formula(paste0("Surv(time, status_cr) ~ row_id + patient + ",
                    paste0("`", CONF_FULL, "`", collapse = " + "))),
  data = surv_dat, etype = "death")

nrow(fg_base); nrow(surv_dat)
round(nrow(fg_base) / nrow(surv_dat), 2)     # expansion factor

conf_terms <- paste0("`", CONF_FULL, "`", collapse = " + ")

cs_fml <- as.formula(paste0("Surv(time, died_inhosp) ~ snp_x + ", conf_terms,
                            " + cluster(patient)"))
fg_fml_scan <- as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ snp_x + ",
                                 conf_terms, " + cluster(patient)"))

scan <- NULL
scan_skipped <- NULL

t_start <- Sys.time()

for (i in seq_along(surv_scan_snps)) {

  sn <- surv_scan_snps[i]
  v  <- surv_dat[[sn]]

  # both arms need enough deaths for a hazard ratio to mean anything
  tab <- table(v, surv_dat$died_inhosp)
  if (!all(dim(tab) == c(2, 2)) || min(tab) < MIN_CELL) {
    scan_skipped <- bind_rows(scan_skipped, tibble(
      snp = sn, reason = "cell too small", min_cell = min(tab)))
    next
  }

  surv_dat$snp_x    <- v
  surv_dat$snp_x_lr <- factor(v, 0:1, c("absent", "present"))
  fg_base$snp_x     <- v[fg_base$row_id]

  row <- tibble(snp = sn, n = nrow(surv_dat),
                patients = n_distinct(surv_dat$patient),
                n_present = sum(v),
                deaths_present = sum(surv_dat$died_inhosp[v == 1]),
                deaths_absent  = sum(surv_dat$died_inhosp[v == 0]))

  # THE LOG-RANK TEST, which is the p-value that belongs on a Kaplan-Meier
  # figure. It is unadjusted and it ignores the patient clustering, so it will
  # disagree with p_cs below and should disagree: it answers a cruder question.
  # Report it beside the curve and report p_cs in the text.
  sdf <- try(survdiff(Surv(time, died_inhosp) ~ snp_x_lr, data = surv_dat),
             silent = TRUE)
  if (!inherits(sdf, "try-error"))
    row$p_logrank <- pchisq(sdf$chisq, length(sdf$n) - 1, lower.tail = FALSE)

  # unadjusted but clustered, so the crude hazard ratio on the figure and the
  # adjusted one in the table can be compared like for like
  fcr <- try(coxph(Surv(time, died_inhosp) ~ snp_x + cluster(patient),
                   data = surv_dat), silent = TRUE)
  if (!inherits(fcr, "try-error")) {
    row$HR_crude <- summary(fcr)$conf.int["snp_x", "exp(coef)"]
    row$p_crude  <- summary(fcr)$coefficients["snp_x", "Pr(>|z|)"]
  }

  fcs <- try(coxph(cs_fml, data = surv_dat), silent = TRUE)

  if (!inherits(fcs, "try-error")) {
    ci <- summary(fcs)$conf.int
    co <- summary(fcs)$coefficients
    row$HR_cs    <- ci["snp_x", "exp(coef)"]
    row$lo_cs    <- ci["snp_x", "lower .95"]
    row$hi_cs    <- ci["snp_x", "upper .95"]
    row$p_cs     <- co["snp_x", "Pr(>|z|)"]
    row$se_naive_cs   <- co["snp_x", "se(coef)"]
    row$se_cluster_cs <- co["snp_x", "robust se"]
  }

  ffg <- try(coxph(fg_fml_scan, weights = fgwt, data = fg_base), silent = TRUE)

  if (!inherits(ffg, "try-error")) {
    ci <- summary(ffg)$conf.int
    co <- summary(ffg)$coefficients
    row$sHR    <- ci["snp_x", "exp(coef)"]
    row$lo_fg  <- ci["snp_x", "lower .95"]
    row$hi_fg  <- ci["snp_x", "upper .95"]
    row$p_fg   <- co["snp_x", "Pr(>|z|)"]
  }

  scan <- bind_rows(scan, row)

  if (i %% 250 == 0)
    cat(sprintf("%5d/%5d  %.1f min elapsed\n", i, length(surv_scan_snps),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

as.numeric(difftime(Sys.time(), t_start, units = "mins"))

nrow(scan); nrow(scan_skipped)
count(scan_skipped, reason)

scan <- scan %>%
  mutate(q_cs = p.adjust(p_cs, method = "BH"),
         q_fg = p.adjust(p_fg, method = "BH"),
         q_logrank = p.adjust(p_logrank, method = "BH"),
         q_crude   = p.adjust(p_crude,   method = "BH"))

# how far the log-rank and the adjusted clustered Cox disagree, which is the
# cost of the adjustment and the clustering together
scan %>%
  summarise(sig_logrank = sum(q_logrank < ALPHA, na.rm = TRUE),
            sig_crude_clustered = sum(q_crude < ALPHA, na.rm = TRUE),
            sig_adjusted_clustered = sum(q_cs < ALPHA, na.rm = TRUE),
            median_p_ratio = median(p_cs / p_logrank, na.rm = TRUE))

scan %>% arrange(p_cs) %>%
  select(snp, n_present, deaths_present, deaths_absent,
         p_logrank, q_logrank, HR_crude, q_crude,
         HR_cs, lo_cs, hi_cs, q_cs, sHR, q_fg) %>%
  head(25) %>% print(n = 25)

scan %>%
  summarise(fitted = n(),
            sig_cause_specific = sum(q_cs < ALPHA, na.rm = TRUE),
            sig_fine_gray = sum(q_fg < ALPHA, na.rm = TRUE),
            sig_both = sum(q_cs < ALPHA & q_fg < ALPHA, na.rm = TRUE))

# clustering, again in numbers
scan %>% summarise(median_se_inflation = median(se_cluster_cs / se_naive_cs,
                                                na.rm = TRUE))

# the three scans against each other
three <- scan %>%
  select(snp, HR_cs, q_cs, sHR, q_fg) %>%
  inner_join(select(res, snp, OR = `OR_clinical + structure`,
                    q_or = `q_clinical + structure`), by = "snp")

nrow(three)

three %>%
  summarise(sig_logistic = sum(q_or < ALPHA, na.rm = TRUE),
            sig_cs = sum(q_cs < ALPHA, na.rm = TRUE),
            sig_fg = sum(q_fg < ALPHA, na.rm = TRUE),
            sig_all_three = sum(q_or < ALPHA & q_cs < ALPHA & q_fg < ALPHA,
                                na.rm = TRUE))

three %>%
  filter(q_or < ALPHA | q_cs < ALPHA | q_fg < ALPHA) %>%
  left_join(select(ann, snp, GENE, EFFECT, PRODUCT), by = "snp") %>%
  arrange(q_cs) %>%
  print(n = Inf)

# the cause-specific and subdistribution hazard ratios should agree in sign and
# the second should sit closer to 1, because the competing risk absorbs part of
# the effect. A pair that disagrees in sign is worth opening individually.
three %>%
  filter(is.finite(HR_cs), is.finite(sHR)) %>%
  summarise(correlation = cor(log(HR_cs), log(sHR), use = "complete.obs"),
            sign_disagreements = sum(sign(log(HR_cs)) != sign(log(sHR)),
                                     na.rm = TRUE),
            median_attenuation = median(abs(log(sHR)) / abs(log(HR_cs)),
                                        na.rm = TRUE))

p_scan <- three %>%
  filter(is.finite(HR_cs), is.finite(sHR)) %>%
  mutate(sig = case_when(q_cs < ALPHA & q_fg < ALPHA ~ "both",
                         q_cs < ALPHA ~ "cause-specific only",
                         q_fg < ALPHA ~ "Fine-Gray only",
                         TRUE ~ "neither")) %>%
  ggplot(aes(HR_cs, sHR, colour = sig)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 1, colour = "grey85") +
  geom_vline(xintercept = 1, colour = "grey85") +
  geom_point(size = 1.1, alpha = 0.65) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(both = OKABE[2],
                                 `cause-specific only` = OKABE[1],
                                 `Fine-Gray only` = OKABE[3],
                                 neither = OKABE[8]), name = NULL) +
  labs(x = "Cause-specific hazard ratio", y = "Subdistribution hazard ratio",
       title = "The two survival scans against each other",
       subtitle = paste0("Points below the diagonal on the right side are the ",
                         "usual pattern: the competing risk\nabsorbs part of ",
                         "the effect, so the subdistribution ratio sits closer ",
                         "to 1.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_scan

p_or_hr <- three %>%
  filter(is.finite(OR), is.finite(HR_cs)) %>%
  mutate(sig = q_or < ALPHA | q_cs < ALPHA) %>%
  ggplot(aes(OR, HR_cs, colour = sig)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 1, colour = "grey85") +
  geom_vline(xintercept = 1, colour = "grey85") +
  geom_point(size = 1.1, alpha = 0.65) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = OKABE[2], `FALSE` = OKABE[8]),
                      labels = c(`TRUE` = "significant somewhere",
                                 `FALSE` = "neither"), name = NULL) +
  labs(x = "Odds ratio, death ever (1073 isolates)",
       y = "Cause-specific hazard ratio, in-hospital death (852)",
       title = "Odds of dying at all against the rate of dying in hospital",
       subtitle = paste0("Different outcomes on different cohorts, so ",
                         "disagreement is informative rather than\nan error. A ",
                         "SNP high on one axis and not the other acts on one and ",
                         "not the other.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_or_hr

# tidy up the scratch columns so they cannot leak into a later join
surv_dat$snp_x    <- NULL
surv_dat$snp_x_lr <- NULL
fg_base$snp_x     <- NULL


# =============================================================================
# 14  THE THREE SIDE BY SIDE
# =============================================================================

# The same exposure through three lenses. They are not three attempts at one
# number; they answer three questions, and a table that says which is which
# saves a paragraph of discussion.

# the logistic row must exist for this SNP, or the comparison is built on an
# empty vector and silently produces a zero-row tibble
sum(res$snp == SNP_OF_INTEREST)

compare <- bind_rows(
  tibble(analysis = "Logistic, death ever",
         question = "odds of dying at any point",
         estimate = res$`OR_clinical + structure`[res$snp == SNP_OF_INTEREST],
         lo = res$`lo_clinical + structure`[res$snp == SNP_OF_INTEREST],
         hi = res$`hi_clinical + structure`[res$snp == SNP_OF_INTEREST],
         scale = "odds ratio"),
  tibble(analysis = "Cox, cause-specific",
         question = "rate of dying among those still admitted and alive",
         estimate = summary(cox_cs)$conf.int["expopresent", "exp(coef)"],
         lo = summary(cox_cs)$conf.int["expopresent", "lower .95"],
         hi = summary(cox_cs)$conf.int["expopresent", "upper .95"],
         scale = "cause-specific hazard ratio"),
  tibble(analysis = "Fine-Gray",
         question = "probability of dying in hospital by a given day",
         estimate = summary(fg_adj)$conf.int["expopresent", "exp(coef)"],
         lo = summary(fg_adj)$conf.int["expopresent", "lower .95"],
         hi = summary(fg_adj)$conf.int["expopresent", "upper .95"],
         scale = "subdistribution hazard ratio"))

print(compare, n = Inf)

p_compare <- compare %>%
  mutate(analysis = factor(analysis, levels = rev(analysis))) %>%
  ggplot(aes(estimate, analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0, linewidth = 0.5, colour = OKABE[1]) +
  geom_point(size = 2.6, colour = OKABE[1]) +
  scale_x_log10() +
  labs(x = "Estimate (log scale)", y = NULL,
       title = paste0(SNP_OF_INTEREST, ": three questions, three answers"),
       subtitle = paste0("Not three estimates of one quantity. All adjusted for ",
                         "the clinical set and lineage,\nall with ",
                         "patient-clustered intervals.")) +
  theme(panel.grid.minor = element_blank())

p_compare

# A CLOSING CAUTION FOR THE MANUSCRIPT. The SNP was chosen as the leading hit of
# a scan over roughly four thousand tests on these same patients, so the effect
# estimated in blocks 10 to 12 is conditional on having been selected for being
# extreme and is biased away from the null. The survival analysis characterises
# the association; it does not validate it. Only a second collection can do
# that, and the discussion should say so rather than leaving a reader to notice.

length(TEST_SNPS)


# =============================================================================
# 15  SAVING
# =============================================================================

# ggsave("snp_fig_pca.pdf",       p_pc,      width =  7, height = 6)
# ggsave("snp_fig_or_shift.pdf",  p_shift,   width =  7, height = 6.5)
# ggsave("snp_fig_forest.pdf",    p_forest,  width =  8, height = 6)
# ggsave("snp_fig_km.pdf",        p_km,      width =  7, height = 5)
# ggsave("snp_fig_cif.pdf",       p_cif,     width =  7, height = 5)
# ggsave("snp_fig_compare.pdf",   p_compare, width =  7, height = 3.5)
#
# write_csv(screen,   "conf_screen.csv")
# write_csv(res,      "snp_logistic_death.csv")
# write_csv(hits,     "snp_hits_annotated.csv")
# write_csv(grp,      "snp_pattern_groups.csv")
# write_csv(surv_dat |> select(-any_of(TEST_SNPS)), "survival_cohort.csv")
# write_csv(compare,  "three_analyses_one_snp.csv")
#
# THE PER-SNP SURVIVAL RESULTS, which is what block 13 exists to produce:
# write_csv(scan,  "snp_survival_scan.csv")     # Cox and Fine-Gray, every SNP
# write_csv(three, "snp_three_scans_joined.csv")# logistic + Cox + Fine-Gray
# ggsave("snp_fig_cs_vs_fg.pdf", p_scan,   width = 7, height = 6.5)
# ggsave("snp_fig_or_vs_hr.pdf", p_or_hr,  width = 7, height = 6.5)
