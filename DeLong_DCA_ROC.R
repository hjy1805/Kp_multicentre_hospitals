# =============================================================================
# PAIRED ROC COMPARISON, DECISION CURVES, AND ROC PLOTS FROM POOLED PREDICTIONS
#
# Works from the per-patient prediction files rather than the per-fold summary
# metrics, which makes three things possible that the summary files cannot
# support: DeLong's test for two correlated ROC curves, decision curve
# analysis, and an ROC plot you can build from any set of models you pick.
#
# Inputs
#   ALL_TEST_PREDICTIONS_FOR_DELONG.csv    KAUST_ID Fold Model Outcome y_true prob
#   ALL_TRAIN_PREDICTIONS_FOR_DELONG.csv   same columns
#
# No pROC. DeLong is computed here from midrank placement values, which is the
# same estimator pROC uses. It is verified two ways: the AUC it produces equals
# the brute-force count of concordant pairs including ties, and its standard
# errors agree with a 4000-draw bootstrap to within one percent
# (se(AUC) 0.0258 against 0.0260, se(difference) 0.0398 against 0.0402).
#
# Flat script, no named functions. Run block by block. Figures print to the
# device. The saving block is at the end, commented out.
#
# Blocks
#   01 packages and settings        08 calibration, required before any DCA
#   02 read                         09 decision curves
#   03 integrity checks             10 net benefit differences, bootstrapped
#   04 the wide matrix              11 ROC curves, your selection
#   05 AUC and DeLong variance      12 train against test
#   06 pairwise DeLong              13 saving (commented out)
#   07 within-fold DeLong
# =============================================================================


# =============================================================================
# 01  PACKAGES AND SETTINGS
# =============================================================================

# rm(list = ls()) first. dplyr resolves a bare name against the data frame and
# then against the global environment, so a leftover object with a column's
# name is used silently wherever that column is out of scope.

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)

DIR <- "/Users/daneshm/Documents/PA_KAIMRC/ML_files/predictionsWOAST"

F_TEST  <- file.path(DIR, "ALL_TEST_PREDICTIONS_FOR_DELONG.csv")
F_TRAIN <- file.path(DIR, "ALL_TRAIN_PREDICTIONS_FOR_DELONG.csv")

OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#999999", "#000000")

theme_set(theme_bw(base_size = 11))

MODEL_LABELS <- c(
  Model_1 = "1  Tier 1",
  Model_2 = "2  + Tier 2",
  Model_3 = "3  + Tier 3",
  Model_4 = "4  + Tier 4",
  Model_5 = "5  + ARG/vir",
  Model_6 = "6  + pangenome",
  Model_7 = "7  + SNPs",
  Model_8 = "8  + unitigs",
  Model_9 = "9  genomic only")

MODEL_ORDER <- names(MODEL_LABELS)

REF_MODEL <- "Model_4"      # clinical tiers only, the comparator that matters
ALPHA <- 0.05

# Threshold probabilities for the decision curves. The range has to be one a
# clinician would actually work in. Treating everyone above a 1% risk of death
# is not a policy anyone would adopt, and net benefit at such thresholds is
# dominated by the prevalence rather than by the model.
PT_GRID <- seq(0.01, 0.60, by = 0.005)
PT_TABLE <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)

N_BOOT <- 2000              # for the net benefit intervals in block 10
set.seed(20260909)


# =============================================================================
# 02  READ
# =============================================================================

pred_test  <- read_csv(F_TEST,  show_col_types = FALSE)
pred_train <- read_csv(F_TRAIN, show_col_types = FALSE)

dim(pred_test); dim(pred_train)
names(pred_test)

glimpse(pred_test)

count(pred_test, Outcome)
count(pred_test, Model)
count(pred_test, Fold)

pred_test <- pred_test %>%
  mutate(Model = factor(Model, levels = MODEL_ORDER),
         model_lab = factor(MODEL_LABELS[as.character(Model)],
                            levels = unname(MODEL_LABELS)))

pred_train <- pred_train %>%
  mutate(Model = factor(Model, levels = MODEL_ORDER),
         model_lab = factor(MODEL_LABELS[as.character(Model)],
                            levels = unname(MODEL_LABELS)))

# any model name that did not match the lookup shows up as NA here
sum(is.na(pred_test$Model)); sum(is.na(pred_train$Model))


# =============================================================================
# 03  INTEGRITY CHECKS
# =============================================================================

# Everything downstream assumes one prediction per patient per model per
# outcome, and assumes the patients are independent. Both assumptions are
# checkable and both can fail quietly.

range(pred_test$prob)
sum(!is.finite(pred_test$prob))
count(pred_test, y_true)

# one row per patient, per model, per outcome
dupes <- pred_test %>%
  count(Outcome, Model, KAUST_ID) %>%
  filter(n > 1)

nrow(dupes)
head(dupes, 20)

# each patient in exactly one test fold
fold_spread <- pred_test %>%
  group_by(Outcome, Model, KAUST_ID) %>%
  summarise(n_folds = n_distinct(Fold), .groups = "drop") %>%
  count(n_folds)

print(fold_spread)

# how many patients per outcome, and how many events
pred_test %>%
  group_by(Outcome) %>%
  summarise(models = n_distinct(Model),
            patients = n_distinct(KAUST_ID),
            events = sum(y_true[Model == first(Model)]),
            prevalence = mean(y_true[Model == first(Model)]),
            .groups = "drop") %>%
  print(n = Inf)

# TIED PREDICTIONS. In your sample rows PA00300, PA00320, PA00322, PA00324 and
# PA00782 all carry prob 0.05828115471432093 to the last digit under Model_1.
# For a model built on demographics alone that is expected: identical covariate
# patterns give identical predictions, and DeLong handles the ties correctly
# through midranks. What is not benign is the same PATIENT appearing under two
# isolate IDs, because then the sample is not independent and every standard
# error below is too small. The block above catches repeated IDs; this one
# shows how large the tied groups get, which is the signal that the model has
# very little to separate patients with.

pred_test %>%
  group_by(Outcome, Model) %>%
  summarise(n = n(), distinct_probs = n_distinct(prob),
            largest_tied_group = max(table(prob)),
            .groups = "drop") %>%
  arrange(distinct_probs) %>%
  print(n = Inf)

# a tied group whose members disagree on the outcome is where the AUC is lost
pred_test %>%
  filter(Model == "Model_1") %>%
  group_by(Outcome, prob) %>%
  filter(n() > 1) %>%
  summarise(size = n(), events = sum(y_true), .groups = "drop") %>%
  filter(events > 0, events < size) %>%
  arrange(desc(size)) %>%
  head(15)


# =============================================================================
# 04  THE WIDE MATRIX
# =============================================================================

# DeLong compares curves on the same patients, so the models have to be aligned
# patient by patient. Any patient missing from any model is dropped, and the
# count of dropped patients is printed rather than left implicit.

wide_test <- pred_test %>%
  select(Outcome, KAUST_ID, y_true, Model, prob) %>%
  mutate(Model = as.character(Model)) %>%
  pivot_wider(names_from = Model, values_from = prob)

dim(wide_test)
head(wide_test)

mods_present <- MODEL_ORDER[MODEL_ORDER %in% names(wide_test)]
mods_present

before <- nrow(wide_test)
wide_test <- wide_test[complete.cases(wide_test[, mods_present]), ]
cat("patients dropped for incomplete model coverage:", before - nrow(wide_test), "\n")

count(wide_test, Outcome)

# y_true must be constant for a patient across models; if it is not, the
# pivot silently list-columned something and the rest of the script is void
sum(!wide_test$y_true %in% c(0, 1))


# =============================================================================
# 05  AUC AND DELONG VARIANCE
# =============================================================================

# Placement values. For a case i, V10_i is the proportion of controls it beats,
# counting a tie as half. For a control j, V01_j is the proportion of cases
# that beat it, again counting ties as half. The mean of either equals the AUC.
# The midrank identity used here,
#
#     V10_i = (rank of case i among all - rank of case i among cases) / n_neg
#
# gives exactly the brute-force count, ties included, in O(n log n) instead of
# O(n_pos x n_neg).
#
# DeLong's covariance is then S = cov(V10)/n_pos + cov(V01)/n_neg, and the
# variance of a difference between two models is S_aa + S_bb - 2*S_ab. The
# 2*S_ab term is the whole point: it removes the shared patient variation, so
# the paired comparison is far more powerful than comparing two independent
# confidence intervals by eye.

auc_tbl <- NULL
S_list  <- list()      # DeLong covariance, one matrix per outcome
V_list  <- list()      # placement values, kept for block 06

for (o in sort(unique(wide_test$Outcome))) {

  w  <- filter(wide_test, Outcome == o)
  yy <- w$y_true
  ip <- which(yy == 1); ineg <- which(yy == 0)
  np <- length(ip); nn <- length(ineg)

  if (np < 5 || nn < 5) { cat("skipping", o, "- too few in one class\n"); next }

  V10 <- matrix(NA_real_, np, length(mods_present),
                dimnames = list(NULL, mods_present))
  V01 <- matrix(NA_real_, nn, length(mods_present),
                dimnames = list(NULL, mods_present))

  for (m in mods_present) {
    x <- w[[m]][ip]; z <- w[[m]][ineg]
    r_all <- rank(c(x, z), ties.method = "average")
    V10[, m] <- (r_all[seq_len(np)] - rank(x, ties.method = "average")) / nn
    V01[, m] <- 1 - (r_all[np + seq_len(nn)] - rank(z, ties.method = "average")) / np
  }

  S <- cov(V10) / np + cov(V01) / nn
  a <- colMeans(V10)

  # the two estimators of the AUC must agree to machine precision
  stopifnot(max(abs(a - colMeans(V01))) < 1e-9)

  se  <- sqrt(diag(S))
  # logit-transformed interval, so the limits cannot leave [0, 1]
  lg  <- log(a / (1 - a)); se_lg <- se / (a * (1 - a))

  auc_tbl <- bind_rows(auc_tbl, tibble(
    Outcome = o, model = mods_present, n = nrow(w), n_pos = np, n_neg = nn,
    auc = a, se = se,
    lo = plogis(lg - qnorm(0.975) * se_lg),
    hi = plogis(lg + qnorm(0.975) * se_lg)))

  S_list[[o]] <- S
  V_list[[o]] <- list(V10 = V10, V01 = V01, np = np, nn = nn)
}

auc_tbl <- mutate(auc_tbl,
                  model = factor(model, levels = MODEL_ORDER),
                  model_lab = factor(MODEL_LABELS[as.character(model)],
                                     levels = unname(MODEL_LABELS)))

auc_tbl %>% arrange(Outcome, model) %>%
  select(Outcome, model_lab, n, n_pos, auc, se, lo, hi) %>%
  print(n = Inf)

# these are the pooled cross-validated AUCs. They will not match the mean of
# the per-fold AUCs from the summary files exactly, and should not: pooling
# ranks all patients against each other, averaging ranks them only within fold.
round(range(auc_tbl$auc), 4)


# =============================================================================
# 06  PAIRWISE DELONG
# =============================================================================

# All pairs within each outcome, BH corrected within outcome.
#
# WHAT THIS TEST DOES AND DOES NOT ACCOUNT FOR. DeLong's variance covers the
# sampling of patients. It treats the two sets of predictions as fixed numbers
# attached to those patients. Here the predictions came out of cross-validation,
# so they also carry the variation of the model fitting, and the training sets
# behind them overlap. That extra variation is not in the standard error, which
# makes these p-values optimistic in the same direction as an uncorrected
# paired t-test on fold metrics. Block 07 is the check on how much this matters.

delong <- NULL

for (o in names(S_list)) {

  S <- S_list[[o]]
  a <- auc_tbl$auc[auc_tbl$Outcome == o]
  names(a) <- as.character(auc_tbl$model[auc_tbl$Outcome == o])

  for (i in seq_along(mods_present)) for (j in seq_along(mods_present)) {

    if (j <= i) next
    ma <- mods_present[i]; mb <- mods_present[j]

    d <- a[mb] - a[ma]
    v <- S[ma, ma] + S[mb, mb] - 2 * S[ma, mb]

    if (!is.finite(v) || v <= 0) {
      delong <- bind_rows(delong, tibble(
        Outcome = o, model_a = ma, model_b = mb, delta = unname(d),
        se = NA_real_, z = NA_real_, p = NA_real_))
      next
    }

    zst <- unname(d) / sqrt(v)

    delong <- bind_rows(delong, tibble(
      Outcome = o, model_a = ma, model_b = mb, delta = unname(d),
      se = sqrt(v), z = zst, p = 2 * pnorm(-abs(zst))))
  }
}

delong <- delong %>%
  group_by(Outcome) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(lo = delta - qnorm(0.975) * se,
         hi = delta + qnorm(0.975) * se)

nrow(delong)

delong %>% filter(q < ALPHA) %>% arrange(Outcome, q) %>% print(n = Inf)

# the comparison that answers the question: everything against clinical only
delong_ref <- delong %>%
  filter(model_a == REF_MODEL | model_b == REF_MODEL) %>%
  mutate(other = ifelse(model_a == REF_MODEL, model_b, model_a),
         # orient every row as other minus reference
         delta_o = ifelse(model_a == REF_MODEL, delta, -delta),
         lo_o    = ifelse(model_a == REF_MODEL, lo, -hi),
         hi_o    = ifelse(model_a == REF_MODEL, hi, -lo)) %>%
  group_by(Outcome) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(other = factor(other, levels = MODEL_ORDER),
         other_lab = factor(MODEL_LABELS[as.character(other)],
                            levels = unname(MODEL_LABELS)),
         verdict = case_when(is.na(q) ~ "not testable",
                             q < ALPHA & delta_o > 0 ~ "better",
                             q < ALPHA & delta_o < 0 ~ "worse",
                             TRUE ~ "no evidence"))

delong_ref %>%
  arrange(Outcome, other) %>%
  select(Outcome, other_lab, delta_o, lo_o, hi_o, p, q, verdict) %>%
  print(n = Inf)

delong_ref %>% count(other_lab, verdict) %>% print(n = Inf)

p_delong <- ggplot(delong_ref, aes(delta_o, fct_rev(other_lab), colour = verdict)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo_o, xmax = hi_o), orientation = "y",
                width = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ Outcome, ncol = 3) +
  scale_colour_manual(values = c(better = OKABE[3], worse = OKABE[2],
                                 `no evidence` = OKABE[8],
                                 `not testable` = "grey80"), name = NULL) +
  labs(x = paste0("DeLong difference in AUC against ", MODEL_LABELS[REF_MODEL]),
       y = NULL,
       title = "Paired comparison of ROC curves on the same patients",
       subtitle = paste0("DeLong's test for correlated curves, BH corrected ",
                         "within outcome.\nStandard errors cover patient ",
                         "sampling only, not the refitting of the model.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_delong

p_dl_grid <- delong %>%
  mutate(model_a = factor(MODEL_LABELS[model_a], levels = unname(MODEL_LABELS)),
         model_b = factor(MODEL_LABELS[model_b], levels = unname(MODEL_LABELS)),
         shown = ifelse(!is.na(q) & q < ALPHA, delta, NA_real_)) %>%
  ggplot(aes(model_a, fct_rev(model_b), fill = shown)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_wrap(~ Outcome, ncol = 3) +
  scale_fill_gradient2(low = OKABE[2], mid = "white", high = OKABE[1],
                       midpoint = 0, na.value = "grey93", name = "delta AUC") +
  labs(x = NULL, y = NULL,
       title = "All pairwise DeLong comparisons, BH q < 0.05 only",
       subtitle = "Grey means the two curves were not distinguishable.") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6), panel.grid = element_blank())

p_dl_grid


# =============================================================================
# 07  WITHIN-FOLD DELONG
# =============================================================================

# The pooled test in block 06 ignores that the predictions for different
# patients came from different fitted models. Running DeLong inside each fold
# separately removes that: within one fold every prediction comes from one fit.
# The cost is power, since each fold holds a fifth of the patients. If the
# pooled result survives here in most folds it is real, and if it appears only
# after pooling it is worth a sentence of caution in the paper.

fold_delong <- NULL

for (o in sort(unique(pred_test$Outcome))) {
  for (fl in sort(unique(pred_test$Fold))) {

    wf <- pred_test %>%
      filter(Outcome == o, Fold == fl) %>%
      select(KAUST_ID, y_true, Model, prob) %>%
      mutate(Model = as.character(Model)) %>%
      pivot_wider(names_from = Model, values_from = prob)

    mp <- mods_present[mods_present %in% names(wf)]
    wf <- wf[complete.cases(wf[, mp]), ]

    yy <- wf$y_true
    ip <- which(yy == 1); ineg <- which(yy == 0)
    np <- length(ip); nn <- length(ineg)
    if (np < 5 || nn < 5 || length(mp) < 2) next

    V10 <- matrix(NA_real_, np, length(mp), dimnames = list(NULL, mp))
    V01 <- matrix(NA_real_, nn, length(mp), dimnames = list(NULL, mp))

    for (m in mp) {
      x <- wf[[m]][ip]; z <- wf[[m]][ineg]
      r_all <- rank(c(x, z), ties.method = "average")
      V10[, m] <- (r_all[seq_len(np)] - rank(x, ties.method = "average")) / nn
      V01[, m] <- 1 - (r_all[np + seq_len(nn)] -
                         rank(z, ties.method = "average")) / np
    }

    S <- cov(V10) / np + cov(V01) / nn
    a <- colMeans(V10)

    for (m in mp) {
      if (m == REF_MODEL) next
      d <- a[m] - a[REF_MODEL]
      v <- S[m, m] + S[REF_MODEL, REF_MODEL] - 2 * S[m, REF_MODEL]
      zz <- if (is.finite(v) && v > 0) unname(d) / sqrt(v) else NA_real_
      fold_delong <- bind_rows(fold_delong, tibble(
        Outcome = o, Fold = fl, model = m, n = nrow(wf), n_pos = np,
        delta = unname(d), z = zz,
        p = if (is.na(zz)) NA_real_ else 2 * pnorm(-abs(zz))))
    }
  }
}

fold_delong %>% arrange(Outcome, model, Fold) %>% print(n = 60)

fold_summary <- fold_delong %>%
  group_by(Outcome, model) %>%
  summarise(folds = n(),
            folds_positive = sum(delta > 0, na.rm = TRUE),
            folds_p05 = sum(p < ALPHA, na.rm = TRUE),
            median_delta = median(delta, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(model = factor(model, levels = MODEL_ORDER)) %>%
  arrange(Outcome, model)

print(fold_summary, n = Inf)

# side by side with the pooled result
fold_summary %>%
  left_join(delong_ref %>%
              transmute(Outcome, model = other, pooled_delta = delta_o,
                        pooled_q = q, pooled_verdict = verdict),
            by = c("Outcome", "model")) %>%
  select(Outcome, model, median_delta, folds_positive, folds, pooled_delta,
         pooled_q, pooled_verdict) %>%
  print(n = Inf)


# =============================================================================
# 08  CALIBRATION, REQUIRED BEFORE ANY DCA
# =============================================================================

# Decision curve analysis reads the predicted probability as a probability. It
# compares it against a threshold that encodes how a clinician trades a missed
# death against an unnecessary intervention. If the probabilities are not
# calibrated the comparison is meaningless, and elastic net probabilities from
# a penalised fit are routinely shrunk towards the prevalence. So this block
# comes first, and its output decides whether block 09 can be believed.
#
# Calibration slope: regress the outcome on the linear predictor. A slope of 1
# is right. Below 1 means the predictions are too extreme, above 1 too flat.
# Calibration in the large: the intercept with the linear predictor as an
# offset. Zero is right; a negative value means systematic over-prediction.

calib <- NULL

for (o in sort(unique(wide_test$Outcome))) {
  w <- filter(wide_test, Outcome == o)
  for (m in mods_present) {

    p <- pmin(pmax(w[[m]], 1e-8), 1 - 1e-8)
    lp <- log(p / (1 - p))
    if (sd(lp) == 0) next

    f1 <- glm(w$y_true ~ lp, family = binomial())
    f0 <- glm(w$y_true ~ offset(lp), family = binomial())

    calib <- bind_rows(calib, tibble(
      Outcome = o, model = m,
      slope = unname(coef(f1)[2]),
      slope_se = unname(coef(summary(f1))[2, 2]),
      intercept_large = unname(coef(f0)[1]),
      mean_pred = mean(p), observed = mean(w$y_true),
      brier = mean((p - w$y_true)^2)))
  }
}

calib <- calib %>%
  mutate(model = factor(model, levels = MODEL_ORDER),
         model_lab = factor(MODEL_LABELS[as.character(model)],
                            levels = unname(MODEL_LABELS)),
         slope_lo = slope - qnorm(0.975) * slope_se,
         slope_hi = slope + qnorm(0.975) * slope_se,
         flag = case_when(slope_hi < 1 ~ "too extreme",
                          slope_lo > 1 ~ "too flat",
                          TRUE ~ "consistent with 1"))

calib %>% arrange(Outcome, model) %>%
  select(Outcome, model_lab, slope, slope_lo, slope_hi, intercept_large,
         mean_pred, observed, brier, flag) %>%
  print(n = Inf)

count(calib, flag)

p_calib <- ggplot(calib, aes(slope, fct_rev(model_lab), colour = flag)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = slope_lo, xmax = slope_hi), orientation = "y",
                width = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ Outcome, ncol = 3) +
  scale_colour_manual(values = c(`too extreme` = OKABE[2],
                                 `too flat` = OKABE[5],
                                 `consistent with 1` = OKABE[3]), name = NULL) +
  labs(x = "Calibration slope", y = NULL,
       title = "Calibration of the pooled cross-validated probabilities",
       subtitle = paste0("A slope away from 1 means the decision curves in ",
                         "block 09 are reading the probabilities\nas something ",
                         "they are not. Recalibrate before drawing conclusions ",
                         "from net benefit.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_calib

# the calibration plot itself, for one outcome, by decile of predicted risk
o <- "Total_Death"

p_calib_curve <- wide_test %>%
  filter(Outcome == o) %>%
  select(y_true, all_of(mods_present)) %>%
  pivot_longer(all_of(mods_present), names_to = "model", values_to = "p") %>%
  group_by(model) %>%
  mutate(bin = ntile(p, 10)) %>%
  group_by(model, bin) %>%
  summarise(pred = mean(p), obs = mean(y_true), n = n(),
            se = sqrt(obs * (1 - obs) / n), .groups = "drop") %>%
  mutate(model_lab = factor(MODEL_LABELS[model], levels = unname(MODEL_LABELS))) %>%
  ggplot(aes(pred, obs)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin = pmax(obs - 1.96 * se, 0),
                    ymax = pmin(obs + 1.96 * se, 1)),
                width = 0, linewidth = 0.35, colour = "grey60") +
  geom_point(size = 1.6, colour = OKABE[1]) +
  geom_line(linewidth = 0.4, colour = OKABE[1]) +
  facet_wrap(~ model_lab, ncol = 3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted risk in decile", y = "Observed proportion",
       title = paste0(o, ": calibration by decile of predicted risk"),
       subtitle = "Points below the diagonal are over-prediction.") +
  theme(panel.grid.minor = element_blank())

p_calib_curve


# =============================================================================
# 09  DECISION CURVES
# =============================================================================

# Net benefit at threshold pt, on the scale of true positives per patient:
#
#     NB = TP/n - (FP/n) * pt/(1 - pt)
#
# The weight pt/(1-pt) is the exchange rate the threshold implies. A clinician
# who would act at a 20% risk is saying one missed death is worth four
# unnecessary interventions, and the arithmetic charges false positives at
# exactly that rate. The two reference strategies are treat everyone,
#
#     NB_all = prevalence - (1 - prevalence) * pt/(1 - pt)
#
# which crosses zero exactly at pt equal to the prevalence, and treat nobody,
# which is zero everywhere. A model is only useful over the range of thresholds
# where its curve sits above both.

dca <- NULL

for (o in sort(unique(wide_test$Outcome))) {

  w  <- filter(wide_test, Outcome == o)
  yy <- w$y_true
  n  <- length(yy); prev <- mean(yy)

  for (pt in PT_GRID) {

    wgt <- pt / (1 - pt)

    dca <- bind_rows(dca, tibble(
      Outcome = o, pt = pt, strategy = "Treat all",
      nb = prev - (1 - prev) * wgt))

    dca <- bind_rows(dca, tibble(
      Outcome = o, pt = pt, strategy = "Treat none", nb = 0))

    for (m in mods_present) {
      flag <- w[[m]] >= pt
      tp <- sum(flag & yy == 1); fp <- sum(flag & yy == 0)
      dca <- bind_rows(dca, tibble(
        Outcome = o, pt = pt, strategy = m,
        nb = tp / n - (fp / n) * wgt))
    }
  }
}

dca <- dca %>%
  mutate(strategy_lab = ifelse(strategy %in% names(MODEL_LABELS),
                               MODEL_LABELS[strategy], strategy),
         strategy_lab = factor(strategy_lab,
                               levels = c(unname(MODEL_LABELS[mods_present]),
                                          "Treat all", "Treat none")))

# named, so the mapping survives plotting a subset of the models
DCA_COLS <- c(setNames(OKABE[seq_along(mods_present)],
                       unname(MODEL_LABELS[mods_present])),
              `Treat all` = "grey30", `Treat none` = "grey60")
DCA_COLS

nrow(dca)

# net benefit at the thresholds worth tabulating
dca %>%
  filter(pt %in% PT_TABLE) %>%
  select(Outcome, pt, strategy_lab, nb) %>%
  pivot_wider(names_from = pt, values_from = nb) %>%
  arrange(Outcome, strategy_lab) %>%
  print(n = Inf)

# where each model beats both reference strategies, which is the only claim
# a decision curve supports
useful_range <- dca %>%
  filter(strategy %in% mods_present) %>%
  left_join(dca %>% filter(strategy == "Treat all") %>%
              select(Outcome, pt, nb_all = nb), by = c("Outcome", "pt")) %>%
  mutate(useful = nb > nb_all & nb > 0) %>%
  group_by(Outcome, strategy) %>%
  summarise(pt_min = if (any(useful)) min(pt[useful]) else NA_real_,
            pt_max = if (any(useful)) max(pt[useful]) else NA_real_,
            width  = sum(useful) * (PT_GRID[2] - PT_GRID[1]),
            .groups = "drop") %>%
  mutate(strategy = factor(strategy, levels = MODEL_ORDER)) %>%
  arrange(Outcome, strategy)

print(useful_range, n = Inf)

o <- "Total_Death"

p_dca <- dca %>%
  filter(Outcome == o) %>%
  ggplot(aes(pt, nb, colour = strategy_lab, linetype = strategy_lab %in%
               c("Treat all", "Treat none"))) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = DCA_COLS, name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"),
                        guide = "none") +
  coord_cartesian(ylim = c(-0.05, NA)) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = paste0(o, ": decision curves"),
       subtitle = paste0("A model earns its place only where it sits above ",
                         "both dashed lines.\nRead block 08 first: net benefit ",
                         "is only meaningful if the probabilities are calibrated.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2))

p_dca

# all outcomes at once, models only, with the two references
p_dca_all <- dca %>%
  ggplot(aes(pt, nb, colour = strategy_lab,
             linetype = strategy_lab %in% c("Treat all", "Treat none"))) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Outcome, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = DCA_COLS, name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"),
                        guide = "none") +
  coord_cartesian(ylim = c(-0.05, NA)) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "Decision curves across outcomes") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2))

p_dca_all


# =============================================================================
# 10  NET BENEFIT DIFFERENCES, BOOTSTRAPPED
# =============================================================================

# A decision curve with no uncertainty on it invites the reader to take a gap
# of 0.002 seriously. The bootstrap resamples patients, recomputes both curves
# on each draw, and gives a percentile interval on the difference. It is
# paired, because both models are recomputed on the same resampled patients.
#
# Also reported: net interventions avoided per 100 patients, which is the
# difference in net benefit divided by the threshold odds. It is the same
# information in units a clinician can act on.

o <- "Total_Death"

w  <- filter(wide_test, Outcome == o)
yy <- w$y_true
n  <- length(yy)
cat("bootstrapping", o, "with n =", n, "and", sum(yy), "events\n")

nb_diff <- NULL

for (m in mods_present) {

  if (m == REF_MODEL) next

  for (pt in PT_TABLE) {

    wgt <- pt / (1 - pt)

    fa <- w[[m]] >= pt; fb <- w[[REF_MODEL]] >= pt
    nb_a <- sum(fa & yy == 1)/n - sum(fa & yy == 0)/n * wgt
    nb_b <- sum(fb & yy == 1)/n - sum(fb & yy == 0)/n * wgt

    boot <- numeric(N_BOOT)
    for (b in seq_len(N_BOOT)) {
      idx <- sample.int(n, n, replace = TRUE)
      ys <- yy[idx]; pa <- w[[m]][idx]; pb <- w[[REF_MODEL]][idx]
      ga <- pa >= pt; gb <- pb >= pt
      boot[b] <- (sum(ga & ys == 1)/n - sum(ga & ys == 0)/n * wgt) -
                 (sum(gb & ys == 1)/n - sum(gb & ys == 0)/n * wgt)
    }

    nb_diff <- bind_rows(nb_diff, tibble(
      Outcome = o, model = m, pt = pt,
      nb_model = nb_a, nb_ref = nb_b, diff = nb_a - nb_b,
      lo = unname(quantile(boot, 0.025)),
      hi = unname(quantile(boot, 0.975)),
      p_boot = 2 * min(mean(boot <= 0), mean(boot >= 0)),
      avoided_per_100 = (nb_a - nb_b) / wgt * 100))
  }
}

nb_diff <- nb_diff %>%
  mutate(model = factor(model, levels = MODEL_ORDER),
         model_lab = factor(MODEL_LABELS[as.character(model)],
                            levels = unname(MODEL_LABELS)),
         crosses_zero = lo <= 0 & hi >= 0)

nb_diff %>%
  select(model_lab, pt, nb_model, nb_ref, diff, lo, hi, avoided_per_100,
         crosses_zero) %>%
  print(n = Inf)

p_nb_diff <- ggplot(nb_diff, aes(diff, fct_rev(model_lab),
                                 colour = !crosses_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0, linewidth = 0.5) +
  geom_point(size = 1.9) +
  facet_wrap(~ pt, ncol = 4, labeller = label_both) +
  scale_colour_manual(values = c(`TRUE` = OKABE[2], `FALSE` = OKABE[1]),
                      labels = c(`TRUE` = "interval excludes zero",
                                 `FALSE` = "interval covers zero"),
                      name = NULL) +
  labs(x = paste0("Net benefit against ", MODEL_LABELS[REF_MODEL]), y = NULL,
       title = paste0(o, ": change in net benefit, bootstrap intervals"),
       subtitle = paste0(N_BOOT, " resamples of patients, paired across models.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_nb_diff


# =============================================================================
# 11  ROC CURVES, YOUR SELECTION
# =============================================================================

# Set the three lines below and rerun the block. SEL_MODELS takes any subset in
# any order; SEL_SPLIT switches the whole plot between held-out and training
# predictions, which is the quickest way to see overfitting as a shape rather
# than as a number.

sort(unique(pred_test$Outcome))
mods_present

SEL_OUTCOME <- "Total_Death"
SEL_MODELS  <- c("Model_1", "Model_4", "Model_5", "Model_8", "Model_9")
SEL_SPLIT   <- "test"        # "test" or "train"

src <- if (SEL_SPLIT == "train") pred_train else pred_test

if (!SEL_OUTCOME %in% src$Outcome)
  stop("outcome '", SEL_OUTCOME, "' is not in the ", SEL_SPLIT, " file")
if (!all(SEL_MODELS %in% as.character(src$Model)))
  stop("not in the file: ",
       paste(setdiff(SEL_MODELS, as.character(src$Model)), collapse = ", "))

roc_pts <- NULL
roc_lab <- NULL

for (m in SEL_MODELS) {

  s <- src %>% filter(Outcome == SEL_OUTCOME, Model == m) %>%
    arrange(desc(prob))

  yy <- s$y_true; pp <- s$prob
  np <- sum(yy == 1); nn <- sum(yy == 0)

  # step only at distinct thresholds. Within a tied group the ROC is a
  # straight diagonal, not a staircase, which is why this block draws the
  # curve with geom_line rather than geom_step. On simulated data with heavy
  # ties the staircase encloses 0.885 while the AUC is 0.837; the straight
  # segments reproduce the AUC exactly. Your Model_1 has large tied groups,
  # so this is not a hypothetical difference.
  keep <- c(which(diff(pp) != 0), length(pp))
  tpr <- c(0, cumsum(yy == 1)[keep] / np, 1)
  fpr <- c(0, cumsum(yy == 0)[keep] / nn, 1)

  x <- pp[yy == 1]; z <- pp[yy == 0]
  r_all <- rank(c(x, z), ties.method = "average")
  v10 <- (r_all[seq_len(np)] - rank(x, ties.method = "average")) / nn
  v01 <- 1 - (r_all[np + seq_len(nn)] - rank(z, ties.method = "average")) / np
  a  <- mean(v10)
  se <- sqrt(var(v10) / np + var(v01) / nn)
  lg <- log(a / (1 - a)); se_lg <- se / (a * (1 - a))

  roc_pts <- bind_rows(roc_pts, tibble(model = m, fpr = fpr, tpr = tpr))
  roc_lab <- bind_rows(roc_lab, tibble(
    model = m, auc = a, se = se,
    lo = plogis(lg - qnorm(0.975) * se_lg),
    hi = plogis(lg + qnorm(0.975) * se_lg)))
}

roc_lab <- roc_lab %>%
  mutate(label = sprintf("%s  AUC %.3f (%.3f-%.3f)",
                         MODEL_LABELS[model], auc, lo, hi))

print(roc_lab)

roc_pts <- roc_pts %>%
  left_join(select(roc_lab, model, label), by = "model") %>%
  mutate(label = factor(label, levels = roc_lab$label))

p_roc_sel <- ggplot(roc_pts, aes(fpr, tpr, colour = label)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.6) +
  coord_equal() +
  scale_colour_manual(values = OKABE, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "1 - specificity", y = "Sensitivity",
       title = paste0(SEL_OUTCOME, ": ROC curves, ", SEL_SPLIT, " predictions"),
       subtitle = "Intervals on the AUC are DeLong.") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = "grey80"),
        legend.text = element_text(size = 7),
        panel.grid.minor = element_blank())

p_roc_sel

# the pairwise DeLong results for exactly the models you plotted
delong %>%
  filter(Outcome == SEL_OUTCOME,
         model_a %in% SEL_MODELS, model_b %in% SEL_MODELS) %>%
  mutate(a = MODEL_LABELS[model_a], b = MODEL_LABELS[model_b]) %>%
  select(a, b, delta, se, p, q) %>%
  arrange(q) %>%
  print(n = Inf)


# =============================================================================
# 12  TRAIN AGAINST TEST
# =============================================================================

# The same pooled AUC computed on the training predictions. The difference is
# the optimism, now measured on patients rather than on fold summaries.

auc_train <- NULL

wide_train <- pred_train %>%
  select(Outcome, KAUST_ID, Fold, y_true, Model, prob) %>%
  mutate(Model = as.character(Model)) %>%
  pivot_wider(names_from = Model, values_from = prob)

# a patient appears in every training fold but one, so this frame is longer
# than the test one and has to be summarised within fold before pooling
dim(wide_train)

for (o in sort(unique(wide_train$Outcome))) {
  for (fl in sort(unique(wide_train$Fold))) {

    w <- filter(wide_train, Outcome == o, Fold == fl)
    mp <- mods_present[mods_present %in% names(w)]
    w <- w[complete.cases(w[, mp]), ]

    yy <- w$y_true
    ip <- which(yy == 1); ineg <- which(yy == 0)
    np <- length(ip); nn <- length(ineg)
    if (np < 5 || nn < 5) next

    for (m in mp) {
      x <- w[[m]][ip]; z <- w[[m]][ineg]
      r_all <- rank(c(x, z), ties.method = "average")
      v10 <- (r_all[seq_len(np)] - rank(x, ties.method = "average")) / nn
      auc_train <- bind_rows(auc_train, tibble(
        Outcome = o, Fold = fl, model = m, auc = mean(v10), n = nrow(w)))
    }
  }
}

auc_train_pooled <- auc_train %>%
  group_by(Outcome, model) %>%
  summarise(auc_train = mean(auc), .groups = "drop")

optimism <- auc_tbl %>%
  transmute(Outcome, model = as.character(model), auc_test = auc) %>%
  left_join(auc_train_pooled, by = c("Outcome", "model")) %>%
  mutate(gap = auc_train - auc_test,
         model = factor(model, levels = MODEL_ORDER),
         model_lab = factor(MODEL_LABELS[as.character(model)],
                            levels = unname(MODEL_LABELS))) %>%
  arrange(Outcome, model)

print(optimism, n = Inf)

p_opt <- optimism %>%
  select(Outcome, model_lab, Train = auc_train, Test = auc_test) %>%
  pivot_longer(c(Train, Test), names_to = "split", values_to = "auc") %>%
  mutate(split = factor(split, levels = c("Train", "Test"))) %>%
  ggplot(aes(model_lab, auc, colour = split, group = split)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.9) +
  facet_wrap(~ Outcome, ncol = 3) +
  scale_colour_manual(values = c(Train = OKABE[5], Test = OKABE[1]), name = NULL) +
  labs(x = NULL, y = "Pooled AUC",
       title = "Training against held-out discrimination, patient level",
       subtitle = "Training AUCs are averaged over the folds a patient was trained in.") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank())

p_opt


# =============================================================================
# 13  SAVING
# =============================================================================

# ggsave("fig_delong_vs_ref.pdf", p_delong,      width = 11, height = 8)
# ggsave("fig_delong_grid.pdf",   p_dl_grid,     width = 11, height = 9)
# ggsave("fig_calibration.pdf",   p_calib,       width = 11, height = 8)
# ggsave("fig_calib_curve.pdf",   p_calib_curve, width =  9, height = 8)
# ggsave("fig_dca.pdf",           p_dca,         width =  8, height = 6)
# ggsave("fig_dca_all.pdf",       p_dca_all,     width = 11, height = 8)
# ggsave("fig_nb_diff.pdf",       p_nb_diff,     width = 11, height = 6)
# ggsave("fig_roc_selected.pdf",  p_roc_sel,     width =  7, height = 7)
# ggsave("fig_optimism_auc.pdf",  p_opt,         width = 11, height = 8)
#
# write_csv(auc_tbl,      "pooled_auc_delong.csv")
# write_csv(delong,       "delong_pairwise.csv")
# write_csv(fold_delong,  "delong_within_fold.csv")
# write_csv(calib,        "calibration.csv")
# write_csv(dca,          "decision_curves.csv")
# write_csv(nb_diff,      "net_benefit_differences.csv")
# write_csv(optimism,     "auc_train_vs_test.csv")
