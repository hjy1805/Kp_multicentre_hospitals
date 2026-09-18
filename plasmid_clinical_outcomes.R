# ============================================================
# CANONICAL hvKp VIRULENCE PLASMIDS AND CLINICAL SEVERITY
#
# Reviewer 2, major comment 4: test whether presence of the
# canonical hypervirulence plasmids is linked to increased
# clinical severity.
#
#   AY378100  pLVPK        IncFIB(K)/IncHI1B(Mar)
#   AP006726  pK2044       IncFIB(K)/IncHI1B(Mar), NTUH-K2044
#
# Plasmid presence is called from breadth of coverage
# (callable bases) after read mapping, at two thresholds.
#
# OUTCOMES
#   Mortality       binary, logistic, odds ratio
#   ICU_Admission   binary, logistic, odds ratio
#   LOS             continuous, log1p linear model, % change
#
# WHAT THIS SCRIPT ADDS TO A PLAIN glm() LOOP
#   1. Reports how many isolates failed to match between the
#      callability reports and the analysis table, instead of
#      letting match() turn them silently into NA.
#   2. Prints the full 2x2 table per test and refuses to
#      report an effect estimate from a table with an empty
#      cell, where the odds ratio is not identifiable.
#   3. Fits an adjusted model alongside the crude one, so the
#      association is not read off a model that ignores
#      resistance status.
#   4. Applies Benjamini-Hochberg across the family of tests.
#   5. Checks how many isolates sit between the two
#      thresholds, which is where a breadth-based call is
#      least stable.
#   6. Reports the agreement between the two plasmids, which
#      are largely homologous and therefore not independent
#      tests.
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)


# ------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------

BASE_DIR <- "/Users/daneshm/Documents/Kp_KAIMRC/revision"

PLASMIDS <- list(
  AP006726 = list(
    label = "pK2044 (AP006726)",
    file = file.path(
      BASE_DIR, "callable_bases_report_plasmid_AP006726.tsv"
    )
  ),
  AY378100 = list(
    label = "pLVPK (AY378100)",
    file = file.path(
      BASE_DIR, "callable_bases_report_plasmid_AY378100.tsv"
    )
  )
)

CUTOFFS <- c(90, 95)

# Column names in the callability reports
ID_COL <- "Sample"
CALLABLE_COL <- "Callable_percent"

# Key linking the analysis table to the callability reports
ANALYSIS_ID_COL <- "Name"

BINARY_OUTCOMES <- c("Mortality", "ICU_Admission")
LOS_OUTCOME <- "LOS"

# Covariates for the adjusted model. Leave empty to fit the
# crude model only. These must be columns of
# kleborate_encoded_df. Set this to the same confounder set
# used elsewhere in the manuscript, e.g. the ESBL/carbapenemase
# indicator plus the demographic and comorbidity variables, so
# that this analysis is comparable with the others.
ADJUST_VARS <- character(0)
# ADJUST_VARS <- c("ESBL_CP_positive", "Age", "Sex", "Charlson")

# Minimum count in any cell of the 2x2 before an odds ratio is
# reported rather than suppressed
MIN_CELL <- 1

OUT_RESULTS <- file.path(
  BASE_DIR, "plasmid_clinical_outcomes_90_95.csv"
)
OUT_COUNTS <- file.path(
  BASE_DIR, "plasmid_clinical_outcomes_counts.csv"
)


# ------------------------------------------------------------
# 2. READ THE CALLABILITY REPORTS
# ------------------------------------------------------------

read_callability <- function(spec, name) {
  d <- read_tsv(spec$file, show_col_types = FALSE)

  missing_cols <- setdiff(c(ID_COL, CALLABLE_COL), names(d))
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s) in ", basename(spec$file), ": ",
      paste(missing_cols, collapse = ", "),
      "\nColumns present: ", paste(names(d), collapse = ", ")
    )
  }

  d %>%
    transmute(
      sample = .data[[ID_COL]],
      callable = as.numeric(.data[[CALLABLE_COL]])
    ) %>%
    distinct(sample, .keep_all = TRUE)
}

callability <- imap(PLASMIDS, read_callability)

cat("\n=== CALLABILITY DISTRIBUTIONS ===\n")
for (p in names(callability)) {
  cat("\n", PLASMIDS[[p]]$label, "  n = ",
      nrow(callability[[p]]), "\n", sep = "")
  print(summary(callability[[p]]$callable))
}


# ------------------------------------------------------------
# 3. ATTACH CALLABILITY TO THE ANALYSIS TABLE
#
# match() returns NA for anything it cannot find. Left
# unchecked, an identifier-format mismatch produces a column
# of NAs, drop_na() removes those rows, and the models run on
# a silently truncated dataset. Report the overlap explicitly.
# ------------------------------------------------------------

stopifnot(exists("kleborate_encoded_df"))
stopifnot(ANALYSIS_ID_COL %in% names(kleborate_encoded_df))

analysis_ids <- kleborate_encoded_df[[ANALYSIS_ID_COL]]

cat("\n=== IDENTIFIER MATCHING ===\n")
for (p in names(callability)) {
  hits <- analysis_ids %in% callability[[p]]$sample

  cat(
    PLASMIDS[[p]]$label, ": ",
    sum(hits), " of ", length(hits),
    " isolates matched (", sum(!hits), " unmatched)\n",
    sep = ""
  )

  if (any(!hits)) {
    cat("  first unmatched: ",
        paste(utils::head(analysis_ids[!hits], 5),
              collapse = ", "),
        "\n", sep = "")
  }

  kleborate_encoded_df[[paste0(p, "_callable")]] <-
    callability[[p]]$callable[
      match(analysis_ids, callability[[p]]$sample)
    ]
}


# ------------------------------------------------------------
# 4. CALL PRESENCE AT EACH THRESHOLD
#
# NA callability stays NA rather than becoming 0, so an
# isolate that was never mapped is dropped from the model
# instead of being counted as plasmid-negative.
# ------------------------------------------------------------

for (p in names(PLASMIDS)) {
  cal <- kleborate_encoded_df[[paste0(p, "_callable")]]
  for (k in CUTOFFS) {
    kleborate_encoded_df[[paste0(p, "_", k)]] <-
      ifelse(is.na(cal), NA_integer_, as.integer(cal >= k))
  }
}

plasmid_vars <- as.vector(outer(
  names(PLASMIDS), CUTOFFS,
  function(p, k) paste0(p, "_", k)
))

cat("\n=== PREVALENCE BY THRESHOLD ===\n")
for (v in plasmid_vars) {
  cat("\n", v, "\n", sep = "")
  print(table(kleborate_encoded_df[[v]], useNA = "ifany"))
}


# ------------------------------------------------------------
# 5. STABILITY OF THE CALL BETWEEN THE TWO THRESHOLDS
#
# Isolates whose breadth falls between 90% and 95% change
# class depending on the threshold. If that grey zone holds a
# large share of the collection, neither threshold is a
# defensible single answer and the result should be reported
# at both, as a sensitivity analysis.
# ------------------------------------------------------------

cat("\n=== GREY ZONE BETWEEN THRESHOLDS ===\n")
grey_zone <- map_dfr(names(PLASMIDS), function(p) {
  cal <- kleborate_encoded_df[[paste0(p, "_callable")]]
  cal <- cal[!is.na(cal)]

  tibble(
    Plasmid = PLASMIDS[[p]]$label,
    n_called = length(cal),
    below_90 = sum(cal < min(CUTOFFS)),
    grey_90_95 = sum(cal >= min(CUTOFFS) & cal < max(CUTOFFS)),
    at_or_above_95 = sum(cal >= max(CUTOFFS)),
    pct_grey = round(
      100 * sum(cal >= min(CUTOFFS) & cal < max(CUTOFFS)) /
        length(cal),
      1
    )
  )
})
print(as.data.frame(grey_zone))


# ------------------------------------------------------------
# 6. AGREEMENT BETWEEN THE TWO PLASMIDS
#
# pLVPK and pK2044 share extensive sequence homology, so
# these are not two independent tests of two independent
# exposures. Quantify the overlap before interpreting them
# as separate findings.
# ------------------------------------------------------------

cat("\n=== AGREEMENT BETWEEN PLASMIDS ===\n")
for (k in CUTOFFS) {
  a <- kleborate_encoded_df[[paste0(names(PLASMIDS)[1], "_", k)]]
  b <- kleborate_encoded_df[[paste0(names(PLASMIDS)[2], "_", k)]]
  ok <- !is.na(a) & !is.na(b)

  cat("\nThreshold ", k, "%  (n = ", sum(ok), ")\n", sep = "")
  print(table(
    a[ok], b[ok],
    dnn = c(names(PLASMIDS)[1], names(PLASMIDS)[2])
  ))
  cat("  agreement: ",
      round(100 * mean(a[ok] == b[ok]), 1), "%\n", sep = "")
}


# ------------------------------------------------------------
# 7. CROSS-TABULATION AND SEPARATION CHECK
#
# A 2x2 with an empty cell gives an odds ratio of 0 or Inf
# with a confidence interval spanning everything. glm() will
# return a large coefficient and a large standard error
# without warning, which reads as a null result when it is
# actually an unidentifiable one.
# ------------------------------------------------------------

cross_tab <- function(data, predictor, outcome) {
  dat <- data %>%
    select(all_of(c(predictor, outcome))) %>%
    drop_na()

  tibble(
    Predictor = predictor,
    Outcome = outcome,
    N = nrow(dat),
    present_event = sum(dat[[predictor]] == 1 &
                          dat[[outcome]] == 1),
    present_no_event = sum(dat[[predictor]] == 1 &
                             dat[[outcome]] == 0),
    absent_event = sum(dat[[predictor]] == 0 &
                         dat[[outcome]] == 1),
    absent_no_event = sum(dat[[predictor]] == 0 &
                            dat[[outcome]] == 0)
  ) %>%
    mutate(
      min_cell = pmin(
        present_event, present_no_event,
        absent_event, absent_no_event
      ),
      estimable = min_cell >= MIN_CELL
    )
}

counts <- expand_grid(
  predictor = plasmid_vars,
  outcome = BINARY_OUTCOMES
) %>%
  pmap_dfr(function(predictor, outcome) {
    cross_tab(kleborate_encoded_df, predictor, outcome)
  })

cat("\n=== 2x2 COUNTS ===\n")
print(as.data.frame(counts))

if (any(!counts$estimable)) {
  cat("\nNOT ESTIMABLE (empty cell), suppressed below:\n")
  print(as.data.frame(
    counts %>% filter(!estimable) %>%
      select(Predictor, Outcome, N, min_cell)
  ))
}


# ------------------------------------------------------------
# 8. MODEL FITTING
#
# Each predictor is fitted twice: crude, and adjusted for
# ADJUST_VARS if any are supplied. Only the plasmid term is
# retained from the adjusted fit.
# ------------------------------------------------------------

build_formula <- function(predictor, response, covariates) {
  reformulate(c(predictor, covariates), response = response)
}

fit_logistic <- function(data, predictor, outcome,
                         covariates = character(0)) {
  vars <- c(predictor, outcome, covariates)

  dat <- data %>%
    select(all_of(vars)) %>%
    drop_na()

  if (n_distinct(dat[[predictor]]) < 2 ||
      n_distinct(dat[[outcome]]) < 2) {
    return(NULL)
  }

  fit <- glm(
    build_formula(predictor, outcome, covariates),
    data = dat,
    family = binomial()
  )

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == predictor) %>%
    mutate(
      Predictor = predictor,
      Outcome = outcome,
      Model = if (length(covariates) == 0) "Crude" else "Adjusted",
      Scale = "Odds ratio",
      N = nrow(dat),
      Present = sum(dat[[predictor]] == 1),
      Absent = sum(dat[[predictor]] == 0)
    )
}

# Length of stay on the log1p scale, back-transformed to a
# percentage change. log1p rather than log because LOS can be
# zero for same-day discharge.
fit_los <- function(data, predictor,
                    covariates = character(0)) {
  vars <- c(predictor, LOS_OUTCOME, covariates)

  dat <- data %>%
    select(all_of(vars)) %>%
    drop_na() %>%
    filter(.data[[LOS_OUTCOME]] >= 0)

  if (n_distinct(dat[[predictor]]) < 2 || nrow(dat) < 3) {
    return(NULL)
  }

  fit <- lm(
    build_formula(
      predictor, paste0("log1p(", LOS_OUTCOME, ")"), covariates
    ),
    data = dat
  )

  tidy(fit, conf.int = TRUE) %>%
    filter(term == predictor) %>%
    mutate(
      log_estimate = estimate,
      log_conf.low = conf.low,
      log_conf.high = conf.high,
      estimate = 100 * (exp(log_estimate) - 1),
      conf.low = 100 * (exp(log_conf.low) - 1),
      conf.high = 100 * (exp(log_conf.high) - 1),
      Predictor = predictor,
      Outcome = "Length of stay",
      Model = if (length(covariates) == 0) "Crude" else "Adjusted",
      Scale = "% change in LOS",
      N = nrow(dat),
      Present = sum(dat[[predictor]] == 1),
      Absent = sum(dat[[predictor]] == 0)
    )
}


# ------------------------------------------------------------
# 9. RUN EVERYTHING
# ------------------------------------------------------------

covariate_sets <- list(crude = character(0))
if (length(ADJUST_VARS) > 0) {
  missing_adj <- setdiff(ADJUST_VARS, names(kleborate_encoded_df))
  if (length(missing_adj) > 0) {
    stop(
      "ADJUST_VARS not found in kleborate_encoded_df: ",
      paste(missing_adj, collapse = ", ")
    )
  }
  covariate_sets$adjusted <- ADJUST_VARS
}

results_raw <- map_dfr(covariate_sets, function(covs) {
  binary <- expand_grid(
    predictor = plasmid_vars,
    outcome = BINARY_OUTCOMES
  ) %>%
    pmap_dfr(function(predictor, outcome) {
      fit_logistic(
        kleborate_encoded_df, predictor, outcome, covs
      )
    })

  los <- map_dfr(plasmid_vars, function(v) {
    fit_los(kleborate_encoded_df, v, covs)
  })

  bind_rows(binary, los)
})


# ------------------------------------------------------------
# 10. ANNOTATE, SUPPRESS UNIDENTIFIABLE ESTIMATES, CORRECT
#
# Correction is applied within each model type, across the
# whole family of plasmid x threshold x outcome tests. The
# two thresholds are nested rather than independent, so this
# is conservative in one direction and anti-conservative in
# the other; report it as a sensitivity analysis, not as
# twice as many findings.
# ------------------------------------------------------------

results <- results_raw %>%
  left_join(
    counts %>% select(Predictor, Outcome, min_cell, estimable),
    by = c("Predictor", "Outcome")
  ) %>%
  mutate(
    estimable = ifelse(is.na(estimable), TRUE, estimable),
    Plasmid = case_when(
      grepl("AP006726", Predictor) ~ PLASMIDS$AP006726$label,
      grepl("AY378100", Predictor) ~ PLASMIDS$AY378100$label
    ),
    Cutoff = paste0(sub(".*_", "", Predictor), "%"),
    across(
      c(estimate, conf.low, conf.high, p.value),
      ~ ifelse(estimable, .x, NA_real_)
    )
  ) %>%
  group_by(Model) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(Model, Plasmid, Outcome, Cutoff) %>%
  select(
    Plasmid, Cutoff, Outcome, Model, Scale,
    N, Present, Absent, min_cell, estimable,
    estimate, conf.low, conf.high, p.value, p.adj
  )

cat("\n=== RESULTS ===\n")
print(as.data.frame(results))


# ------------------------------------------------------------
# 11. SAVE
# ------------------------------------------------------------

write_csv(results, OUT_RESULTS)
write_csv(counts, OUT_COUNTS)

cat("\nWritten:\n  ", OUT_RESULTS, "\n  ", OUT_COUNTS, "\n", sep = "")
