# =============================================================================
# TEMPORAL VALIDATION: TRAIN ON THE PAST, PREDICT THE NEXT YEAR
#
# The same file layout as the regional and lineage folders, and the same
# random-control design, but the held-out unit is a calendar year and the
# training set is only the years before it. That difference matters more than it
# looks. Holding out a region or a clone is a symmetric experiment: any of them
# could have been the held-out one. Holding out a year is not symmetric, because
# the training data can only come from the past, and there is less of it in 2023
# than in 2025.
#
# This is also the only split in the project that matches how the model would
# actually be used. A model deployed in 2026 is fitted on everything up to 2025
# and asked about patients who have not arrived yet. Cross-validation never asks
# that question, so the numbers here, not the fold numbers, are the ones to
# quote as expected performance. Block 12 makes that comparison explicit.
#
# File names carry five fields separated by "|"
#   split | year | scheme | model | outcome.csv
#   e.g. test|2025|nonrandom|Model_8|Total_Death.csv
#
# TWO CONFOUNDS THIS DESIGN CARRIES, both checkable and both checked below.
#
#   TRAINING SIZE. The 2023 model saw the least data and the 2025 model the
#   most, so year is confounded with training size, and in opposite directions
#   across the sequence: a large penalty in 2023 may be a small training set
#   rather than drift, while a large penalty in 2025 cannot be, because that
#   model had the most data available. A penalty that persists into the last
#   year is the strong version of the finding. Block 05.
#
#   THE CONTROL'S TRAINING SET. The random control is a control for the test set
#   only if it also trained on the same number of patients. If the control was
#   fitted on a random majority of the whole collection while the temporal model
#   was fitted on two years, the gap between them measures training size at
#   least as much as it measures drift. This is the check that did not matter in
#   the regional folder and matters here. Block 04.
#
# WHAT THREE YEARS CAN SUPPORT. The sign test floor at n = 3 is p = 0.25, so
# three years all pointing one way cannot be significant. For the paired t with
# the correlation correction, the mean penalty has to exceed roughly three times
# its own spread across years to reach 0.05. Report the three numbers and their
# direction over time, not a pooled p-value.
#
# Flat script, no named functions. Run block by block. Figures print to the
# device. The saving block is at the end, commented out.
#
# Blocks
#   01 packages and settings        08 is the penalty growing or shrinking
#   02 find and parse the files     09 does genomics drift faster
#   03 read                         10 train against test
#   04 the split-matching check     11 outcome and case-mix drift
#   05 year composition and growth  12 the deployment estimate
#   06 the temporal penalty         13 the summary figure
#   07 how far toward chance        14 saving (commented out)
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

DIR <- "/Users/daneshm/Documents/PA_KAIMRC/ML_files/Temporal"

MODEL_LABELS <- c(
  Model_1 = "1  Tier 1",
  Model_2 = "2  + Tier 2",
  Model_3 = "3  + Tier 3",
  Model_4 = "4  Clinical only",
  Model_5 = "5  + ARG/vir",
  Model_6 = "6  + pangenome",
  Model_7 = "7  + SNPs",
  Model_8 = "8  Clinical + genomic",
  Model_9 = "9  Genomic only")

MODEL_ORDER <- names(MODEL_LABELS)

REF_MODEL   <- "Model_4"
TEST_METRIC <- "ROC_AUC"
ALPHA <- 0.05

SCHEME_GEO  <- "nonrandom"   # trained on earlier years only
SCHEME_CTRL <- "random"      # the size-matched control

# If you also ran the fold pipeline on the same cohort, put its pooled ROC AUC
# per model here and block 12 will quantify how optimistic it was. Leave as NULL
# to skip that comparison.
CV_REFERENCE <- NULL
# CV_REFERENCE <- c(Model_4 = 0.000, Model_8 = 0.000, Model_9 = 0.000)

OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#999999", "#000000")

theme_set(theme_bw(base_size = 11))

NB_CORRECT  <- TRUE
NB_RHO_MULT <- 1


# =============================================================================
# 02  FIND AND PARSE THE FILE NAMES
# =============================================================================

files_all <- list.files(DIR, pattern = "\\.csv$", full.names = TRUE)
length(files_all)

n_pipes <- str_count(basename(files_all), fixed("|"))
table(n_pipes)

basename(files_all)[n_pipes != 4]     # anything excluded, look at it

files <- files_all[n_pipes == 4]
length(files)

parts <- str_split_fixed(basename(files), fixed("|"), 5)
head(parts)

file_index <- tibble(
  path    = files,
  split   = parts[, 1],
  year    = parts[, 2],
  scheme  = parts[, 3],
  model   = parts[, 4],
  outcome = str_remove(parts[, 5], "\\.csv$"))

count(file_index, split)
count(file_index, year)
count(file_index, scheme)
count(file_index, model)
count(file_index, outcome)

file_index %>% count(year, scheme, model) %>% print(n = Inf)
file_index %>% count(year, scheme, model, outcome) %>% count(n)

YEARS  <- sort(unique(file_index$year))          # ordered, unlike regions
MODELS <- MODEL_ORDER[MODEL_ORDER %in% file_index$model]
YEARS; MODELS

length(YEARS)   # the sample size for every pooled test in this script

# the ordering is the one thing this design has that the others do not
as.integer(YEARS)


# =============================================================================
# 03  READ
# =============================================================================

raw <- NULL

for (i in seq_len(nrow(file_index))) {
  d <- read_csv(file_index$path[i], show_col_types = FALSE)
  d$row_in_file <- seq_len(nrow(d))
  d$split   <- file_index$split[i]
  d$year    <- file_index$year[i]
  d$scheme  <- file_index$scheme[i]
  d$model   <- file_index$model[i]
  d$outcome <- file_index$outcome[i]
  raw <- bind_rows(raw, d)
}

dim(raw); names(raw)

raw <- select(raw, -any_of(c("...1", "X1")))
glimpse(raw)

# rows per file. One means the random control is a single draw with its own
# noise; several mean it averages, which matters here because three years is
# very little to work with.
raw %>% count(split, year, scheme, model) %>% count(n)

count(raw, Dataset)

dat <- raw %>%
  mutate(model = factor(model, levels = MODEL_ORDER),
         model_lab = factor(MODEL_LABELS[as.character(model)],
                            levels = unname(MODEL_LABELS)),
         year = factor(year, levels = YEARS),
         year_num = as.integer(as.character(year)),
         prevalence = (TP + FN) / N)

test  <- filter(dat, split == "test")
train <- filter(dat, split == "train")

nrow(test); nrow(train)

METRICS <- c("ROC_AUC", "PR_AUC", "Sensitivity", "Specificity",
             "Balanced_Accuracy", "Precision", "F1", "MCC")

max(abs(test$Sensitivity - test$Recall), na.rm = TRUE)   # expect 0


# =============================================================================
# 04  THE SPLIT-MATCHING CHECK
# =============================================================================

# Two matchings to check here, not one. The test sets must be the same size, as
# in the other folders. The TRAINING sets must also be the same size, which did
# not arise before: in the regional design the control trained on three regions
# worth of patients just as the geographic model did, but a temporal model for
# 2023 may have trained on two years while its random control trained on a
# majority of the whole collection. If the training sizes differ, the gap in
# block 06 is partly a learning-curve effect and cannot be reported as drift.

match_test <- test %>%
  group_by(year, model, scheme) %>%
  summarise(N = mean(N), events = mean(TP + FN), prev = mean(prevalence),
            .groups = "drop") %>%
  pivot_wider(names_from = scheme, values_from = c(N, events, prev)) %>%
  mutate(dN     = .data[[paste0("N_", SCHEME_GEO)]] -
                  .data[[paste0("N_", SCHEME_CTRL)]],
         dEvent = .data[[paste0("events_", SCHEME_GEO)]] -
                  .data[[paste0("events_", SCHEME_CTRL)]],
         dPrev  = .data[[paste0("prev_", SCHEME_GEO)]] -
                  .data[[paste0("prev_", SCHEME_CTRL)]])

print(match_test, n = Inf)

match_train <- train %>%
  group_by(year, model, scheme) %>%
  summarise(N = mean(N), events = mean(TP + FN), .groups = "drop") %>%
  pivot_wider(names_from = scheme, values_from = c(N, events)) %>%
  mutate(dN_train     = .data[[paste0("N_", SCHEME_GEO)]] -
                        .data[[paste0("N_", SCHEME_CTRL)]],
         dEvent_train = .data[[paste0("events_", SCHEME_GEO)]] -
                        .data[[paste0("events_", SCHEME_CTRL)]],
         ratio_train  = .data[[paste0("N_", SCHEME_GEO)]] /
                        .data[[paste0("N_", SCHEME_CTRL)]])

print(match_train, n = Inf)

# the four numbers that decide whether block 06 measures drift or measures
# having less data to learn from
match_test  %>% summarise(max_abs_dN_test = max(abs(dN)),
                          max_abs_dEvent_test = max(abs(dEvent)))
match_train %>% summarise(max_abs_dN_train = max(abs(dN_train)),
                          min_ratio_train = min(ratio_train),
                          max_ratio_train = max(ratio_train))

# ratio_train well below 1 for the early years means the temporal model was
# handicapped by having less to learn from, and the penalty for those years has
# to be reported as drift plus learning curve rather than as drift.


# =============================================================================
# 05  YEAR COMPOSITION AND TRAINING SET GROWTH
# =============================================================================

comp <- test %>%
  filter(scheme == SCHEME_GEO) %>%
  group_by(year, year_num) %>%
  summarise(N_test = mean(N), events_test = mean(TP + FN),
            prev_test = mean(prevalence), .groups = "drop") %>%
  left_join(train %>%
              filter(scheme == SCHEME_GEO) %>%
              group_by(year) %>%
              summarise(N_train = mean(N), events_train = mean(TP + FN),
                        prev_train = mean(prevalence), .groups = "drop"),
            by = "year") %>%
  mutate(rho = N_test / N_train,
         cumulative = N_train + N_test) %>%
  arrange(year_num)

print(comp, n = Inf)

# THE GROWTH. N_train should rise across the years, because each model gets one
# more year of history. How steeply it rises tells you how much of any trend in
# block 08 could be a learning curve rather than drift.
comp %>% select(year, N_train, N_test, events_test, prev_test, rho) %>%
  print(n = Inf)

diff(comp$N_train)          # how much history each year added

# how far back the record goes, inferred: the first held-out year trained on
# everything before it, so N_train for that year is the size of the pre-period
comp$N_train[1]

# events in the held-out year governs how noisy that year's AUC is. Below about
# 25 deaths an AUC interval spans most of the useful range.
comp %>% select(year, events_test) %>% arrange(events_test) %>% print(n = Inf)

RHO <- mean(comp$rho, na.rm = TRUE) * NB_RHO_MULT
RHO

p_comp <- comp %>%
  select(year, Training = N_train, `Held-out year` = N_test) %>%
  pivot_longer(-year, names_to = "part", values_to = "N") %>%
  ggplot(aes(year, N, fill = part)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(Training = OKABE[8],
                               `Held-out year` = OKABE[2]), name = NULL) +
  labs(x = NULL, y = "Patients",
       title = "How much history each year's model had",
       subtitle = paste0("The training bar grows across the sequence, so a ",
                         "penalty that persists in the last\nyear cannot be ",
                         "explained by having too little to learn from.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_comp


# =============================================================================
# 06  THE TEMPORAL PENALTY
# =============================================================================

# penalty = control minus temporal. Positive means the model did worse on the
# next year than on the same number of patients removed at random from across
# the whole period, and that difference is what predicting forward cost.

pen <- test %>%
  select(year, year_num, scheme, model, model_lab, all_of(METRICS)) %>%
  pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value") %>%
  group_by(year, year_num, model, model_lab, metric, scheme) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = scheme, values_from = value) %>%
  mutate(geo = .data[[SCHEME_GEO]], ctrl = .data[[SCHEME_CTRL]],
         penalty = ctrl - geo)

pen %>%
  filter(metric == TEST_METRIC) %>%
  arrange(model, year_num) %>%
  select(model_lab, year, ctrl, geo, penalty) %>%
  print(n = Inf)

pen_sum <- pen %>%
  filter(metric == TEST_METRIC) %>%
  group_by(model, model_lab) %>%
  summarise(k_years = n(), mean_ctrl = mean(ctrl), mean_geo = mean(geo),
            mean_penalty = mean(penalty), sd_penalty = sd(penalty),
            min_penalty = min(penalty), max_penalty = max(penalty),
            years_worse = sum(penalty > 0),
            worst_year = as.character(year)[which.max(penalty)],
            last_year_penalty =
              penalty[which.max(year_num)],
            .groups = "drop")

print(pen_sum, n = Inf)

# last_year_penalty is the column to read first. It is the penalty for the model
# with the most history behind it, so it is the one least confounded with
# training size and the best estimate of what deploying now would cost.

pen_test <- NULL

for (m in MODELS) {

  v <- pen$penalty[pen$metric == TEST_METRIC & pen$model == m]
  v <- v[!is.na(v)]
  kk <- length(v)
  if (kk < 2) next

  s <- sd(v)
  se_use <- if (NB_CORRECT) s * sqrt(1 / kk + RHO) else s / sqrt(kk)
  tstat  <- if (se_use > 0) mean(v) / se_use else NA_real_
  mult   <- qt(0.975, kk - 1)

  pen_test <- bind_rows(pen_test, tibble(
    model = m, k = kk, mean_penalty = mean(v), sd_penalty = s,
    lo = mean(v) - mult * se_use, hi = mean(v) + mult * se_use,
    p = if (is.na(tstat)) NA_real_ else 2 * pt(-abs(tstat), df = kk - 1),
    years_up = sum(v > 0), years_down = sum(v < 0),
    mean_over_sd = mean(v) / s,
    threshold_ratio = mult * sqrt(1 / kk + RHO)))
}

pen_test <- pen_test %>%
  mutate(model_lab = factor(MODEL_LABELS[model], levels = unname(MODEL_LABELS)))

print(pen_test, n = Inf)

# mean_over_sd has to exceed threshold_ratio for p < 0.05. At three years the
# threshold is around three, which is why the p column is not the headline.
binom.test(length(YEARS), length(YEARS), 0.5)$p.value   # sign test floor, 0.25

p_paired <- pen %>%
  filter(metric == TEST_METRIC) %>%
  select(year, model_lab, Control = ctrl, `Next year` = geo) %>%
  pivot_longer(c(Control, `Next year`), names_to = "scheme",
               values_to = "auc") %>%
  mutate(scheme = factor(scheme, levels = c("Control", "Next year"))) %>%
  ggplot(aes(scheme, auc, group = year, colour = year)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ model_lab, nrow = 1) +
  scale_colour_manual(values = OKABE, name = NULL) +
  scale_y_continuous(limits = c(0.35, 1)) +
  labs(x = NULL, y = TEST_METRIC,
       title = "The same patients removed at random, and removed as a future year",
       subtitle = "A line sloping down is a model that does not hold up over time.") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_paired


# =============================================================================
# 07  HOW FAR TOWARD CHANCE
# =============================================================================

# The penalty in AUC points is hard to compare across models, because one
# starting at 0.62 cannot lose what one starting at 0.85 can. Expressed as a
# share of the control model's distance above chance,
#
#     collapse = (ctrl - geo) / (ctrl - 0.5)

collapse <- pen %>%
  filter(metric == TEST_METRIC) %>%
  mutate(above_chance = ctrl - 0.5,
         collapse = ifelse(above_chance > 0, penalty / above_chance, NA_real_))

collapse %>%
  arrange(model, year_num) %>%
  select(model_lab, year, ctrl, geo, collapse) %>%
  print(n = Inf)

collapse %>%
  group_by(model_lab) %>%
  summarise(median_collapse = median(collapse, na.rm = TRUE),
            worst_collapse = max(collapse, na.rm = TRUE),
            years_at_or_below_chance = sum(geo <= 0.5), .groups = "drop") %>%
  print(n = Inf)

p_collapse <- ggplot(collapse, aes(collapse, fct_rev(model_lab), colour = year)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = OKABE, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Share of the control model's advantage over chance that was lost",
       y = NULL,
       title = "How much of the signal survives into the next year",
       subtitle = paste0("The dashed line at 100% is a model that fell to ",
                         "chance on the year it had not seen.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_collapse


# =============================================================================
# 08  IS THE PENALTY GROWING OR SHRINKING
# =============================================================================

# The one thing this design has that the regional and lineage ones do not: the
# held-out units are ordered, so the penalty has a direction in time. Three
# readings, and they mean different things.
#
#   PENALTY SHRINKING   Later models have more history and lose less. The drift
#                       is being absorbed by accumulating data, which is an
#                       argument for periodic refitting and for nothing worse.
#   PENALTY FLAT        A stable cost of predicting forward, unaffected by how
#                       much history is available. This is the number to quote
#                       as the cost of deployment.
#   PENALTY GROWING     The population is moving faster than the training data
#                       accumulates. This is the finding that would matter most
#                       and the one that most needs more than three years
#                       before it is claimed.
#
# Three points. This is a direction, not a trend test, and the script says so
# rather than fitting a line and reporting its p-value.

trend <- pen %>%
  filter(metric == TEST_METRIC) %>%
  arrange(model, year_num) %>%
  group_by(model, model_lab) %>%
  summarise(years = paste(year, collapse = " -> "),
            penalties = paste(sprintf("%+.3f", penalty), collapse = " -> "),
            first_penalty = penalty[1],
            last_penalty = penalty[n()],
            change = penalty[n()] - penalty[1],
            monotone = all(diff(penalty) > 0) | all(diff(penalty) < 0),
            direction = case_when(all(diff(penalty) > 0) ~ "growing",
                                  all(diff(penalty) < 0) ~ "shrinking",
                                  TRUE ~ "not monotone"),
            .groups = "drop")

print(trend, n = Inf)

# the competing explanation, printed next to it rather than joined onto it:
# training size grew over the same sequence, so a shrinking penalty may be the
# extra data rather than the drift settling down
comp %>% select(year, N_train, N_test) %>% print(n = Inf)

trend %>% select(model_lab, direction, first_penalty, last_penalty, change) %>%
  print(n = Inf)

p_trend <- pen %>%
  filter(metric == TEST_METRIC) %>%
  ggplot(aes(year_num, penalty, colour = model_lab, group = model_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = comp$year_num) +
  scale_colour_manual(values = OKABE, name = NULL) +
  labs(x = "Held-out year", y = paste0("Penalty in ", TEST_METRIC),
       title = "Does predicting forward get harder or easier",
       subtitle = paste0("Training data grows from left to right, so a flat or ",
                         "rising line is drift outpacing\nthe data. Three ",
                         "points, so read the direction and not a slope.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_trend


# =============================================================================
# 09  DOES GENOMICS DRIFT FASTER
# =============================================================================

# Clinical predictors have little to drift: age and comorbidity mean the same
# thing in 2025 as in 2023, so Model_4's penalty is the floor, the cost that any
# model pays for a differently composed later year. Anything Model_8 or Model_9
# loses beyond that floor is the genomic features going out of date, which for
# this organism means the lineage mix changed, new resistance determinants
# arrived, or both.
#
# This connects directly to the lineage folder. If the genomic features are
# largely clone markers, and the circulating clones turn over between years,
# then the temporal penalty and the leave-one-lineage-out penalty are two
# measurements of the same weakness, and they should agree.

pen_wide <- pen %>%
  filter(metric == TEST_METRIC) %>%
  select(year, year_num, model, penalty) %>%
  mutate(model = as.character(model)) %>%
  pivot_wider(names_from = model, values_from = penalty) %>%
  arrange(year_num)

print(pen_wide, n = Inf)

drift <- NULL

for (m in MODELS) {

  if (m == REF_MODEL || !m %in% names(pen_wide)) next

  v <- pen_wide[[m]] - pen_wide[[REF_MODEL]]
  v <- v[!is.na(v)]
  kk <- length(v)
  if (kk < 2) next

  s <- sd(v)
  se_use <- if (NB_CORRECT) s * sqrt(1 / kk + RHO) else s / sqrt(kk)
  tstat  <- if (se_use > 0) mean(v) / se_use else NA_real_
  mult   <- qt(0.975, kk - 1)

  drift <- bind_rows(drift, tibble(
    model = m, k = kk, extra_penalty = mean(v), sd = s,
    lo = mean(v) - mult * se_use, hi = mean(v) + mult * se_use,
    p = if (is.na(tstat)) NA_real_ else 2 * pt(-abs(tstat), df = kk - 1),
    years_worse_than_clinical = sum(v > 0),
    min_extra = min(v), max_extra = max(v)))
}

drift <- drift %>%
  mutate(model_lab = factor(MODEL_LABELS[model], levels = unname(MODEL_LABELS)))

print(drift, n = Inf)

# per-year detail, the more honest display at three years
pen_wide %>%
  pivot_longer(any_of(setdiff(MODELS, REF_MODEL)), names_to = "model",
               values_to = "pen_model") %>%
  left_join(transmute(pen_wide, year, pen_ref = .data[[REF_MODEL]]),
            by = "year") %>%
  mutate(extra = pen_model - pen_ref, model_lab = MODEL_LABELS[model]) %>%
  select(year, model_lab, pen_ref, pen_model, extra) %>%
  arrange(model_lab, year) %>%
  print(n = Inf)

p_drift <- ggplot(drift, aes(extra_penalty, fct_rev(model_lab))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0, linewidth = 0.5, colour = "grey50") +
  geom_point(size = 2.6, colour = OKABE[2]) +
  labs(x = paste0("Loss beyond what the clinical model loses (", TEST_METRIC, ")"),
       y = NULL,
       title = "What the genomic features cost when the year is new",
       subtitle = paste0("Positive means the genomic model ages faster than ",
                         "clinical data alone.\nCompare with the leave-one-",
                         "lineage-out result: they measure related things.")) +
  theme(panel.grid.minor = element_blank())

p_drift


# =============================================================================
# 10  TRAIN AGAINST TEST
# =============================================================================

gap <- bind_rows(train %>% mutate(split = "train"),
                 test  %>% mutate(split = "test")) %>%
  select(split, year, year_num, scheme, model, model_lab, all_of(METRICS)) %>%
  pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value") %>%
  group_by(year, year_num, scheme, model, model_lab, metric, split) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = split, values_from = value) %>%
  mutate(gap = train - test)

gap %>%
  filter(metric == TEST_METRIC) %>%
  arrange(model, scheme, year_num) %>%
  select(model_lab, scheme, year, train, test, gap) %>%
  print(n = Inf)

gap %>%
  filter(metric == TEST_METRIC) %>%
  group_by(model_lab, scheme) %>%
  summarise(mean_train = mean(train), mean_test = mean(test),
            mean_gap = mean(gap), .groups = "drop") %>%
  pivot_wider(names_from = scheme,
              values_from = c(mean_train, mean_test, mean_gap)) %>%
  print(n = Inf)

p_gap <- gap %>%
  filter(metric == TEST_METRIC) %>%
  select(year, scheme, model_lab, Train = train, Test = test) %>%
  pivot_longer(c(Train, Test), names_to = "split", values_to = "auc") %>%
  mutate(split = factor(split, levels = c("Train", "Test"))) %>%
  ggplot(aes(split, auc, group = interaction(year, scheme),
             colour = year, linetype = scheme)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_wrap(~ model_lab, nrow = 1) +
  scale_colour_manual(values = OKABE, name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "dotted"),
                                          c(SCHEME_GEO, SCHEME_CTRL)),
                        name = NULL) +
  labs(x = NULL, y = TEST_METRIC,
       title = "Training against held-out performance, by year and scheme",
       subtitle = "Solid lines predict forward, dotted hold out at random.") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_gap


# =============================================================================
# 11  OUTCOME AND CASE-MIX DRIFT
# =============================================================================

# Whether patients died at a different rate in the held-out year than in the
# years before it. Discrimination is invariant to a base rate shift, but
# calibration is not, and neither is PR AUC, precision, or anything read at a
# fixed threshold. A model whose ROC AUC survives the year can still be
# unusable without recalibration, and this is the block that says so.

shift <- comp %>%
  select(year, year_num, prev_test, prev_train, N_test, N_train) %>%
  mutate(shift = prev_test - prev_train, ratio = prev_test / prev_train) %>%
  arrange(year_num)

print(shift, n = Inf)

# is the death rate moving in one direction across the period
shift$prev_test
diff(shift$prev_test)

pen_shift <- pen %>%
  filter(metric == TEST_METRIC) %>%
  left_join(select(shift, year, shift, ratio), by = "year")

pen_shift %>%
  select(model_lab, year, shift, penalty) %>%
  arrange(model_lab, year) %>%
  print(n = Inf)

p_shift <- ggplot(pen_shift, aes(shift, penalty, colour = year)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 2.8) +
  facet_wrap(~ model_lab, nrow = 1) +
  scale_colour_manual(values = OKABE, name = NULL) +
  labs(x = "Held-out year death rate minus training death rate",
       y = paste0("Penalty in ", TEST_METRIC),
       title = "Does the penalty follow the change in how often patients died",
       subtitle = paste0("Three points per panel. A visible slope points at ",
                         "recalibration rather than\nat the model needing to be ",
                         "refitted.")) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_shift


# =============================================================================
# 12  THE DEPLOYMENT ESTIMATE
# =============================================================================

# The number to put in the abstract. Cross-validation asks how well the model
# predicts a patient drawn from the same period, which is not a question anyone
# will ever face. The last held-out year answers the question that deployment
# actually poses: fitted on everything up to then, asked about the year after.
# It is also the least confounded of the three, having the most history behind
# it.

deployment <- pen %>%
  filter(metric == TEST_METRIC) %>%
  group_by(model, model_lab) %>%
  slice_max(year_num, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(model, model_lab, year, forward_auc = geo, control_auc = ctrl,
         penalty)

print(deployment, n = Inf)

if (!is.null(CV_REFERENCE)) {
  deployment %>%
    mutate(cv_auc = unname(CV_REFERENCE[as.character(model)]),
           optimism_of_cv = cv_auc - forward_auc,
           pct_of_cv_advantage_lost =
             100 * (cv_auc - forward_auc) / (cv_auc - 0.5)) %>%
    print(n = Inf)
}

# the sentence this supports: fitted on data to <year-1> and evaluated on
# <year>, the model achieved <forward_auc>, against <control_auc> when the same
# number of patients were held out at random from across the period.


# =============================================================================
# 13  THE SUMMARY FIGURE
# =============================================================================

p_summary <- pen %>%
  filter(metric == TEST_METRIC) %>%
  ggplot(aes(y = fct_rev(year))) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = ctrl, xend = geo, yend = fct_rev(year)),
               colour = "grey60", linewidth = 0.6,
               arrow = arrow(length = unit(0.07, "inches"), type = "closed")) +
  geom_point(aes(x = ctrl), size = 2.3, colour = OKABE[8]) +
  geom_point(aes(x = geo),  size = 2.7, colour = OKABE[2]) +
  facet_wrap(~ model_lab, nrow = 1) +
  labs(x = TEST_METRIC, y = NULL,
       title = "From a random split to the following year",
       subtitle = paste0("Grey point is the size-matched random control, orange ",
                         "the year predicted forward.\nThe arrow is what time ",
                         "cost. Dashed line is chance.")) +
  theme(panel.grid.minor = element_blank())

p_summary

# the table for the paper
pen %>%
  filter(metric %in% c("ROC_AUC", "PR_AUC", "Balanced_Accuracy", "MCC")) %>%
  mutate(cell = sprintf("%.3f / %.3f (%+.3f)", ctrl, geo, -penalty)) %>%
  select(model_lab, year, metric, cell) %>%
  pivot_wider(names_from = metric, values_from = cell) %>%
  arrange(model_lab, year) %>%
  print(n = Inf)


# =============================================================================
# 14  SAVING
# =============================================================================

# ggsave("tmp_fig_composition.pdf", p_comp,     width =  7, height = 4)
# ggsave("tmp_fig_paired.pdf",      p_paired,   width = 10, height = 5)
# ggsave("tmp_fig_collapse.pdf",    p_collapse, width =  9, height = 5)
# ggsave("tmp_fig_trend.pdf",       p_trend,    width =  8, height = 5)
# ggsave("tmp_fig_drift.pdf",       p_drift,    width =  8, height = 4)
# ggsave("tmp_fig_train_test.pdf",  p_gap,      width = 10, height = 5)
# ggsave("tmp_fig_shift.pdf",       p_shift,    width = 10, height = 4.5)
# ggsave("tmp_fig_summary.pdf",     p_summary,  width = 11, height = 4)
#
# write_csv(match_test,  "tmp_split_matching_test.csv")
# write_csv(match_train, "tmp_split_matching_train.csv")
# write_csv(comp,        "tmp_year_composition.csv")
# write_csv(pen,         "tmp_temporal_penalty.csv")
# write_csv(pen_test,    "tmp_penalty_test.csv")
# write_csv(trend,       "tmp_penalty_direction.csv")
# write_csv(drift,       "tmp_genomic_drift.csv")
# write_csv(deployment,  "tmp_deployment_estimate.csv")
