# ============================================================
# ICU ADMISSION
# One adjusted model per resistance / virulence determinant
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)

MIN_POSITIVE <- 10
MIN_NEGATIVE <- 10

# Stability limits for the adjusted OR confidence interval
MIN_CI <- 0.1
MAX_CI <- 10

COVARIATES <- c("AGE", "CHARLSON", "BMI", "GENDER", "SOURCE", "RGN")


# ============================================================
# 1. LOAD DATA
# ============================================================

kleborate_encoded_df <- read_csv(
  "/Users/daneshm/Documents/Kp_KAIMRC/revision/Kleborate_OD_calculation.csv",
  show_col_types = FALSE
)


# ============================================================
# 2. DETERMINANTS AND LABELS
# ============================================================

determinants <- c(
  # Hypervirulence
  "hv",

  # Colistin resistance
  "mcr",
  "Col_mutations",
  "PmrB",
  "MgrB",

  # Other resistance determinants
  "Tet_acquired",
  "MLS_acquired",
  "Flq_acquired",
  "Flq_mutations",
  "AGly_acquired",
  "SHV_mutations",
  "Bla_ESBL_acquired",
  "OXY",
  "CTX.M",

  # Porins
  "Omp_mutations",
  "OmpK36TD",
  "OmpK36GD",
  "OmpK36.",
  "OmpK35.",

  # Carbapenemases
  "Bla_Carb_acquired",
  "KPC.2",
  "OXA.181",
  "NDM.1",
  "NDM.5",
  "OXA.48",
  "OXA.232"
)

# Every determinant labels as itself except these
LABEL_OVERRIDES <- c(
  "hv" = "5-marker hvKp",
  "CTX.M" = "CTX-M",
  "OmpK36." = "OmpK36",
  "OmpK35." = "OmpK35",
  "KPC.2" = "KPC-2",
  "OXA.181" = "OXA-181",
  "NDM.1" = "NDM-1",
  "NDM.5" = "NDM-5",
  "OXA.48" = "OXA-48",
  "OXA.232" = "OXA-232"
)

gene_labels <- setNames(determinants, determinants)
gene_labels[names(LABEL_OVERRIDES)] <- LABEL_OVERRIDES


# ============================================================
# 3. PREPARE ICU DATA
# ============================================================

icu_data <- kleborate_encoded_df %>%
  mutate(
    across(c(ICU_Admission, AGE, BMI, CHARLSON), as.numeric),
    across(c(GENDER, SOURCE, RGN), factor)
  )


# ============================================================
# 4. RESULT ROW HELPER
# ============================================================

result_row <- function(gene,
                       n_obs,
                       n_positive,
                       n_negative,
                       estimate = NA_real_,
                       conf.low = NA_real_,
                       conf.high = NA_real_,
                       p.value = NA_real_,
                       status) {
  data.frame(
    Gene = gene,
    N = n_obs,
    n_positive = n_positive,
    n_negative = n_negative,
    estimate = estimate,
    conf.low = conf.low,
    conf.high = conf.high,
    p.value = p.value,
    Status = status
  )
}


# ============================================================
# 5. ADJUSTED ICU MODEL FOR ONE DETERMINANT
# ============================================================

run_icu_model <- function(gene) {
  model_df <- icu_data %>%
    select(ICU_Admission, all_of(gene), all_of(COVARIATES)) %>%
    filter(complete.cases(.))

  model_df[[gene]] <- as.numeric(model_df[[gene]])

  n_positive <- sum(model_df[[gene]] == 1, na.rm = TRUE)
  n_negative <- sum(model_df[[gene]] == 0, na.rm = TRUE)

  exclude <- function(status) {
    result_row(
      gene, nrow(model_df), n_positive, n_negative, status = status
    )
  }

  # Rare determinants
  if (n_positive < MIN_POSITIVE || n_negative < MIN_NEGATIVE) {
    return(exclude("Excluded: rare"))
  }

  # Outcome must vary
  if (length(unique(model_df$ICU_Admission)) < 2) {
    return(exclude("Excluded: no outcome variation"))
  }

  fit <- tryCatch(
    suppressWarnings(
      glm(
        reformulate(c(gene, COVARIATES), response = "ICU_Admission"),
        data = model_df,
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(exclude("Excluded: model failed"))
  }

  coef_table <- summary(fit)$coefficients

  if (!gene %in% rownames(coef_table)) {
    return(exclude("Excluded: coefficient unavailable"))
  }

  beta <- coef_table[gene, "Estimate"]
  se <- coef_table[gene, "Std. Error"]

  estimate <- exp(beta)
  conf_low <- exp(beta - 1.96 * se)
  conf_high <- exp(beta + 1.96 * se)

  # || short-circuits, so the bound checks never see an NA
  unstable <- !all(is.finite(c(estimate, conf_low, conf_high))) ||
    conf_low < MIN_CI ||
    conf_high > MAX_CI

  result_row(
    gene,
    nrow(model_df),
    n_positive,
    n_negative,
    estimate = estimate,
    conf.low = conf_low,
    conf.high = conf_high,
    p.value = coef_table[gene, "Pr(>|z|)"],
    status = if (unstable) "Excluded: unstable" else "Included"
  )
}


# ============================================================
# 6. RUN ALL ICU MODELS
# ============================================================

icu_all <- map_dfr(determinants, run_icu_model)

SUMMARY_COLS <- c(
  "Gene", "n_positive", "n_negative",
  "estimate", "conf.low", "conf.high", "Status"
)

cat("\nALL DETERMINANTS:\n")
print(select(icu_all, all_of(SUMMARY_COLS)))

cat("\nEXCLUDED DETERMINANTS:\n")
icu_all %>%
  filter(Status != "Included") %>%
  select(all_of(SUMMARY_COLS)) %>%
  print()


# ============================================================
# 7. BH CORRECTION AND PLOT ORDER
# ============================================================

included_order <- determinants[determinants %in% icu_all$Gene[
  icu_all$Status == "Included"
]]

plot_order <- gene_labels[included_order]

icu_results <- icu_all %>%
  filter(Status == "Included") %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    Significant = p_adj < 0.055,
    Label = factor(gene_labels[Gene], levels = rev(plot_order))
  )

cat("\nFINAL ICU RESULTS:\n")
icu_results %>%
  select(
    Gene, n_positive, N, estimate, conf.low, conf.high,
    p.value, p_adj, Significant
  ) %>%
  print()


# ============================================================
# 8. ICU FOREST PLOT
# ============================================================

SIG_COLOURS <- c("FALSE" = "#4C78A8", "TRUE" = "red")

p_icu <- ggplot(icu_results, aes(x = estimate, y = Label)) +
  # OR = 1
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high, colour = Significant),
    height = 0.15,
    linewidth = 0.7
  ) +
  geom_point(aes(colour = Significant), size = 2.7) +
  # BH-significant results = red
  scale_colour_manual(values = SIG_COLOURS) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = label_number(accuracy = 0.01)
  ) +
  coord_cartesian(xlim = c(0.2, 8)) +
  labs(
    title = "ICU\nadmission",
    x = "Odds Ratio",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


# ============================================================
# 9. SHOW ICU PLOT
# ============================================================

p_icu
