# ============================================================
# ICU
# ROC WITH 95% CI RIBBONS + DELONG + DECISION CURVE
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)
library(pROC)

if (!requireNamespace("dcurves", quietly = TRUE)) {
  install.packages("dcurves")
}
library(dcurves)


# ============================================================
# 1. PATHS AND MODELS
# ============================================================

base <- "/Users/daneshm/Documents/Kp_KAIMRC"
output_dir <- file.path(base, "ICU_ROC_DCA_DeLong")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Display name -> postfix of the XGBoost run, in plotting order
MODEL_POSTFIX <- c(
  "Clinical" = "_1_ICU",
  "Clinical + Genomic" = "_5_ICU",
  "Genomic" = "_6_ICU"
)

MODEL_ORDER <- names(MODEL_POSTFIX)

# Column names used inside `comparison` (syntactic, no backticks)
SCORE_COLS <- c(
  "Clinical" = "Clinical",
  "Clinical + Genomic" = "Clinical_Genomic",
  "Genomic" = "Genomic"
)


# ============================================================
# 2. READ OOF PREDICTIONS
# ============================================================

read_oof <- function(postfix, score_col) {
  path <- file.path(
    base,
    paste0("preds_xgb", postfix),
    paste0("preds_test_oof", postfix, ".csv")
  )

  d <- read_csv(path, show_col_types = FALSE) %>%
    select(Name, index, fold, y_true, y_score)

  names(d)[names(d) == "y_score"] <- score_col
  d
}

clinical <- read_oof(MODEL_POSTFIX[["Clinical"]], "Clinical")
combined <- read_oof(
  MODEL_POSTFIX[["Clinical + Genomic"]], "Clinical_Genomic"
)
Genomic <- read_oof(MODEL_POSTFIX[["Genomic"]], "Genomic")


# ============================================================
# 3. MERGE
# ============================================================

join_keys <- c("Name", "index", "fold")

comparison <- clinical %>%
  inner_join(
    combined %>% rename(y_true_combined = y_true),
    by = join_keys
  ) %>%
  inner_join(
    Genomic %>% rename(y_true_Genomic = y_true),
    by = join_keys
  )


# ============================================================
# 4. CHECKS
# ============================================================

stopifnot(nrow(comparison) == nrow(clinical))
stopifnot(all(comparison$y_true == comparison$y_true_combined))
stopifnot(all(comparison$y_true == comparison$y_true_Genomic))
stopifnot(!anyNA(comparison[, unname(SCORE_COLS)]))

write_csv(
  comparison,
  file.path(output_dir, "OOF_predictions_all_ICU_models.csv")
)


# ============================================================
# 5. ROC OBJECTS, AUC AND 95% CI
# ============================================================

make_roc <- function(scores) {
  roc(
    comparison$y_true,
    scores,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

rocs <- lapply(SCORE_COLS, function(col) make_roc(comparison[[col]]))
names(rocs) <- MODEL_ORDER

auc_value <- vapply(rocs, function(r) as.numeric(auc(r)), numeric(1))

auc_ci <- vapply(
  rocs,
  function(r) as.numeric(ci.auc(r, method = "delong"))[c(1, 3)],
  numeric(2)
)

auc_table <- tibble(
  Model = MODEL_ORDER,
  AUC = auc_value[MODEL_ORDER],
  CI_low = auc_ci[1, MODEL_ORDER],
  CI_high = auc_ci[2, MODEL_ORDER]
)

print(auc_table)

write_csv(auc_table, file.path(output_dir, "ICU_AUC_95CI.csv"))


# ============================================================
# 6. PAIRED DELONG TESTS
# ============================================================

COMPARISONS <- list(
  c("Clinical", "Clinical + Genomic"),
  c("Clinical", "Genomic"),
  c("Genomic", "Clinical + Genomic")
)

delong_table <- bind_rows(lapply(COMPARISONS, function(pair) {
  m1 <- pair[1]
  m2 <- pair[2]

  test <- roc.test(
    rocs[[m1]],
    rocs[[m2]],
    method = "delong",
    paired = TRUE
  )

  tibble(
    Comparison = paste(m1, "vs", m2),
    AUC_model_1 = auc_value[[m1]],
    AUC_model_2 = auc_value[[m2]],
    AUC_difference = auc_value[[m2]] - auc_value[[m1]],
    Z = as.numeric(test$statistic),
    P_value = test$p.value
  )
})) %>%
  mutate(P_BH = p.adjust(P_value, method = "BH"))

print(delong_table, n = Inf)

write_csv(delong_table, file.path(output_dir, "ICU_DeLong_tests.csv"))


# ============================================================
# 7. ROC 95% CI RIBBONS
# Bootstrap sensitivity CI at fixed specificity values
# ============================================================

specificity_grid <- seq(0, 1, length.out = 101)

get_roc_ci <- function(roc_object, model_name) {
  ci_matrix <- as.matrix(
    ci.se(
      roc_object,
      specificities = specificity_grid,
      boot.n = 2000,
      conf.level = 0.95,
      progress = "none"
    )
  )

  tibble(
    FPR = 1 - specificity_grid,
    lower = ci_matrix[, 1],
    middle = ci_matrix[, 2],
    upper = ci_matrix[, 3],
    Model = model_name
  ) %>%
    arrange(FPR)
}

set.seed(12345)

roc_ci_df <- bind_rows(
  lapply(MODEL_ORDER, function(m) get_roc_ci(rocs[[m]], m))
)


# ============================================================
# 8. ROC CURVE DATA
# ============================================================

roc_to_df <- function(roc_object, model_name) {
  tibble(
    FPR = 1 - roc_object$specificities,
    TPR = roc_object$sensitivities,
    Model = model_name
  ) %>%
    arrange(FPR)
}

roc_df <- bind_rows(
  lapply(MODEL_ORDER, function(m) roc_to_df(rocs[[m]], m))
)


# ============================================================
# 9. LEGEND LABELS WITH AUC + 95% CI
# ============================================================

model_labels <- sprintf(
  "%s: AUC %.3f (95%% CI %.3f–%.3f)",
  MODEL_ORDER,
  auc_value[MODEL_ORDER],
  auc_ci[1, MODEL_ORDER],
  auc_ci[2, MODEL_ORDER]
)

model_labels <- sprintf(
  MODEL_ORDER
)
names(model_labels) <- MODEL_ORDER


# ============================================================
# 10. SHARED PLOT THEME
# ============================================================

plot_theme <- list(
  theme_classic(base_size = 14),
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )
)

save_plot <- function(plot, stem, width = 9, height = 7) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 600
  )
}


# ============================================================
# 11. ROC PLOT
# CURVES + 95% CI RIBBONS + LEGEND
# ============================================================
model_colors <- c(
  "Clinical" = "red",
  "Genomic" = "#B8860B",          # dark golden yellow
  "Clinical + Genomic" = "orange"
)

p_roc <- ggplot() +
  geom_ribbon(
    data = roc_ci_df,
    aes(x = FPR, ymin = lower, ymax = upper, fill = Model),
    alpha = 0.15,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = roc_df,
    aes(x = FPR, y = TPR, colour = Model),
    linewidth = 1.2
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.8
  ) +
  # Keep curve and ribbon colours consistent
  scale_fill_manual(
    values = model_colors,
    limits = MODEL_ORDER
  ) +
  scale_colour_manual(
    values = model_colors,
    limits = MODEL_ORDER,
    labels = model_labels
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(x = "1 − Specificity", y = "Sensitivity", colour = NULL) +
  coord_equal() +
  plot_theme

print(p_roc)

save_plot(p_roc, "ROC_ICU_95CI_ribbon")


# ============================================================
# 12. DECISION CURVE ANALYSIS
# ============================================================

dca_data <- comparison %>%
  transmute(
    outcome = as.numeric(y_true),
    Clinical = Clinical,
    `Clinical + Genomic` = Clinical_Genomic,
    Genomic = Genomic
  )

# Scores must be probabilities
stopifnot(
  all(vapply(
    dca_data[MODEL_ORDER],
    function(x) all(x >= 0 & x <= 1),
    logical(1)
  ))
)

thresholds_dca <- seq(0.01, 0.50, by = 0.01)

dca_result <- dca(
  outcome ~ Clinical + `Clinical + Genomic` + Genomic,
  data = dca_data,
  thresholds = thresholds_dca
)


# ============================================================
# 13. DECISION CURVE PLOT
# ============================================================

p_dca <- plot(dca_result, smooth = FALSE) +
  labs(
    x = "Threshold probability",
    y = "Net benefit",
    colour = NULL
  ) +
  plot_theme

print(p_dca)

save_plot(p_dca, "Decision_Curve_ICU_final")

write_csv(
  as_tibble(dca_result),
  file.path(output_dir, "ICU_Decision_Curve_net_benefit.csv")
)


# ============================================================
# 14. FINISHED
# ============================================================

cat("\n========================================\n")
cat("ICU ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "\nROC with 95% CI ribbon:\n",
  file.path(output_dir, "ROC_ICU_95CI_ribbon.pdf"),
  "\n"
)
cat(
  "\nDecision curve:\n",
  file.path(output_dir, "Decision_Curve_ICU_final.pdf"),
  "\n"
)

cat("\nAUC results:\n")
print(auc_table, n = Inf)

cat("\nDeLong results:\n")
print(delong_table, n = Inf)




# ============================================================
# 13. DECISION CURVE PLOT
# ============================================================

dca_colors <- c(
  "Treat All" = "black",
  "Treat None" = "grey50",
  "Clinical" = "red",
  "Clinical + Genomic" = "orange",
  "Genomic" = "#B8860B"
)

p_dca <- plot(dca_result, smooth = FALSE)

# Standardise model names generated internally by dcurves
p_dca$data <- p_dca$data %>%
  mutate(
    label = case_when(
      label %in% c(
        "Clinical + genomics",
        "Clinical + Genomics",
        "Clinical + genomic"
      ) ~ "Clinical + Genomic",
      
      label %in% c(
        "Genomics",
        "genomics",
        "genomic"
      ) ~ "Genomic",
      
      TRUE ~ label
    )
  )

p_dca <- p_dca +
  scale_colour_manual(
    values = dca_colors,
    breaks = c(
      "Treat All",
      "Treat None",
      "Clinical",
      "Clinical + Genomic",
      "Genomic"
    )
  ) +
  labs(
    x = "Threshold probability",
    y = "Net benefit",
    colour = NULL
  ) +
  plot_theme

print(p_dca)

save_plot(p_dca, "Decision_Curve_ICU_final")
