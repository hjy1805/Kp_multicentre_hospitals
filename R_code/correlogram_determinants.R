# ============================================================
# PAIRWISE CORRELATION OF RESISTANCE / VIRULENCE DETERMINANTS
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)
library(purrr)

MIN_COUNT <- 10


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
  "hv",
  "mcr",
  "Col_mutations",
  "PmrB",
  "MgrB",
  "Tet_acquired",
  "MLS_acquired",
  "Flq_acquired",
  "Flq_mutations",
  "AGly_acquired",
  "SHV_mutations",
  "Bla_ESBL_acquired",
  "OXY",
  "CTX.M",
  "Omp_mutations",
  "OmpK36TD",
  "OmpK36GD",
  "OmpK36.",
  "OmpK35.",
  "Bla_Carb_acquired",
  "NDM.1",
  "NDM.5",
  "OXA.48",
  "OXA.232"
)

# Every determinant labels as itself except these
LABEL_OVERRIDES <- c(
  "hv" = "hvKp",
  "CTX.M" = "CTX-M",
  "OmpK36." = "OmpK36",
  "OmpK35." = "OmpK35",
  "NDM.1" = "NDM-1",
  "NDM.5" = "NDM-5",
  "OXA.48" = "OXA-48",
  "OXA.232" = "OXA-232"
)

gene_labels <- setNames(determinants, determinants)
gene_labels[names(LABEL_OVERRIDES)] <- LABEL_OVERRIDES


# ============================================================
# 3. REMOVE RARE FEATURES
# Both classes must be seen at least MIN_COUNT times
# ============================================================

cor_data <- kleborate_encoded_df %>%
  select(all_of(determinants)) %>%
  mutate(across(everything(), as.numeric))

feature_counts <- data.frame(
  Gene = names(cor_data),
  n_positive = sapply(cor_data, function(x) sum(x == 1, na.rm = TRUE)),
  n_negative = sapply(cor_data, function(x) sum(x == 0, na.rm = TRUE)),
  row.names = NULL
)

eligible_genes <- feature_counts %>%
  filter(n_positive >= MIN_COUNT, n_negative >= MIN_COUNT) %>%
  pull(Gene)

cat("\nGenes included in correlogram:\n")
print(eligible_genes)

cat("\nGenes excluded because rare:\n")
feature_counts %>%
  filter(n_positive < MIN_COUNT | n_negative < MIN_COUNT) %>%
  print()

cor_data <- cor_data %>%
  select(all_of(eligible_genes))


# ============================================================
# 4. CORRELATION + P VALUE FOR ONE PAIR
# ============================================================

cor_test_pair <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  # Cannot calculate a correlation without variation
  if (length(x) < 3 ||
      length(unique(x)) < 2 ||
      length(unique(y)) < 2) {
    return(c(correlation = NA_real_, p.value = NA_real_))
  }

  test <- suppressWarnings(cor.test(x, y, method = "pearson"))

  c(
    correlation = unname(test$estimate),
    p.value = test$p.value
  )
}


# ============================================================
# 5. PAIRWISE CORRELATIONS (LOWER TRIANGLE ONLY)
# ============================================================

cor_results <- expand.grid(
  Gene1 = eligible_genes,
  Gene2 = eligible_genes,
  stringsAsFactors = FALSE
) %>%
  mutate(
    position1 = match(Gene1, eligible_genes),
    position2 = match(Gene2, eligible_genes)
  ) %>%
  filter(position1 <= position2)

pair_results <- map2_dfr(
  cor_results$Gene1,
  cor_results$Gene2,
  function(g1, g2) {
    result <- cor_test_pair(cor_data[[g1]], cor_data[[g2]])
    data.frame(
      correlation = unname(result["correlation"]),
      p.value = unname(result["p.value"])
    )
  }
)

cor_results <- bind_cols(cor_results, pair_results)


# ============================================================
# 6. SIGNIFICANCE AND PLOT LABELS
# ============================================================

plot_order <- gene_labels[eligible_genes]

cor_results <- cor_results %>%
  mutate(
    Significant = case_when(
      Gene1 == Gene2 ~ TRUE,
      !is.na(p.value) & p.value < 0.05 ~ TRUE,
      TRUE ~ FALSE
    ),
    Label1 = factor(gene_labels[Gene1], levels = plot_order),
    Label2 = factor(gene_labels[Gene2], levels = rev(plot_order))
  )


# ============================================================
# 7. CORRELOGRAM
# ============================================================

p_corr <- ggplot(cor_results, aes(x = Label1, y = Label2)) +
  # Grey background for non-significant correlations
  geom_tile(fill = "grey85", colour = "white", linewidth = 0.3) +
  # Colour only significant correlations
  geom_tile(
    data = filter(cor_results, Significant),
    aes(fill = correlation),
    colour = "white",
    linewidth = 0.3
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  labs(title = "Pairwise correlation", x = NULL, y = NULL) +
  coord_fixed() +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  )


# ============================================================
# 8. SHOW CORRELOGRAM
# ============================================================

p_corr

#Methods 

#Captions

#Mantel test 
#plasmid analysis 

#results ML+survival improve 







