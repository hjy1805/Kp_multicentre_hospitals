# ============================================================
# CUMULATIVE MORTALITY BY PATHOTYPE AT 30, 60 AND 90 DAYS
# Proportions with exact binomial 95% CIs.
# Convergent hvKp is excluded from the figure (n = 5).
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

base <- "/Users/daneshm/Documents/Kp_KAIMRC"
revision_dir <- file.path(base, "revision")

PATHOTYPE_LEVELS <- c(
  "ESBL/CP-negative non-hvKp",
  "hvKp only",
  "ESBL(+)/CP(+) only",
  "Convergent hvKp"
)

CUTOFFS <- c(30, 60, 90)

PATHOTYPE_COLOURS <- c(
  "ESBL(+)/CP(+) only" = "#D62728",
  "hvKp only" = "#4C78A8",
  "ESBL/CP-negative non-hvKp" = "#BDBDBD"
)


# ============================================================
# 1. GENOTYPE DATA AND PATHOTYPE DEFINITION
# ============================================================

genotype_df <- read_csv(
  file.path(revision_dir, "Combined_df_withST.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    # Complete five-marker hv profile
    hv_status = iuc_bin == 1 &
      iro_bin == 1 &
      peg_bin == 1 &
      rmpA_bin == 1 &
      rmpA2_bin == 1,
    # Kleborate resistance score >= 1 = ESBL and/or carbapenemase
    AMR_status = resistance_score_total >= 1,
    pathotype = case_when(
      !AMR_status & !hv_status ~ "ESBL/CP-negative non-hvKp",
      AMR_status & !hv_status ~ "ESBL(+)/CP(+) only",
      !AMR_status & hv_status ~ "hvKp only",
      AMR_status & hv_status ~ "Convergent hvKp"
    )
  )


# ============================================================
# 2. PATIENT DATA (INPATIENT VISITS ONLY)
# ============================================================

patient_df <- read_csv(
  file.path(revision_dir, "metadata_death_df_1Sep.csv"),
  show_col_types = FALSE
)

table(patient_df$PT_TYPE)

patient_df <- filter(patient_df, PT_TYPE == "IP visit")

patient_df$pathotype <- genotype_df$pathotype[
  match(patient_df$KAUST_ID, genotype_df$KAUST_ID_total)
]


# ============================================================
# 3. DAYS FROM COLLECTION TO DEATH
# ============================================================

parse_date <- function(x) as.Date(x, format = "%d/%m/%Y")

# 1 if death is recorded within `cutoff` days of collection
died_within <- function(days, cutoff) {
  as.integer(!is.na(days) & days >= 0 & days <= cutoff)
}

plot_data <- patient_df %>%
  mutate(
    COLL_DT_Main = parse_date(COLL_DT_Main),
    DTH_DT_Main = parse_date(DTH_DT_Main),
    days_to_death = as.numeric(DTH_DT_Main - COLL_DT_Main),
    pathotype = factor(pathotype, levels = PATHOTYPE_LEVELS),
    death_30 = died_within(days_to_death, 30),
    death_60 = died_within(days_to_death, 60),
    death_90 = died_within(days_to_death, 90)
  )


# ============================================================
# 4. LONG FORMAT
# ============================================================

mortality_long <- plot_data %>%
  select(pathotype, death_30, death_60, death_90) %>%
  pivot_longer(
    cols = starts_with("death_"),
    names_to = "Time",
    values_to = "Death"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = paste0("death_", CUTOFFS),
      labels = as.character(CUTOFFS)
    )
  )


# ============================================================
# 5. MORTALITY RATE + EXACT 95% BINOMIAL CI
# ============================================================

mortality_summary <- mortality_long %>%
  group_by(pathotype, Time) %>%
  summarise(
    N = n(),
    Deaths = sum(Death),
    Mortality = Deaths / N,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(ci = list(binom.test(Deaths, N)$conf.int)) %>%
  ungroup() %>%
  mutate(
    lower = vapply(ci, function(x) x[1], numeric(1)),
    upper = vapply(ci, function(x) x[2], numeric(1))
  ) %>%
  select(-ci) %>%
  filter(pathotype != "Convergent hvKp") %>%
  droplevels()

print(mortality_summary)


# ============================================================
# 6. FIGURE
# ============================================================

dodge <- position_dodge(width = 0.8)

p_mortality <- ggplot(
  mortality_summary,
  aes(x = Time, y = Mortality, fill = pathotype)
) +
  geom_col(position = dodge, width = 0.7) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = dodge,
    width = 0.15,
    linewidth = 0.6
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(values = PATHOTYPE_COLOURS) +
  labs(
    x = "Days since Collection",
    y = "Mortality Rate",
    fill = "Type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 13, hjust = 1),
    axis.title.x = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 15)
  )

p_mortality

# ============================================================
# 5B. STATISTICAL COMPARISON OF MORTALITY BETWEEN PATHOTYPES
#     Separate tests at 30, 60 and 90 days
# ============================================================

# Exclude convergent hvKp, consistent with the figure
mortality_test_data <- mortality_long %>%
  filter(pathotype != "Convergent hvKp") %>%
  droplevels()


# ------------------------------------------------------------
# Global Fisher's exact test at each time point
# ------------------------------------------------------------

global_tests <- lapply(levels(mortality_test_data$Time), function(t) {
  
  dat_t <- mortality_test_data %>%
    filter(Time == t)
  
  tab <- table(dat_t$pathotype, dat_t$Death)
  
  test <- fisher.test(tab)
  
  data.frame(
    Time = t,
    P_value = test$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    P_adj_BH = p.adjust(P_value, method = "BH")
  )

print(global_tests)


# ------------------------------------------------------------
# Show contingency tables as well
# ------------------------------------------------------------

for (t in levels(mortality_test_data$Time)) {
  
  cat("\n====================================\n")
  cat("Mortality at", t, "days\n")
  cat("====================================\n")
  
  tab <- mortality_test_data %>%
    filter(Time == t) %>%
    with(table(pathotype, Death))
  
  print(tab)
  
  cat("\nFisher's exact test:\n")
  print(fisher.test(tab))
}


# ============================================================
# 5C. PAIRWISE FISHER TESTS
# ============================================================

pathotypes <- levels(mortality_test_data$pathotype)

comparisons <- combn(pathotypes, 2, simplify = FALSE)

pairwise_results <- lapply(levels(mortality_test_data$Time), function(t) {
  
  dat_t <- mortality_test_data %>%
    filter(Time == t)
  
  res <- lapply(comparisons, function(comp) {
    
    dat_pair <- dat_t %>%
      filter(pathotype %in% comp) %>%
      droplevels()
    
    tab <- table(dat_pair$pathotype, dat_pair$Death)
    
    test <- fisher.test(tab)
    
    data.frame(
      Time = t,
      Group1 = comp[1],
      Group2 = comp[2],
      P_value = test$p.value
    )
  }) %>%
    bind_rows()
  
  # Adjust the three pairwise comparisons within each time point
  res %>%
    mutate(
      P_adj_BH = p.adjust(P_value, method = "BH")
    )
  
}) %>%
  bind_rows()

print(pairwise_results)

# ============================================================
# 7. STATISTICAL TESTS OF MORTALITY BETWEEN PATHOTYPES
#    Fisher's exact test at 30, 60 and 90 days
# ============================================================

# Exclude convergent hvKp, consistent with the figure
test_data <- mortality_long %>%
  filter(pathotype != "Convergent hvKp") %>%
  droplevels()


# ============================================================
# 7A. GLOBAL FISHER'S EXACT TEST
#     Tests whether mortality differs across ANY of the
#     three pathotypes at each time point
# ============================================================

global_results <- lapply(levels(test_data$Time), function(tt) {
  
  dat <- test_data %>%
    filter(Time == tt)
  
  tab <- table(dat$pathotype, dat$Death)
  
  ft <- fisher.test(tab)
  
  data.frame(
    Time = tt,
    P_value = ft$p.value
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    P_adj_BH = p.adjust(P_value, method = "BH"),
    
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

cat("\n\nGLOBAL FISHER TESTS\n")
print(global_results)


# ============================================================
# 7B. PAIRWISE FISHER'S EXACT TESTS
#     All three pairwise comparisons at each time point
# ============================================================

pathotype_names <- levels(test_data$pathotype)

pair_list <- combn(
  pathotype_names,
  2,
  simplify = FALSE
)


pairwise_results <- lapply(levels(test_data$Time), function(tt) {
  
  dat_time <- test_data %>%
    filter(Time == tt)
  
  tmp <- lapply(pair_list, function(pair) {
    
    dat_pair <- dat_time %>%
      filter(pathotype %in% pair) %>%
      droplevels()
    
    tab <- table(
      dat_pair$pathotype,
      dat_pair$Death
    )
    
    ft <- fisher.test(tab)
    
    data.frame(
      Time = tt,
      Group1 = pair[1],
      Group2 = pair[2],
      P_value = ft$p.value
    )
  })
  
  tmp <- bind_rows(tmp)
  
  # BH correction separately within each time point
  tmp <- tmp %>%
    mutate(
      P_adj_BH = p.adjust(P_value, method = "BH"),
      
      Significance = case_when(
        P_adj_BH < 0.001 ~ "***",
        P_adj_BH < 0.01  ~ "**",
        P_adj_BH < 0.05  ~ "*",
        TRUE             ~ "ns"
      )
    )
  
  return(tmp)
})

pairwise_results <- bind_rows(pairwise_results)

cat("\n\nPAIRWISE FISHER TESTS\n")
print(pairwise_results)


# ============================================================
# 7C. OPTIONAL: FORMAT P VALUES FOR EASY READING
# ============================================================

pairwise_results_formatted <- pairwise_results %>%
  mutate(
    P_value = signif(P_value, 4),
    P_adj_BH = signif(P_adj_BH, 4)
  )

print(pairwise_results_formatted)


# ============================================================
# 7D. ONLY SIGNIFICANT PAIRWISE RESULTS
# ============================================================

significant_results <- pairwise_results %>%
  filter(P_adj_BH < 0.05)

cat("\n\nSIGNIFICANT PAIRWISE COMPARISONS\n")
print(significant_results)
#Update SNP numbr 

