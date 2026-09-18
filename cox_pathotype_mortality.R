# ============================================================
# IN-HOSPITAL MORTALITY BY PATHOTYPE
# Cox models, each exposure against the Non-AMR/Non-hvKp
# reference, clustered on patient ID.
# ============================================================

library(readr)
library(dplyr)
library(survival)
library(broom)
library(forestmodel)

base <- "/Users/daneshm/Documents/Kp_KAIMRC"
revision_dir <- file.path(base, "revision")

REFERENCE_PATHOTYPE <- "ESBL/CP-negative non-hvKp"

COX_COVARIATES <- c(
  "AGE", "SOURCE", "BMI", "CHARLSON", "GENDER", "RGN", "NUM_ADM"
)

CONFOUNDER_COLS <- COX_COVARIATES


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

table(genotype_df$pathotype)


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
# 3. CLINICAL CONFOUNDERS
# ============================================================

Confounders <- read_csv(
  file.path(base, "files_kleborate", "kleborate_results_totall.csv"),
  show_col_types = FALSE
)

# Negative ages are treated as missing, so those isolates drop
# out of the models at the na.omit() step below
Confounders$AGE <- ifelse(Confounders$AGE < 0, NA, Confounders$AGE)

confounder_idx <- match(patient_df$KAUST_ID, Confounders$strain)
patient_df[CONFOUNDER_COLS] <- Confounders[
  confounder_idx, CONFOUNDER_COLS
]


# ============================================================
# 4. SURVIVAL DATA
#
# Follow-up starts at isolate collection and ends at death
# (event) or discharge alive (censored).
# ============================================================

parse_date <- function(x) as.Date(x, format = "%d/%m/%Y")

cox_data <- patient_df %>%
  mutate(
    collection_date = parse_date(COLL_DT_Main),
    discharge_date = parse_date(DS_DT_Main),
    death_date = parse_date(DTH_DT_Main),
    event = ifelse(HS_DTH_STATUS_Main == "Deceased", 1, 0),
    end_date = as.Date(
      ifelse(event == 1, death_date, discharge_date),
      origin = "1970-01-01"
    ),
    time = as.numeric(end_date - collection_date)
  ) %>%
  filter(
    !is.na(ID),
    !is.na(time),
    time >= 0,
    !is.na(pathotype)
  )


# ============================================================
# 5. BASIC CHECKS
# ============================================================

table(cox_data$pathotype, useNA = "ifany")

cat("Total isolates:", nrow(cox_data), "\n")
cat("Unique patients:", n_distinct(cox_data$ID), "\n")

cox_data %>%
  count(ID) %>%
  summarise(
    patients = n(),
    patients_with_multiple_isolates = sum(n > 1),
    max_isolates_per_patient = max(n)
  ) %>%
  print()


# ============================================================
# 6. BUILD ONE EXPOSURE-VERSUS-REFERENCE DATASET
# ============================================================

build_cox_data <- function(exposure_level) {
  cox_data %>%
    filter(
      pathotype %in% c(REFERENCE_PATHOTYPE, exposure_level),
      # Sparse category: one observation and zero deaths
      SOURCE != "Others"
    ) %>%
    mutate(
      Exposure = factor(
        pathotype,
        levels = c(REFERENCE_PATHOTYPE, exposure_level),
        labels = c("ESBL/CP-negative non-hvKp", exposure_level)
      ),
      SOURCE = relevel(factor(SOURCE), ref = "Blood"),
      GENDER = relevel(factor(GENDER), ref = "F"),
      RGN = relevel(factor(RGN), ref = "Central")
    ) %>%
    select(ID, time, event, Exposure, all_of(COX_COVARIATES)) %>%
    na.omit() %>%
    droplevels()
}


# ============================================================
# 7. FIT UNADJUSTED AND ADJUSTED COX MODELS
# ============================================================

run_cox_analysis <- function(exposure_level, label) {
  dat <- build_cox_data(exposure_level)

  cat("\n============================================\n")
  cat(label, "\n")
  cat("============================================\n")
  cat("Isolates:", nrow(dat), "\n")
  cat("Patients:", n_distinct(dat$ID), "\n")
  cat("Deaths:", sum(dat$event), "\n\n")

  print(table(dat$Exposure))
  print(table(dat$Exposure, dat$event))
  print(table(dat$SOURCE, dat$event))

  unadjusted <- coxph(
    Surv(time, event) ~ Exposure,
    data = dat,
    cluster = ID
  )

  adjusted <- coxph(
    reformulate(
      c("Exposure", COX_COVARIATES),
      response = "Surv(time, event)"
    ),
    data = dat,
    cluster = ID
  )

  cat("\n--- Unadjusted ---\n")
  print(summary(unadjusted))

  cat("\n--- Adjusted ---\n")
  print(summary(adjusted))

  results <- tidy(adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      HR_95CI = sprintf(
        "%.2f (%.2f–%.2f)", estimate, conf.low, conf.high
      ),
      p = if_else(
        p.value < 0.001,
        "<0.001",
        sprintf("%.3f", p.value)
      )
    )

  # Proportional hazards assumption
  ph <- cox.zph(adjusted)
  cat("\n--- Proportional hazards ---\n")
  print(ph)

  list(
    data = dat,
    unadjusted = unadjusted,
    adjusted = adjusted,
    results = results,
    ph = ph,
    exposure = filter(results, grepl("^Exposure", term)) %>%
      select(term, HR_95CI, p)
  )
}


# ============================================================
# 8. RUN BOTH COMPARISONS
# ============================================================

hv <- run_cox_analysis("hvKp only", "hvKp only vs Non-AMR/Non-hvKp")

amr <- run_cox_analysis(
  "ESBL(+)/CP(+) only",
  "ESBL(+)/CP(+) only vs Non-AMR/Non-hvKp"
)


# ============================================================
# 9. DIAGNOSTIC AND FOREST PLOTS
# ============================================================

plot(hv$ph)
plot(amr$ph)

forest_model(hv$adjusted)
forest_model(amr$adjusted)


# ============================================================
# 10. EXPOSURE RESULTS
# ============================================================

cat("\nhvKp exposure:\n")
print(hv$exposure)

cat("\nESBL/CP exposure:\n")
print(amr$exposure)
