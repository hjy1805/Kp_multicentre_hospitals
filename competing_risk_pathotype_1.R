# ============================================================
# COMPETING-RISK ANALYSIS BY PATHOTYPE
# In-hospital death (event 1) with discharge alive (event 2)
# as the competing event, each exposure against the
# Non-AMR/Non-hvKp reference group.
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)
library(cmprsk)
library(survminer)

base <- "/Users/daneshm/Documents/Kp_KAIMRC"
revision_dir <- file.path(base, "revision")

REFERENCE_PATHOTYPE <- "Non-AMR/Non-hvKp"
REFERENCE_LABEL <- "Other"

# cuminc() names its output elements paste(group, event), and
# ggcompetingrisks() splits those on whitespace. A group label
# containing a space is therefore torn in half and its second
# word is drawn as a phantom event. Group labels must have no
# spaces; these are used for modelling, PLOT_TITLES for display.
GROUP_LABELS <- c(
  "hvKp only" = "hvKp",
  "ESBL(+)/CP(+) only" = "ESBL_CP"
)

EVENT_LABELS <- c(
  "1" = "In-hospital death",
  "2" = "Discharged alive"
)

FOLLOW_UP_DAYS <- 90
REPORT_TIMES <- c(30, 60, 90)


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
      !AMR_status & !hv_status ~ "Non-AMR/Non-hvKp",
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
# 3. COMPETING-RISK DATA
#
# status_cr: 0 = censored, 1 = in-hospital death,
#            2 = discharged alive
# Follow-up starts at K. pneumoniae isolation.
# ============================================================

parse_date <- function(x) as.Date(x, format = "%d/%m/%Y")

cr_data <- patient_df %>%
  mutate(
    collection_date = parse_date(COLL_DT_Main),
    discharge_date = parse_date(DS_DT_Main),
    death_date = parse_date(DTH_DT_Main),
    status_cr = case_when(
      HS_DTH_STATUS_Main == "Deceased" & !is.na(death_date) ~ 1,
      HS_DTH_STATUS_Main == "Alive" & !is.na(discharge_date) ~ 2,
      TRUE ~ 0
    ),
    end_date = case_when(
      status_cr == 1 ~ death_date,
      status_cr == 2 ~ discharge_date,
      TRUE ~ as.Date(NA)
    ),
    time_cr = as.numeric(end_date - collection_date)
  ) %>%
  filter(
    !is.na(collection_date),
    !is.na(time_cr),
    time_cr >= 0,
    !is.na(pathotype)
  )


# ============================================================
# 4. CHECK EVENTS
# ============================================================

cat("\nStatus codes (0 censored, 1 death, 2 discharged):\n")
print(table(cr_data$status_cr))

cat("\nPathotype by status:\n")
print(table(cr_data$pathotype, cr_data$status_cr))


# ============================================================
# 5. ONE EXPOSURE VERSUS THE REFERENCE GROUP
# ============================================================

# Turn a cuminc object into a tidy data frame of step curves
cif_to_df <- function(cif) {
  curves <- cif[names(cif) != "Tests"]

  bind_rows(lapply(names(curves), function(nm) {
    tokens <- strsplit(nm, " ")[[1]]
    tibble(
      group = paste(tokens[-length(tokens)], collapse = " "),
      event = tokens[length(tokens)],
      time = curves[[nm]]$time,
      est = curves[[nm]]$est,
      se = sqrt(curves[[nm]]$var)
    )
  })) %>%
    mutate(
      lower = pmax(0, est - 1.96 * se),
      upper = pmin(1, est + 1.96 * se)
    )
}

run_competing_risk <- function(exposure_level) {
  exposure_label <- GROUP_LABELS[[exposure_level]]

  dat <- cr_data %>%
    filter(pathotype %in% c(REFERENCE_PATHOTYPE, exposure_level)) %>%
    mutate(
      group = factor(
        pathotype,
        levels = c(REFERENCE_PATHOTYPE, exposure_level),
        labels = c(REFERENCE_LABEL, exposure_label)
      )
    )

  cat("\n============================================\n")
  cat(exposure_level, "vs", REFERENCE_PATHOTYPE, "\n")
  cat("============================================\n")
  print(table(dat$group, dat$status_cr))

  cif <- cuminc(
    ftime = dat$time_cr,
    fstatus = dat$status_cr,
    group = dat$group,
    cencode = 0
  )

  cat("\n--- Cumulative incidence at", REPORT_TIMES, "days ---\n")
  print(timepoints(cif, times = REPORT_TIMES))

  cat("\n--- Gray's tests ---\n")
  print(cif$Tests)

  curves <- cif_to_df(cif) %>%
    mutate(
      group = factor(group, levels = c(REFERENCE_LABEL, exposure_label)),
      Event = factor(EVENT_LABELS[event], levels = EVENT_LABELS)
    )

  # Default two-panel view, both events per group
  panel_plot <- ggcompetingrisks(
    cif,
    multiple_panels = TRUE,
    conf.int = TRUE
  ) +
    coord_cartesian(xlim = c(0, FOLLOW_UP_DAYS))

  # Publication figure: death only, both groups on one panel
  death_plot <- curves %>%
    filter(event == "1") %>%
    ggplot(aes(x = time, y = est, colour = group, fill = group)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.15,
      colour = NA
    ) +
    geom_step(linewidth = 1) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_cartesian(xlim = c(0, FOLLOW_UP_DAYS)) +
    labs(
      x = "Days since isolation",
      y = "Cumulative incidence of in-hospital death",
      colour = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_blank(),
      panel.grid.minor = element_blank()
    )

  list(
    data = dat,
    cif = cif,
    curves = curves,
    panel_plot = panel_plot,
    death_plot = death_plot
  )
}


# ============================================================
# 6. hvKp ONLY
# ============================================================

hv <- run_competing_risk("hvKp only")
print(hv$panel_plot)
print(hv$death_plot)


# ============================================================
# 7. ESBL(+)/CP(+) ONLY
# ============================================================

amr <- run_competing_risk("ESBL(+)/CP(+) only")
print(amr$panel_plot)
print(amr$death_plot)
