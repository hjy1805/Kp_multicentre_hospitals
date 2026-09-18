# =============================================================================
# K-MER THRESHOLD CALIBRATION AND TRANSMISSION SENSITIVITY ANALYSIS
#
# Part A  Calibrate k-mer distance against the <=20 core-genome SNP reference
#         within ST11, ST14, ST45, ST147, ST307 and ST2096. Uses the original
#         isolate set, and additionally reports the calibration with
#         same-patient pairs excluded.
# Part B  Reduce the isolate set to one sample per patient before any
#         transmission-network reconstruction.
# Part C  Rebuild seqTrack networks across k-mer thresholds and recompute the
#         pathotype, acquisition and epidemiological associations.
# Part D  Assemble the supplementary sensitivity table and the network figure.
#
# Required objects:
#   kmers                 square pairwise k-mer distance matrix, isolate IDs in
#                         the column names
#   metadata              must contain KAUST_ID, a patient identifier (MRN or
#                         PT_NO) and Collection_Date_Corrected
#   kleborate             Kleborate output, isolate column named "strain"
#   pathotype_meta        pathotype assignment, isolate column named
#                         "KAUST_ID_total"
#   visits_at_collection  admission records, one row per isolate
#
# Optional:
#   dates_for()           function mapping isolate IDs to collection dates.
#                         Only used if metadata dates cannot be parsed.
# =============================================================================


# =============================================================================
# 00  PACKAGES AND SETTINGS
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(adegenet)
  library(igraph)
  library(gt)
})

DATA_DIR  <- "/Users/daneshm/Documents/Kp_KAIMRC"
CLONE_DIR <- file.path(DATA_DIR, "revision", "transmission", "clones")
OUT_DIR   <- file.path(DATA_DIR, "revision", "transmission",
                       "threshold_sensitivity")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STS <- c("ST11", "ST14", "ST45", "ST147", "ST307", "ST2096")

SNP_THRESHOLD <- 20

# --- k-mer distance scale ----------------------------------------------------
#
# KMER_THRESHOLDS are applied to `kmers / KMER_SCALE`. The raw matrix in this
# project holds k-mer difference counts in the hundreds of thousands, which is
# why the earlier per-patient collapsing code used a cut-off of 250000. The
# thresholds 2 to 8 only mean anything once the matrix has been divided down to
# that scale, and the guard in section 03 refuses to run if the two disagree.
#
# Set KMER_SCALE to 1 if `kmers` is already on the threshold scale.
KMER_SCALE      <- 1
KMER_THRESHOLDS <- c(2, 3, 4, 5, 6, 8)
PRIMARY_KMER    <- 4

stopifnot(PRIMARY_KMER %in% KMER_THRESHOLDS)

# --- patient deduplication ---------------------------------------------------
#
# "patient"         one isolate per patient, the earliest by collection date.
#                   This is the strict reading of one sample per patient.
# "patient_strain"  one isolate per patient per within-patient genomic cluster,
#                   so a patient carrying two unrelated strains contributes two
#                   independent isolates. This reproduces the earlier
#                   per-patient collapsing loop.
DEDUP_MODE <- "patient"

# Only used when DEDUP_MODE is "patient_strain".
WITHIN_PATIENT_KMER <- PRIMARY_KMER

stopifnot(DEDUP_MODE %in% c("patient", "patient_strain"))

N_GRID <- 400

# Colour nodes of the network figure by this isolate-level variable.
NETWORK_COLOUR_BY <- "pathotype"


# =============================================================================
# 01  REQUIRED OBJECTS
# =============================================================================

required_objects <- c(
  "kmers", "metadata", "kleborate", "pathotype_meta", "visits_at_collection"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  stop("Missing objects: ", paste(missing_objects, collapse = ", "),
       call. = FALSE)
}

for (col in c("KAUST_ID", "Collection_Date_Corrected")) {
  if (!col %in% names(metadata)) {
    stop("metadata must contain a ", col, " column.", call. = FALSE)
  }
}

patient_id_col <- intersect(c("MRN", "PT_NO"), names(metadata))

if (length(patient_id_col) == 0) {
  stop("metadata must contain MRN or PT_NO.", call. = FALSE)
}

patient_id_col <- patient_id_col[1]
message("Patient identifier: metadata$", patient_id_col)


# =============================================================================
# 02  PREPARE THE K-MER MATRIX
# =============================================================================

#' Accepts a matrix whose row names are positional (1, 2, 3 ...) rather than
#' isolate IDs, which is how the distance matrix arrives from the k-mer tool.
#' Row names are only copied from the columns once the matrix has been shown to
#' be symmetric, because that is the condition under which the two orderings
#' must agree.
prepare_kmer_matrix <- function(m) {

  m <- as.matrix(m)

  if (nrow(m) != ncol(m)) {
    stop("kmers is not square.", call. = FALSE)
  }

  if (is.null(colnames(m))) {
    stop("kmers has no isolate IDs in its column names.", call. = FALSE)
  }

  if (anyDuplicated(colnames(m)) > 0) {
    stop("kmers contains duplicated isolate IDs.", call. = FALSE)
  }

  if (!isTRUE(all.equal(unname(m), unname(t(m)), check.attributes = FALSE))) {
    stop(
      "kmers is not symmetric, so row identity cannot be inferred from the ",
      "column names and pairwise lookups would be ambiguous.",
      call. = FALSE
    )
  }

  positional <- is.null(rownames(m)) ||
    identical(rownames(m), as.character(seq_len(nrow(m)))) ||
    !any(rownames(m) %in% colnames(m))

  if (positional) {
    message("Row names look positional. Copying isolate IDs from the columns.")
    rownames(m) <- colnames(m)
  } else if (!identical(rownames(m), colnames(m))) {
    if (setequal(rownames(m), colnames(m))) {
      message("Row and column IDs agree but are ordered differently. ",
              "Reordering rows to match the columns.")
      m <- m[colnames(m), , drop = FALSE]
    } else {
      stop("Row and column isolate IDs in kmers are not the same set.",
           call. = FALSE)
    }
  }

  storage.mode(m) <- "double"
  m
}

kmers <- prepare_kmer_matrix(kmers)
message("Original k-mer matrix: ", nrow(kmers), " isolates")


# =============================================================================
# 03  SCALE THE MATRIX AND VERIFY THE THRESHOLDS ARE ON THAT SCALE
# =============================================================================

kmers_scaled <- kmers / KMER_SCALE

off_diag <- kmers_scaled[upper.tri(kmers_scaled)]
off_diag <- off_diag[is.finite(off_diag)]

if (length(off_diag) == 0) {
  stop("No finite off-diagonal k-mer distances.", call. = FALSE)
}

kmer_quantiles <- quantile(
  off_diag,
  probs = c(0, 0.0001, 0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 1)
)

message("\nScaled k-mer distance distribution (KMER_SCALE = ", KMER_SCALE, "):")
print(kmer_quantiles)

# Where the primary threshold sits in the observed distribution. A transmission
# threshold should select a small tail. Sitting near the median means almost
# every pair is called linked, and sitting below the minimum means none is.
primary_percentile <- 100 * mean(off_diag <= PRIMARY_KMER)

message(
  "Primary threshold <=", PRIMARY_KMER, " selects ",
  signif(primary_percentile, 3), "% of all isolate pairs."
)

if (PRIMARY_KMER >= max(off_diag)) {
  stop(
    "The primary threshold is at or above the largest observed distance (",
    signif(max(off_diag), 4), "). Every pair would be called linked. Set ",
    "KMER_SCALE to the divisor that brings the matrix onto the threshold ",
    "scale, or set KMER_THRESHOLDS to the scale of the matrix.",
    call. = FALSE
  )
}

if (PRIMARY_KMER < min(off_diag)) {
  stop(
    "The primary threshold is below the smallest observed distance (",
    signif(min(off_diag), 4), "). No pair would be called linked.",
    call. = FALSE
  )
}

if (primary_percentile > 5) {
  warning(
    "The primary threshold links ", signif(primary_percentile, 3),
    "% of all pairs, which is high for a transmission threshold. Check that ",
    "KMER_SCALE is correct.",
    call. = FALSE
  )
}


# =============================================================================
# 04  COLLECTION DATES
# =============================================================================

#' Collection dates arrive as text in several layouts. Sorting the raw strings
#' would order 12/03/2019 before 02/04/2020 correctly by accident and
#' 1/5/2019 after 10/1/2018 incorrectly, so the representative isolate for a
#' patient could silently be the wrong one. Parsing is therefore explicit and
#' the format that resolves the most values wins.
parse_collection_date <- function(x) {

  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))

  x <- as.character(x)

  formats <- c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%d-%m-%Y", "%d.%m.%Y",
               "%Y/%m/%d", "%d-%b-%Y", "%d %b %Y")

  parsed <- map(formats, ~ suppressWarnings(as.Date(x, format = .x)))
  n_ok   <- map_int(parsed, ~ sum(!is.na(.x)))

  if (max(n_ok) == 0) {
    stop("Collection_Date_Corrected could not be parsed in any known format.",
         call. = FALSE)
  }

  best <- which.max(n_ok)

  message(
    "Collection dates parsed with format ", formats[best], ": ",
    n_ok[best], " of ", length(x), " resolved."
  )

  parsed[[best]]
}


# =============================================================================
# 05  PART A HELPERS
# =============================================================================

get_snp_dist <- function(st, dir = CLONE_DIR) {

  file <- file.path(dir, paste0(st, ".filtered_polymorphic_sites.fasta"))

  if (!file.exists(file)) {
    stop("Cannot find alignment: ", file, call. = FALSE)
  }

  x <- read.dna(file, format = "fasta")

  message(st, ": ", nrow(x), " isolates, ", ncol(x), " polymorphic sites")

  d <- dist.dna(x, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)

  idx <- which(upper.tri(d), arr.ind = TRUE)

  tibble(
    ST       = st,
    isolate1 = rownames(d)[idx[, 1]],
    isolate2 = colnames(d)[idx[, 2]],
    SNP      = as.numeric(d[idx])
  )
}


add_kmer_distance <- function(pairs, m) {

  ok   <- pairs$isolate1 %in% rownames(m) & pairs$isolate2 %in% colnames(m)
  kmer <- rep(NA_real_, nrow(pairs))

  if (any(ok)) {
    kmer[ok] <- m[cbind(pairs$isolate1[ok], pairs$isolate2[ok])]
  }

  mutate(pairs, kmer = kmer)
}


evaluate_thresholds <- function(dat, thresholds) {

  truth <- dat$SNP <= SNP_THRESHOLD

  safe_ratio <- function(num, den) if (den > 0) num / den else NA_real_

  map_dfr(thresholds, function(t) {

    pred <- dat$kmer <= t

    TP <- sum(pred & truth, na.rm = TRUE)
    TN <- sum(!pred & !truth, na.rm = TRUE)
    FP <- sum(pred & !truth, na.rm = TRUE)
    FN <- sum(!pred & truth, na.rm = TRUE)

    sensitivity <- safe_ratio(TP, TP + FN)
    specificity <- safe_ratio(TN, TN + FP)

    tibble(
      threshold = t,
      TP = TP, FP = FP, TN = TN, FN = FN,
      sensitivity       = sensitivity,
      specificity       = specificity,
      PPV               = safe_ratio(TP, TP + FP),
      NPV               = safe_ratio(TN, TN + FN),
      F1                = safe_ratio(2 * TP, 2 * TP + FP + FN),
      balanced_accuracy = (sensitivity + specificity) / 2,
      Youden            = sensitivity + specificity - 1
    )
  })
}


best_threshold <- function(perf, rule = c("youden", "spec"), min_spec = 0.95) {

  rule <- match.arg(rule)

  if (rule == "youden") {
    perf %>%
      filter(!is.na(Youden)) %>%
      arrange(desc(Youden), desc(specificity), threshold) %>%
      slice(1)
  } else {
    perf %>%
      filter(!is.na(specificity), specificity >= min_spec) %>%
      arrange(desc(sensitivity), threshold) %>%
      slice(1)
  }
}


#' Run the whole calibration on one pair set. Returning a list keeps the
#' all-pairs analysis and the same-patient-excluded analysis on identical code.
calibrate <- function(pairs, label) {

  grid <- pairs$kmer %>%
    quantile(probs = seq(0, 1, length.out = N_GRID), na.rm = TRUE) %>%
    c(KMER_THRESHOLDS) %>%
    unique() %>%
    sort()

  by_st <- pairs %>%
    group_split(ST) %>%
    map_dfr(function(x) {
      evaluate_thresholds(x, grid) %>% mutate(ST = unique(x$ST), .before = 1)
    })

  overall <- evaluate_thresholds(pairs, grid)

  youden <- by_st %>%
    group_by(ST) %>%
    group_modify(~ best_threshold(.x, "youden")) %>%
    ungroup()

  s95 <- by_st %>%
    group_by(ST) %>%
    group_modify(~ best_threshold(.x, "spec", min_spec = 0.95)) %>%
    ungroup()

  s99 <- by_st %>%
    group_by(ST) %>%
    group_modify(~ best_threshold(.x, "spec", min_spec = 0.99)) %>%
    ungroup()

  dropped_95 <- setdiff(unique(pairs$ST), s95$ST)
  if (length(dropped_95) > 0) {
    message(
      "  [", label, "] no threshold reached 95% specificity in: ",
      paste(dropped_95, collapse = ", ")
    )
  }

  list(
    label       = label,
    pairs       = pairs,
    by_st       = by_st,
    overall     = overall,
    youden      = youden,
    spec_95     = s95,
    spec_99     = s99,
    candidates  = by_st %>%
      filter(threshold %in% KMER_THRESHOLDS) %>%
      select(ST, threshold, sensitivity, specificity, PPV, NPV, F1,
             balanced_accuracy, TP, FP, TN, FN) %>%
      arrange(ST, threshold),
    overall_candidates = overall %>% filter(threshold %in% KMER_THRESHOLDS)
  )
}


# =============================================================================
# 06  PART A  SNP TO K-MER CALIBRATION
# =============================================================================

message("\n============================================================")
message("PART A: SNP-to-k-mer calibration")
message("============================================================")

snp_pairs <- map_dfr(STS, get_snp_dist)
message("Total pairwise SNP comparisons: ", nrow(snp_pairs))

snp_isolates <- unique(c(snp_pairs$isolate1, snp_pairs$isolate2))
missing_kmer <- setdiff(snp_isolates, colnames(kmers_scaled))

message(
  "Alignment isolates: ", length(snp_isolates),
  " | matched: ", sum(snp_isolates %in% colnames(kmers_scaled)),
  " | unmatched: ", length(missing_kmer)
)

if (length(missing_kmer) > 0) {
  message("Unmatched IDs (first 20): ",
          paste(head(missing_kmer, 20), collapse = ", "))
}

# Patient identity is attached here so that the calibration can be repeated
# with same-patient pairs removed. Repeat isolates from one patient are
# near-identical by construction and would otherwise dominate the linked class
# and flatter the apparent performance of any threshold.
patient_of <- metadata %>%
  transmute(
    KAUST_ID,
    patient_id = as.character(.data[[patient_id_col]])
  ) %>%
  filter(!is.na(KAUST_ID)) %>%
  distinct(KAUST_ID, .keep_all = TRUE)

calib <- snp_pairs %>%
  add_kmer_distance(kmers_scaled) %>%
  filter(!is.na(SNP), !is.na(kmer)) %>%
  left_join(rename(patient_of, patient_1 = patient_id),
            by = c("isolate1" = "KAUST_ID")) %>%
  left_join(rename(patient_of, patient_2 = patient_id),
            by = c("isolate2" = "KAUST_ID")) %>%
  mutate(
    SNP_link     = SNP <= SNP_THRESHOLD,
    same_patient = as.integer(
      !is.na(patient_1) & !is.na(patient_2) & patient_1 == patient_2
    )
  )

if (nrow(calib) == 0) {
  stop("No SNP pairs matched the k-mer matrix.", call. = FALSE)
}

n_same_patient_calib <- sum(calib$same_patient == 1)

message(
  "Pairs with both distances: ", nrow(calib),
  " | same-patient pairs: ", n_same_patient_calib,
  " (", round(100 * n_same_patient_calib / nrow(calib), 1), "%)"
)

calib_between <- filter(calib, same_patient == 0)

# ---- per-lineage summary ----------------------------------------------------

summary_ST <- calib %>%
  group_by(ST) %>%
  summarise(
    isolates          = n_distinct(c(isolate1, isolate2)),
    pairs             = n(),
    same_patient_pairs = sum(same_patient),
    SNP_le_ref        = sum(SNP_link),
    SNP_gt_ref        = sum(!SNP_link),
    proportion_link   = mean(SNP_link),
    median_SNP        = median(SNP),
    max_SNP           = max(SNP),
    median_kmer       = median(kmer),
    max_kmer          = max(kmer),
    Spearman_rho      = cor(SNP, kmer, method = "spearman",
                            use = "complete.obs"),
    .groups = "drop"
  )

print(summary_ST)

# ---- k-mer distribution in windows around the SNP reference -----------------

boundary_summary <- map_dfr(c(0, 1, 2, 3, 5, 8), function(w) {

  temp <- filter(calib,
                 SNP >= SNP_THRESHOLD - w,
                 SNP <= SNP_THRESHOLD + w)

  if (nrow(temp) == 0) return(NULL)

  temp %>%
    group_by(ST) %>%
    summarise(
      window      = paste0(SNP_THRESHOLD, " +/- ", w),
      lower_SNP   = SNP_THRESHOLD - w,
      upper_SNP   = SNP_THRESHOLD + w,
      n_pairs     = n(),
      median_kmer = median(kmer),
      q25         = unname(quantile(kmer, 0.25)),
      q75         = unname(quantile(kmer, 0.75)),
      q90         = unname(quantile(kmer, 0.90)),
      q95         = unname(quantile(kmer, 0.95)),
      .groups = "drop"
    )
})

print(boundary_summary)

# ---- calibration, all pairs and between-patient pairs -----------------------

cal_all     <- calibrate(calib,         "all pairs")
cal_between <- calibrate(calib_between, "between-patient pairs")

performance         <- cal_all$by_st
overall_performance <- cal_all$overall
best_youden         <- cal_all$youden
spec_95             <- cal_all$spec_95
spec_99             <- cal_all$spec_99
candidate_performance <- cal_all$candidates
overall_candidates    <- cal_all$overall_candidates

message("\nOptimal threshold by Youden index (all pairs):")
print(best_youden)

message("\nOptimal threshold by Youden index (between-patient pairs only):")
print(cal_between$youden)

message("\nPerformance at the configured thresholds (all pairs):")
print(candidate_performance)

# ---- misclassified pairs at the primary threshold ---------------------------

false_positives <- calib %>%
  filter(kmer <= PRIMARY_KMER, !SNP_link) %>%
  arrange(ST, desc(SNP), kmer)

false_negatives <- calib %>%
  filter(SNP_link, kmer > PRIMARY_KMER) %>%
  arrange(ST, SNP, kmer)

message(
  "\nAt scaled k-mer <= ", PRIMARY_KMER, ": ",
  nrow(false_positives), " false positives, ",
  nrow(false_negatives), " false negatives"
)

# ---- reviewer table ---------------------------------------------------------

build_reviewer_table <- function(cal, summ) {
  summ %>%
    select(ST, isolates, pairs, same_patient_pairs,
           SNP_le_ref, SNP_gt_ref, Spearman_rho) %>%
    left_join(
      cal$youden %>%
        select(ST,
               optimal_kmer              = threshold,
               optimal_sensitivity       = sensitivity,
               optimal_specificity       = specificity,
               optimal_PPV               = PPV,
               optimal_balanced_accuracy = balanced_accuracy),
      by = "ST"
    ) %>%
    left_join(
      cal$candidates %>%
        filter(threshold == PRIMARY_KMER) %>%
        select(ST,
               primary_sensitivity       = sensitivity,
               primary_specificity       = specificity,
               primary_PPV               = PPV,
               primary_balanced_accuracy = balanced_accuracy,
               primary_FP                = FP,
               primary_FN                = FN),
      by = "ST"
    )
}

reviewer_table <- build_reviewer_table(cal_all, summary_ST)

reviewer_table_between <- build_reviewer_table(
  cal_between,
  calib_between %>%
    group_by(ST) %>%
    summarise(
      isolates = n_distinct(c(isolate1, isolate2)),
      pairs = n(),
      same_patient_pairs = 0L,
      SNP_le_ref = sum(SNP_link),
      SNP_gt_ref = sum(!SNP_link),
      Spearman_rho = cor(SNP, kmer, method = "spearman",
                         use = "complete.obs"),
      .groups = "drop"
    )
)

message("\nReviewer calibration table (all pairs):")
print(reviewer_table)

# ---- calibration figures ----------------------------------------------------

p_scatter <- ggplot(calib, aes(SNP, kmer)) +
  geom_point(aes(colour = factor(same_patient)), alpha = 0.25, size = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.5,
              colour = "black") +
  geom_vline(xintercept = SNP_THRESHOLD, linetype = "dashed") +
  geom_hline(yintercept = PRIMARY_KMER, linetype = "dotted") +
  scale_colour_manual(
    values = c(`0` = "grey40", `1` = "firebrick"),
    labels = c(`0` = "Between patients", `1` = "Same patient"),
    name   = NULL
  ) +
  facet_wrap(~ ST, scales = "free") +
  labs(x = "Pairwise SNP distance", y = "Scaled pairwise k-mer distance") +
  theme_bw() +
  theme(legend.position = "bottom")

p_roc <- performance %>%
  filter(threshold <= max(KMER_THRESHOLDS) * 3) %>%
  select(ST, threshold, sensitivity, specificity) %>%
  pivot_longer(c(sensitivity, specificity),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(threshold, value, linetype = metric)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = PRIMARY_KMER, linetype = "dotted") +
  facet_wrap(~ ST) +
  labs(x = "Scaled k-mer distance threshold",
       y = "Classification performance", linetype = NULL) +
  theme_bw()


# =============================================================================
# 07  PART B  ONE SAMPLE PER PATIENT
# =============================================================================

message("\n============================================================")
message("PART B: patient deduplication (mode: ", DEDUP_MODE, ")")
message("============================================================")

ids <- colnames(kmers_scaled)

md_raw <- metadata %>%
  transmute(
    KAUST_ID,
    patient_id      = as.character(.data[[patient_id_col]]),
    collection_date = parse_collection_date(Collection_Date_Corrected)
  ) %>%
  filter(KAUST_ID %in% ids)

if (anyDuplicated(md_raw$KAUST_ID) > 0) {
  dup_ids <- unique(md_raw$KAUST_ID[duplicated(md_raw$KAUST_ID)])
  conflicting <- md_raw %>%
    filter(KAUST_ID %in% dup_ids) %>%
    group_by(KAUST_ID) %>%
    filter(n_distinct(patient_id) > 1) %>%
    ungroup()
  warning(
    length(dup_ids), " isolates have more than one metadata row",
    if (nrow(conflicting) > 0) {
      paste0(", and ", n_distinct(conflicting$KAUST_ID),
             " of them carry conflicting patient IDs")
    } else "",
    ". Keeping the first row for each.",
    call. = FALSE
  )
}

md_dedup <- distinct(md_raw, KAUST_ID, .keep_all = TRUE)

ids_without_metadata <- setdiff(ids, md_dedup$KAUST_ID)
ids_missing_patient  <- md_dedup$KAUST_ID[is.na(md_dedup$patient_id)]

n_missing_date <- sum(
  is.na(md_dedup$collection_date) & !is.na(md_dedup$patient_id)
)

message("Isolates in k-mer matrix: ", length(ids))
message("Represented in metadata: ", nrow(md_dedup))
message("With a known patient ID: ", sum(!is.na(md_dedup$patient_id)))
message("Missing a patient ID: ", length(ids_missing_patient))
message("Absent from metadata: ", length(ids_without_metadata))
message("Known patient but unparsed collection date: ", n_missing_date)

patient_sampling <- md_dedup %>%
  filter(!is.na(patient_id)) %>%
  count(patient_id, name = "n_original_isolates")

n_known_patients  <- nrow(patient_sampling)
n_repeat_patients <- sum(patient_sampling$n_original_isolates > 1)

message("Known patients: ", n_known_patients)
message("Patients with more than one isolate: ", n_repeat_patients)


#' Earliest isolate per patient. Isolates whose date could not be parsed sort
#' last within a patient, so a dated isolate is always preferred over an
#' undated one.
representatives_by_patient <- function(md) {
  md %>%
    filter(!is.na(patient_id)) %>%
    arrange(patient_id, is.na(collection_date), collection_date, KAUST_ID) %>%
    group_by(patient_id) %>%
    slice(1L) %>%
    ungroup() %>%
    mutate(within_patient_cluster = 1L)
}


#' Earliest isolate per patient per within-patient genomic cluster. Isolates
#' from one patient that sit further apart than WITHIN_PATIENT_KMER are treated
#' as separate acquisitions and each contributes a representative.
representatives_by_patient_strain <- function(md, m, cutoff) {

  md <- filter(md, !is.na(patient_id))

  md %>%
    group_split(patient_id) %>%
    map_dfr(function(p) {

      s <- intersect(p$KAUST_ID, colnames(m))

      if (length(s) <= 1) {
        return(mutate(p, within_patient_cluster = 1L))
      }

      sub <- m[s, s, drop = FALSE]

      adj <- (sub <= cutoff) & !is.na(sub)
      diag(adj) <- FALSE

      g    <- graph_from_adjacency_matrix(adj, mode = "undirected")
      comp <- components(g)$membership

      p %>%
        left_join(
          tibble(KAUST_ID = names(comp),
                 within_patient_cluster = as.integer(comp)),
          by = "KAUST_ID"
        ) %>%
        mutate(
          within_patient_cluster = coalesce(within_patient_cluster, 0L)
        ) %>%
        arrange(within_patient_cluster, is.na(collection_date),
                collection_date, KAUST_ID) %>%
        group_by(within_patient_cluster) %>%
        slice(1L) %>%
        ungroup()
    })
}

patient_representatives <- if (DEDUP_MODE == "patient") {
  representatives_by_patient(md_dedup)
} else {
  representatives_by_patient_strain(md_dedup, kmers_scaled,
                                    WITHIN_PATIENT_KMER)
}

# Isolates whose patient identity is unknown cannot be checked for repetition,
# so they are retained. If there are many of them the one-sample-per-patient
# claim is weaker than it looks, hence the warning.
unresolved_isolates <- unique(c(ids_missing_patient, ids_without_metadata))

if (length(unresolved_isolates) > 0.05 * length(ids)) {
  warning(
    length(unresolved_isolates), " isolates (",
    round(100 * length(unresolved_isolates) / length(ids), 1),
    "%) have no patient identifier and are retained without a repetition ",
    "check.",
    call. = FALSE
  )
}

to_keep <- unique(c(patient_representatives$KAUST_ID, unresolved_isolates))
to_keep <- ids[ids %in% to_keep]          # preserve original ordering

removed_isolates <- setdiff(ids, to_keep)

kmers_dedup <- kmers_scaled[to_keep, to_keep, drop = FALSE]

stopifnot(
  nrow(kmers_dedup) == ncol(kmers_dedup),
  identical(rownames(kmers_dedup), colnames(kmers_dedup))
)

# ---- verify the deduplication -----------------------------------------------

max_per_patient <- if (DEDUP_MODE == "patient") 1L else NA_integer_

dedup_check <- tibble(KAUST_ID = colnames(kmers_dedup)) %>%
  left_join(md_dedup, by = "KAUST_ID") %>%
  filter(!is.na(patient_id)) %>%
  count(patient_id, name = "n_retained")

if (DEDUP_MODE == "patient" && any(dedup_check$n_retained > 1)) {
  stop("Deduplication failed: more than one isolate remains for a patient.",
       call. = FALSE)
}

dedup_audit <- md_dedup %>%
  left_join(patient_sampling, by = "patient_id") %>%
  left_join(
    patient_representatives %>%
      select(KAUST_ID, within_patient_cluster),
    by = "KAUST_ID"
  ) %>%
  mutate(retained_for_transmission = KAUST_ID %in% to_keep) %>%
  arrange(patient_id, collection_date, KAUST_ID)

message("\nIsolates retained for transmission analysis: ", ncol(kmers_dedup))
message("Repeat isolates removed: ", length(removed_isolates))
message("Unresolved isolates retained: ", length(unresolved_isolates))

if (nrow(dedup_check) > 0) {
  message("Maximum isolates per known patient after filtering: ",
          max(dedup_check$n_retained))
}

# ---- collection dates for the retained isolates -----------------------------

dedup_dates <- tibble(KAUST_ID = colnames(kmers_dedup)) %>%
  left_join(select(md_dedup, KAUST_ID, collection_date), by = "KAUST_ID") %>%
  pull(collection_date)

if (anyNA(dedup_dates)) {

  if (exists("dates_for") && is.function(dates_for)) {
    message("Filling ", sum(is.na(dedup_dates)),
            " missing collection dates from dates_for().")
    fallback <- as.Date(dates_for(colnames(kmers_dedup)))
    dedup_dates[is.na(dedup_dates)] <- fallback[is.na(dedup_dates)]
  }

  if (anyNA(dedup_dates)) {
    stop(
      sum(is.na(dedup_dates)), " retained isolates have no collection date. ",
      "seqTrack cannot order them, so these must be resolved or excluded ",
      "before the network is built.",
      call. = FALSE
    )
  }
}


# =============================================================================
# 08  TRANSMISSION-NETWORK HELPERS
# =============================================================================

#' seqTrack returns positional indices in its id and ances columns and orients
#' every edge from ancestor to descendant by date, so the adjacency it produces
#' is asymmetric and roughly half the edges fall below the diagonal. Both the
#' upper-triangle pair extraction and igraph's undirected conversion read only
#' the upper triangle, so without symmetrising here a real transmission edge is
#' dropped from the network and then counted in the unlinked comparison group.
build_network <- function(dmat, dates, max_distance) {

  ids <- colnames(dmat)

  if (length(dates) != length(ids)) {
    stop("Collection dates and isolate IDs differ in length.", call. = FALSE)
  }

  if (anyNA(dates)) {
    stop("seqTrack cannot accept missing collection dates.", call. = FALSE)
  }

  st <- seqTrack(dmat, x.names = ids, x.dates = dates)

  edges <- as_tibble(st) %>%
    filter(!is.na(ances), !is.na(weight), weight <= max_distance)

  adj <- matrix(0L, length(ids), length(ids), dimnames = list(ids, ids))

  if (nrow(edges) > 0) {
    adj[cbind(edges$ances, edges$id)] <- 1L
  }

  adj <- pmax(adj, t(adj))
  diag(adj) <- 0L

  linked  <- rowSums(adj) > 0
  adj_sub <- adj[linked, linked, drop = FALSE]

  list(
    edges     = edges,
    adjacency = adj_sub,
    graph     = graph_from_adjacency_matrix(adj_sub, mode = "undirected",
                                            diag = FALSE)
  )
}


extract_pairs <- function(adj, value = 1L) {

  idx <- which(adj == value & upper.tri(adj), arr.ind = TRUE)

  tibble(
    KAUST_ID1 = rownames(adj)[idx[, 1]],
    KAUST_ID2 = colnames(adj)[idx[, 2]]
  )
}


eq_flag <- function(a, b) {
  as.integer(!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
}


annotate_pairs <- function(pairs, visits) {

  vis <- visits %>%
    select(
      KAUST_ID,
      MRN     = PT_NO,
      FCLT_NM,
      DS_DT,
      ADS_DT,
      ADS_WD  = ADS_WD_DEPT_CD,
      DS_WD   = DS_WD_DEPT_CD
    ) %>%
    mutate(MRN = as.character(MRN))

  if (anyDuplicated(vis$KAUST_ID) > 0) {
    warning(
      "visits_at_collection has ", sum(duplicated(vis$KAUST_ID)),
      " duplicated isolate IDs. A many-to-many join would multiply the pair ",
      "counts every odds ratio is computed from, so the first record for ",
      "each isolate is kept.",
      call. = FALSE
    )
    vis <- distinct(vis, KAUST_ID, .keep_all = TRUE)
  }

  suffixed <- function(df, s) rename_with(df, ~ paste0(.x, s), -KAUST_ID)

  pairs %>%
    left_join(suffixed(vis, "_1"), by = c("KAUST_ID1" = "KAUST_ID")) %>%
    left_join(suffixed(vis, "_2"), by = c("KAUST_ID2" = "KAUST_ID")) %>%
    mutate(
      same_patient    = eq_flag(MRN_1, MRN_2),
      same_region     = eq_flag(FCLT_NM_1, FCLT_NM_2),
      same_ward_adm   = eq_flag(ADS_WD_1, ADS_WD_2),
      same_ward_dis   = eq_flag(DS_WD_1, DS_WD_2),
      same_ward_cross = as.integer(
        eq_flag(ADS_WD_1, DS_WD_2) == 1 | eq_flag(ADS_WD_2, DS_WD_1) == 1
      ),
      same_ward_any = as.integer(
        same_ward_adm == 1 | same_ward_dis == 1 | same_ward_cross == 1
      ),
      start1  = pmin(ADS_DT_1, DS_DT_1),
      end1    = pmax(ADS_DT_1, DS_DT_1),
      start2  = pmin(ADS_DT_2, DS_DT_2),
      end2    = pmax(ADS_DT_2, DS_DT_2),
      overlap = as.integer(
        !is.na(start1) & !is.na(end1) & !is.na(start2) & !is.na(end2) &
          pmax(start1, start2) <= pmin(end1, end2)
      ),
      overlap_same_ward = as.integer(overlap == 1 & same_ward_any == 1)
    )
}


fisher_or <- function(df, group_col, flag_col,
                      group_levels = c("linked", "unlinked")) {

  tab <- table(
    factor(df[[group_col]], levels = group_levels),
    factor(df[[flag_col]],  levels = c(1, 0))
  )

  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    return(tibble(variable = flag_col, OR = NA_real_,
                  CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_))
  }

  ft <- fisher.test(tab)

  tibble(
    variable = flag_col,
    OR       = unname(ft$estimate),
    CI_low   = ft$conf.int[1],
    CI_high  = ft$conf.int[2],
    p_value  = ft$p.value
  )
}


enrichment_by_level <- function(df, group_col, level_col, focus) {

  d <- df %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[level_col]]))

  g <- as.character(d[[group_col]])
  v <- as.character(d[[level_col]])

  map_dfr(sort(unique(v)), function(lv) {

    tab <- matrix(
      c(sum(g == focus & v == lv),
        sum(g != focus & v == lv),
        sum(g == focus & v != lv),
        sum(g != focus & v != lv)),
      nrow = 2
    )

    base <- tibble(
      level = lv, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p_value = NA_real_,
      a = tab[1, 1], b = tab[2, 1], c = tab[1, 2], d = tab[2, 2]
    )

    if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) return(base)

    ft <- fisher.test(tab)

    base %>%
      mutate(
        OR      = unname(ft$estimate),
        CI_low  = ft$conf.int[1],
        CI_high = ft$conf.int[2],
        p_value = ft$p.value
      )
  })
}


EPI_PREDICTORS <- c(
  "same_ward_any", "same_ward_dis", "same_ward_adm",
  "overlap_same_ward", "overlap"
)


run_threshold_analysis <- function(threshold) {

  message("\n-- scaled k-mer threshold <= ", threshold, " --")

  net <- build_network(kmers_dedup, dedup_dates, threshold)

  if (vcount(net$graph) == 0) {
    warning("No linked isolates at threshold ", threshold, call. = FALSE)
    return(NULL)
  }

  comp <- components(net$graph)$membership

  isolate_table <- tibble(KAUST_ID = colnames(kmers_dedup)) %>%
    mutate(
      cluster_status = factor(
        if_else(KAUST_ID %in% names(comp), "Cluster", "None"),
        levels = c("None", "Cluster")
      ),
      component = as.integer(comp[match(KAUST_ID, names(comp))])
    ) %>%
    left_join(distinct(kleborate, strain, .keep_all = TRUE),
              by = c("KAUST_ID" = "strain")) %>%
    left_join(
      metadata %>%
        select(KAUST_ID, Acquisition) %>%
        distinct(KAUST_ID, .keep_all = TRUE),
      by = "KAUST_ID"
    ) %>%
    left_join(
      distinct(pathotype_meta, KAUST_ID_total, .keep_all = TRUE),
      by = c("KAUST_ID" = "KAUST_ID_total")
    )

  or_pathotype <- isolate_table %>%
    filter(!is.na(pathotype), pathotype != "Convergent hvKp") %>%
    enrichment_by_level("cluster_status", "pathotype", focus = "Cluster") %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
           threshold = threshold, analysis = "Pathotype")

  or_acquisition <- isolate_table %>%
    filter(!is.na(Acquisition), Acquisition != "other") %>%
    enrichment_by_level("cluster_status", "Acquisition", focus = "Cluster") %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
           threshold = threshold, analysis = "Acquisition")

  # Linked pairs are those joined by a seqTrack edge. Unlinked pairs are the
  # remaining pairs among isolates that entered the network, so the contrast is
  # edge versus no edge within the clustered set rather than clustered versus
  # unclustered isolates.
  annotated <- bind_rows(
    annotate_pairs(extract_pairs(net$adjacency, 1L),
                   visits_at_collection) %>% mutate(pair_type = "linked"),
    annotate_pairs(extract_pairs(net$adjacency, 0L),
                   visits_at_collection) %>% mutate(pair_type = "unlinked")
  ) %>%
    filter(!is.na(MRN_1), !is.na(MRN_2))

  n_same_patient_pairs <- sum(annotated$same_patient == 1, na.rm = TRUE)

  if (n_same_patient_pairs > 0) {
    warning(
      n_same_patient_pairs, " same-patient pairs remain at threshold ",
      threshold, ". Check that metadata$", patient_id_col,
      " and visits_at_collection$PT_NO use the same identifier.",
      call. = FALSE
    )
  }

  annotated_epi <- filter(annotated, same_patient == 0)

  or_epi <- map_dfr(EPI_PREDICTORS,
                    ~ fisher_or(annotated_epi, "pair_type", .x)) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
           threshold = threshold, analysis = "Epidemiology")

  summary_tbl <- tibble(
    threshold                   = threshold,
    n_edges                     = nrow(net$edges),
    n_clustered_isolates        = sum(isolate_table$cluster_status == "Cluster"),
    n_unclustered_isolates      = sum(isolate_table$cluster_status == "None"),
    n_components                = length(unique(comp)),
    largest_component           = max(table(comp)),
    linked_pairs_with_visits    = sum(annotated_epi$pair_type == "linked"),
    unlinked_pairs_with_visits  = sum(annotated_epi$pair_type == "unlinked"),
    same_patient_pairs_detected = n_same_patient_pairs
  )

  message(
    "edges=", summary_tbl$n_edges,
    " | clustered=", summary_tbl$n_clustered_isolates,
    " | components=", summary_tbl$n_components,
    " | largest=", summary_tbl$largest_component,
    " | same-patient pairs=", n_same_patient_pairs
  )

  list(
    threshold       = threshold,
    network         = net,
    isolate_table   = isolate_table,
    annotated_pairs = annotated_epi,
    pathotype       = or_pathotype,
    acquisition     = or_acquisition,
    epidemiology    = or_epi,
    summary         = summary_tbl
  )
}


# =============================================================================
# 09  PART C  THRESHOLD SENSITIVITY
# =============================================================================

message("\n============================================================")
message("PART C: transmission threshold sensitivity")
message("============================================================")

threshold_results <- KMER_THRESHOLDS %>%
  set_names(paste0("kmer_", KMER_THRESHOLDS)) %>%
  map(run_threshold_analysis) %>%
  compact()

if (length(threshold_results) == 0) {
  stop("No threshold produced a network.", call. = FALSE)
}

network_sensitivity      <- map_dfr(threshold_results, "summary")
pathotype_sensitivity    <- map_dfr(threshold_results, "pathotype")
acquisition_sensitivity  <- map_dfr(threshold_results, "acquisition")
epidemiology_sensitivity <- map_dfr(threshold_results, "epidemiology")

message("\nNetwork size across thresholds:")
print(network_sensitivity)

nosocomial_sensitivity <- acquisition_sensitivity %>%
  filter(level == "nosocomial") %>%
  select(threshold, OR, CI_low, CI_high, p_value, p_adj_BH)

esbl_sensitivity <- pathotype_sensitivity %>%
  filter(level == "ESBL(+)/CP(+) only") %>%
  select(threshold, OR, CI_low, CI_high, p_value, p_adj_BH)

message("\nNosocomial acquisition:")
print(nosocomial_sensitivity)

message("\nESBL(+)/CP(+) only:")
print(esbl_sensitivity)

message("\nEpidemiological associations:")
print(
  epidemiology_sensitivity %>%
    select(threshold, variable, OR, CI_low, CI_high, p_value, p_adj_BH) %>%
    arrange(variable, threshold)
)


# =============================================================================
# 10  SENSITIVITY FIGURES
# =============================================================================

or_plot <- function(dat, facet_var = NULL, title = NULL) {

  p <- ggplot(dat, aes(factor(threshold), OR)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15) +
    geom_point(size = 2.5) +
    scale_y_log10() +
    labs(x = "Scaled k-mer distance threshold", y = "Odds ratio",
         title = title) +
    theme_bw()

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(vars(.data[[facet_var]]), scales = "free_y")
  }

  p
}

p_nosocomial <- or_plot(nosocomial_sensitivity,
                        title = "Nosocomial acquisition")
p_pathotype  <- or_plot(pathotype_sensitivity, "level",
                        title = "Pathotype associations")
p_epi        <- or_plot(epidemiology_sensitivity, "variable",
                        title = "Epidemiological associations")

p_cluster_size <- ggplot(network_sensitivity,
                         aes(threshold, n_clustered_isolates)) +
  geom_line() +
  geom_point(size = 3) +
  geom_vline(xintercept = PRIMARY_KMER, linetype = "dashed") +
  scale_x_continuous(breaks = KMER_THRESHOLDS) +
  labs(x = "Scaled k-mer distance threshold",
       y = "Clustered isolates (one per patient)",
       title = "Network size sensitivity") +
  theme_bw()


# =============================================================================
# 11  NETWORK FIGURE AT THE PRIMARY THRESHOLD
# =============================================================================

primary_result <- threshold_results[[paste0("kmer_", PRIMARY_KMER)]]

plot_transmission_network <- function(result, colour_by = NETWORK_COLOUR_BY,
                                      file = NULL) {

  g <- result$network$graph

  if (!colour_by %in% names(result$isolate_table)) {
    warning("Colour variable '", colour_by,
            "' is not in the isolate table. Drawing an uncoloured network.",
            call. = FALSE)
    values <- rep("unknown", vcount(g))
  } else {
    lookup <- result$isolate_table[[colour_by]]
    names(lookup) <- result$isolate_table$KAUST_ID
    values <- as.character(lookup[V(g)$name])
  }

  values[is.na(values)] <- "unknown"

  lev <- sort(unique(values))
  pal <- setNames(
    grDevices::hcl.colors(max(length(lev), 2), palette = "Dark 3")[
      seq_along(lev)
    ],
    lev
  )

  draw <- function() {
    op <- par(mar = c(0, 0, 2, 0))
    on.exit(par(op), add = TRUE)
    set.seed(1)
    plot(
      g,
      layout       = layout_with_fr(g),
      vertex.size  = 5,
      vertex.label = NA,
      vertex.color = pal[values],
      vertex.frame.color = "grey30",
      edge.width   = 1,
      edge.color   = "gray70",
      main = paste0("Hospital transmission network, scaled k-mer <= ",
                    result$threshold)
    )
    legend("bottomleft", legend = names(pal), pt.bg = pal,
           pch = 21, pt.cex = 1.4, bty = "n", cex = 0.8)
  }

  if (!is.null(file)) {
    pdf(file, width = 8, height = 8)
    draw()
    dev.off()
  }

  draw()
  invisible(pal)
}

plot_transmission_network(
  primary_result,
  file = file.path(OUT_DIR, "transmission_network_primary_threshold.pdf")
)


# =============================================================================
# 12  PART D  PUBLICATION SUPPLEMENTARY TABLE
# =============================================================================

message("\n============================================================")
message("PART D: publication table")
message("============================================================")

fmt_or <- function(or, lo, hi) {
  if_else(is.na(or), "—", sprintf("%.2f (%.2f–%.2f)", or, lo, hi))
}

EPI_LABELS <- c(
  same_ward_any     = "Same ward (any)",
  same_ward_adm     = "Same admission ward",
  same_ward_dis     = "Same discharge ward",
  overlap           = "Overlapping hospital stay",
  overlap_same_ward = "Overlapping stay and same ward"
)

CATEGORY_ORDER <- c(
  "ESBL/CP-negative non-hvKp",
  "hvKp only",
  "ESBL(+)/CP(+) only",
  "Nosocomial acquisition",
  unname(EPI_LABELS)
)

publication_long <- bind_rows(

  pathotype_sensitivity %>%
    transmute(Analysis = "Pathotype",
              Category = as.character(level),
              threshold,
              result = fmt_or(OR, CI_low, CI_high)),

  acquisition_sensitivity %>%
    filter(level == "nosocomial") %>%
    transmute(Analysis = "Acquisition",
              Category = "Nosocomial acquisition",
              threshold,
              result = fmt_or(OR, CI_low, CI_high)),

  epidemiology_sensitivity %>%
    transmute(Analysis = "Epidemiological overlap",
              Category = unname(EPI_LABELS[variable]),
              threshold,
              result = fmt_or(OR, CI_low, CI_high))
)

unmapped <- setdiff(publication_long$Category, CATEGORY_ORDER)

if (length(unmapped) > 0) {
  warning("Categories absent from CATEGORY_ORDER and dropped: ",
          paste(unmapped, collapse = ", "), call. = FALSE)
}

threshold_keys   <- paste0("k", KMER_THRESHOLDS)
threshold_labels <- set_names(paste0("≤", KMER_THRESHOLDS),
                              threshold_keys)
primary_key      <- paste0("k", PRIMARY_KMER)

publication_table <- publication_long %>%
  filter(Category %in% CATEGORY_ORDER) %>%
  mutate(
    Analysis  = factor(Analysis, levels = c("Pathotype", "Acquisition",
                                            "Epidemiological overlap")),
    Category  = factor(Category, levels = CATEGORY_ORDER),
    threshold = paste0("k", threshold)
  ) %>%
  arrange(Analysis, Category) %>%
  mutate(across(c(Analysis, Category), as.character)) %>%
  pivot_wider(names_from = threshold, values_from = result) %>%
  select(Analysis, Category, all_of(threshold_keys))

publication_gt <- publication_table %>%
  gt(groupname_col = "Analysis")

# cols_label() takes the labels as named arguments. Splicing them with do.call
# avoids the deprecated .list argument, which warns on current gt.
publication_gt <- do.call(
  cols_label,
  c(list(publication_gt), list(Category = ""), as.list(threshold_labels))
)

publication_gt <- publication_gt %>%
  tab_spanner(label = "Scaled k-mer distance threshold",
              columns = all_of(threshold_keys)) %>%
  tab_header(title = md(paste0(
    "**Sensitivity of associations with putative genomic ",
    "transmission-cluster membership to alternative k-mer distance ",
    "thresholds**"
  ))) %>%
  tab_footnote(
    footnote  = "Primary scaled k-mer distance threshold.",
    locations = cells_column_labels(columns = all_of(primary_key))
  ) %>%
  tab_source_note(source_note = md(paste0(
    "Values are odds ratios (95% confidence intervals) for membership in a ",
    "putative genomic transmission cluster. Transmission-network analyses ",
    "were restricted to one isolate per patient, retaining the earliest ",
    "collected isolate where a patient contributed more than one. The ",
    "primary analysis used a scaled k-mer distance threshold of ≤",
    PRIMARY_KMER, ". Thresholds of ≤",
    paste(setdiff(KMER_THRESHOLDS, PRIMARY_KMER), collapse = ", ≤"),
    " were evaluated in sensitivity analyses. Em dashes mark comparisons ",
    "with no informative two-by-two table."
  ))) %>%
  cols_align("left", columns = Category) %>%
  cols_align("center", columns = all_of(threshold_keys)) %>%
  cols_width(Category ~ px(240), all_of(threshold_keys) ~ px(130)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = all_of(primary_key))) %>%
  tab_options(
    table.font.names                  = "Arial",
    table.font.size                   = px(12),
    heading.title.font.size           = px(14),
    heading.title.font.weight         = "bold",
    column_labels.font.weight         = "bold",
    row_group.font.weight             = "bold",
    table.border.top.width            = px(2),
    table.border.bottom.width         = px(2),
    column_labels.border.top.width    = px(1),
    column_labels.border.bottom.width = px(1),
    row_group.border.top.width        = px(1),
    data_row.padding                  = px(6),
    source_notes.font.size            = px(10)
  )

print(publication_gt)


# =============================================================================
# 13  CLUSTER MEMBERSHIP AT THE PRIMARY THRESHOLD
# =============================================================================

primary_patient_table <- primary_result$isolate_table %>%
  select(KAUST_ID, cluster_status, component) %>%
  left_join(select(md_dedup, KAUST_ID, patient_id), by = "KAUST_ID")

n_clustered_patients <- primary_patient_table %>%
  filter(!is.na(patient_id), cluster_status == "Cluster") %>%
  distinct(patient_id) %>%
  nrow()

n_total_patients <- primary_patient_table %>%
  filter(!is.na(patient_id)) %>%
  distinct(patient_id) %>%
  nrow()

clustered_percentage <- 100 * n_clustered_patients / n_total_patients

message(
  "\nPrimary threshold <= ", PRIMARY_KMER, ": ",
  n_clustered_patients, " of ", n_total_patients,
  " patients clustered (", round(clustered_percentage, 1), "%)"
)


# =============================================================================
# 14  FINAL QC
# =============================================================================

message("\n============================================================")
message("FINAL QC")
message("============================================================")

qc <- function(condition, pass_msg, fail_msg, fatal = TRUE) {
  if (isTRUE(condition)) {
    message("PASS: ", pass_msg)
  } else if (fatal) {
    stop("QC FAILED: ", fail_msg, call. = FALSE)
  } else {
    warning("QC WARNING: ", fail_msg, call. = FALSE)
  }
}

qc(
  DEDUP_MODE != "patient" || all(dedup_check$n_retained <= 1),
  "one isolate per known patient.",
  "more than one isolate remains for a known patient."
)

qc(
  identical(rownames(kmers_dedup), colnames(kmers_dedup)),
  "k-mer row and column IDs match.",
  "k-mer row and column IDs differ."
)

qc(
  !anyNA(dedup_dates),
  "every retained isolate has a collection date.",
  "some retained isolates have no collection date."
)

qc(
  sum(primary_result$annotated_pairs$same_patient == 1, na.rm = TRUE) == 0,
  "no same-patient pairs in the primary epidemiological analysis.",
  "same-patient pairs remain in the primary epidemiological analysis.",
  fatal = FALSE
)

qc(
  primary_percentile <= 5,
  paste0("the primary threshold links ", signif(primary_percentile, 3),
         "% of pairs, consistent with a transmission threshold."),
  paste0("the primary threshold links ", signif(primary_percentile, 3),
         "% of pairs, which suggests KMER_SCALE is wrong."),
  fatal = FALSE
)


# =============================================================================
# 15  WRITE OUTPUTS
# =============================================================================

csv_outputs <- list(
  patient_deduplication_audit        = dedup_audit,
  patient_representatives            = patient_representatives,
  removed_repeat_isolates            = tibble(KAUST_ID = removed_isolates),
  calibration_pairwise_distances     = calib,
  calibration_summary_by_ST          = summary_ST,
  calibration_kmer_around_SNP_ref    = boundary_summary,
  calibration_performance_by_ST      = performance,
  calibration_performance_overall    = overall_performance,
  calibration_best_youden_by_ST      = best_youden,
  calibration_spec95_by_ST           = spec_95,
  calibration_spec99_by_ST           = spec_99,
  calibration_candidate_thresholds   = candidate_performance,
  calibration_overall_candidates     = overall_candidates,
  calibration_false_positives        = false_positives,
  calibration_false_negatives        = false_negatives,
  calibration_reviewer_table         = reviewer_table,
  calibration_reviewer_table_between = reviewer_table_between,
  calibration_performance_between    = cal_between$by_st,
  calibration_candidates_between     = cal_between$candidates,
  network_summary_by_threshold       = network_sensitivity,
  acquisition_OR_by_threshold        = acquisition_sensitivity,
  pathotype_OR_by_threshold          = pathotype_sensitivity,
  epidemiology_OR_by_threshold       = epidemiology_sensitivity,
  nosocomial_OR_by_threshold         = nosocomial_sensitivity,
  ESBL_CP_OR_by_threshold            = esbl_sensitivity,
  supplementary_table_wide           = publication_table,
  primary_threshold_patient_table    = primary_patient_table
)

iwalk(csv_outputs, function(x, nm) {
  write_csv(x, file.path(OUT_DIR, paste0(nm, ".csv")))
})

plot_outputs <- list(
  SNP_vs_kmer_by_ST           = list(p_scatter,      10, 7),
  kmer_threshold_performance  = list(p_roc,          10, 7),
  nosocomial_OR_sensitivity   = list(p_nosocomial,    7, 5),
  pathotype_OR_sensitivity    = list(p_pathotype,    10, 5),
  epidemiology_OR_sensitivity = list(p_epi,          11, 7),
  network_size_sensitivity    = list(p_cluster_size,  7, 5)
)

iwalk(plot_outputs, function(spec, nm) {
  ggsave(file.path(OUT_DIR, paste0(nm, ".pdf")),
         spec[[1]], width = spec[[2]], height = spec[[3]])
  ggsave(file.path(OUT_DIR, paste0(nm, ".png")),
         spec[[1]], width = spec[[2]], height = spec[[3]], dpi = 600)
})

for (ext in c("html", "rtf")) {
  gtsave(publication_gt,
         file.path(OUT_DIR,
                   paste0("Table_Sx_kmer_threshold_sensitivity.", ext)))
}


# =============================================================================
# 16  SUMMARY
# =============================================================================

message("\n============================================================")
message("ANALYSIS COMPLETE")
message("============================================================")
message("Isolates used for genomic calibration: ", ncol(kmers_scaled))
message("Same-patient pairs in the calibration set: ", n_same_patient_calib)
message("Isolates after patient deduplication: ", ncol(kmers_dedup))
message("Known patients in the transmission analysis: ", n_total_patients)
message("Repeat isolates removed: ", length(removed_isolates))
message("Deduplication mode: ", DEDUP_MODE)
message("KMER_SCALE: ", KMER_SCALE)
message("Primary scaled k-mer threshold: <= ", PRIMARY_KMER)
message("Clustered patients at the primary threshold: ",
        n_clustered_patients, "/", n_total_patients,
        " (", round(clustered_percentage, 1), "%)")
message("Outputs written to: ", OUT_DIR)
message("============================================================")
