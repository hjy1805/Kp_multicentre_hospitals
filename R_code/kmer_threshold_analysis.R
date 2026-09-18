# =============================================================================
# K-MER DISTANCE THRESHOLD: CALIBRATION AND SENSITIVITY ANALYSIS
#
# Part A  Calibrate k-mer distance against the <=20 SNP reference threshold
#         within each major lineage (ST11, ST14, ST45, ST147, ST307, ST2096).
# Part B  Rebuild the seqTrack transmission network at alternative k-mer
#         thresholds and recompute every association reported in Figure B.
# Part C  Assemble the publication supplementary table.
#
# Objects that must already exist in the environment before sourcing:
#   kmers                 square pairwise k-mer distance matrix, isolate IDs
#                         in dimnames
#   kmers_dedup           the same matrix after de-duplication, used for the
#                         network analysis
#   metadata              isolate metadata, must contain KAUST_ID
#   kleborate             Kleborate output, isolate column named "strain"
#   pathotype_meta        pathotype assignment, isolate column named
#                         "KAUST_ID_total"
#   visits_at_collection  one admission record per isolate
#   dates_for()           function mapping isolate IDs to collection dates
# =============================================================================


# =============================================================================
# 00  CONFIGURATION
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(adegenet)
  library(igraph)
  library(gt)
})

# Root of the project tree. Everything else is derived from it, so this is the
# only path that changes between machines.
DATA_DIR <- "/Users/daneshm/Documents/Kp_KAIMRC"

CLONE_DIR <- file.path(DATA_DIR, "revision", "transmission", "clones")
OUT_DIR   <- file.path(DATA_DIR, "revision", "transmission",
                       "threshold_sensitivity")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STS <- c("ST11", "ST14", "ST45", "ST147", "ST307", "ST2096")

# Reference definition of a genomic link, from the literature.
SNP_THRESHOLD <- 20

# K-mer thresholds carried through the sensitivity analysis. PRIMARY_KMER must
# be one of them. These are on the same scale as the values in `kmers`, which
# check_kmer_scale() verifies before anything else runs.
KMER_THRESHOLDS <- c(2, 3, 4, 5, 6, 8)
PRIMARY_KMER    <- 4

stopifnot(PRIMARY_KMER %in% KMER_THRESHOLDS)

# Number of points in the k-mer threshold sweep used for the ROC-style
# calibration in Part A. The grid is taken from observed quantiles, so the
# analysis works whether k-mer distances are counts or fractions.
N_GRID <- 400


# =============================================================================
# 01  PRECONDITIONS
# =============================================================================

require_objects <- function(...) {
  needed  <- c(...)
  missing <- needed[!map_lgl(needed, exists)]
  if (length(missing) > 0) {
    stop(
      "Missing from the environment: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_objects(
  "kmers", "kmers_dedup", "metadata", "kleborate",
  "pathotype_meta", "visits_at_collection", "dates_for"
)


#' Fail early if the configured thresholds are on a different scale from the
#' data. Mash-style distances sit near 0.001 to 0.01 while raw k-mer difference
#' counts run into the hundreds. Applying integer thresholds to fractional
#' distances silently classifies every pair as linked, which is the failure
#' this guard exists to prevent.
check_kmer_scale <- function(m, thresholds) {

  off <- m[upper.tri(m)]
  off <- off[is.finite(off)]

  rng <- range(off)
  med <- median(off)

  message(
    "K-mer distances: min ", signif(rng[1], 3),
    ", median ", signif(med, 3),
    ", max ", signif(rng[2], 3)
  )

  if (max(thresholds) >= rng[2]) {
    stop(
      "Every configured k-mer threshold is at or above the largest observed ",
      "distance (", signif(rng[2], 3), "). Sensitivity would be 1 and ",
      "specificity 0 at every threshold. Set KMER_THRESHOLDS to the scale of ",
      "the data.",
      call. = FALSE
    )
  }

  if (min(thresholds) <= rng[1]) {
    warning(
      "The smallest configured threshold is at or below the smallest ",
      "observed distance, so it will link no pairs.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# =============================================================================
# 02  PREPARE THE K-MER MATRIX
# =============================================================================

prepare_kmer_matrix <- function(m) {

  m <- as.matrix(m)

  if (nrow(m) != ncol(m)) {
    stop("kmers is not square.", call. = FALSE)
  }

  if (is.null(colnames(m))) {
    stop("kmers has no isolate IDs in its column names.", call. = FALSE)
  }

  # The matrix as supplied carries positional row names and isolate IDs in the
  # columns. Copying the column names onto the rows is only safe because the
  # matrix is symmetric, which is asserted below.
  if (is.null(rownames(m)) || !identical(rownames(m), colnames(m))) {
    rownames(m) <- colnames(m)
  }

  if (!isTRUE(all.equal(m, t(m), check.attributes = FALSE))) {
    stop(
      "kmers is not symmetric. Row names cannot be inferred from column ",
      "names and pairwise lookups would be ambiguous.",
      call. = FALSE
    )
  }

  if (any(duplicated(colnames(m)))) {
    stop("kmers has duplicated isolate IDs.", call. = FALSE)
  }

  m
}

kmers <- prepare_kmer_matrix(kmers)

check_kmer_scale(kmers, KMER_THRESHOLDS)

message("K-mer matrix: ", nrow(kmers), " isolates")


# =============================================================================
# 03  PART A HELPERS
# =============================================================================

#' Pairwise SNP distances from a lineage polymorphic-site alignment.
get_snp_dist <- function(st, dir = CLONE_DIR) {

  file <- file.path(dir, paste0(st, ".filtered_polymorphic_sites.fasta"))

  if (!file.exists(file)) {
    stop("Cannot find alignment: ", file, call. = FALSE)
  }

  x <- read.dna(file, format = "fasta")

  message(
    "  ", st, ": ", nrow(x), " isolates, ", ncol(x), " polymorphic sites"
  )

  # model = "N" returns the raw count of differing sites.
  d <- dist.dna(x, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)

  idx <- which(upper.tri(d), arr.ind = TRUE)

  tibble(
    ST       = st,
    isolate1 = rownames(d)[idx[, 1]],
    isolate2 = colnames(d)[idx[, 2]],
    SNP      = as.numeric(d[idx])
  )
}


#' Vectorised k-mer lookup for a table of isolate pairs. The original
#' map2_dbl() version made one matrix subscript call per pair, which dominates
#' runtime once the pair table passes a few tens of thousands of rows.
add_kmer_distance <- function(pairs, m) {

  ok <- pairs$isolate1 %in% rownames(m) & pairs$isolate2 %in% colnames(m)

  kmer <- rep(NA_real_, nrow(pairs))

  if (any(ok)) {
    kmer[ok] <- m[cbind(pairs$isolate1[ok], pairs$isolate2[ok])]
  }

  mutate(pairs, kmer = kmer)
}


#' Classification performance of a set of k-mer thresholds against the SNP
#' reference, for one lineage or for the pooled data.
evaluate_thresholds <- function(dat, thresholds) {

  truth <- dat$SNP <= SNP_THRESHOLD

  map_dfr(thresholds, function(t) {

    pred <- dat$kmer <= t

    TP <- sum(pred & truth, na.rm = TRUE)
    TN <- sum(!pred & !truth, na.rm = TRUE)
    FP <- sum(pred & !truth, na.rm = TRUE)
    FN <- sum(!pred & truth, na.rm = TRUE)

    safe_ratio <- function(num, den) if (den > 0) num / den else NA_real_

    sens <- safe_ratio(TP, TP + FN)
    spec <- safe_ratio(TN, TN + FP)

    tibble(
      threshold         = t,
      TP = TP, FP = FP, TN = TN, FN = FN,
      sensitivity       = sens,
      specificity       = spec,
      PPV               = safe_ratio(TP, TP + FP),
      NPV               = safe_ratio(TN, TN + FN),
      F1                = safe_ratio(2 * TP, 2 * TP + FP + FN),
      balanced_accuracy = (sens + spec) / 2,
      Youden            = sens + spec - 1
    )
  })
}


#' Best threshold under a stated rule.
#'   rule = "youden"  maximise sensitivity + specificity - 1
#'   rule = "spec"    among thresholds meeting min_spec, take the highest
#'                    sensitivity, breaking ties towards the smaller threshold
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


# =============================================================================
# 04  PART A  SNP TO K-MER CALIBRATION
# =============================================================================

message("\n== Part A: SNP to k-mer calibration ==")

snp_pairs <- map_dfr(STS, get_snp_dist)

message("Pairwise SNP comparisons: ", nrow(snp_pairs))

# ---- isolate ID matching ----------------------------------------------------

snp_isolates <- unique(c(snp_pairs$isolate1, snp_pairs$isolate2))
missing_kmer <- setdiff(snp_isolates, colnames(kmers))

message(
  "Isolates in alignments: ", length(snp_isolates),
  " | matched to k-mer matrix: ",
  sum(snp_isolates %in% colnames(kmers)),
  " | unmatched: ", length(missing_kmer)
)

if (length(missing_kmer) > 0) {
  message(
    "Unmatched IDs (first 20): ",
    paste(head(missing_kmer, 20), collapse = ", ")
  )
  # FASTA headers frequently carry a suffix the k-mer matrix does not. If the
  # counts above look wrong, normalise IDs here rather than downstream.
}

# ---- attach k-mer distances -------------------------------------------------

calib <- snp_pairs %>%
  add_kmer_distance(kmers) %>%
  filter(!is.na(SNP), !is.na(kmer)) %>%
  mutate(SNP_link = SNP <= SNP_THRESHOLD)

message(
  "Pairs with both distances: ", nrow(calib),
  " of ", nrow(snp_pairs)
)

if (nrow(calib) == 0) {
  stop("No pairs matched between the alignments and the k-mer matrix.",
       call. = FALSE)
}

# ---- per-lineage summary ----------------------------------------------------

summary_ST <- calib %>%
  group_by(ST) %>%
  summarise(
    isolates        = n_distinct(c(isolate1, isolate2)),
    pairs           = n(),
    SNP_le_ref      = sum(SNP_link),
    SNP_gt_ref      = sum(!SNP_link),
    proportion_link = mean(SNP_link),
    median_SNP      = median(SNP),
    max_SNP         = max(SNP),
    median_kmer     = median(kmer),
    max_kmer        = max(kmer),
    Spearman_rho    = cor(SNP, kmer, method = "spearman",
                          use = "complete.obs"),
    .groups = "drop"
  )

print(summary_ST)

# ---- k-mer distribution in windows around the SNP reference -----------------

boundary_summary <- map_dfr(c(0, 1, 2, 3, 5, 8), function(w) {

  temp <- filter(
    calib,
    SNP >= SNP_THRESHOLD - w,
    SNP <= SNP_THRESHOLD + w
  )

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

# ---- threshold sweep --------------------------------------------------------

# The grid comes from the observed distance distribution rather than from a
# hard-coded integer sequence, so the sweep is correct on any distance scale.
kmer_grid <- calib$kmer %>%
  quantile(probs = seq(0, 1, length.out = N_GRID), na.rm = TRUE) %>%
  c(KMER_THRESHOLDS) %>%
  unique() %>%
  sort()

performance <- calib %>%
  group_split(ST) %>%
  map_dfr(function(x) {
    evaluate_thresholds(x, kmer_grid) %>%
      mutate(ST = unique(x$ST), .before = 1)
  })

overall_performance <- evaluate_thresholds(calib, kmer_grid)

# ---- optimal and conservative thresholds ------------------------------------

best_youden <- performance %>%
  group_by(ST) %>%
  group_modify(~ best_threshold(.x, "youden")) %>%
  ungroup()

spec_95 <- performance %>%
  group_by(ST) %>%
  group_modify(~ best_threshold(.x, "spec", min_spec = 0.95)) %>%
  ungroup()

spec_99 <- performance %>%
  group_by(ST) %>%
  group_modify(~ best_threshold(.x, "spec", min_spec = 0.99)) %>%
  ungroup()

overall_best_youden <- best_threshold(overall_performance, "youden")
overall_spec_95     <- best_threshold(overall_performance, "spec", 0.95)
overall_spec_99     <- best_threshold(overall_performance, "spec", 0.99)

message("\nOptimal threshold by Youden index, per lineage:")
print(best_youden)

message("\nThreshold at >=95% specificity, per lineage:")
print(spec_95)

# ---- performance at the thresholds actually used ----------------------------

candidate_performance <- performance %>%
  filter(threshold %in% KMER_THRESHOLDS) %>%
  select(ST, threshold, sensitivity, specificity, PPV, NPV, F1,
         balanced_accuracy, TP, FP, TN, FN) %>%
  arrange(ST, threshold)

overall_candidates <- overall_performance %>%
  filter(threshold %in% KMER_THRESHOLDS)

message("\nPerformance at the configured thresholds:")
print(candidate_performance)

# ---- misclassified pairs at the primary threshold ---------------------------

false_positives <- calib %>%
  filter(kmer <= PRIMARY_KMER, !SNP_link) %>%
  arrange(ST, desc(SNP), kmer)

false_negatives <- calib %>%
  filter(SNP_link, kmer > PRIMARY_KMER) %>%
  arrange(ST, SNP, kmer)

message(
  "\nAt k-mer <= ", PRIMARY_KMER, ": ",
  nrow(false_positives), " false positives, ",
  nrow(false_negatives), " false negatives"
)

# ---- reviewer summary table -------------------------------------------------

reviewer_table <- summary_ST %>%
  select(ST, isolates, pairs, SNP_le_ref, SNP_gt_ref, Spearman_rho) %>%
  left_join(
    best_youden %>%
      select(ST,
             optimal_kmer              = threshold,
             optimal_sensitivity       = sensitivity,
             optimal_specificity       = specificity,
             optimal_PPV               = PPV,
             optimal_balanced_accuracy = balanced_accuracy),
    by = "ST"
  ) %>%
  left_join(
    candidate_performance %>%
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

message("\nReviewer summary table:")
print(reviewer_table)

# ---- figures ----------------------------------------------------------------

p_scatter <- ggplot(calib, aes(SNP, kmer)) +
  geom_point(alpha = 0.20, size = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.5) +
  geom_vline(xintercept = SNP_THRESHOLD, linetype = "dashed") +
  geom_hline(yintercept = PRIMARY_KMER, linetype = "dotted") +
  facet_wrap(~ ST, scales = "free") +
  labs(x = "Pairwise SNP distance", y = "Pairwise k-mer distance") +
  theme_bw()

p_roc <- performance %>%
  filter(threshold <= max(KMER_THRESHOLDS) * 3) %>%
  select(ST, threshold, sensitivity, specificity) %>%
  pivot_longer(c(sensitivity, specificity),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(threshold, value, linetype = metric)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = PRIMARY_KMER, linetype = "dotted") +
  facet_wrap(~ ST) +
  labs(x = "k-mer distance threshold",
       y = "Classification performance",
       linetype = NULL) +
  theme_bw()


# =============================================================================
# 05  PART B HELPERS
# =============================================================================

#' Build the seqTrack network at one k-mer threshold.
#'
#' seqTrack returns integer positional indices in its `id` and `ances` columns,
#' and orients every edge from ancestor to descendant by collection date. The
#' resulting adjacency is therefore asymmetric, with roughly half the edges
#' below the diagonal. Both the pair extraction and igraph's undirected
#' conversion read only the upper triangle, so the matrix is symmetrised here.
#' Without this step a genuine transmission edge stored below the diagonal is
#' dropped from the network and then counted in the unlinked comparison group.
build_network <- function(dmat, dates, max_distance) {

  ids <- colnames(dmat)
  stopifnot(length(dates) == length(ids))

  st <- seqTrack(dmat, x.names = ids, x.dates = dates)

  edges <- st %>%
    as_tibble() %>%
    filter(!is.na(ances), !is.na(weight), weight <= max_distance)

  adj <- matrix(0L, length(ids), length(ids), dimnames = list(ids, ids))

  if (nrow(edges) > 0) {
    adj[cbind(edges$ances, edges$id)] <- 1L
  }

  adj <- pmax(adj, t(adj))          # symmetrise
  diag(adj) <- 0L

  linked  <- rowSums(adj) > 0
  adj_sub <- adj[linked, linked, drop = FALSE]

  list(
    edges     = edges,
    adjacency = adj_sub,
    graph     = graph_from_adjacency_matrix(
      adj_sub, mode = "undirected", diag = FALSE
    )
  )
}


#' Upper-triangle pairs of a symmetric adjacency matrix carrying `value`.
extract_pairs <- function(adj, value = 1L) {

  idx <- which(adj == value & upper.tri(adj), arr.ind = TRUE)

  tibble(
    KAUST_ID1 = rownames(adj)[idx[, 1]],
    KAUST_ID2 = colnames(adj)[idx[, 2]]
  )
}


eq_flag <- function(a, b) as.integer(!is.na(a) & !is.na(b) & a == b)


#' Annotate isolate pairs with ward, admission-window and patient overlap.
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
    )

  # One admission record per isolate. Without this a many-to-many join
  # silently multiplies the pair counts that every odds ratio below is
  # computed from.
  if (anyDuplicated(vis$KAUST_ID) > 0) {
    warning(
      "visits_at_collection has ", sum(duplicated(vis$KAUST_ID)),
      " duplicated isolate IDs. Keeping the first record for each.",
      call. = FALSE
    )
    vis <- distinct(vis, KAUST_ID, .keep_all = TRUE)
  }

  suffixed <- function(df, s) rename_with(df, ~ paste0(.x, s), -KAUST_ID)

  pairs %>%
    left_join(suffixed(vis, "_1"), by = c("KAUST_ID1" = "KAUST_ID")) %>%
    left_join(suffixed(vis, "_2"), by = c("KAUST_ID2" = "KAUST_ID")) %>%
    mutate(
      same_patient   = eq_flag(MRN_1, MRN_2),
      same_region    = eq_flag(FCLT_NM_1, FCLT_NM_2),
      same_ward_adm  = eq_flag(ADS_WD_1, ADS_WD_2),
      same_ward_dis  = eq_flag(DS_WD_1, DS_WD_2),
      same_ward_cross = as.integer(
        eq_flag(ADS_WD_1, DS_WD_2) == 1 | eq_flag(ADS_WD_2, DS_WD_1) == 1
      ),
      same_ward_any = as.integer(
        same_ward_adm == 1 | same_ward_dis == 1 | same_ward_cross == 1
      ),
      start1 = pmin(ADS_DT_1, DS_DT_1),
      end1   = pmax(ADS_DT_1, DS_DT_1),
      start2 = pmin(ADS_DT_2, DS_DT_2),
      end2   = pmax(ADS_DT_2, DS_DT_2),
      overlap = as.integer(
        !is.na(start1) & !is.na(end1) & !is.na(start2) & !is.na(end2) &
          pmax(start1, start2) <= pmin(end1, end2)
      ),
      overlap_same_ward = as.integer(overlap == 1 & same_ward_any == 1)
    )
}


#' Fisher odds ratio for a binary flag between linked and unlinked pairs.
#' Both variables are coerced to factors with explicit levels so the 2x2 table
#' orientation, and therefore the direction of the odds ratio, cannot depend on
#' which categories happen to be present.
fisher_or <- function(df, group_col, flag_col,
                      group_levels = c("linked", "unlinked")) {

  tab <- table(
    factor(df[[group_col]], levels = group_levels),
    factor(df[[flag_col]],  levels = c(1, 0))
  )

  empty <- tibble(
    variable = flag_col, OR = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_
  )

  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) return(empty)

  ft <- fisher.test(tab)

  tibble(
    variable = flag_col,
    OR       = unname(ft$estimate),
    CI_low   = ft$conf.int[1],
    CI_high  = ft$conf.int[2],
    p_value  = ft$p.value
  )
}


#' One-versus-rest enrichment of each level of `level_col` in the focus group.
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

    if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
      return(tibble(
        level = lv, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        p_value = NA_real_,
        a = tab[1, 1], b = tab[2, 1], c = tab[1, 2], d = tab[2, 2]
      ))
    }

    ft <- fisher.test(tab)

    tibble(
      level   = lv,
      OR      = unname(ft$estimate),
      CI_low  = ft$conf.int[1],
      CI_high = ft$conf.int[2],
      p_value = ft$p.value,
      a = tab[1, 1], b = tab[2, 1], c = tab[1, 2], d = tab[2, 2]
    )
  })
}


EPI_PREDICTORS <- c(
  "same_ward_any", "same_ward_dis", "same_ward_adm",
  "overlap_same_ward", "overlap"
)


#' Rerun the whole Figure B analysis at one k-mer threshold.
run_threshold_analysis <- function(threshold) {

  message("\n-- k-mer threshold <= ", threshold, " --")

  net <- build_network(
    kmers_dedup,
    dates_for(colnames(kmers_dedup)),
    threshold
  )

  if (vcount(net$graph) == 0) {
    warning("No linked isolates at threshold ", threshold, call. = FALSE)
    return(NULL)
  }

  comp <- components(net$graph)$membership

  comp_lookup <- tibble(
    KAUST_ID = names(comp),
    comp     = as.integer(comp)
  )

  # ---- isolate level ----

  isolate_table <- tibble(KAUST_ID = colnames(kmers_dedup)) %>%
    mutate(
      cluster_status = factor(
        if_else(KAUST_ID %in% names(comp), "Cluster", "None"),
        levels = c("None", "Cluster")
      )
    ) %>%
    left_join(kleborate, by = c("KAUST_ID" = "strain")) %>%
    left_join(select(metadata, KAUST_ID, Acquisition), by = "KAUST_ID") %>%
    left_join(pathotype_meta, by = c("KAUST_ID" = "KAUST_ID_total"))

  or_pathotype <- isolate_table %>%
    filter(!is.na(pathotype), pathotype != "Convergent hvKp") %>%
    enrichment_by_level("cluster_status", "pathotype", focus = "Cluster") %>%
    mutate(
      p_adj_BH  = p.adjust(p_value, method = "BH"),
      threshold = threshold,
      analysis  = "Pathotype"
    )

  or_acquisition <- isolate_table %>%
    filter(!is.na(Acquisition), Acquisition != "other") %>%
    enrichment_by_level("cluster_status", "Acquisition", focus = "Cluster") %>%
    mutate(
      p_adj_BH  = p.adjust(p_value, method = "BH"),
      threshold = threshold,
      analysis  = "Acquisition"
    )

  # ---- pair level ----
  #
  # The comparison is between isolate pairs joined by a seqTrack edge and pairs
  # of network members that are not directly joined. Both groups are drawn from
  # the linked subgraph, so the contrast is edge versus no edge among clustered
  # isolates rather than clustered versus unclustered isolates.

  annotated <- bind_rows(
    annotate_pairs(extract_pairs(net$adjacency, 1L), visits_at_collection) %>%
      mutate(pair_type = "linked"),
    annotate_pairs(extract_pairs(net$adjacency, 0L), visits_at_collection) %>%
      mutate(pair_type = "unlinked")
  ) %>%
    filter(!is.na(MRN_1), !is.na(MRN_2))

  or_epi <- map_dfr(
    EPI_PREDICTORS,
    ~ fisher_or(annotated, "pair_type", .x)
  ) %>%
    mutate(
      p_adj_BH  = p.adjust(p_value, method = "BH"),
      threshold = threshold,
      analysis  = "Epidemiology"
    )

  # ---- network summary ----

  summary_tbl <- tibble(
    threshold                  = threshold,
    n_edges                    = nrow(net$edges),
    n_clustered_isolates       = sum(isolate_table$cluster_status == "Cluster"),
    n_unclustered_isolates     = sum(isolate_table$cluster_status == "None"),
    n_components               = length(unique(comp)),
    largest_component          = max(table(comp)),
    linked_pairs_with_visits   = sum(annotated$pair_type == "linked"),
    unlinked_pairs_with_visits = sum(annotated$pair_type == "unlinked")
  )

  message(
    "   edges ", summary_tbl$n_edges,
    " | clustered ", summary_tbl$n_clustered_isolates,
    " | components ", summary_tbl$n_components,
    " | largest ", summary_tbl$largest_component
  )

  list(
    threshold       = threshold,
    network         = net,
    isolate_table   = isolate_table,
    annotated_pairs = annotated,
    pathotype       = or_pathotype,
    acquisition     = or_acquisition,
    epidemiology    = or_epi,
    summary         = summary_tbl
  )
}


# =============================================================================
# 06  PART B  THRESHOLD SENSITIVITY OF THE NETWORK ANALYSIS
# =============================================================================

message("\n== Part B: network sensitivity across k-mer thresholds ==")

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

message("\nESBL(+)/CP(+):")
print(esbl_sensitivity)

message("\nEpidemiological overlap:")
print(
  epidemiology_sensitivity %>%
    select(threshold, variable, OR, CI_low, CI_high, p_value, p_adj_BH) %>%
    arrange(variable, threshold)
)

# ---- figures ----------------------------------------------------------------

or_plot <- function(dat, facet_var = NULL, title = NULL) {

  p <- ggplot(dat, aes(factor(threshold), OR)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15) +
    geom_point(size = 2.5) +
    scale_y_log10() +
    labs(x = "k-mer distance threshold", y = "Odds ratio", title = title) +
    theme_bw()

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(vars(.data[[facet_var]]), scales = "free_y")
  }

  p
}

p_nosocomial   <- or_plot(nosocomial_sensitivity,
                          title = "Nosocomial acquisition")
p_pathotype    <- or_plot(pathotype_sensitivity, "level",
                          title = "Pathotype associations")
p_epi          <- or_plot(epidemiology_sensitivity, "variable",
                          title = "Epidemiological associations")

p_cluster_size <- ggplot(network_sensitivity,
                         aes(threshold, n_clustered_isolates)) +
  geom_line() +
  geom_point(size = 3) +
  geom_vline(xintercept = PRIMARY_KMER, linetype = "dashed") +
  scale_x_continuous(breaks = KMER_THRESHOLDS) +
  labs(x = "k-mer distance threshold",
       y = "Number of clustered isolates",
       title = "Network size sensitivity") +
  theme_bw()


# =============================================================================
# 07  PART C  PUBLICATION SUPPLEMENTARY TABLE
# =============================================================================

message("\n== Part C: publication table ==")

fmt_or <- function(or, lo, hi) {
  if_else(
    is.na(or),
    "—",
    sprintf("%.2f (%.2f–%.2f)", or, lo, hi)
  )
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
  EPI_LABELS[["same_ward_any"]],
  EPI_LABELS[["same_ward_adm"]],
  EPI_LABELS[["same_ward_dis"]],
  EPI_LABELS[["overlap"]],
  EPI_LABELS[["overlap_same_ward"]]
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
  warning(
    "Categories absent from CATEGORY_ORDER and dropped from the table: ",
    paste(unmapped, collapse = ", "),
    call. = FALSE
  )
}

# Column keys are generated from the thresholds actually analysed, so the table
# follows KMER_THRESHOLDS instead of a hard-coded k2..k8 list.
threshold_keys <- paste0("k", KMER_THRESHOLDS)

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

# Labels state the thresholds as configured. The earlier draft printed these as
# 0.002 to 0.008 while the analysis ran at 2 to 8, so the column headings and
# the numbers underneath described different analyses.
threshold_labels <- set_names(
  paste0("≤", KMER_THRESHOLDS),
  threshold_keys
)
primary_key <- paste0("k", PRIMARY_KMER)

publication_gt <- publication_table %>%
  gt(groupname_col = "Analysis") %>%
  cols_label(.list = c(list(Category = ""), as.list(threshold_labels))) %>%
  tab_spanner(label = "k-mer distance threshold",
              columns = all_of(threshold_keys)) %>%
  tab_header(title = md(paste0(
    "**Sensitivity of associations with putative genomic ",
    "transmission-cluster membership to alternative k-mer distance ",
    "thresholds**"
  ))) %>%
  tab_footnote(
    footnote  = "Primary k-mer distance threshold.",
    locations = cells_column_labels(columns = all_of(primary_key))
  ) %>%
  tab_source_note(source_note = md(paste0(
    "Values are odds ratios (95% confidence intervals) for membership in a ",
    "putative genomic transmission cluster. The primary analysis used a ",
    "k-mer distance threshold of ≤", PRIMARY_KMER, ". Thresholds of ",
    paste(setdiff(KMER_THRESHOLDS, PRIMARY_KMER), collapse = ", "),
    " were evaluated as sensitivity analyses. Em dashes mark comparisons ",
    "with no informative 2x2 table."
  ))) %>%
  cols_align("left", columns = Category) %>%
  cols_align("center", columns = all_of(threshold_keys)) %>%
  cols_width(
    Category ~ px(240),
    all_of(threshold_keys) ~ px(130)
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = all_of(primary_key))
  ) %>%
  tab_options(
    table.font.names             = "Arial",
    table.font.size              = px(12),
    heading.title.font.size      = px(14),
    heading.title.font.weight    = "bold",
    column_labels.font.weight    = "bold",
    row_group.font.weight        = "bold",
    table.border.top.width       = px(2),
    table.border.bottom.width    = px(2),
    column_labels.border.top.width    = px(1),
    column_labels.border.bottom.width = px(1),
    row_group.border.top.width   = px(1),
    data_row.padding             = px(6),
    source_notes.font.size       = px(10)
  )

print(publication_gt)


# =============================================================================
# 08  CLUSTERED PATIENTS AT THE PRIMARY THRESHOLD
# =============================================================================

primary_result <- threshold_results[[paste0("kmer_", PRIMARY_KMER)]]

patient_id_col <- intersect(c("MRN", "PT_NO"), names(metadata))

if (length(patient_id_col) == 0) {

  warning(
    "metadata has no MRN or PT_NO column, so clustered patients cannot be ",
    "counted.",
    call. = FALSE
  )
  clustered_patients <- tibble()

} else {

  clustered_patients <- primary_result$isolate_table %>%
    select(KAUST_ID, cluster_status) %>%
    left_join(
      metadata %>%
        select(KAUST_ID, patient_id = all_of(patient_id_col[1])),
      by = "KAUST_ID"
    ) %>%
    filter(!is.na(patient_id))

  n_clustered <- clustered_patients %>%
    filter(cluster_status == "Cluster") %>%
    distinct(patient_id) %>%
    nrow()

  # The denominator is derived rather than hard-coded, so it cannot drift away
  # from the cohort the manuscript describes.
  n_total <- n_distinct(clustered_patients$patient_id)

  message(
    "\nClustered patients: ", n_clustered,
    " of ", n_total,
    " (", round(100 * n_clustered / n_total, 1), "%)"
  )
}


# =============================================================================
# 09  WRITE OUTPUTS
# =============================================================================

csv_outputs <- list(
  calibration_pairwise_distances   = calib,
  calibration_summary_by_ST        = summary_ST,
  calibration_kmer_around_SNP_ref  = boundary_summary,
  calibration_performance_by_ST    = performance,
  calibration_performance_overall  = overall_performance,
  calibration_best_youden_by_ST    = best_youden,
  calibration_spec95_by_ST         = spec_95,
  calibration_spec99_by_ST         = spec_99,
  calibration_candidate_thresholds = candidate_performance,
  calibration_false_positives      = false_positives,
  calibration_false_negatives      = false_negatives,
  calibration_reviewer_table       = reviewer_table,
  network_summary_by_threshold     = network_sensitivity,
  acquisition_OR_by_threshold      = acquisition_sensitivity,
  pathotype_OR_by_threshold        = pathotype_sensitivity,
  epidemiology_OR_by_threshold     = epidemiology_sensitivity,
  nosocomial_OR_by_threshold       = nosocomial_sensitivity,
  ESBL_CP_OR_by_threshold          = esbl_sensitivity,
  supplementary_table_wide         = publication_table
)

iwalk(csv_outputs, function(x, nm) {
  write_csv(x, file.path(OUT_DIR, paste0(nm, ".csv")))
})

plot_outputs <- list(
  SNP_vs_kmer_by_ST              = list(p_scatter, 10, 7),
  kmer_threshold_performance     = list(p_roc, 10, 7),
  nosocomial_OR_sensitivity      = list(p_nosocomial, 7, 5),
  pathotype_OR_sensitivity       = list(p_pathotype, 10, 5),
  epidemiology_OR_sensitivity    = list(p_epi, 11, 7),
  network_size_sensitivity       = list(p_cluster_size, 7, 5)
)

iwalk(plot_outputs, function(spec, nm) {
  ggsave(file.path(OUT_DIR, paste0(nm, ".pdf")),
         spec[[1]], width = spec[[2]], height = spec[[3]])
  ggsave(file.path(OUT_DIR, paste0(nm, ".png")),
         spec[[1]], width = spec[[2]], height = spec[[3]], dpi = 600)
})

for (ext in c("html", "rtf")) {
  gtsave(
    publication_gt,
    file.path(OUT_DIR, paste0("Table_Sx_kmer_threshold_sensitivity.", ext))
  )
}

message("\nAll outputs written to: ", OUT_DIR)
