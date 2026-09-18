# ============================================================
# PUTATIVE TRANSMISSION NETWORK
# seqTrack ancestries from summed k-mer distances, plotted
# four times on one fixed layout with different node colours.
# ============================================================

library(readr)
library(dplyr)
library(adegenet)
library(igraph)

base <- "/Users/daneshm/Documents/Kp_KAIMRC"
revision_dir <- file.path(base, "revision")

# Ancestries above this k-mer distance are discarded
MAX_LINK_DISTANCE <- 4

STUDY_YEARS <- as.character(2018:2024)

TOP_STS <- c("ST2096", "ST14", "ST147", "ST307", "ST45", "ST11")

SEED <- 12345

REGION_COLOURS <- c(
  "Central" = "red",
  "Eastern" = "blue",
  "Western" = "green",
  "Madinah" = "yellow"
)

PATHOTYPE_COLOURS <- c(
  "ESBL/CP-negative non-hvKp" = "gray70",
  "hvKp only" = "steelblue",
  "ESBL(+)/CP(+) only" = "red",
  "Convergent hvKp" = "forestgreen"
)

YEAR_COLOURS <- c(
  "2018" = "#c6dbef",
  "2019" = "#9ecae1",
  "2020" = "#6baed6",
  "2021" = "#4292c6",
  "2022" = "#2171b5",
  "2023" = "#08519c",
  "2024" = "#08306b"
)

ST_COLOURS <- c(
  "ST2096" = "#0000FF",
  "ST14" = "#8B0000",
  "ST147" = "#009E73",
  "ST307" = "#F39C12",
  "ST45" = "#56B4E9",
  "ST11" = "#CC79A7",
  "Others" = "#8FA5A8"
)


# ============================================================
# 1. K-MER DISTANCE MATRIX
# ============================================================

kmers <- read_csv(
  file.path(base, "distancematrix", "summed_kmer_distance.csv"),
  show_col_types = FALSE
)

kmers <- round(1000 * kmers)


# ============================================================
# 2. METADATA, PATHOTYPE AND ST GROUP
# ============================================================

metadata_df <- read_csv(
  file.path(revision_dir, "Combined_df_withST.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    ST = sub("-.*$", "", ST_total),
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
    ),
    Date_total = as.Date(Date_total, format = "%d/%m/%Y"),
    Year = format(Date_total, "%Y"),
    ST_group = if_else(ST %in% TOP_STS, ST, "Others")
  ) %>%
  filter(!is.na(Date_total), Year %in% STUDY_YEARS)


# ============================================================
# 3. SUBSET THE DISTANCE MATRIX TO THESE ISOLATES
# ============================================================

isolate_idx <- match(metadata_df$KAUST_ID_total, colnames(kmers))

stopifnot(!anyNA(isolate_idx))

kmers_short <- as.matrix(kmers[isolate_idx, isolate_idx])
isolate_names <- colnames(kmers_short)


# ============================================================
# 4. SEQTRACK ANCESTRIES
# ============================================================

collection_dates <- as.POSIXct(metadata_df$Date_total)

outbreak_result <- seqTrack(
  kmers_short,
  x.names = isolate_names,
  x.dates = collection_dates
)

# Drop links that are too distant to be plausible transmission
outbreak_result <- outbreak_result %>%
  mutate(
    ances = if_else(
      weight >= MAX_LINK_DISTANCE & !is.na(weight),
      NA_real_,
      ances
    )
  ) %>%
  tidyr::drop_na()

dim(outbreak_result)


# ============================================================
# 5. ADJACENCY MATRIX
# ============================================================

n_isolates <- nrow(kmers_short)
adjmat <- matrix(0, nrow = n_isolates, ncol = n_isolates)

for (i in seq_len(nrow(outbreak_result))) {
  if (!is.na(outbreak_result$ances[i])) {
    adjmat[outbreak_result$ances[i], outbreak_result$id[i]] <- 1
  }
}

# Keep only isolates that appear in at least one link
connected <- (rowSums(adjmat) + colSums(adjmat)) != 0

adjmat_binary <- adjmat[connected, connected]
dimnames(adjmat_binary) <- list(
  isolate_names[connected],
  isolate_names[connected]
)


# ============================================================
# 6. GRAPH AND FIXED LAYOUT
#
# The layout is computed once so that all four panels below
# show the same network with the same node positions.
# ============================================================

net <- graph_from_adjacency_matrix(
  adjmat_binary,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# `area` and `repulserad` were removed from layout_with_fr() in
# igraph 0.8.0; niter and start.temp are the remaining controls
set.seed(SEED)
layout_fixed <- layout_with_fr(
  net,
  niter = 1000,
  start.temp = sqrt(vcount(net)) / 10
)


# ============================================================
# 7. PLOT ONE COLOURING OF THE NETWORK
# ============================================================

plot_network <- function(values, palette, title, legend_title) {
  node_colours <- palette[
    values[match(V(net)$name, metadata_df$KAUST_ID_total)]
  ]

  plot(
    net,
    layout = layout_fixed,
    vertex.size = 4,
    vertex.label = NA,
    vertex.color = node_colours,
    edge.width = E(net)$weight,
    edge.color = "gray70",
    main = title
  )

  legend(
    "topleft",
    legend = names(palette),
    fill = palette,
    title = legend_title,
    bty = "n",
    cex = 0.8
  )
}


# ============================================================
# 8. FOUR COLOURINGS
# ============================================================

plot_network(
  metadata_df$RGN_total,
  REGION_COLOURS,
  "Region",
  "Region"
)

plot_network(
  metadata_df$pathotype,
  PATHOTYPE_COLOURS,
  "Pathotype",
  "Pathotype"
)

plot_network(
  metadata_df$Year,
  YEAR_COLOURS,
  "Collection year",
  "Year"
)

plot_network(
  metadata_df$ST_group,
  ST_COLOURS,
  "Sequence type",
  "ST"
)
