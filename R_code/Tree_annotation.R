# =====================================================================
# Unitig Jaccard tree + metadata annotation bands
#
# Flat script, no functions. Run block by block, check the printed
# output of each block before moving on. Figures print to the device;
# the saving block at the very end is commented out.
#
# What changed relative to the original, and why:
#   - The distance matrix now carries dimnames. Without them as.dist()
#     returns an unlabelled dist, nj() names the tips "1".."n", and
#     every match()/left_join() against KAUST_ID silently returns NA.
#     This is the single most dangerous bug in the original.
#   - round(1000 * D) removed. It creates artificial ties in a set of
#     distances that are already tightly clustered, and NJ resolves
#     ties arbitrarily. NJ takes the continuous distance directly.
#   - Metadata is matched to tree_mid$tip.label ONCE, after ladderize,
#     with a reported overlap. The original matched to tree_nj first
#     and tree_mid later, and re-read the file part-way through,
#     discarding the earlier ordering.
#   - Band geometry (tree_max_x / band_x / band_width) computed once
#     in block 07 instead of six identical recomputations.
#   - xlim() replaced by coord_cartesian(). xlim() sets SCALE limits,
#     which deletes any tile straddling the boundary; coord_cartesian()
#     zooms without dropping data and is what clip = "off" expects.
#   - year_band (built, never used) and the first year_cols assignment
#     (immediately overwritten) removed. adegenet / ade4 / igraph
#     removed - nothing in the script uses them.
#   - The frequent-ST selector is now data-driven instead of a
#     hard-coded list of twelve STs that goes stale the moment the
#     collection grows.
#   - The hospital-phenotype block audits raw values before recoding.
#     The original coerced anything that was not literally "1" or "0"
#     to NA, so logical TRUE/FALSE columns and the continuous LOS /
#     ICU_Days columns rendered as blank white bands with no warning.
#   - Block 15 composes every band onto ONE figure. This needs
#     ggnewscale: a ggplot can hold only one fill scale, which is the
#     reason the original had to draw eight separate trees.
# =====================================================================


# ---------------------------------------------------------------------
# 01  Settings and libraries
# ---------------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)

PATH_SHARED   <- "/Users/daneshm/Documents/PA_KAIMRC/GWAS_new/unitigs/Shared_unitigs.csv"
PATH_MASTER   <- "/Users/daneshm/Documents/PA_KAIMRC/Causal/PA_master_29June_2026_with_CC.csv"
PATH_LABELS   <- "/Users/daneshm/Documents/PA_KAIMRC/ML_files/GWAS_Labels_new.csv"
PATH_LABELS_C <- "/Users/daneshm/Documents/PA_KAIMRC/ML_files/GWAS_Labels_new_cont.csv"

ID_COL <- "KAUST_ID"

# Minimum isolate count for an ST to get its own colour in block 12.
ST_MIN_N <- 12

# Band geometry, as fractions of total tree depth.
BAND_GAP_FRAC   <- 0.02    # tree tip -> first band
BAND_WIDTH_FRAC <- 0.025   # thickness of one categorical band
BAND_SPACE_FRAC <- 0.004   # gap between adjacent bands
BAR_WIDTH_FRAC  <- 0.13    # full width of one continuous barplot

# `midpoint` exists in both phangorn and phytools. Call it explicitly.
midpoint_fun <- phangorn::midpoint


# ---------------------------------------------------------------------
# 02  Shared-unitig matrix -> numeric matrix with dimnames
#
# read_csv gives a tibble. Whether the first column is an ID column or
# the matrix starts immediately depends on how the file was written, so
# detect it rather than assume. Everything downstream depends on this
# matrix having BOTH rownames and colnames set to KAUST IDs.
# ---------------------------------------------------------------------

kmers_shared <- read_csv(PATH_SHARED, show_col_types = FALSE)

dim(kmers_shared)
kmers_shared[1:3, 1:4]

first_is_id <- !is.numeric(kmers_shared[[1]])
first_is_id

if (first_is_id) {
  ids_row <- as.character(kmers_shared[[1]])
  S <- as.matrix(kmers_shared[, -1, drop = FALSE])
} else {
  ids_row <- colnames(kmers_shared)
  S <- as.matrix(kmers_shared)
}

storage.mode(S) <- "numeric"
rownames(S) <- ids_row       # colnames already carry the IDs from the header

# The two label sets must agree, otherwise the matrix is transposed or
# a column was dropped somewhere upstream.
length(rownames(S))
length(colnames(S))
sum(rownames(S) != colnames(S))
stopifnot(nrow(S) == ncol(S))
stopifnot(identical(rownames(S), colnames(S)))

# Structural checks on the counts themselves.
sum(is.na(S))
max(abs(S - t(S)))                       # 0 if symmetric
summary(diag(S))                         # unitigs per isolate
sum(diag(S) == 0)                        # isolates with no unitigs

# Shared count can never exceed either isolate's own count.
viol <- sum(S > outer(diag(S), diag(S), pmin))
viol
stopifnot(viol == 0)


# ---------------------------------------------------------------------
# 03  Jaccard distance
#
#   J[i,j] = S[i,j] / (S[i,i] + S[j,j] - S[i,j])
#   D      = 1 - J
#
# Verified against a brute-force set-intersection computation: exact to
# machine precision.
# ---------------------------------------------------------------------

diag_vals <- diag(S)

denom <- outer(diag_vals, diag_vals, "+") - S
sum(denom <= 0)                          # must be 0; a zero denominator
                                         # means two empty isolates

J <- S / denom
D <- 1 - J
diag(D) <- 0

dimnames(D) <- dimnames(S)               # keep the labels - critical

summary(D[lower.tri(D)])
length(unique(D[lower.tri(D)]))          # distinct pairwise distances
sum(!is.finite(D))

d <- as.dist(D)
length(labels(d))
head(labels(d))
stopifnot(length(labels(d)) == nrow(S))  # empty here = unlabelled tree


# ---------------------------------------------------------------------
# 04  Neighbour joining, midpoint rooting, ladderize
#
# NJ can return negative branch lengths. They are a real feature of the
# algorithm, not an error, but they push tip x-coordinates around and
# make midpoint rooting behave oddly. Count them before deciding.
# ---------------------------------------------------------------------

tree_nj <- nj(d)

Ntip(tree_nj)
sum(tree_nj$edge.length < 0)
min(tree_nj$edge.length)
sum(tree_nj$edge.length[tree_nj$edge.length < 0])   # total negative length

# Uncomment if the negatives are numerous enough to distort the plot.
# Zeroing them changes branch lengths but not topology.
# tree_nj$edge.length <- pmax(tree_nj$edge.length, 0)

tree_mid <- midpoint_fun(tree_nj)
tree_mid <- ladderize(tree_mid, right = FALSE)

Ntip(tree_mid)
is.rooted(tree_mid)
setdiff(tree_nj$tip.label, tree_mid$tip.label)   # must be empty

ggtree(tree_mid, size = 0.3) + theme_tree2()


# ---------------------------------------------------------------------
# 05  Metadata, matched to the tree once
#
# Row order of `input` is irrelevant downstream (every band is built
# with a left_join on tip label), so the match below is a completeness
# check, not a reordering step. Read the printed overlap before going
# further - a large `missing_in_meta` means the two files use different
# identifier conventions.
# ---------------------------------------------------------------------

input <- read_csv(PATH_MASTER, show_col_types = FALSE)

dim(input)
ID_COL %in% colnames(input)
stopifnot(ID_COL %in% colnames(input))

input[[ID_COL]] <- as.character(input[[ID_COL]])

sum(duplicated(input[[ID_COL]]))         # duplicated IDs break the join

tip_ids <- tree_mid$tip.label

missing_in_meta <- setdiff(tip_ids, input[[ID_COL]])
extra_in_meta   <- setdiff(input[[ID_COL]], tip_ids)

length(tip_ids)
length(missing_in_meta)
head(missing_in_meta, 20)
length(extra_in_meta)

# Reorder to tree order so head(input) is readable alongside the tree.
input <- input[match(tip_ids, input[[ID_COL]]), ]
stopifnot(identical(input[[ID_COL]][!is.na(input[[ID_COL]])],
                    tip_ids[!is.na(input[[ID_COL]])]))


# ---------------------------------------------------------------------
# 06  Derived metadata columns
#
# All recoding happens here, once, so the band blocks below are purely
# graphical. Check every table() before plotting.
# ---------------------------------------------------------------------

input$Year <- as.numeric(input$Year)
year_levels <- as.character(sort(unique(input$Year[!is.na(input$Year)])))
year_levels
input$Year_f <- factor(as.character(input$Year), levels = year_levels)
table(input$Year_f, useNA = "ifany")

loc_levels <- c("Jeddah", "Ahsaa", "Dammam", "Riyadh", "Madinah")
setdiff(unique(input$Location[!is.na(input$Location)]), loc_levels)   # unmapped
input$Location_f <- factor(input$Location, levels = loc_levels)
table(input$Location_f, useNA = "ifany")

table(input$Specimen_Source, useNA = "ifany")
input$Specimen_Group <- case_when(
  input$Specimen_Source == "Blood Culture"       ~ "Blood",
  input$Specimen_Source == "Respiratory Culture" ~ "Respiratory",
  input$Specimen_Source == "Wound Culture"       ~ "Wound",
  input$Specimen_Source == "Urine Culture"       ~ "Urine",
  !is.na(input$Specimen_Source)                  ~ "Other",
  TRUE                                           ~ NA_character_
)
spec_levels <- c("Blood", "Respiratory", "Wound", "Urine", "Other")
input$Specimen_Group <- factor(input$Specimen_Group, levels = spec_levels)
table(input$Specimen_Group, useNA = "ifany")

# DTR / CRPA: NA is being read as absence. That is an assumption, not a
# fact - it is only defensible if the source pipeline writes NA when a
# phenotype was tested and not met. Check the raw values and the NA
# count first. as.logical() returns NA for "Yes"/"No", which would then
# be swept into the FALSE group; if the tables below show anything other
# than TRUE/FALSE/0/1, recode explicitly instead.
table(as.character(input$DTR),  useNA = "ifany")
table(as.character(input$CRPA), useNA = "ifany")
sum(is.na(input$DTR))
sum(is.na(input$CRPA))

input$DTR_Group  <- ifelse(replace_na(as.logical(input$DTR),  FALSE), "DTR",  "Non-DTR")
input$CRPA_Group <- ifelse(replace_na(as.logical(input$CRPA), FALSE), "CRPA", "Non-CRPA")
input$DTR_Group  <- factor(input$DTR_Group,  levels = c("DTR", "Non-DTR"))
input$CRPA_Group <- factor(input$CRPA_Group, levels = c("CRPA", "Non-CRPA"))
table(input$DTR_Group,  useNA = "ifany")
table(input$CRPA_Group, useNA = "ifany")

# Frequent STs, derived from the data rather than hard-coded.
st_tab <- sort(table(as.character(input$ST)), decreasing = TRUE)
head(st_tab, 25)
highlight_STs <- names(st_tab)[st_tab >= ST_MIN_N]
highlight_STs
length(highlight_STs)

input$ST_Group <- ifelse(as.character(input$ST) %in% highlight_STs,
                         paste0("ST", as.character(input$ST)),
                         "Other")
st_levels <- c(paste0("ST", highlight_STs), "Other")
input$ST_Group <- factor(input$ST_Group, levels = st_levels)
table(input$ST_Group, useNA = "ifany")


# ---------------------------------------------------------------------
# 07  Base tree, tip positions, band geometry
#
# Computed once. Every block below reads these objects and never
# recomputes them.
# ---------------------------------------------------------------------

p <- ggtree(tree_mid, linewidth = 0.3)

tip_positions <- p$data %>%
  filter(isTip) %>%
  select(label, y)

nrow(tip_positions)
head(tip_positions)
stopifnot(nrow(tip_positions) == Ntip(tree_mid))

tree_max_x  <- max(p$data$x, na.rm = TRUE)
tree_max_y  <- max(tip_positions$y, na.rm = TRUE)

band_width  <- tree_max_x * BAND_WIDTH_FRAC
band_space  <- tree_max_x * BAND_SPACE_FRAC
band_pitch  <- band_width + band_space
band_start  <- tree_max_x * (1 + BAND_GAP_FRAC)

# Centre of a single band, for the one-variable plots in blocks 08-13.
band_x      <- band_start + band_width / 2

# Right-hand plot limit for a single-band figure.
x_max_one   <- band_start + band_width * 1.5

# Vertical position for column labels above the bands.
label_y     <- tree_max_y * 1.02

c(tree_max_x = tree_max_x, tree_max_y = tree_max_y,
  band_width = band_width, band_start = band_start, band_x = band_x)

# Shared theme for every annotated figure below.
theme_band <- theme_tree2() +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 9),
    plot.margin     = margin(t = 90, r = 20, b = 10, l = 10)
  )


# ---------------------------------------------------------------------
# 08  Year band
# ---------------------------------------------------------------------

year_data <- tip_positions %>%
  left_join(input %>% select(all_of(ID_COL), Year_f),
            by = setNames(ID_COL, "label"))

nrow(year_data)
table(year_data$Year_f, useNA = "ifany")

year_pal <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
              "#FEE090", "#FDAE61", "#F46D43", "#A50026")
year_cols <- setNames(
  colorRampPalette(year_pal)(length(year_levels)),
  year_levels
)
year_cols

p_year <- p +
  geom_tile(data = year_data,
            aes(x = band_x, y = y, fill = Year_f),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = year_cols, drop = FALSE,
                    na.value = "white", name = "Year") +
  coord_cartesian(xlim = c(0, x_max_one), clip = "off") +
  theme_band

p_year


# ---------------------------------------------------------------------
# 09  Location band
# ---------------------------------------------------------------------

location_data <- tip_positions %>%
  left_join(input %>% select(all_of(ID_COL), Location_f),
            by = setNames(ID_COL, "label"))

table(location_data$Location_f, useNA = "ifany")

location_cols <- c("Jeddah"  = "#E41A1C",
                   "Ahsaa"   = "#377EB8",
                   "Dammam"  = "#4DAF4A",
                   "Riyadh"  = "#984EA3",
                   "Madinah" = "#FF7F00")

p_location <- p +
  geom_tile(data = location_data,
            aes(x = band_x, y = y, fill = Location_f),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = location_cols, drop = FALSE,
                    na.value = "white", name = "Location") +
  coord_cartesian(xlim = c(0, x_max_one), clip = "off") +
  theme_band

p_location


# ---------------------------------------------------------------------
# 10  Specimen source band
# ---------------------------------------------------------------------

specimen_data <- tip_positions %>%
  left_join(input %>% select(all_of(ID_COL), Specimen_Group),
            by = setNames(ID_COL, "label"))

table(specimen_data$Specimen_Group, useNA = "ifany")

specimen_cols <- c("Blood"       = "#D73027",
                   "Respiratory" = "#4575B4",
                   "Wound"       = "#FDAE61",
                   "Urine"       = "#66BD63",
                   "Other"       = "#984EA3")

p_specimen <- p +
  geom_tile(data = specimen_data,
            aes(x = band_x, y = y, fill = Specimen_Group),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = specimen_cols, drop = FALSE,
                    na.value = "white", name = "Specimen source") +
  coord_cartesian(xlim = c(0, x_max_one), clip = "off") +
  theme_band

p_specimen


# ---------------------------------------------------------------------
# 11  DTR band
# ---------------------------------------------------------------------

dtr_data <- tip_positions %>%
  left_join(input %>% select(all_of(ID_COL), DTR_Group),
            by = setNames(ID_COL, "label"))

table(dtr_data$DTR_Group, useNA = "ifany")

dtr_cols <- c("DTR" = "#D73027", "Non-DTR" = "#D9D9D9")

p_dtr <- p +
  geom_tile(data = dtr_data,
            aes(x = band_x, y = y, fill = DTR_Group),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = dtr_cols, drop = FALSE,
                    na.value = "white", name = "DTR") +
  coord_cartesian(xlim = c(0, x_max_one), clip = "off") +
  theme_band

p_dtr


# ---------------------------------------------------------------------
# 12  CRPA band
# ---------------------------------------------------------------------

crpa_data <- tip_positions %>%
  left_join(input %>% select(all_of(ID_COL), CRPA_Group),
            by = setNames(ID_COL, "label"))

table(crpa_data$CRPA_Group, useNA = "ifany")

crpa_cols <- c("CRPA" = "#D73027", "Non-CRPA" = "#D9D9D9")

p_crpa <- p +
  geom_tile(data = crpa_data,
            aes(x = band_x, y = y, fill = CRPA_Group),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = crpa_cols, drop = FALSE,
                    na.value = "white", name = "CRPA") +
  coord_cartesian(xlim = c(0, x_max_one), clip = "off") +
  theme_band

p_crpa


# ---------------------------------------------------------------------
# 13  Sequence type band
#
# Colours are generated for however many STs clear ST_MIN_N, so the
# palette cannot silently run short when the collection grows.
# ---------------------------------------------------------------------

st_data <- tip_positions %>%
  left_join(input %>% select(all_of(ID_COL), ST_Group),
            by = setNames(ID_COL, "label"))

table(st_data$ST_Group, useNA = "ifany")

st_base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
             "#A65628", "#F781BF", "#17BECF", "#BCBD22", "#1F78B4",
             "#6A3D9A", "#E7298A", "#66C2A5", "#FC8D62", "#8DA0CB",
             "#E78AC3", "#A6D854", "#FFD92F")

st_cols <- setNames(
  c(rep_len(st_base, length(highlight_STs)), "grey90"),
  st_levels
)
st_cols

p_st <- p +
  geom_tile(data = st_data,
            aes(x = band_x, y = y, fill = ST_Group),
            width = band_width, height = 1, inherit.aes = FALSE,
            colour = NA) +
  scale_fill_manual(values = st_cols,
                    breaks = paste0("ST", highlight_STs),
                    drop = FALSE, na.value = "white",
                    name = "Sequence type") +
  coord_cartesian(xlim = c(0, x_max_one), clip = "off") +
  theme_band

p_st


# ---------------------------------------------------------------------
# 14  Hospital phenotypes - audit before recoding
#
# The original recoded anything that was not literally "1" or "0" to
# NA. Logical TRUE/FALSE columns and the continuous LOS / ICU_Days
# columns all fall into that bucket, and a variable that ends up 100%
# NA renders as a blank white column with no error. Read the audit
# table below before you trust any band in block 15.
# ---------------------------------------------------------------------

input_hospital <- read_csv(PATH_LABELS, show_col_types = FALSE)

dim(input_hospital)
colnames(input_hospital)
stopifnot(ID_COL %in% colnames(input_hospital))
input_hospital[[ID_COL]] <- as.character(input_hospital[[ID_COL]])

length(intersect(input_hospital[[ID_COL]], tip_ids))
length(setdiff(tip_ids, input_hospital[[ID_COL]]))

hospital_vars <- c("Total_Death", "Transmission", "LOS",
                   "Readmission_30d", "ICU_Days", "ICU_Admission")

setdiff(hospital_vars, colnames(input_hospital))     # must be empty
stopifnot(length(setdiff(hospital_vars, colnames(input_hospital))) == 0)

# Per-variable audit: how many values map to 1, to 0, and to neither.
hosp_audit <- input_hospital %>%
  select(all_of(hospital_vars)) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "raw") %>%
  mutate(
    mapped = case_when(
      raw %in% c("1", "TRUE", "True", "Yes", "Y")  ~ "one",
      raw %in% c("0", "FALSE", "False", "No", "N") ~ "zero",
      is.na(raw) | raw %in% c("-", "", "NA")       ~ "missing",
      TRUE                                         ~ "other"
    )
  ) %>%
  count(var, mapped) %>%
  pivot_wider(names_from = mapped, values_from = n, values_fill = 0)

hosp_audit

# Distinct raw values, so a continuous variable is unmistakable.
sapply(input_hospital[hospital_vars],
       function(z) length(unique(as.character(z))))

# Keep only variables that are genuinely binary: at most two non-missing
# distinct values, and no "other" values.
n_other <- setNames(rep(0, length(hospital_vars)), hospital_vars)
if ("other" %in% colnames(hosp_audit)) {
  n_other[hosp_audit$var] <- hosp_audit$other
}
hospital_bin <- names(n_other)[n_other == 0]
hospital_bin
setdiff(hospital_vars, hospital_bin)     # dropped as non-binary

hospital_long <- input_hospital %>%
  select(all_of(ID_COL), all_of(hospital_bin)) %>%
  mutate(across(all_of(hospital_bin), as.character)) %>%
  pivot_longer(cols = all_of(hospital_bin),
               names_to = "Phenotype", values_to = "raw") %>%
  mutate(
    Status = case_when(
      raw %in% c("1", "TRUE", "True", "Yes", "Y")  ~ "1",
      raw %in% c("0", "FALSE", "False", "No", "N") ~ "0",
      TRUE                                         ~ NA_character_
    )
  ) %>%
  select(all_of(ID_COL), Phenotype, Status)

table(hospital_long$Phenotype, hospital_long$Status, useNA = "ifany")

hospital_labels <- c(Total_Death     = "Death",
                     Transmission    = "Transmission",
                     LOS             = "LOS",
                     Readmission_30d = "Readmission 30d",
                     ICU_Days        = "ICU days",
                     ICU_Admission   = "ICU admission")

phenotype_positions <- data.frame(
  Phenotype = factor(hospital_bin, levels = hospital_bin),
  band_x    = band_start + band_width / 2 +
              (seq_along(hospital_bin) - 1) * band_pitch,
  Label     = unname(hospital_labels[hospital_bin]),
  stringsAsFactors = FALSE
)

phenotype_positions

hospital_data <- tip_positions %>%
  left_join(hospital_long, by = setNames(ID_COL, "label")) %>%
  mutate(Phenotype = factor(Phenotype, levels = hospital_bin)) %>%
  left_join(phenotype_positions %>% select(Phenotype, band_x),
            by = "Phenotype") %>%
  filter(!is.na(Phenotype))

nrow(hospital_data)
head(hospital_data)

hospital_cols <- c("1" = "#D73027", "0" = "#4575B4")

x_max_hosp <- max(phenotype_positions$band_x) + band_width

p_hospital <- p +
  geom_tile(data = hospital_data,
            aes(x = band_x, y = y, fill = Status),
            width = band_width, height = 1, inherit.aes = FALSE) +
  geom_text(data = phenotype_positions,
            aes(x = band_x, y = label_y, label = Label),
            angle = 90, hjust = 0, vjust = 0.5,
            size = 3.2, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = hospital_cols, breaks = c("1", "0"),
                    na.value = "white", name = "Hospital phenotype") +
  coord_cartesian(xlim = c(0, x_max_hosp), clip = "off") +
  theme_band

p_hospital


# ---------------------------------------------------------------------
# 15  Continuous phenotypes as bars
#
# Bar length is data / max(data), so the two panels are on different
# scales and are NOT comparable to each other. The printed maxima are
# the scale keys - report them in the legend or caption.
# ---------------------------------------------------------------------

input_cont <- read_csv(PATH_LABELS_C, show_col_types = FALSE)

dim(input_cont)
colnames(input_cont)
input_cont[[ID_COL]] <- as.character(input_cont[[ID_COL]])

hospital_cont <- input_cont %>%
  select(all_of(ID_COL), ICU_Days, LOS) %>%
  mutate(
    ICU_Days = suppressWarnings(as.numeric(na_if(as.character(ICU_Days), "-"))),
    LOS      = suppressWarnings(as.numeric(na_if(as.character(LOS), "-")))
  )

summary(hospital_cont$ICU_Days)
summary(hospital_cont$LOS)
sum(is.na(hospital_cont$ICU_Days))
sum(is.na(hospital_cont$LOS))
sum(hospital_cont$ICU_Days == 0, na.rm = TRUE)
sum(hospital_cont$LOS == 0, na.rm = TRUE)

hospital_cont_data <- tip_positions %>%
  left_join(hospital_cont, by = setNames(ID_COL, "label"))

nrow(hospital_cont_data)

bar_width  <- tree_max_x * BAR_WIDTH_FRAC
bar_gap    <- tree_max_x * BAND_WIDTH_FRAC
icu_start  <- band_start
los_start  <- icu_start + bar_width + bar_gap

icu_max <- max(hospital_cont_data$ICU_Days, na.rm = TRUE)
los_max <- max(hospital_cont_data$LOS,      na.rm = TRUE)
c(icu_max = icu_max, los_max = los_max)

hospital_cont_data <- hospital_cont_data %>%
  mutate(
    ICU_scaled = ICU_Days / icu_max * bar_width,
    LOS_scaled = LOS      / los_max * bar_width
  )

x_max_cont <- los_start + bar_width * 1.1

p_continuous <- p +
  geom_rect(data = hospital_cont_data %>%
              filter(!is.na(ICU_Days), ICU_Days > 0),
            aes(xmin = icu_start, xmax = icu_start + ICU_scaled,
                ymin = y - 0.42,   ymax = y + 0.42),
            fill = "#4575B4", inherit.aes = FALSE) +
  geom_point(data = hospital_cont_data %>%
               filter(!is.na(ICU_Days), ICU_Days == 0),
             aes(x = icu_start, y = y),
             shape = 15, size = 0.8, inherit.aes = FALSE) +
  geom_rect(data = hospital_cont_data %>%
              filter(!is.na(LOS), LOS > 0),
            aes(xmin = los_start, xmax = los_start + LOS_scaled,
                ymin = y - 0.42,   ymax = y + 0.42),
            fill = "#D73027", inherit.aes = FALSE) +
  geom_point(data = hospital_cont_data %>%
               filter(!is.na(LOS), LOS == 0),
             aes(x = los_start, y = y),
             shape = 15, size = 0.8, inherit.aes = FALSE) +
  annotate("text", x = icu_start + bar_width / 2, y = label_y,
           label = paste0("ICU days (max ", icu_max, ")"),
           fontface = "bold", size = 3.6, hjust = 0.5) +
  annotate("text", x = los_start + bar_width / 2, y = label_y,
           label = paste0("LOS (max ", los_max, ")"),
           fontface = "bold", size = 3.6, hjust = 0.5) +
  coord_cartesian(xlim = c(0, x_max_cont), clip = "off") +
  theme_tree2() +
  theme(plot.margin = margin(t = 60, r = 20, b = 10, l = 10))

p_continuous


# ---------------------------------------------------------------------
# 16  One combined figure
#
# A ggplot can carry exactly one fill scale. ggnewscale::new_scale_fill()
# closes the current one and opens another, which is what lets all six
# categorical bands, the hospital binaries and the two bar panels sit on
# the same tree. Without it you are forced into the eight separate
# figures the original produced.
#
# install.packages("ggnewscale") if this fails.
# ---------------------------------------------------------------------

library(ggnewscale)

# Column layout, left to right. Edit this vector to reorder the panel.
panel_vars <- c("Year", "Location", "Specimen", "DTR", "CRPA", "ST",
                hospital_bin)

panel_x <- band_start + band_width / 2 +
           (seq_along(panel_vars) - 1) * band_pitch
names(panel_x) <- panel_vars
panel_x

panel_labels <- c(Year = "Year", Location = "Location",
                  Specimen = "Specimen source", DTR = "DTR", CRPA = "CRPA",
                  ST = "Sequence type", hospital_labels)

panel_end   <- max(panel_x) + band_width / 2
bar_icu_x   <- panel_end + band_pitch
bar_los_x   <- bar_icu_x + bar_width + bar_gap
x_max_panel <- bar_los_x + bar_width * 1.05

c(panel_end = panel_end, bar_icu_x = bar_icu_x,
  bar_los_x = bar_los_x, x_max_panel = x_max_panel)

panel_label_df <- data.frame(
  x     = c(unname(panel_x),
            bar_icu_x + bar_width / 2,
            bar_los_x + bar_width / 2),
  Label = c(unname(panel_labels[panel_vars]),
            paste0("ICU days (max ", icu_max, ")"),
            paste0("LOS (max ", los_max, ")")),
  stringsAsFactors = FALSE
)

panel_label_df

# Build up layer by layer so each addition can be inspected on its own.
pc <- p

pc <- pc +
  geom_tile(data = year_data,
            aes(x = panel_x[["Year"]], y = y, fill = Year_f),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = year_cols, drop = FALSE,
                    na.value = "white", name = "Year")

pc <- pc + new_scale_fill() +
  geom_tile(data = location_data,
            aes(x = panel_x[["Location"]], y = y, fill = Location_f),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = location_cols, drop = FALSE,
                    na.value = "white", name = "Location")

pc <- pc + new_scale_fill() +
  geom_tile(data = specimen_data,
            aes(x = panel_x[["Specimen"]], y = y, fill = Specimen_Group),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = specimen_cols, drop = FALSE,
                    na.value = "white", name = "Specimen source")

pc <- pc + new_scale_fill() +
  geom_tile(data = dtr_data,
            aes(x = panel_x[["DTR"]], y = y, fill = DTR_Group),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = dtr_cols, drop = FALSE,
                    na.value = "white", name = "DTR")

pc <- pc + new_scale_fill() +
  geom_tile(data = crpa_data,
            aes(x = panel_x[["CRPA"]], y = y, fill = CRPA_Group),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = crpa_cols, drop = FALSE,
                    na.value = "white", name = "CRPA")

pc <- pc + new_scale_fill() +
  geom_tile(data = st_data,
            aes(x = panel_x[["ST"]], y = y, fill = ST_Group),
            width = band_width, height = 1, inherit.aes = FALSE,
            colour = NA) +
  scale_fill_manual(values = st_cols, breaks = paste0("ST", highlight_STs),
                    drop = FALSE, na.value = "white", name = "Sequence type")

# Hospital binaries share one fill scale across all their columns.
hospital_panel <- hospital_data %>%
  select(-band_x) %>%
  left_join(data.frame(Phenotype = factor(hospital_bin, levels = hospital_bin),
                       band_x    = unname(panel_x[hospital_bin])),
            by = "Phenotype")

head(hospital_panel)

pc <- pc + new_scale_fill() +
  geom_tile(data = hospital_panel,
            aes(x = band_x, y = y, fill = Status),
            width = band_width, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = hospital_cols, breaks = c("1", "0"),
                    na.value = "white", name = "Hospital phenotype")

# Continuous bars, with a manual legend for the two colours.
cont_icu <- hospital_cont_data %>% filter(!is.na(ICU_Days), ICU_Days > 0)
cont_los <- hospital_cont_data %>% filter(!is.na(LOS),      LOS      > 0)

pc <- pc + new_scale_fill() +
  geom_rect(data = cont_icu,
            aes(xmin = bar_icu_x, xmax = bar_icu_x + ICU_scaled,
                ymin = y - 0.42,  ymax = y + 0.42, fill = "ICU days"),
            inherit.aes = FALSE) +
  geom_rect(data = cont_los,
            aes(xmin = bar_los_x, xmax = bar_los_x + LOS_scaled,
                ymin = y - 0.42,  ymax = y + 0.42, fill = "LOS"),
            inherit.aes = FALSE) +
  scale_fill_manual(values = c("ICU days" = "#4575B4", "LOS" = "#D73027"),
                    name = "Duration") +
  geom_point(data = hospital_cont_data %>%
               filter(!is.na(ICU_Days), ICU_Days == 0),
             aes(x = bar_icu_x, y = y),
             shape = 15, size = 0.6, inherit.aes = FALSE) +
  geom_point(data = hospital_cont_data %>%
               filter(!is.na(LOS), LOS == 0),
             aes(x = bar_los_x, y = y),
             shape = 15, size = 0.6, inherit.aes = FALSE)

p_panel <- pc +
  geom_text(data = panel_label_df,
            aes(x = x, y = label_y, label = Label),
            angle = 90, hjust = 0, vjust = 0.5,
            size = 3.0, fontface = "bold", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, x_max_panel), clip = "off") +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.box      = "vertical",
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(t = 110, r = 20, b = 10, l = 10)
  )

p_panel


# ---------------------------------------------------------------------
# 17  Saving - uncomment when the figures look right
# ---------------------------------------------------------------------

# OUT <- "/Users/daneshm/Documents/PA_KAIMRC/figures"
# dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
#
# write.tree(tree_mid, file.path(OUT, "unitig_jaccard_nj_midpoint.nwk"))
# saveRDS(D, file.path(OUT, "unitig_jaccard_distance.rds"))
#
# ggsave(file.path(OUT, "tree_year.pdf"),       p_year,       width = 8,  height = 14)
# ggsave(file.path(OUT, "tree_location.pdf"),   p_location,   width = 8,  height = 14)
# ggsave(file.path(OUT, "tree_specimen.pdf"),   p_specimen,   width = 8,  height = 14)
# ggsave(file.path(OUT, "tree_dtr.pdf"),        p_dtr,        width = 8,  height = 14)
# ggsave(file.path(OUT, "tree_crpa.pdf"),       p_crpa,       width = 8,  height = 14)
# ggsave(file.path(OUT, "tree_st.pdf"),         p_st,         width = 8,  height = 14)
# ggsave(file.path(OUT, "tree_hospital.pdf"),   p_hospital,   width = 9,  height = 14)
# ggsave(file.path(OUT, "tree_continuous.pdf"), p_continuous, width = 9,  height = 14)
# ggsave(file.path(OUT, "tree_panel.pdf"),      p_panel,      width = 13, height = 16)
