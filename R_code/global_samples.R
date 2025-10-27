# ===========================
# Load Required Libraries
# ===========================
library(tidyverse)
library(ape)
library(treedataverse)
library(ggnewscale)
library(rhierbaps)
library(phytools)
library(ggsci)
library(aplot)
library(ggimage)
library(rsvg)
library(cutpointr)
library(ggsignif)

# ===========================
# Load and Merge Metadata
# ===========================
external_meta <- read_tsv("./all_external_meta_clean.tsv")
external_kleborate <- read_tsv("./kleborate_output.tsv")

# Merge metadata with Kleborate outputs
external_meta_kleborate <- external_meta %>%
  left_join(external_kleborate, by = c("Run" = "strain"))

# Count number of isolates per SNP cluster and ST
summaries_st_pds <- external_meta_kleborate %>%
  group_by(SNP_cluster, ST) %>%
  summarise(count = n())

# Export metadata subsets for selected SNP clusters
clusters <- c("PDS000060581.72", "PDS000006578.129", "PDS000091501.207")
for (cl in clusters) {
  df <- external_meta_kleborate %>% filter(SNP_cluster == cl)
  write_tsv(df, paste0("./", cl, "_meta.tsv"))
}

# =======================================
# Load PDG and PDS Accession Mappings
# =======================================
big_meta <- read_tsv("../../Kp_ML/PDG000000012.2086.metadata.tsv") %>%
  distinct(biosample_acc, .keep_all = TRUE)

pds_acc <- read_tsv("../../Kp_ML/PDS.tsv") %>%
  distinct(BioSample, .keep_all = TRUE)

# Map between BioSample, Isolate, and Run
big_meta <- big_meta %>%
  left_join(pds_acc, by = c("biosample_acc" = "BioSample")) %>%
  select(biosample_acc, Isolate, Run) %>%
  distinct(biosample_acc, .keep_all = TRUE) %>%
  distinct(Isolate, .keep_all = TRUE)

# Create a named vector mapping isolate → biosample
name_map <- setNames(big_meta$biosample_acc, big_meta$Isolate)

# Load GenBank accession list
all_genbank <- read_csv("../../Kp_ML/all_Genbank.csv")

# =====================================================
# Compute External–Internal SNP Distance (Example 1)
# =====================================================
# --- Function to calculate and format distance matrices ---
compute_distance_matrix <- function(tree_path, name_map) {
  tree <- read.tree(tree_path)
  dist_mat <- round(cophenetic.phylo(tree))
  
  # Clean row/column names
  rownames(dist_mat) <- gsub("^\\s*|'+$", "", sub(".*,(.*)", "\\1", rownames(dist_mat)))
  colnames(dist_mat) <- gsub("^\\s*|'+$", "", sub(".*,(.*)", "\\1", colnames(dist_mat)))
  
  # Replace isolate names with BioSample accession
  rownames(dist_mat) <- name_map[rownames(dist_mat)]
  colnames(dist_mat) <- name_map[colnames(dist_mat)]
  
  return(dist_mat)
}

# --- Run for specific cluster ---
PDS000060581.72_dist <- compute_distance_matrix("./PDS000060581.72.newick", name_map)

# Separate internal and external isolates
internal_list <- intersect(all_genbank$Genbank, rownames(PDS000060581.72_dist))
external_list <- setdiff(rownames(PDS000060581.72_dist), internal_list)

# Compute external–internal distance pairs
EI_matrix <- PDS000060581.72_dist[external_list, internal_list]
EI_matrix_long <- EI_matrix %>%
  as.data.frame() %>%
  rownames_to_column("external") %>%
  pivot_longer(cols = -external, names_to = "internal", values_to = "distance")

# Get minimum SNP distance per external isolate
EI_matrix_long_min <- EI_matrix_long %>%
  group_by(external) %>%
  summarise(min_dist = min(distance, na.rm = TRUE))

# =======================================
#  Map Countries to Global Regions
# =======================================
cleaned_metadata <- read_tsv("./all_external_meta_clean.tsv")

region_map <- c(
  "China" = "East Asia", "Japan" = "East Asia", "South Korea" = "East Asia",
  "India" = "South Asia", "Nepal" = "South Asia", "Pakistan" = "South Asia",
  "Bangladesh" = "South Asia", "Thailand" = "Southeast Asia",
  "Cambodia" = "Southeast Asia", "Singapore" = "Southeast Asia",
  "Saudi Arabia" = "Middle East", "Oman" = "Middle East", "Qatar" = "Middle East",
  "Kuwait" = "Middle East", "Lebanon" = "Middle East", "Israel" = "Middle East",
  "Armenia" = "Middle East", "Ireland" = "Europe", "Italy" = "Europe",
  "Norway" = "Europe", "Turkey" = "Europe", "Spain" = "Europe", "France" = "Europe",
  "Finland" = "Europe", "Germany" = "Europe", "Netherlands" = "Europe",
  "Switzerland" = "Europe", "Slovakia" = "Europe", "Slovenia" = "Europe",
  "United Kingdom" = "Europe", "Denmark" = "Europe", "Sweden" = "Europe",
  "Belgium" = "Europe", "Portugal" = "Europe", "Bulgaria" = "Europe",
  "Greece" = "Europe", "Croatia" = "Europe", "Serbia" = "Europe",
  "Luxembourg" = "Europe", "Romania" = "Europe", "Hungary" = "Europe",
  "Czech Republic" = "Europe", "USA" = "North America", "Canada" = "North America",
  "Mexico" = "North America", "Guatemala" = "North America",
  "Guadeloupe" = "North America", "Paraguay" = "South America",
  "Colombia" = "South America", "Peru" = "South America", "Chile" = "South America",
  "Egypt" = "Africa", "Senegal" = "Africa", "South Africa" = "Africa",
  "Kenya" = "Africa", "Ghana" = "Africa", "Tunisia" = "Africa",
  "Australia" = "Oceania", "New Zealand" = "Oceania", "Fiji" = "Oceania",
  "not collected" = "Unknown"
)

cleaned_metadata$region <- region_map[cleaned_metadata$country]

# Join with distance data
EI_matrix_long_min_region <- EI_matrix_long_min %>%
  left_join(cleaned_metadata, by = c("external" = "BioSample")) %>%
  filter(!is.na(SNP_cluster))

# =======================================
# Cumulative SNP Distance by Region
# =======================================
EI_cumulative <- expand.grid(
  distance = 1:50,
  region = unique(EI_matrix_long_min_region$region)
) %>%
  rowwise() %>%
  mutate(freq = sum(EI_matrix_long_min_region$min_dist <= distance &
                      EI_matrix_long_min_region$region == region)) %>%
  ungroup()

EI_cumulative$region <- factor(
  EI_cumulative$region,
  levels = c("Africa", "Europe", "Oceania", "South Asia", "Southeast Asia",
             "Middle East", "North America", "South America")
)

# Define region colors
region_colors <- c(
  "Africa" = "#a6cee3", "Europe" = "#1f78b4", "Oceania" = "#6a3d9a",
  "East Asia" = "#fb9a99", "South Asia" = "#fdbf6f",
  "Southeast Asia" = "#ff7f00", "Middle East" = "#e31a1c",
  "North America" = "#b15928", "South America" = "#33a02c"
)

# Plot cumulative distribution of SNP distances by region
p1 <- ggplot(EI_cumulative, aes(x = distance, y = freq, fill = region)) +
  geom_col(position = "stack", width = 1) +
  labs(x = "SNP Distance", y = "Frequency", fill = "Region") +
  scale_x_continuous(breaks = seq(0, 50, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 150, 25), expand = c(0, 0)) +
  scale_fill_manual(values = region_colors) +
  annotate("rect", fill = "#767B91", alpha = 0.4,
           xmin = 0, xmax = 20, ymin = 0, ymax = Inf) +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 15, face = "bold", color = "black"),
    legend.text = element_text(size = 18, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 20)
  )

############################################################
# Combine Kleborate Output with External Metadata
############################################################

# Read Kleborate output
cleaned_kleborate <- read_tsv("./kleborate_output_need_filter.tsv")

# Merge with external metadata, exclude unknown or GenBank strains
external_metadata_kleborate <- cleaned_metadata %>%
  left_join(cleaned_kleborate, by = c("Run" = "strain")) %>%
  filter(region != "Unknown") %>%
  filter(!BioSample %in% all_genbank$Genbank)

# Save merged file
write_tsv(external_metadata_kleborate, "./external_metadata_kleborate.tsv")

############################################################
# Combine Kleborate Output with Internal Metadata
############################################################

# Merge with SNP cluster info
all_pds <- all_genbank %>%
  left_join(pds_acc, by = c("Genbank" = "BioSample")) %>%
  select(sample, Genbank, `SNP cluster`)

# Read and preprocess internal clinical metadata
internal_metadata <- read_tsv("./Kp_clincal_metadata_18May.tsv") %>%
  separate(Collection_Date, into = c("collect_year", "collect_month", "collect_day"), sep = "/", remove = FALSE)

# Join with PDS cluster and Kleborate report
internal_kleb <- read_tsv("./Kp_kleborate_report_23Jun.tsv")

internal_metadata_kleborate <- internal_metadata %>%
  left_join(all_pds, by = c("strain" = "sample")) %>%
  left_join(internal_kleb, by = c("strain" = "strain"))

# Save merged file
write_tsv(internal_metadata_kleborate, "./internal_metadata_kleborate.tsv")

############################################################
# Data Preparation for Plotting
############################################################

cleaned_metadata_kleborate_plot <- read_csv("../../Kp_ML/external_internal_metadata_kleborate.csv") %>%
  filter(!is.na(collect_year)) %>%
  filter(collect_year != "2007")

# Define region order for plotting
cleaned_metadata_kleborate_plot$region <- factor(
  cleaned_metadata_kleborate_plot$region,
  levels = c(
    "Africa", "Europe", "Oceania", "East Asia",
    "South Asia", "Southeast Asia", "Middle East",
    "North America", "South America"
  )
)

############################################################
# Boxplots of Resistance and Virulence Scores by Year
############################################################

# ---- Resistance over time ----
p2 <- ggplot(cleaned_metadata_kleborate_plot, aes(x = factor(collect_year), y = resistance_score)) +
  geom_boxplot() +
  labs(x = "Year", y = "Resistance Score") +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 20)
  )

# ---- Virulence over time ----
p3 <- ggplot(cleaned_metadata_kleborate_plot, aes(x = factor(collect_year), y = virulence_score)) +
  geom_boxplot() +
  labs(x = "Year", y = "Virulence Score") +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 20)
  )

############################################################
# Regional Distribution of Resistance and Virulence
############################################################

# ---- Resistance by region ----
p4 <- ggplot(cleaned_metadata_kleborate_plot, aes(x = factor(collect_year), y = resistance_score, fill = region)) +
  geom_boxplot() +
  labs(x = "Year", y = "Resistance Score") +
  scale_fill_manual(values = region_colors, name = "Region") +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 20)
  )

# ---- Virulence by region ----
p5 <- ggplot(cleaned_metadata_kleborate_plot, aes(x = factor(collect_year), y = virulence_score, fill = region)) +
  geom_boxplot() +
  labs(x = "Year", y = "Virulence Score") +
  scale_fill_manual(values = region_colors, name = "Region") +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 20)
  )

############################################################
# Dot Plots with Trend Lines
############################################################

# ---- Resistance scatter + regression ----
p6 <- ggplot(cleaned_metadata_kleborate_plot, aes(x = factor(collect_year), y = resistance_score, color = region)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Resistance Score") +
  scale_color_manual(values = region_colors, name = "Region") +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black")
  )

# ---- Correlation between year and resistance ----
cor.test(cleaned_metadata_kleborate_plot$collect_year, cleaned_metadata_kleborate_plot$resistance_score)
summary(lm(resistance_score ~ collect_year, data = cleaned_metadata_kleborate_plot))

# ---- Virulence scatter + regression ----
p7 <- ggplot(cleaned_metadata_kleborate_plot, aes(x = factor(collect_year), y = virulence_score, color = region)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Virulence Score") +
  scale_color_manual(values = region_colors, name = "Region") +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black")
  )

# ---- Correlation between year and virulence ----
cor.test(cleaned_metadata_kleborate_plot$collect_year, cleaned_metadata_kleborate_plot$virulence_score)
summary(lm(virulence_score ~ collect_year, data = cleaned_metadata_kleborate_plot))

############################################################
# AMR Gene Class Correlation with Year
############################################################

# Count number of acquired AMR classes per isolate
cleaned_metadata_kleborate_mutate <- cleaned_metadata_kleborate_plot %>%
  mutate(across(
    ends_with("_acquired"),
    .fns = list(num = ~ case_when(.x == "-" ~ 0L, TRUE ~ str_count(.x, ";") + 1L)),
    .names = "{.col}_num"
  ))

# Compute correlations between gene counts and year
acquired_num_cols <- cleaned_metadata_kleborate_mutate %>%
  select(ends_with("_acquired_num")) %>%
  colnames()

cor_result_amr_class <- map_dfr(acquired_num_cols, function(col) {
  test <- cor.test(cleaned_metadata_kleborate_mutate[[col]], cleaned_metadata_kleborate_mutate$collect_year, method = "pearson")
  tibble(variable = col, correlation = test$estimate, p_value = test$p.value)
}) %>%
  filter(!is.na(p_value) & p_value < 0.05)

############################################################
# Temporal Distribution of Carbapenemase Genes
############################################################

carb_gene_df <- cleaned_metadata_kleborate_mutate %>%
  select(sample, Bla_Carb_acquired, collect_year) %>%
  filter(Bla_Carb_acquired != "-") %>%
  separate_rows(Bla_Carb_acquired, sep = ";") %>%
  mutate(
    Bla_Carb_acquired = Bla_Carb_acquired %>%
      str_remove_all("[*?^]") %>%     # remove noisy symbols
      str_remove("^bla") %>%          # remove 'bla' prefix
      str_trim()
  ) %>%
  group_by(collect_year, Bla_Carb_acquired) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collect_year) %>%
  mutate(total = sum(count), proportion = count / total)

# Define key carbapenemase genes
key_genes <- c("KPC-2", "KPC-3", "NDM-1", "NDM-5", "OXA-48", "OXA-232", "VIM-1")

# Group less frequent genes as "Other"
carb_gene_df_grouped <- carb_gene_df %>%
  mutate(Gene_group = if_else(Bla_Carb_acquired %in% key_genes, Bla_Carb_acquired, "Other")) %>%
  group_by(collect_year, Gene_group) %>%
  summarise(proportion = sum(proportion), .groups = "drop") %>%
  mutate(Gene_group = factor(Gene_group, levels = c(key_genes, "Other")))

# ---- Plot carbapenemase gene proportions ----
p8 <- ggplot(carb_gene_df_grouped, aes(x = collect_year, y = proportion, fill = Gene_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "Year", y = "Proportion of Carbapenemase Genes", fill = "Gene") +
  scale_x_continuous(breaks = seq(2011, 2025, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "KPC-2" = "#1f78b4",
    "KPC-3" = "#33a02c",
    "NDM-1" = "#e31a1c",
    "NDM-5" = "#ff7f00",
    "OXA-48" = "#6a3d9a",
    "OXA-232" = "#b15928",
    "VIM-1" = "#a6cee3",
    "Other" = "grey70"
  )) +
  theme_bw() +
  theme(
    text = element_text(size = 15, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black")
  )

############################################################
# Focused Analysis: ST2096 and ST147 Clones
############################################################

# ---- ST2096 ----
st2096_meta <- cleaned_metadata_kleborate_plot %>% filter(ST == "ST2096")

p9 <- ggplot(st2096_meta, aes(x = factor(collect_year), y = resistance_score)) +
  geom_boxplot() +
  labs(x = "Year", y = "Resistance Score", title = "ST2096") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold", color = "black"))

p10 <- ggplot(st2096_meta, aes(x = factor(collect_year), y = virulence_score)) +
  geom_boxplot() +
  labs(x = "Year", y = "Virulence Score", title = "ST2096") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold", color = "black"))

# ---- ST147 ----
st147_meta <- cleaned_metadata_kleborate_plot %>% filter(ST == "ST147")

p11 <- ggplot(st147_meta, aes(x = factor(collect_year), y = resistance_score)) +
  geom_boxplot() +
  labs(x = "Year", y = "Resistance Score", title = "ST147") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold", color = "black"))

p12 <- ggplot(st147_meta, aes(x = factor(collect_year), y = virulence_score)) +
  geom_boxplot() +
  labs(x = "Year", y = "Virulence Score", title = "ST147") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold", color = "black"))

# ---- ST2096 resistance regression ----
p13 <- ggplot(st2096_meta, aes(x = factor(collect_year), y = resistance_score, color = region)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Resistance Score", title = "ST2096") +
  scale_color_manual(values = region_colors, name = "Region") +
  theme_bw() +
  theme(text = element_text(size = 15, face = "bold", color = "black"))

# ---- ST2096 correlation ----
cor.test(st2096_meta$collect_year, st2096_meta$resistance_score)
summary(lm(resistance_score ~ collect_year, data = st2096_meta))


# ST2096: Virulence score over collection year
p14 <- ggplot(st2096_meta, aes(x = factor(collect_year), y = virulence_score)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(color = region)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Virulence Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    plot.margin  = margin(t = 10, r = 0, b = 10, l = 20)
  ) +
  scale_color_manual(values = region_colors, name = "Region")

# Correlation and linear regression for virulence vs year
cor.test(st2096_meta$collect_year, st2096_meta$virulence_score)
model <- lm(virulence_score ~ collect_year, data = st2096_meta)
summary(model)


# ST147: Resistance score over collection year
p15 <- ggplot(st147_meta, aes(x = factor(collect_year), y = resistance_score)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(color = region)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Resistance Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    plot.margin  = margin(t = 10, r = 0, b = 10, l = 20)
  ) +
  scale_color_manual(values = region_colors, name = "Region")

cor.test(st147_meta$collect_year, st147_meta$resistance_score)
model <- lm(resistance_score ~ collect_year, data = st147_meta)
summary(model)


# ST147: Virulence score over collection year
p16 <- ggplot(st147_meta, aes(x = factor(collect_year), y = virulence_score)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(color = region)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Virulence Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    plot.margin  = margin(t = 10, r = 0, b = 10, l = 20)
  ) +
  scale_color_manual(values = region_colors, name = "Region")

cor.test(st147_meta$collect_year, st147_meta$virulence_score)
model <- lm(virulence_score ~ collect_year, data = st147_meta)
summary(model)



### === AMR GENE CORRELATION WITH YEAR === ###

# Count number of acquired AMR genes for each class
cleaned_metadata_kleborate_mutate <- st2096_meta %>%
  mutate(across(
    ends_with("_acquired"),
    .fns = list(num = ~ case_when(.x == "-" ~ 0L, TRUE ~ str_count(.x, ";") + 1L)),
    .names = "{.col}_num"
  ))

# Select numeric AMR columns
acquired_num_cols <- cleaned_metadata_kleborate_mutate %>%
  select(ends_with("_acquired_num")) %>%
  colnames()

# Compute Pearson correlations between AMR gene count and collection year
cor_result_amr_class <- map_dfr(acquired_num_cols, function(col) {
  test <- cor.test(cleaned_metadata_kleborate_mutate[[col]],
                   cleaned_metadata_kleborate_mutate$collect_year,
                   method = "pearson")
  tibble(
    variable    = col,
    correlation = test$estimate,
    p_value     = test$p.value
  )
}) %>%
  filter(!is.na(p_value), p_value < 0.05)



### === TEMPORAL DISTRIBUTION OF SPECIFIC GENES === ###

# ESBL gene dynamics in ST2096
st2096_esbl_gene_df <- cleaned_metadata_kleborate_mutate %>%
  select(sample, Bla_ESBL_acquired, collect_year) %>%
  filter(Bla_ESBL_acquired != "-") %>%
  separate_rows(Bla_ESBL_acquired, sep = ";") %>%
  mutate(Bla_ESBL_acquired = Bla_ESBL_acquired %>%
           str_remove_all("[*?^]") %>%
           str_remove("^bla") %>%
           str_trim()) %>%
  group_by(collect_year, Bla_ESBL_acquired) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collect_year) %>%
  mutate(total = sum(count),
         proportion = count / total)

# Tetracycline gene dynamics in ST2096
st2096_tet_gene_df <- cleaned_metadata_kleborate_mutate %>%
  select(sample, Tet_acquired, collect_year) %>%
  filter(Tet_acquired != "-") %>%
  separate_rows(Tet_acquired, sep = ";") %>%
  mutate(Tet_acquired = Tet_acquired %>%
           str_remove_all("[*?^]") %>%
           str_remove("^bla") %>%
           str_trim()) %>%
  group_by(collect_year, Tet_acquired) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collect_year) %>%
  mutate(total = sum(count),
         proportion = count / total)

# Proportion plot for carbapenemase genes
p8 <- ggplot(st2096_esbl_gene_df, aes(x = collect_year, y = proportion, fill = Gene_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "Year", y = "Proportion of Carb Genes", fill = "Gene") +
  theme_bw() +
  scale_x_continuous(breaks = seq(2011, 2025, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),
    text = element_text(size = 15, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 20)
  ) +
  scale_fill_manual(values = c(
    "KPC-2"   = "#1f78b4",
    "KPC-3"   = "#33a02c",
    "NDM-1"   = "#e31a1c",
    "NDM-5"   = "#ff7f00",
    "OXA-48"  = "#6a3d9a",
    "OXA-232" = "#b15928",
    "VIM-1"   = "#a6cee3",
    "Other"   = "grey70"
  ))

# Order gene groups for consistent legend
carb_gene_df_grouped <- carb_gene_df_grouped %>%
  mutate(Gene_group = factor(Gene_group,
                             levels = c("KPC-2", "KPC-3", "NDM-1", "NDM-5",
                                        "OXA-48", "OXA-232", "VIM-1", "Other")))



# Load required metadata files
city_info <- read_tsv("./internal_metadata_with_city.tsv") %>%
  select(strain, City_Name)

snp_cluster_info <- read_csv("./external_internal_genbank.csv") %>%
  left_join(pds_acc, by = c("Genbank" = "BioSample")) %>%
  select(sample, Genbank, type, `SNP cluster`)


# --------------------------------------------------------------
# Clean and standardize date fields
# --------------------------------------------------------------
st2096_meta$collect_day <- as.numeric(st2096_meta$collect_day)

# Replace missing month/day with default values for date construction
st2096_meta <- st2096_meta %>%
  mutate(
    collect_month = if_else(is.na(collect_month), 6L, collect_month),  # Default month = June
    collect_day = case_when(
      is.na(collect_day) & is.na(collect_month) ~ 30L,  # If both missing, assume 30th June
      is.na(collect_day) ~ 15L,                        # If only day missing, assume 15th
      TRUE ~ collect_day
    )
  )

# Collect column names from Klebsiella dataset (excluding 'strain')
kleb_col <- setdiff(colnames(internal_kleb), "strain")

# Create proper collection date column
st2096_meta <- st2096_meta %>%
  mutate(collection_date = make_date(
    year = collect_year,
    month = collect_month,
    day = collect_day
  ))


# --------------------------------------------------------------
# Split metadata into internal and external datasets
# --------------------------------------------------------------
st2096_meta_external <- st2096_meta %>%
  filter(type == "external") %>%
  select(all_of(c(
    "sample", "type", "collect_year", "collect_month", "collect_day",
    "collection_date", "country", "region", kleb_col
  )))

st2096_meta_internal <- st2096_meta %>%
  filter(type == "internal") %>%
  left_join(city_info, by = c("sample" = "strain")) %>%
  select(all_of(c(
    "sample", "type", "collect_year", "collect_month", "collect_day",
    "collection_date", "country", "region", "City_Name", kleb_col
  )))

# Save processed metadata
write_tsv(st2096_meta_external, "./st2096_meta_external.tsv")
write_tsv(st2096_meta_internal, "./st2096_meta_internal.tsv")


# --------------------------------------------------------------
# Add SNP cluster grouping
# --------------------------------------------------------------
st2096_meta <- st2096_meta %>%
  left_join(snp_cluster_info, by = c("sample" = "sample", "type" = "type")) %>%
  mutate(
    grouped_snp_cluster = if_else(
      `SNP cluster` %in% c("PDS000060581.72", "PDS000166904.9", "PDS000166905.6"),
      `SNP cluster`,
      "Others"
    )
  )


# --------------------------------------------------------------
# Build Neighbor-Joining tree from SNP alignment
# --------------------------------------------------------------
st2096_dna <- read.dna("./st2096.snp_sites.aln", format = "fasta")
st2096_distdna <- dist.dna(st2096_dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
st2096_njtree <- nj(st2096_distdna) %>% midpoint.root()

plot(st2096_njtree)


# --------------------------------------------------------------
# Annotate tree with metadata
# --------------------------------------------------------------
st2096_njtree <- left_join(st2096_njtree, st2096_meta, by = c("label" = "sample"))

# Base tree plot colored by SNP cluster
p9 <- ggtree(st2096_njtree) +
  geom_tippoint(aes(color = factor(grouped_snp_cluster)), size = 1.5, alpha = 1) +
  theme_tree2(legend.position = 'right') +
  scale_colour_discrete("SNP cluster")
p9


# --------------------------------------------------------------
# Layer 1: Region heatmap
# --------------------------------------------------------------
tree_continent <- as.data.frame(st2096_meta$region)
colnames(tree_continent) <- "region"
rownames(tree_continent) <- st2096_meta$sample

p10 <- gheatmap(
  p9, tree_continent,
  offset = 50, width = 0.1, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(values = region_colors, name = "Region")
p10


# --------------------------------------------------------------
# Layer 2: Country heatmap
# --------------------------------------------------------------
tree_country <- as.data.frame(st2096_meta$country)
rownames(tree_country) <- st2096_meta$sample
length(unique(tree_country$country))

p11 <- p10 + new_scale_fill()

p12 <- gheatmap(
  p11, tree_country,
  offset = 150, width = 0.1, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(
    values = c(
      "#DD6E42", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF",
      "#800000", "#008000", "#000080", "#808000", "#800080", "#008080",
      "#B1740F", "#FFD07B", "#FFA500", "#FFFF80", "#80FF00", "#80FFFF",
      "#FF80C0", "#FF0080", "#8000FF"
    ),
    name = "Country"
  )
p12


# --------------------------------------------------------------
# Layer 3: Collection year heatmap
# --------------------------------------------------------------
tree_year <- as.data.frame(as.character(st2096_meta$collect_year))
colnames(tree_year) <- "collect_year"
rownames(tree_year) <- st2096_meta$sample

p13 <- p12 + new_scale_fill()

p14 <- gheatmap(
  p13, tree_year,
  offset = 250, width = 0.1, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(
    values = c(
      "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6",
      "#2171b5", "#08519c", "#08306b", "#041f4a", "#021024"
    ),
    name = "Year"
  )
p14


# --------------------------------------------------------------
# Layer 4: Carbapenemase gene presence/absence heatmap
# --------------------------------------------------------------
st2096_meta_pivot <- st2096_meta %>%
  mutate(Bla_Carb_acquired = na_if(Bla_Carb_acquired, "-")) %>%
  separate_rows(Bla_Carb_acquired, sep = ";") %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Bla_Carb_acquired,
    values_from = present,
    values_fill = 0
  ) %>%
  select(-c("NA"))

# Extract carbapenemase genes of interest
tree_carb <- st2096_meta_pivot %>%
  select(c(`OXA-232`, `OXA-48`, `NDM-1`)) %>%
  mutate(across(everything(), as.character))

rownames(tree_carb) <- st2096_meta_pivot$sample

p15 <- p14 + new_scale_fill()

p16 <- gheatmap(
  p15, tree_carb,
  offset = 350, width = 0.3, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(
    breaks = c("0", "1"),
    labels = c("Absence", "Presence"),
    values = c('#C0C0C0', 'red'),
    name = "Carb gene"
  )
p16

# --------------------------------------------------------------
# 1. Clean and standardize collection date information
# --------------------------------------------------------------
st147_meta$collect_day <- as.numeric(st147_meta$collect_day)

# Replace missing month/day values with sensible defaults
st147_meta <- st147_meta %>%
  mutate(
    collect_month = if_else(is.na(collect_month), 6L, collect_month),   # Default to June if missing
    collect_day = case_when(
      is.na(collect_day) & is.na(collect_month) ~ 30L,  # If both NA → assume June 30
      is.na(collect_day) ~ 15L,                        # If only day NA → assume 15th
      TRUE ~ collect_day
    )
  )

# Construct a valid date column
st147_meta <- st147_meta %>%
  mutate(collection_date = make_date(
    year = collect_year,
    month = collect_month,
    day = collect_day
  ))


# --------------------------------------------------------------
# 2. Split internal and external metadata tables
# --------------------------------------------------------------
st147_meta_external <- st147_meta %>%
  filter(type == "external") %>%
  select(all_of(c(
    "sample", "type", "collect_year", "collect_month", "collect_day",
    "collection_date", "country", "region", kleb_col
  )))

st147_meta_internal <- st147_meta %>%
  filter(type == "internal") %>%
  left_join(city_info, by = c("sample" = "strain")) %>%
  select(all_of(c(
    "sample", "type", "collect_year", "collect_month", "collect_day",
    "collection_date", "country", "region", "City_Name", kleb_col
  )))

# Save processed metadata
write_tsv(st147_meta_external, "./st147_meta_external.tsv")
write_tsv(st147_meta_internal, "./st147_meta_internal.tsv")


# --------------------------------------------------------------
# 3. Join SNP cluster data and group clusters of interest
# --------------------------------------------------------------
st147_meta <- st147_meta %>%
  left_join(snp_cluster_info, by = c("sample" = "sample", "type" = "type")) %>%
  mutate(
    grouped_snp_cluster = if_else(
      `SNP cluster` %in% c("PDS000091501.207", "PDS000006578.129"),
      `SNP cluster`,
      "Others"
    )
  )


# --------------------------------------------------------------
# 4. Reshape carbapenemase gene presence/absence data
# --------------------------------------------------------------
st147_meta_pivot <- st147_meta %>%
  mutate(Bla_Carb_acquired = na_if(Bla_Carb_acquired, "-")) %>%
  separate_rows(Bla_Carb_acquired, sep = ";") %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Bla_Carb_acquired,
    values_from = present,
    values_fill = 0
  ) %>%
  select(-c("NA"))


# --------------------------------------------------------------
# 5. Construct Neighbor-Joining phylogenetic tree
# --------------------------------------------------------------
st147_dna <- read.dna("./st147.snp_sites.aln", format = "fasta")
st147_distdna <- dist.dna(st147_dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
st147_njtree <- nj(st147_distdna) %>% midpoint.root()

plot(st147_njtree)  # Quick visualization


# --------------------------------------------------------------
# 6. Annotate tree with metadata
# --------------------------------------------------------------
st147_njtree <- left_join(st147_njtree, st147_meta_pivot, by = c("label" = "sample"))

# Base tree visualization (colored by SNP cluster)
p17 <- ggtree(st147_njtree) +
  geom_tippoint(aes(color = factor(grouped_snp_cluster)), size = 1.5, alpha = 1) +
  theme_tree2(legend.position = 'right') +
  scale_colour_discrete("SNP cluster")
p17


# --------------------------------------------------------------
# 7. Add region-level annotation
# --------------------------------------------------------------
tree_continent <- as.data.frame(st147_meta$region)
colnames(tree_continent) <- "region"
rownames(tree_continent) <- st147_meta$sample

p18 <- gheatmap(
  p17, tree_continent,
  offset = 100, width = 0.1, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(values = region_colors, name = "Region")
p18


# --------------------------------------------------------------
# 8. Add country-level annotation
# --------------------------------------------------------------
tree_country <- as.data.frame(st147_meta$country)
rownames(tree_country) <- st147_meta$sample
length(unique(tree_country$country))  # Check number of countries

p19 <- p18 + new_scale_fill()

p20 <- gheatmap(
  p19, tree_country,
  offset = 800, width = 0.1, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(
    values = c(
      "#DD6E42", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF",
      "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#8AA8A1",
      "#B1740F", "#FFD07B", "#FFA500", "#FFFF80", "#80FF00", "#80FFFF",
      "#FF80C0", "#FF0080", "#8000FF", "#D1B490"
    ),
    name = "Country"
  )
p20


# --------------------------------------------------------------
# 9. Add collection year annotation
# --------------------------------------------------------------
tree_year <- as.data.frame(as.character(st147_meta$collect_year))
colnames(tree_year) <- "collect_year"
rownames(tree_year) <- st147_meta$sample
length(unique(tree_year$collect_year))

p21 <- p20 + new_scale_fill()

p22 <- gheatmap(
  p21, tree_year,
  offset = 1500, width = 0.1, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(
    values = c(
      "#f1f7fc", "#deebf7", "#c6dbef", "#aad2e6", "#9ecae1", "#85bcdb",
      "#6baed6", "#4292c6", "#2171b5", "#1762a9", "#08519c", "#08306b",
      "#041f4a", "#021024"
    ),
    name = "Year"
  )
p22


# --------------------------------------------------------------
# 10. Add carbapenemase gene heatmap (presence/absence)
# --------------------------------------------------------------
tree_carb <- st147_meta_pivot %>%
  select(c(`OXA-232`, `OXA-48`, `OXA-181`, `NDM-1`, `NDM-5`, `KPC-3`)) %>%
  mutate(across(everything(), as.character))

rownames(tree_carb) <- st147_meta_pivot$sample

p23 <- p22 + new_scale_fill()

p24 <- gheatmap(
  p23, tree_carb,
  offset = 2200, width = 0.3, font.size = 3,
  colnames = FALSE, hjust = 0, color = NA
) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(
    breaks = c("0", "1"),
    labels = c("Absence", "Presence"),
    values = c('#C0C0C0', 'red'),
    name = "Carb gene"
  )
p24


# --------------------------------------------------------------
# 11. Prepare metadata subsets for phylodynamic analysis
# --------------------------------------------------------------
# --- ST2096-related clusters (for cross-comparison) ---
pds72 <- st2096_meta %>%
  filter(`SNP cluster` == "PDS000060581.72") %>%
  select(sample, type, collection_date, country)

pds72_internal <- pds72 %>%
  filter(type == "internal") %>%
  left_join(city_info, by = c("sample" = "strain")) %>%
  select(sample, type, collection_date, City_Name)

write_tsv(pds72, "./PDS72_world.tsv")
write_tsv(pds72_internal, "./PDS72_saudi.tsv")

pds9 <- st2096_meta %>%
  filter(`SNP cluster` == "PDS000166904.9")

pds9_internal <- pds9 %>%
  filter(type == "internal") %>%
  left_join(city_info, by = c("sample" = "strain")) %>%
  select(sample, type, collection_date, City_Name)

write_tsv(pds9_internal, "./PDS9_saudi.tsv")


# --- ST147-related clusters for global analysis ---
pds207 <- st147_meta %>%
  filter(`SNP cluster` == "PDS000091501.207") %>%
  select(sample, type, collection_date, country)
write_tsv(pds207, "./PDS207_world.tsv")

pds129 <- st147_meta %>%
  filter(`SNP cluster` == "PDS000006578.129") %>%
  select(sample, type, collection_date, country)
write_tsv(pds129, "./PDS129_world.tsv")






