# Load libraries -----------------------------------------------------------
library(tidyverse)
library(ape)
library(treedataverse)
library(ggnewscale)
# library(rhierbaps)   # optional for hierarchical clustering
library(phytools)
library(fastbaps)
library(skygrowth)

# Load metadata ------------------------------------------------------------
kleborate <- read_tsv("./Kp_kleborate_report_18May.tsv")
meta <- read_tsv("./Kp_clincal_metadata_18May.tsv")

# Merge Kleborate report with clinical metadata
Kp_kleborate_meta <- kleborate %>%
  left_join(meta, by = c("strain" = "strain"))

# -------------------------------------------------------------------------
# ST2096 isolate information
# -------------------------------------------------------------------------
ST2096_info <- Kp_kleborate_meta %>%
  filter(ST == "ST2096", !is.na(Collection_Date)) %>%
  select(strain, Collection_Date, Hospital_Name, Regional_Name)

write_tsv(ST2096_info, "./ST2096_info.tsv")

# Clean date format for ST2096
ST2096_info <- read_tsv("./ST2096_info.tsv") %>%
  mutate(Collection_Date = format(as.Date(Collection_Date, format = "%Y/%m/%d"), "%Y-%m-%d")) %>%
  select(strain, Collection_Date, Regional_Name, City)

write_tsv(ST2096_info, "./ST2096_info.tsv")

# -------------------------------------------------------------------------
# ST147 isolate information
# -------------------------------------------------------------------------
ST147_info <- Kp_kleborate_meta %>%
  filter(ST == "ST147", !is.na(Collection_Date)) %>%
  select(strain, Collection_Date, Hospital_Name, Regional_Name)

write_tsv(ST147_info, "./ST147_info.tsv")

# Clean date format for ST147
ST147_info <- read_tsv("./ST147_info.tsv") %>%
  mutate(Collection_Date = format(as.Date(Collection_Date, format = "%Y/%m/%d"), "%Y-%m-%d")) %>%
  select(strain, Collection_Date, Regional_Name, City)

write_tsv(ST147_info, "./ST147_info.tsv")

# -------------------------------------------------------------------------
# Integrate external SNP cluster and GenBank data
# -------------------------------------------------------------------------
MM_100_ena <- read_csv("../../../../Documents/Kp_ML/MM100_ena.csv")
big_meta <- read_tsv("../../../../Documents/Kp_ML/PDG000000012.2086.metadata.tsv")
pds_acc <- read_tsv("../../../../Documents/Kp_ML/PDS.tsv")
all_genbank <- read_csv("../../../../Documents/Kp_ML/all_Genbank.csv")

# Merge ENA and GenBank with internal metadata
MM_100_ena_genbank <- MM_100_ena %>%
  left_join(big_meta, by = c("ENA_acc" = "Run"))
write_tsv(MM_100_ena_genbank, "./MM_100_ena_genbank.tsv")

# Merge GenBank info for ST2096 and ST147
ST2096_genbank <- ST2096_info %>%
  left_join(all_genbank, by = c("strain" = "sample"))
ST147_genbank <- ST147_info %>%
  left_join(all_genbank, by = c("strain" = "sample"))

# Link SNP cluster information
all_pds <- all_genbank %>%
  left_join(pds_acc, by = c("Genbank" = "BioSample")) %>%
  select(sample, Genbank, `SNP cluster`)

# Identify available SNP clusters
all_pds_value <- na.omit(unique(all_pds$`SNP cluster`))

# Identify external (non-local) BioSamples belonging to known clusters
all_external <- pds_acc %>%
  filter(`SNP cluster` %in% all_pds_value) %>%
  filter(!BioSample %in% all_genbank) %>%
  distinct(BioSample, .keep_all = TRUE)

# Annotate external metadata with global context
all_external_meta <- all_external %>%
  left_join(big_meta, by = c("BioSample" = "biosample_acc")) %>%
  arrange(BioSample, desc(!is.na(Run))) %>%
  distinct(BioSample, .keep_all = TRUE)

# -------------------------------------------------------------------------
# Filter and clean external metadata
# -------------------------------------------------------------------------
all_external_meta_filter <- all_external_meta %>%
  filter(!Run %in% c("NULL"),
         !geo_loc_name %in% c("NULL"),
         !collection_date %in% c("NULL", "missing", "2922-08-23", "2014-11/2016-01", "2019-10/2020-09")) %>%
  separate(geo_loc_name, sep = ":", into = c("country", "city"), remove = FALSE) %>%
  separate(collection_date, sep = "-", into = c("year", "month", "day"), remove = FALSE)

# Keep SNP clusters with at least 10 samples
snp_cluster <- all_external_meta_filter %>%
  group_by(`SNP cluster`) %>%
  summarise(unique_sample_count = n()) %>%
  filter(unique_sample_count >= 10)

# Filter by cluster list
all_external_meta_filter <- all_external_meta_filter %>%
  filter(`SNP cluster` %in% snp_cluster$`SNP cluster`)

# Select essential fields for output
all_external_meta_clean <- all_external_meta_filter %>%
  select(`SNP cluster`, BioSample, Run, bioproject_acc, collection_date,
         year, month, day, geo_loc_name, country, city, host)

write_tsv(all_external_meta_clean, "./all_external_meta_clean.tsv")

# -------------------------------------------------------------------------
# Filter internal STs for phylodynamic analysis
# -------------------------------------------------------------------------
internal_meta <- read_tsv("./internal_metadata_kleborate.tsv") %>%
  mutate(ST_clean = sub("-.*", "", ST)) %>%
  filter(!is.na(ST_clean))

# Summarize ST counts
st_summarise <- internal_meta %>%
  group_by(ST_clean) %>%
  summarise(count = n())

write_tsv(st_summarise, "./Kp_st_count.tsv")

# Keep STs with ≥15 isolates
st_summarise <- st_summarise %>%
  filter(count >= 15)

# Select representative genomes (highest quality)
st_beast_meta <- internal_meta %>%
  filter(ST_clean %in% st_summarise$ST_clean) %>%
  group_by(ST_clean) %>%
  arrange(contig_count, desc(N50)) %>%
  mutate(best_assembly = first(strain)) %>%
  ungroup()

st_beast_mapping <- st_beast_meta %>%
  select(strain, ST_clean, best_assembly)

write_tsv(st_beast_mapping, "./st_beast_mapping.tsv")

# -------------------------------------------------------------------------
# SNP and K-mer distance comparison
# -------------------------------------------------------------------------
kmer_dist <- read_csv("./new_phylodynamic/summed_kmer_distance.csv")
kmer_dist_mat <- as.matrix(kmer_dist)
rownames(kmer_dist_mat) <- colnames(kmer_dist_mat)

# Convert distance matrix to long format
kmer_long <- as.data.frame(as.table(as.matrix(kmer_dist_mat))) %>%
  rename(sample1 = Var1, sample2 = Var2, kmer_distance = Freq)

# Initialize SNP distance dataframe
snp_long <- data.frame(
  sample1 = character(),
  sample2 = character(),
  snp_distance = numeric(),
  stringsAsFactors = FALSE
)

# Define ST groups for analysis
st_group <- c("ST101", "ST11", "ST14", "ST147", "ST17", "ST20", "ST2096",
              "ST219", "ST23", "ST231", "ST268", "ST29", "ST307", "ST35",
              "ST37", "ST39", "ST45", "ST48", "ST661")

# Compute SNP distance matrices for each ST
for (i in st_group) {
  file_path <- paste0("./new_phylodynamic/snp_aln/", i, ".filtered_polymorphic_sites.fasta")
  dna_baps <- read.dna(file_path, format = "fasta")
  distdna_baps <- dist.dna(dna_baps, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
  
  snp_mat <- as.data.frame(as.table(as.matrix(distdna_baps))) %>%
    rename(sample1 = Var1, sample2 = Var2, snp_distance = Freq)
  
  snp_long <- rbind(snp_long, snp_mat)
}

# Combine SNP and K-mer distances
snp_long <- snp_long %>%
  mutate(sample2 = as.character(sample2)) %>%
  left_join(st_beast_mapping, by = c("sample1" = "strain")) %>%
  select(sample1, sample2, ST_clean, snp_distance)

snp_kmer_dist <- snp_long %>%
  left_join(kmer_long, by = c("sample1", "sample2")) %>%
  filter(sample1 < sample2)

# -------------------------------------------------------------------------
# Plot correlation between SNP and K-mer distances
# -------------------------------------------------------------------------
p3 <- ggplot(snp_kmer_dist, aes(x = kmer_distance, y = snp_distance, fill = factor(ST_clean))) +
  geom_point(shape = 21) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 0.6, 0.05)) +
  scale_y_continuous(breaks = seq(0, 52000, 5000)) +
  labs(x = "K-mer distance", y = "SNP distance", fill = "Sequence Type (ST)") +
  theme(
    text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 15, color = "black", face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold"),
    plot.margin = margin(t = 10, r = 25, b = 10, l = 30)
  ) +
  scale_color_manual(values = c(
    "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4",
    "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080", "#E6BEFF",
    "#9A6324", "#FFFAC8", "#800000", "#Aaffc3", "#FFD8B1", "#000075",
    "#808000", "#FF4500", "#808080"
  ), na.translate = FALSE)

p3

# -------------------------------------------------------------------------
# Filter SNP-Kmer distance pairs for different thresholds
# -------------------------------------------------------------------------
snp_kmer_dist_20 <- snp_kmer_dist %>%
  filter(snp_distance <= 20, snp_distance != 0) %>%
  mutate(sample1 = as.character(sample1), sample2 = as.character(sample2)) %>%
  filter(sample1 < sample2)   # keep unique pairs

snp_kmer_dist_50 <- snp_kmer_dist %>%
  filter(snp_distance <= 50, snp_distance != 0) %>%
  mutate(sample1 = as.character(sample1), sample2 = as.character(sample2)) %>%
  filter(sample1 < sample2)

snp_kmer_dist_0.01 <- snp_kmer_dist %>%
  filter(kmer_distance <= 0.01, kmer_distance != 0) %>%
  mutate(sample1 = as.character(sample1), sample2 = as.character(sample2)) %>%
  filter(sample1 < sample2)

# -------------------------------------------------------------------------
# Plot histogram of K-mer distances for SNP <= 20
# -------------------------------------------------------------------------
p4 <- ggplot(snp_kmer_dist_20, aes(x = kmer_distance)) +
  geom_histogram(binwidth = 0.001, fill = "#107E7D", color = "black", alpha = 0.8) +
  theme_bw() +
  labs(x = "K-mer Distance", y = "Count") +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p4

# -------------------------------------------------------------------------
# Function: Remove outliers based on IQR
# -------------------------------------------------------------------------
remove_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  x[x >= lower & x <= upper]
}

# Remove outliers from K-mer distances
snp_kmer_dist_20_rm <- remove_outliers(snp_kmer_dist_20$kmer_distance)

# Calculate mean and median of cleaned data
snp_kmer_dist_20_mean <- mean(snp_kmer_dist_20_rm)
snp_kmer_dist_20_median <- median(snp_kmer_dist_20_rm)

# -------------------------------------------------------------------------
# Hierarchical clustering based on K-mer distances
# -------------------------------------------------------------------------
kmer_dist_obj <- as.dist(kmer_dist_mat)
hc <- hclust(kmer_dist_obj, method = "average")

# Cut tree at chosen height to define clusters
clusters <- cutree(hc, h = 0.00038052)

# Create data frame of sample and cluster assignment
cluster_df <- tibble(
  sample = names(clusters),
  cluster = as.integer(factor(clusters))
)

write_tsv(cluster_df, "./new_phylodynamic/final_cluster_26Aug.tsv")

# -------------------------------------------------------------------------
# Plot histogram of cluster sizes
# -------------------------------------------------------------------------
p5 <- ggplot(cluster_df, aes(x = cluster)) +
  geom_histogram(binwidth = 1, fill = "#107E7D", color = "#107E7D", alpha = 0.8) +
  theme_bw() +
  labs(x = "Cluster No.", y = "Count") +
  scale_x_continuous(breaks = seq(0, 1340, 50), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = seq(0, 50, 5), expand = expansion(mult = c(0, 0.05))) +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p5

# -------------------------------------------------------------------------
# Filter clusters for downstream analysis (≥5 samples)
# -------------------------------------------------------------------------
cluster_summary <- cluster_df %>%
  group_by(cluster) %>%
  summarise(count = n()) %>%
  filter(count >= 5)

cluster_beast <- cluster_df %>%
  filter(cluster %in% cluster_summary$cluster)

write_tsv(cluster_beast, "./new_phylodynamic/cluster_beast_26Aug.tsv")

# Merge cluster info with internal metadata
cluster_beast_meta <- cluster_beast %>%
  left_join(internal_meta, by = c("sample" = "strain")) %>%
  select(sample, cluster, collect_year, collect_month, collect_day, Genbank, ST_clean)

# Annotate with SNP cluster info
big_meta <- read_tsv("./new_phylodynamic/PDG000000012.2086.metadata.tsv")
pds_acc <- read_tsv("./new_phylodynamic/PDS.tsv")

cluster_beast_meta <- cluster_beast_meta %>%
  left_join(pds_acc, by = c("Genbank" = "BioSample")) %>%
  select(sample, cluster, collect_year, collect_month, collect_day, Genbank, `SNP cluster`, ST_clean)

# -------------------------------------------------------------------------
# ST summary with resistance and virulence annotation
# -------------------------------------------------------------------------
st_summary <- internal_meta %>%
  mutate(
    esbl_gene = if_else(Bla_ESBL_acquired == "-", 0, 1),
    carb_gene = if_else(Bla_Carb_acquired == "-", 0, 1),
    mdr = if_else(resistance_score >= 1, 1, 0),
    hv = if_else(virulence_score >= 3, 1, 0),
    mdr_hv = if_else(mdr == 1 & hv == 1, 1, 0)
  ) %>%
  group_by(ST_clean) %>%
  summarise(
    count = n(),
    resistance_max = max(resistance_score, na.rm = TRUE),
    resistance_mean = mean(resistance_score, na.rm = TRUE),
    virulence_max = max(virulence_score, na.rm = TRUE),
    virulence_mean = mean(virulence_score, na.rm = TRUE),
    esbl_ratio = mean(esbl_gene, na.rm = TRUE),
    carb_ratio = mean(carb_gene, na.rm = TRUE),
    mdr_ratio = mean(mdr, na.rm = TRUE),
    hv_ratio = mean(hv, na.rm = TRUE),
    dual_risk_ratio = mean(mdr_hv, na.rm = TRUE)
  ) %>%
  filter(count >= 10)

write_tsv(st_summary, "./st_summary.tsv")

# -------------------------------------------------------------------------
# Filter internal metadata for selected STs and best assemblies
# -------------------------------------------------------------------------
selected_st <- read_tsv("./select_ST.tsv")

internal_meta_select <- internal_meta %>%
  filter(ST_clean %in% selected_st$ST) %>%
  group_by(ST_clean) %>%
  arrange(contig_count, desc(N50)) %>%
  mutate(best_assembly = first(strain)) %>%
  ungroup() %>%
  select(strain, collect_year, collect_month, collect_day, Hospital_Name,
         Genbank, `SNP cluster`, ST_clean, best_assembly) %>%
  filter(!is.na(collect_year))

write_tsv(internal_meta_select, "./internal_meta_select.tsv")

# -------------------------------------------------------------------------
# External samples corresponding to selected SNP clusters
# -------------------------------------------------------------------------
external_samples_genbank <- pds_acc %>%
  filter(`SNP cluster` %in% unique(internal_meta_select$`SNP cluster`)) %>%
  filter(!is.na(`SNP cluster`)) %>%
  filter(!BioSample %in% internal_meta_select$Genbank) %>%
  select(BioSample, `SNP cluster`)

external_samples_meta <- big_meta %>%
  filter(biosample_acc %in% external_samples_genbank$BioSample) %>%
  filter(!is.na(Run), Run != "NULL") %>%
  filter(!collection_date %in% c("NULL", "missing", "2014-11/2016-01", "2015-01-01/2015-08-31", "2019-10/2020-09")) %>%
  filter(!geo_loc_name %in% c("NULL", "not collected")) %>%
  separate(collection_date, into = c("year", "month", "day"), sep = "-", remove = FALSE) %>%
  mutate(
    month = ifelse(is.na(month) & is.na(day), "06", month),
    day = ifelse(is.na(month) & is.na(day), "01", day),
    day = ifelse(!is.na(month) & is.na(day), "15", day)
  ) %>%
  separate(geo_loc_name, into = c("country", "city"), sep = ":", remove = FALSE) %>%
  select(Run, biosample_acc, year, month, day, epi_type, country) %>%
  mutate(collection_date = paste(year, month, day, sep = "-")) %>%
  left_join(external_samples_genbank, by = c("biosample_acc" = "BioSample"))

# Join SNP cluster info and best assembly table
best_assembly_table <- internal_meta_select %>%
  select(ST_clean, best_assembly) %>%
  distinct(ST_clean, .keep_all = TRUE)

external_samples_meta <- external_samples_meta %>%
  left_join(best_assembly_table, by = c("ST_clean" = "ST_clean"))

write_tsv(external_samples_meta, "./external_samples_meta.tsv")

# Remove placeholder ST147 entries
ST147_remove <- external_samples_meta %>%
  filter(ST_clean == "ST147", month == "06", day == "15")
write_tsv(ST147_remove, "./ST147_remove.tsv")

# -------------------------------------------------------------------------
# Add new isolates and merge with internal/external metadata
# -------------------------------------------------------------------------
new_pds <- read_tsv("/Users/huanj0f/Downloads/isolates.tsv")
new_big_meta <- read_tsv("./PDG000000012.2169.metadata.tsv")

internal_meta_add <- internal_meta_select %>%
  filter(ST_clean %in% c("ST17","ST23","ST268","ST35","ST37","ST39","ST45","ST29")) %>%
  left_join(new_pds, by = c("Genbank" = "BioSample")) %>%
  select(strain, collect_year, collect_month, collect_day, Hospital_Name,
         Genbank, `SNP cluster.x`, `SNP cluster.y`, ST_clean, best_assembly)

external_add <- new_pds %>%
  filter(`SNP cluster` %in% internal_meta_add$`SNP cluster.y`) %>%
  filter(!BioSample %in% internal_meta_add$Genbank) %>%
  filter(!is.na(`SNP cluster`))

# -------------------------------------------------------------------------
# Build neighbor-joining tree for ST35
# -------------------------------------------------------------------------
dna <- read.dna("./phylodynamic_st/ST35.snp_sites.aln", format = "fasta")
distdna <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
njtree <- nj(distdna)
njtree <- midpoint.root(njtree)
plot(njtree)

write.tree(njtree, file = "./phylodynamic_st/ST35.tree")
# -------------------------------------------------------------------------
# Load BEAST phylodynamic results
# -------------------------------------------------------------------------
beast_result_table <- read_tsv("./select_ST_beast.tsv")

# Filter out rows with missing age estimates
plot_table <- beast_result_table %>% filter(!is.na(age))

# -------------------------------------------------------------------------
# Common theme for all bar plots
# -------------------------------------------------------------------------
common_theme <- theme_bw() +
  theme(
    text = element_text(size = 25, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", color = "black"),
    axis.text.x = element_text(size = 25, angle = 45, vjust = 1, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(size = 25, color = "black", face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", color = "black"),
    axis.title.y = element_text(size = 25, face = "bold", color = "black"),
    legend.title = element_text(size = 25, face = "bold", color = "black"),
    legend.text = element_text(size = 25, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 35)
  )

# -------------------------------------------------------------------------
# MRCA Age Plot
# -------------------------------------------------------------------------
p1 <- ggplot(plot_table, aes(x = ST, y = age)) +
  # geom_bar(stat="identity", fill="#FF6666") + # optional fill
  geom_errorbar(aes(ymin = age_lower, ymax = age_upper),
                width = 0.4, colour = "black", alpha = 0.9, size = 1) +
  coord_cartesian(ylim = c(1900, 2025)) +
  scale_y_continuous(breaks = seq(1900, 2025, 20)) +
  labs(x = "ST", y = "MRCA age") +
  common_theme
p1

# -------------------------------------------------------------------------
# Clock Rate Plot
# -------------------------------------------------------------------------
p2 <- ggplot(plot_table, aes(x = ST, y = clock_rate)) +
  geom_bar(stat = "identity", fill = "#73AB84") +
  geom_errorbar(aes(ymin = clock_rate_lower, ymax = clock_rate_upper),
                width = 0.4, colour = "black", alpha = 0.9, size = 1,
                position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0.0005, 0.0065)) +
  scale_y_continuous(breaks = seq(0.0005, 0.0065, 0.001)) +
  labs(x = "ST", y = "Clock rate") +
  common_theme
p2

# -------------------------------------------------------------------------
# Migration Rate Plot
# -------------------------------------------------------------------------
p3 <- ggplot(plot_table, aes(x = ST, y = migration_rate)) +
  geom_bar(stat = "identity", fill = "#D4B483") +
  geom_errorbar(aes(ymin = migration_rate_lower, ymax = migration_rate_upper),
                width = 0.4, colour = "black", alpha = 0.9, size = 1,
                position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0.0005, 0.1)) +
  scale_y_continuous(breaks = seq(0.0005, 0.1, 0.01)) +
  labs(x = "ST", y = "Migration rate") +
  common_theme
p3

# -------------------------------------------------------------------------
# Skygrowth Analysis for ST661
# -------------------------------------------------------------------------
# Load the phylogenetic tree
tree <- read.nexus("./phylodynamic_st/ST661_added.tree")

# Fit Ne (effective population size) changes over time
fit <- skygrowth.map(tree,
                     # res = 24*13,  # optional: Ne changes every 2 weeks
                     tau0 = 0.8    # smoothing parameter
)

# Plot skygrowth result
p_sky <- plot(fit) +
  theme_bw() +
  labs(title = "ST661") +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, vjust = 0.5, face = "bold", color = "black"),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 25, face = "bold", color = "black"),
    axis.title.y = element_text(size = 25, face = "bold", color = "black"),
    strip.text.x = element_text(size = 25, face = "bold", color = "black")
  )
p_sky











