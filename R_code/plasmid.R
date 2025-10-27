#----------------------#
# Load Required Libraries
#----------------------#
library(factoextra)
library(fpc)
library(RColorBrewer)
library(tidyverse)
library(Rtsne)

#----------------------#
# Construct Similarity Matrix from kmer
#----------------------#
lines <- readLines("./output.txt")             # Input similarity lines
files <- readLines("./input.list")             # List of tags (sample names)
tags <- as.character(sapply(files, function(x) strsplit(x, " ")[[1]][1]))

sim_mat <- matrix(0, nrow = length(tags), ncol = length(tags))
name <- "chunk"

# Parse similarity data and populate sim_mat
for (i in seq_along(lines)) {
  print(i / length(lines))
  parts <- strsplit(lines[i], "\\|")[[1]]
  elements <- strsplit(trimws(parts[2]), " ")[[1]]
  column_names <- as.character(sapply(elements, function(x) strsplit(x, ":")[[1]][1]))
  
  idx <- match(column_names, tags)
  sim_mat[idx, idx] <- sim_mat[idx, idx] + 1
}

# Export similarity matrix
write.table(sim_mat, file = paste0(name, "output_matrix.txt"), sep = "\t", row.names = FALSE)

#----------------------#
# Compute Distance Matrix
#----------------------#
dist_mat <- matrix(0, nrow = length(tags), ncol = length(tags))
for (i in seq_along(tags)) {
  for (j in seq_along(tags)) {
    denom <- sim_mat[i, i] + sim_mat[j, j]
    dist_mat[i, j] <- ifelse(denom == 0, 1, 1 - 2 * sim_mat[i, j] / denom)
  }
}
rownames(dist_mat) <- tags
colnames(dist_mat) <- tags
write.table(dist_mat, "./dist_matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#----------------------#
# K-means Clustering and PCA Visualization
#----------------------#
pca_res <- prcomp(dist_mat)
k_max <- 15

# Evaluate total within-cluster sum of squares for elbow method
res <- numeric(k_max)
for (i in 2:k_max) {
  kmeans_res <- kmeans(dist_mat, centers = i)
  res[i] <- kmeans_res$tot.withinss
}

# Plot elbow method
fviz_nbclust(dist_mat, kmeans, method = "wss")

# Choose optimal cluster number (manually set to 8)
kmeans_res <- kmeans(dist_mat, centers = 8)

# Define colors for clusters
n_clusters <- length(unique(kmeans_res$cluster))
cluster_colors <- brewer.pal(n_clusters, "Set1")

# PCA visualization
fviz_pca_ind(
  pca_res,
  geom = "point",
  col.ind = as.factor(kmeans_res$cluster),
  palette = cluster_colors
)

# Save cluster assignments
tags_clusters_df <- data.frame(Tag = tags, Cluster = kmeans_res$cluster)
write_tsv(tags_clusters_df, "./plasmid_cluster.tsv")

#----------------------#
# t-SNE Visualization
#----------------------#
dist_obj <- as.dist(dist_mat)
set.seed(42)
tsne_results <- Rtsne(dist_obj, is_distance = TRUE, perplexity = 5)
plot(tsne_results$Y, col = "blue", pch = 19, main = "t-SNE Clustering")

#----------------------#
# Load and Visualize mge-cluster Results
#----------------------#
mge_data <- read.csv("./new-model_results.csv")
mge_data$Standard_Cluster <- as.factor(mge_data$Standard_Cluster)

ggplot(mge_data, aes(x = tsne1D, y = tsne2D, color = Standard_Cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "t-SNE visualization of plasmids from mge-cluster",
    x = "t-SNE Dimension 1", y = "t-SNE Dimension 2", color = "Cluster"
  )

#----------------------#
# Load Annotation Data (CARD, VFDB, Replicon)
#----------------------#
card <- read_tsv("./card_summary.tsv")
vfdb <- read_tsv("./vfdb_summary.tsv")
replicon <- read_tsv("./plasmid_summary.tsv")

# Extract plasmid names
replicon <- replicon %>%
  mutate(plasmid = str_extract(`#FILE`, "(?<=/)[^_]+_[^_]+")) %>%
  select(-`#FILE`) %>%
  select(plasmid, everything())

# Pivot replicon table to long format and group
replicon_pivot <- replicon %>%
  pivot_longer(cols = 3:47, names_to = "replicon", values_to = "value") %>%
  filter(value != ".") %>%
  group_by(plasmid) %>%
  summarise(plasmid_replicons = paste(replicon, collapse = ","), .groups = "drop")

# Repeat for CARD data (AMR genes)
card <- card %>%
  mutate(plasmid = str_extract(`#FILE`, "(?<=/)[^_]+_[^_]+")) %>%
  select(-`#FILE`) %>%
  select(plasmid, everything())

card_pivot <- card %>%
  pivot_longer(cols = 3:89, names_to = "AMG", values_to = "value") %>%
  filter(value != ".") %>%
  group_by(plasmid) %>%
  summarise(AMGs = paste(AMG, collapse = ","), .groups = "drop")

#----------------------#
# Replicon Family Grouping
#----------------------#
replicon_pivot <- replicon_pivot %>%
  mutate(plasmid_group = case_when(
    str_detect(plasmid_replicons, "Col440") ~ "Col440 Family",
    str_detect(plasmid_replicons, "ColKP3_1") ~ "ColKP3 Family",
    str_detect(plasmid_replicons, "ColRNAI_1") ~ "ColRNAI Family",
    str_detect(plasmid_replicons, "ColpVC_1") ~ "ColpVC Family",
    str_detect(plasmid_replicons, "Col\\(MG828\\)_1|Col156_1") ~ "Col(MG828)/Col156 Group",
    str_detect(plasmid_replicons, "FII\\(pBK30683\\)_1|IncFII_1_pKP91") ~ "FII-related (pBK30683/pKP91)",
    str_detect(plasmid_replicons, "IncFIB\\(K\\)_1_Kpn3") ~ "IncFIB(K)_1_Kpn3 Group",
    str_detect(plasmid_replicons, "IncFIA\\(HI1\\)_1_HI1") ~ "IncFIA(HI1) Family",
    str_detect(plasmid_replicons, "IncA/C2_1") ~ "IncA/C2 Group",
    str_detect(plasmid_replicons, "IncB/O/K/Z_1") ~ "IncB/O/K/Z",
    TRUE ~ "Unclassified"
  ))

#----------------------#
# AMR Gene Categorization (ESBL, Carbapenem, Colistin)
#----------------------#
ESBL_genes <- c("CTX-M", "SHV-12", "TEM-150", "OXA-1", "PER", "VEB", "GES", "TEM-141", "TEM-206", "CTX-M-9", "SHV-134", "SHV-187", "TEM-1")
Carbapenem_genes <- c("NDM", "OXA-48", "OXA-181", "OXA-232", "IMP", "VIM", "KPC", "OXA-9")
Colistin_genes <- c("MCR")

card_pivot2 <- card_pivot %>%
  rowwise() %>%
  mutate(
    ESBL = any(str_detect(AMGs, str_c(ESBL_genes, collapse = "|"))),
    Carb = any(str_detect(AMGs, str_c(Carbapenem_genes, collapse = "|"))),
    Colistin = any(str_detect(AMGs, str_c(Colistin_genes, collapse = "|"))),
    class_combo = str_c(
      "ESBL", if_else(ESBL, "+", "-"),
      ",Carb", if_else(Carb, "+", "-"),
      ",Colistin", if_else(Colistin, "+", "-")
    )
  ) %>%
  ungroup()

#----------------------#
# Merge Annotations with mge Data
#----------------------#
mge_data_replicon_iuc_amr <- mge_data %>%
  left_join(replicon_pivot, by = c("Sample_Name" = "plasmid")) %>%
  left_join(card_pivot2, by = c("Sample_Name" = "plasmid"))

write_tsv(mge_data_replicon_iuc_amr, "./mge_data_replicon_iuc_amr.tsv")

#----------------------#
# Visualization (t-SNE overlays)
#----------------------#
# Replicon groups
p1 <- ggplot(mge_data_replicon_iuc_amr, aes(x = tsne1D, y = tsne2D, color = plasmid_group)) +
  geom_point(size = 3) +
  coord_cartesian(xlim = c(-25, 25), ylim = c(-25, 25)) +
  theme_bw() +
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2", color = "Plasmid replicon")

p1

# AMR gene combinations
resistance_colors <- c(
  "ESBL-,Carb-,Colistin-" = "#ffe6e6",
  "ESBL+,Carb-,Colistin-" = "#ff9999",
  "ESBL-,Carb+,Colistin-" = "#ff3333",
  "ESBL+,Carb+,Colistin-" = "#ff3333",
  "ESBL-,Carb-,Colistin+" = "#cc0000",
  "ESBL+,Carb-,Colistin+" = "#cc0000"
)

p3 <- ggplot(mge_data_replicon_iuc_amr, aes(x = tsne1D, y = tsne2D, color = class_combo)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2", color = "AMR genes") +
  scale_color_manual(values = resistance_colors)

p3

###############################################################
# End of Script
###############################################################
















