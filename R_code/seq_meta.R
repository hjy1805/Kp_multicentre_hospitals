# Load required libraries ------------------------------------------------------
library(tidyverse)
library(readxl)
library(ggrepel)

# -----------------------------------------------------------------------------
# Combine sequencing metadata
# -----------------------------------------------------------------------------

# Read and process forward read size file
forward_size <- read_tsv("./forward_size_b4_5.tsv") %>%
  separate(file, into = c("text1", "tag", "text3"), sep = "/", remove = TRUE) %>%
  select(-c(text1, text3)) %>%
  rename(forward_size = size) %>%
  relocate(tag)

# Read and process reverse read size file
reverse_size <- read_tsv("./reverse_size_b4_5.tsv") %>%
  separate(file, into = c("text1", "tag", "text3"), sep = "/", remove = TRUE) %>%
  select(-c(text1, text3)) %>%
  rename(reverse_size = size)

# Read mapping coverage data
map_coverage <- read_tsv("./coverage_b4_5.tsv") %>%
  mutate(coverage = 100 - percent_gap)

# Read Kleborate result table
kleborate <- read_tsv("./kleborate_result_b4_5.txt")

# Read ENA accession metadata
ena <- read_csv("./run_kp_b4_5.csv") %>%
  rename(ENA_acc = id) %>%
  select(ENA_acc, alias)

# Merge all sequencing-related metadata
meta <- forward_size %>%
  left_join(reverse_size, by = "tag") %>%
  left_join(ena, by = c("tag" = "alias")) %>%
  left_join(map_coverage, by = "tag") %>%
  left_join(kleborate, by = c("tag" = "strain"))

# Export merged metadata
write_tsv(meta, "Kp_b4_5_sequencing_meta.tsv")

# -----------------------------------------------------------------------------
# Integrate GenBank accession metadata (Batch B1/B2)
# -----------------------------------------------------------------------------

meta_b1b2 <- read_tsv("kp_b1b2_sequencing_meta.tsv")
genbank_acc <- read_csv("./Genbank_acc_b1b2.csv")

meta_b1b2 <- meta_b1b2 %>%
  left_join(genbank_acc, by = c("tag" = "File"))

write_tsv(meta_b1b2, "kp_b1b2_sequencing_meta_25Jun.tsv")

# -----------------------------------------------------------------------------
# Bubble plot: Resistance vs Virulence by Sequence Type (ST)
# -----------------------------------------------------------------------------

# Summarize internal metadata by ST
plot_bubble <- internal_meta %>%
  group_by(ST_clean) %>%
  summarise(
    count = n(),
    mean_resistance = mean(resistance_score, na.rm = TRUE),
    mean_virulence = mean(virulence_score, na.rm = TRUE)
  )

# Base scatter plot
p1 <- ggplot(plot_bubble, aes(x = mean_virulence, y = mean_resistance))

# Add regions, points, and labels
p2 <- p1 +
  # Annotated background regions
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = Inf, alpha = 0.2, fill = "red") +
  annotate("rect", xmin = 3, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightblue") +
  # Bubble points
  geom_point(aes(size = count), shape = 21, fill = "grey20", alpha = 0.7) +
  # Labels for major STs
  geom_text_repel(
    aes(label = ifelse(count > 15, ST_clean, "")),
    size = 8, box.padding = 1, point.padding = 1.5,
    segment.color = "grey30", segment.size = 0.6, force = 8
  ) +
  # Labels for MDR and hv regions
  geom_text_repel(
    aes(label = ifelse(count <= 15 & mean_resistance > 1, ST_clean, "")),
    size = 8, box.padding = 1, point.padding = 1.5,
    segment.color = "grey30", segment.size = 0.6, force = 8
  ) +
  geom_text_repel(
    aes(label = ifelse(count <= 15 & mean_resistance <= 1 & mean_virulence > 3, ST_clean, "")),
    size = 8, box.padding = 1, point.padding = 1.5,
    segment.color = "grey30", segment.size = 0.6, force = 8
  ) +
  scale_size(range = c(3, 12)) +
  labs(x = "Virulence Score", y = "Resistance Score", size = "Count") +
  theme_bw() +
  theme(
    text = element_text(size = 40, face = "bold", color = "black"),
    axis.text = element_text(size = 40, face = "bold", color = "black"),
    legend.title = element_text(size = 40, face = "bold", color = "black"),
    legend.text = element_text(size = 40, face = "bold", color = "black")
  )

# Save bubble plot
ggsave(filename = "../bubble_plot.svg", plot = p2, device = "svg", width = 20, height = 16, units = "in")

# -----------------------------------------------------------------------------
# Combine metadata with GenBank accessions
# -----------------------------------------------------------------------------

meta <- read_tsv("./Kp_clincal_metadata_21Oct.tsv")
all_genbank <- read_csv("../Kp_ML/all_Genbank.csv")

meta_genbank <- meta %>%
  left_join(all_genbank, by = c("strain" = "sample"))

# Identify isolates without GenBank accession
no_genbank <- meta_genbank %>% filter(is.na(Genbank))

write_tsv(meta_genbank, "./Kp_clincal_metadata_21Oct.tsv")

# -----------------------------------------------------------------------------
# Kleborate annotation and longitudinal trends
# -----------------------------------------------------------------------------

kleborate <- read_tsv("../Kp_ML/Kp_kleborate_report_23Jun.tsv")

meta_kleborate <- meta_genbank %>%
  left_join(kleborate, by = "strain") %>%
  separate(ST, into = c("ST_clean", "st_lv"), sep = "-", remove = FALSE) %>%
  select(-st_lv)

# Define major STs of interest
st_list <- c("ST2096", "ST14", "ST147", "ST307", "ST45", "ST101", "ST11", "ST37", "ST23", "ST35")
meta_kleborate <- meta_kleborate %>%
  mutate(
    ST_focus = ifelse(ST_clean %in% st_list, ST_clean, "Others")
  ) %>%
  separate(Collection_Date, into = c("collect_year", "collect_month", "collect_day"), sep = "/", remove = FALSE)

# -----------------------------------------------------------------------------
# Temporal distribution of major STs
# -----------------------------------------------------------------------------

bar_st_df <- meta_kleborate %>%
  group_by(collect_year, ST_focus) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collect_year) %>%
  mutate(ratio = count / sum(count)) %>%
  filter(!collect_year %in% c(2014, 2015, 2016, 2017), !is.na(collect_year))

bar_st_df$ST_focus <- factor(bar_st_df$ST_focus, levels = c(st_list, "Others"))

ggplot(bar_st_df, aes(x = factor(collect_year), y = ratio, fill = ST_focus)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("#0000FF", "#800000", "#16a085", "#FFA500", "#2980b9",
               "#8e44ad", "#DB7F67", "#DBBEA1", "#A37B73", "#D34F73", "#799496"),
    name = "ST"
  ) +
  labs(x = "Year", y = "Proportion of Isolates") +
  theme_bw() +
  theme(
    text = element_text(size = 25, face = "bold", color = "black"),
    axis.text = element_text(size = 25, face = "bold", color = "black"),
    legend.title = element_text(size = 25, face = "bold", color = "black"),
    legend.text = element_text(size = 25, face = "bold", color = "black")
  )

# -----------------------------------------------------------------------------
# Classification by MDR and hypervirulence profiles
# -----------------------------------------------------------------------------

meta_kleborate <- meta_kleborate %>%
  mutate(
    MDR_VF = case_when(
      resistance_score >= 1 & virulence_score < 3 ~ "MDR",
      resistance_score < 1 & virulence_score >= 3 ~ "hv",
      resistance_score >= 1 & virulence_score >= 3 ~ "MDR-hv",
      TRUE ~ "Others"
    )
  )

bar_mdr_hv_df <- meta_kleborate %>%
  group_by(collect_year, MDR_VF) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collect_year) %>%
  mutate(ratio = count / sum(count)) %>%
  filter(!collect_year %in% c(2014, 2015, 2016, 2017), !is.na(collect_year))

bar_mdr_hv_df$MDR_VF <- factor(bar_mdr_hv_df$MDR_VF, levels = c("MDR", "hv", "MDR-hv", "Others"))

ggplot(bar_mdr_hv_df, aes(x = factor(collect_year), y = ratio, fill = MDR_VF)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("#c82926", "#3d68a4", "#5ca149", "grey"),
    labels = c(
      "MDR" = "Only Multi-drug Resistant",
      "hv" = "Only Hyper-virulent",
      "MDR-hv" = "Hyper-virulent Multi-drug Resistant",
      "Others" = "Others"
    ),
    name = "Strain Category"
  ) +
  labs(x = "Year", y = "Proportion of Isolates") +
  theme_bw() +
  theme(
    text = element_text(size = 25, face = "bold", color = "black"),
    axis.text = element_text(size = 25, face = "bold", color = "black"),
    legend.title = element_text(size = 25, face = "bold", color = "black"),
    legend.text = element_text(size = 25, face = "bold", color = "black")
  )

# -----------------------------------------------------------------------------
# Combine plasmid replicon data with Kleborate annotations
# -----------------------------------------------------------------------------

plasmid_replicon <- read_tsv("./plasmid/plasmid_amr_vf_replicon.tsv")
plasmid_kleborate <- read_tsv("./plasmid/plasmid_kleborate_report.tsv")

plasmid_replicon_kleborate <- plasmid_replicon %>%
  left_join(plasmid_kleborate, by = c("plasmid" = "strain"))

write_tsv(plasmid_replicon_kleborate, "./plasmid_replicon_kleborate.tsv")

############################################################
# End of Script
############################################################

