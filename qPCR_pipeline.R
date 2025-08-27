# ─────────────────────────────────────────────
# qPCR Analysis Script (ΔΔCt and RQ calculation)
# Author: Nicolas Trenk
# ─────────────────────────────────────────────

# ─────────────────────────────────────────────
# USER-DEFINED VARIABLES
# ─────────────────────────────────────────────

# Input files
raw_data_file <- "data/qPCR_raw.xls"       # StepOnePlus raw export file
metadata_file <- "data/sample_metadata.csv" # Sample key / grouping info

# Analysis settings
housekeeping_gene <- "ACTIN"                 # Replace with your housekeeping gene
target_primer <- c("Primer pair 1", "Primer pair 2") # Replace with your primer names to include
baseline_sample_group <- "Baseline sample"           # Group used as ΔΔCt baseline
omit_samples <- c("Sample to omit")  # Optional samples to exclude

# Plot settings
x_axis_label <- "Sample group"
y_axis_label <- "Relative quantification"
output_stats_file <- "results/qPCR_stats.xlsx"

# ─────────────────────────────────────────────
# LOAD PACKAGES
# ─────────────────────────────────────────────
library(dplyr)
library(tidyverse)
library(EnvStats)
library(rstatix)
library(readxl)
library(ggpubr)

# ─────────────────────────────────────────────
# DATA IMPORT & CLEANING
# ─────────────────────────────────────────────

raw_data <- read_csv(raw_data_file)
key <- read_csv(metadata_file)

# Merge raw data with metadata
raw_data_joined <- left_join(
  raw_data,
  key,
  by = c("Sample Name" = "Sample.ID")
)

# Clean up group names (generic example, edit if needed)
raw_data_renamed <- raw_data_joined %>%
  mutate(
    Sample_grouped = recode(Sample_grouped,
      "Old_name" = "New_name",
      "Old_name_2" = "New_name_2",
      .default = Sample_grouped
    )
  )

# Filter dataset (remove NTCs and empty wells)
dat_filtered <- raw_data_renamed %>%
  filter(`Sample Name` != "" & Task != "NTC") %>%
  filter(!(Sample_grouped %in% omit_samples)) %>%
  filter(`Target Name` %in% target_primer | `Target Name` == housekeeping_gene)

# Flag late housekeeping amplification
dat_filtered <- dat_filtered %>%
  mutate(late_housekeeping = ifelse(`Target Name` == housekeeping_gene & Cq > 30, TRUE, FALSE)) %>%
  filter(!late_housekeeping)

# ─────────────────────────────────────────────
# ΔCt CALCULATIONS
# ─────────────────────────────────────────────

housekeeping_values <- dat_filtered %>%
  mutate(housekeep_Cq = as.numeric(Cq)) %>%
  filter(`Target Name` == housekeeping_gene) %>%
  select(Sample_unique, housekeep_Cq)

dat_ΔCT <- dat_filtered %>%
  select(Sample_grouped, Sample_unique, Sample.name, Accession, Timepoint, Experiment, `Target Name`, Cq) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  filter(`Target Name` != housekeeping_gene) %>%
  left_join(housekeeping_values, by = "Sample_unique") %>%
  mutate(ΔCT = Cq - housekeep_Cq)

# ─────────────────────────────────────────────
# ΔΔCt AND RELATIVE QUANTIFICATION
# ─────────────────────────────────────────────

baseline <- dat_ΔCT %>%
  filter(Sample_grouped %in% baseline_sample_group) %>%
  group_by(`Target Name`) %>%
  summarise(ΔCT_baseline = mean(ΔCT, na.rm = TRUE))

dat_ΔΔCT <- dat_ΔCT %>%
  left_join(baseline, by = "Target Name") %>%
  mutate(ΔΔCT = ΔCT - ΔCT_baseline)

dat_RQ <- dat_ΔΔCT %>%
  mutate(RQ = 2^(-ΔΔCT),
         log10RQ = log10(RQ))

# ─────────────────────────────────────────────
# STATISTICAL COMPARISON
# ─────────────────────────────────────────────

# Normality check
dat_RQ %>%
  group_by(Sample_grouped) %>%
  summarise(p_value = shapiro.test(RQ)$p.value)

# Pairwise Wilcoxon tests
stats_results <- dat_RQ %>%
  pairwise_wilcox_test(RQ ~ Sample_grouped, p.adjust.method = "BH")

# Save stats
write.csv(stats_results, output_stats_file, row.names = FALSE)

# ─────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────

dat_RQ <- dat_RQ %>%
  mutate(Sample_grouped = factor(Sample_grouped))

rq_plot <- ggplot(data=dat_RQ, aes(x = Sample_grouped, y = log10(RQ))) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(colour = Accession), width = 0.1) +
  geom_hline(yintercept = log10(1), linetype = "dashed", color = "gray40") +
  labs(x = x_axis_label, y = y_axis_label) +
  stat_compare_means() +
  stat_n_text() +
  theme_pubr()

ggsave("results/qPCR_boxplot.png", rq_plot, width = 6, height = 4)

