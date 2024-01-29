#!/usr/bin/env Rscript

source("src/00_compose_ehi_metadata_file.R")
source("src/00_load_data.R")
source("src/05_count_data.R")


taxonomy_data <-
  counts_filtered %>%
  pivot_longer(-mag_id) %>%
  left_join(mag_data) %>%
  left_join(ehi_colors) %>%
  filter(value > 0) %>%
  select(mag_id, name, value, phylum, colors)

colors <-
  taxonomy_data %>%
  select(phylum, colors) %>%
  distinct()
color_vector <- colors$colors
names(color_vector) <- colors$phylum


ggplot(taxonomy_data, aes(x = name, y = value, fill = phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_color_manual(values = color_vector)
