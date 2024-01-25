library(tidyverse)

source("src/00_load_data.R")

# 5.1  Coverage filtering ----
min_coverage <- 0.3

counts_long <-
  counts %>%
  pivot_longer(-sequence_id, names_to = "sample_id", values_to = "counts") %>%
  rename(mag_id = sequence_id)

coverage_long <-
  coverage %>%
  pivot_longer(-sequence_id, names_to = "sample_id", values_to = "coverage") %>%
  rename(mag_id = sequence_id) %>%
  left_join(mag_data %>% select(mag_id, genome_size)) %>%
  mutate(coverage_ratio = coverage / genome_size)

counts_filtered <-
  coverage_long %>%
  mutate(binary = coverage_ratio > 0.3) %>%
  left_join(counts_long) %>%
  mutate(counts = counts * binary) %>%
  select(mag_id, sample_id, counts) %>%
  pivot_wider(names_from = sample_id, values_from = counts)

# counts_long %>%
#   separate(sep = "\\.", col = "sample_id", into = c("sample", "tissue"), remove = F) %>%
#   group_by(tissue) %>%
#   summarise(counts = sum(counts)) %>%
#   mutate(percent = counts / sum(counts) * 100)


# 5.2 Genome size normalization
read_length <- 150

mag_size_scaled <-
  mag_data %>%
  filter(mag_id %in% counts_filtered$mag_id) %>%
  select(mag_id, genome_size) %>%
  mutate(mag_size_scaled = genome_size / 150) %>%
  pull(mag_size_scaled)


## 5.3 Count table

vertical_tree <-
  tree %>%
  phytools::force.ultrametric(method = "extend", message = FALSE) %>%
  ggtree::ggtree(size = 0.3) %>%
  ggtree::gheatmap(data = mag_to_phylum, offset = 0.85, width = 0.1, colnames = FALSE) +
  scale_fill_manual(
    breaks = phylum_to_color$phylum,
    values = phylum_to_color$colors
  )

vertical_tree <- vertical_tree + ggnewscale::new_scale_fill()

vertical_tree <-
  vertical_tree %>%
  ggtree::gheatmap(
    counts_filtered %>% column_to_rownames("mag_id") %>% log10(),
    offset = 0.04,
    width = 3.5,
    colnames = TRUE,
    colnames_angle = 90,
    font.size = 2,
    colnames_position = "top",
    colnames_offset_y = 9
  ) +
  ggtree::vexpand(.08) +
  coord_cartesian(clip = "off") +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "white") +
  theme(legend.position = "none")


# Add count scale

count_scale <-
  counts_filtered %>%
  column_to_rownames("mag_id") %>%
  max() %>%
  log() %>%
  seq(to = 0, length.out = 5)

count_legend <-
  tibble(
    value = count_scale,
    x = c(
      count_scale %>% exp() %>% round(2)
    )
  ) %>%
  mutate(x = factor(x, levels = x)) %>%
  ggplot(., aes(x = x, y = 0.2)) +
  geom_tile(aes(fill = value, y = 0.2), color = "#CCCCCC") +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "white") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8)
  )

# Arrange both legends
gridExtra::grid.arrange(
  grobs = list(
    vertical_tree,
    count_legend
  ),
  layout_matrix = rbind(
    c(1, 1, 2),
    c(1, 1, 3)
  )
)
