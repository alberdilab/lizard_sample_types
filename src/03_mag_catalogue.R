library(tidyverse)
source("src/00_load_data.R")



circular_tree <-
  tree %>%
  phytools::force.ultrametric(method = "extend", message = FALSE) %>%
  ggtree::ggtree(layout = "circular", size = 0.3, open.angle = 45)

circular_tree <-
  circular_tree %>%
  ggtree::gheatmap(data = mag_to_phylum, offset = 0.85, width = 0.1, colnames = FALSE) +
  scale_fill_manual(
    breaks = phylum_to_color$phylum,
    values = phylum_to_color$colors
  ) +
  ggtree::geom_tiplab2(size = 1, hjust = -0.1) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = margin(0, 0, 0, 0)
  )

# Flush color scale to enable a new color scheme in the next ring
circular_tree <- circular_tree + ggnewscale::new_scale_fill()

# Add completeness ring
circular_tree <-
  circular_tree +
  ggnewscale::new_scale_fill() +
  scale_fill_gradient(low = "#d1f4ba", high = "#f4baba") +
  ggtreeExtra::geom_fruit(
    data = mag_data,
    geom = geom_bar,
    mapping = aes(x = completeness, y = mag_id, fill = contamination),
    offset = 0.55,
    orientation = "y",
    stat = "identity"
  )


# Add genome-size ring
circular_tree <-
  circular_tree +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = "#cccccc") +
  ggtreeExtra::geom_fruit(
    data = mag_data,
    geom = geom_bar,
    mapping = aes(x = genome_size, y = mag_id),
    offset = 0.05,
    orientation = "y",
    stat = "identity"
  )

circular_tree


phyla_data <-
  mag_data %>%
  left_join(ehi_colors) %>%
  arrange(match(mag_id, tree$tip.label)) %>%
  select(phylum, colors) %>%
  unique() %>%
  mutate(phylum = str_remove_all(phylum, "p__") %>% factor())

phyla_data %>%
  ggplot() +
  geom_blank() +
  geom_rect(
    aes(
      xmin = 1:(nrow(phyla_data)) - 0.5,
      xmax = 1:(nrow(phyla)) + 0.5,
      ymin = 0.19,
      ymax = 0.2,
      fill = phylum
    )
  ) +
  scale_fill_manual(
    values = rev(phyla$colors)
  ) +
  geom_text(
    aes(
      x = 1:(nrow(phyla)),
      y = 0.15,
      label = rev(phylum)
    ),
    angle = 90,
    hjust = 0,
    size = 3
  ) +
  theme_void() +
  theme(legend.position = "none")


## 3.2 MAG quality ----

mag_quality <-
  mag_data %>%
  select(mag_id, domain, phylum, completeness, contamination, genome_size) %>%
  mutate(genome_size = round(genome_size / 1e6, 2)) %>%
  arrange(match(mag_id, tree$tip.label))

### 3.2.1 Biplot chart ----

mag_quality_biplot <-
  mag_quality %>%
  ggplot(
    aes(x = completeness, y = contamination, size = genome_size, color = phylum)
  ) +
  geom_point(alpha = 0.7) +
  ylim(c(25, 0)) +
  scale_color_manual(values = colors_alphabetic) +
  labs(y = "Contamination", x = "Completeness") +
  theme_classic() +
  theme(legend.position = "none")


### 3.2.2 Boxplots ----
contamination_bar <-
  mag_quality %>%
  ggplot(aes(y = contamination)) +
  ylim(c(25, 0)) +
  geom_boxplot(colour = "#999999", fill = "#cccccc") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0.40, 0), "inches")
  )


completeness_bar <-
  mag_data %>%
  ggplot(aes(x = completeness)) +
  xlim(c(50, 100)) +
  geom_boxplot(colour = "#999999", fill = "#cccccc") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0.50), "inches")
  )


layout_matrix <- matrix(2, nrow = 11, ncol = 11)
layout_matrix[1, ] <- 1
layout_matrix[, 11] <- 3
layout_matrix[11, 11] <- 4

gridExtra::grid.arrange(
  grobs = list(completeness_bar, mag_quality_biplot, contamination_bar),
  layout_matrix = layout_matrix
)

## 3.3 Functional attributes of MAGs ----

kegg_tree <-
  tree %>%
  phytools::force.ultrametric(method = "extend", message = FALSE) %>%
  ggtree::ggtree(size = 0.3) %>%
  ggtree::gheatmap(mag_to_phylum, offset = 0, width = 0.1, colnames = FALSE) +
  scale_fill_manual(values = colors_alphabetic) +
  ggnewscale::new_scale_fill()

kegg_tree %>%
  ggtree::gheatmap(kegg %>% column_to_rownames("mag_id"), offset = 0.5, width = 3.5, colnames = FALSE) +
  # vexpand(.08) +
  coord_cartesian(clip = "off") +
  scale_fill_gradient(low = "#f4f4f4", high = "steelblue", na.value = "white") +
  theme(legend.position = "none")


## 3.4 Functional ordination of MAGs ----
functional_tsne <-
  Rtsne::Rtsne(
    X = kegg %>% column_to_rownames("mag_id"),
    dims = 2,
    check_duplicates = FALSE
  )

functional_tsne$Y %>%
  as_tibble() %>%
  rename(tsne1 = V1, tsne2 = V2) %>%
  mutate(mag_id = kegg$mag_id) %>%
  left_join(kegg) %>%
  left_join(mag_data) %>%
  select(mag_id, phylum, tsne1, tsne2, completeness) %>%
  ggplot(
    aes(x = tsne1, y = tsne2, color = phylum, size = completeness)
  ) +
  geom_point(shape = 16, alpha = 0.7) +
  scale_color_manual(values = colors_alphabetic) +
  theme_minimal() +
  theme(legend.position = "none")
