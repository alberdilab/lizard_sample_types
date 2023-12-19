source("src/00_compose_ehi_metadata_file.R")
source("src/00_load_data.R")
source("src/04_dna_fractions.R")
source("src/05_count_data.R")

library(hilldiv2)
library(vegan)
library(ggpubr)

#7.1 Data prepraration

present_mags <-
  counts_filtered %>%
  pivot_longer(-mag_id) %>%
  filter(sum(value) > 0, .by = mag_id) %>%
  pull(mag_id) %>%
  unique()

counts_filtered <-
  counts_filtered %>%
  pivot_longer(-mag_id) %>%
  filter(sum(value) > 0, .by = name) %>%
  pivot_wider()

present_mags <- present_mags[present_mags %in% kegg$mag_id]

kegg_filtered <-
  kegg %>%
  filter(mag_id %in% present_mags) %>%
  pivot_longer(-mag_id) %>%
  filter(0 < sum(value), sum(value) < 1, .by = name) %>%
  pivot_wider()

counts_filtered_again <- counts_filtered %>% filter(mag_id %in% present_mags)

## 7.2 Alpha diversity



q0n <- counts_filtered %>%
  column_to_rownames("mag_id") %>%
  hilldiv(q = 0) %>%
  c()

q1n <-
  counts_filtered %>%
  column_to_rownames("mag_id") %>%
  hilldiv(q = 1) %>%
  c()

q1p <-   counts_filtered %>%
  column_to_rownames("mag_id") %>%
  hilldiv(q = 1, tree = tree) %>%
  c()

# dist is already a function!
dist <-
  kegg_filtered %>%
  column_to_rownames("mag_id") %>%
  traits2dist(method = "gower")

q1f <- counts_filtered_again %>%
  column_to_rownames("mag_id") %>%
  hilldiv(q = 1, dist = dist) %>%
  c()


alpha_diversity <-
  tibble(
    sample_id = colnames(counts_filtered)[2:ncol(counts_filtered)],
    richness = q0n,
    neutral = q1n %>% round(3),
    phylogenetic = q1p %>% round(3),
    functional = q1f %>% round(3),
  )

columns <- c(
  "richness", "neutral", "phylogenetic", "functional", "mapped", "total"
)

reads <-
  reads %>%
  pivot_wider()

### Plot it
alpha_diversity <-
  alpha_diversity %>%
  left_join(reads) %>%
  mutate(
    total = mapped_mags + mapped_SceUnd + unmapped_reads + low_quality,
    mapped = mapped_mags / 1e6,
    total = total / 1e6,
  ) %>%
  select(-low_quality, -mapped_SceUnd, -mapped_mags, -unmapped_reads)

alpha_diversity %>%
  pivot_longer(-sample_id, names_to = "metric") %>%
  mutate(
    metric = factor(metric, levels = columns)
  ) %>%
  ggplot(aes(x = value, y = sample_id)) +
  geom_bar(stat = "identity", fill = "#6c9ebc") +
  facet_wrap(~ metric, scales = "free_x", ncol = 6) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_line(linewidth = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### Alpha diversity comparisons
alpha_colors <- c(
  "#e5bd5b", "#6b7398", "#76b183", "#d57d2c", "#2a2d26", "#f9d4cc", "#3c634e",
  "#ea68c3"
)

# I don't know how the samples are distributed. Let's compare cloaca vs faecal
group_n <- 2
sample_data <-
  tibble(sample_id = alpha_diversity$sample_id) %>%
  rowwise() %>%
  mutate(
    tissue = sample_id %>% str_split("\\.") %>% unlist() %>% .[2]
  ) %>%
  ungroup()


 ### Neutral diversity ----
alpha_diversity %>%
  select(sample_id, neutral) %>%
  pivot_longer(-sample_id) %>%
  mutate(name = factor(name, levels = columns)) %>%
  left_join(sample_data) %>%
  ggboxplot(x = "tissue", y = "value", color = "tissue", fill = "tissue") +
  scale_color_manual(values = alpha_colors[c(1:group_n)]) +
  scale_fill_manual(values = paste0(alpha_colors[c(1:group_n)], "50")) +
  stat_compare_means() +
  theme_classic() +
  labs(y = "Neutral Hill Numbers") +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(color = guide_legend(title = "Tissue"), fill = "none")


### Phylogenetic diversity ----
alpha_diversity %>%
  select(sample_id, phylogenetic) %>%
  pivot_longer(-sample_id) %>%
  mutate(name = factor(name, levels = columns)) %>%
  left_join(sample_data) %>%
  ggboxplot(x = "tissue", y = "value", color = "tissue", fill = "tissue") +
  scale_color_manual(values = alpha_colors[c(1:group_n)]) +
  scale_fill_manual(values = paste0(alpha_colors[c(1:group_n)], "50")) +
  stat_compare_means() +
  theme_classic() +
  labs(y = "Phylogenetic Hill Numbers") +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(color = guide_legend(title = "Tissue"), fill = "none")


### Functional diversity ----
alpha_diversity %>%
  select(sample_id, phylogenetic) %>%
  pivot_longer(-sample_id) %>%
  mutate(name = factor(name, levels = columns)) %>%
  left_join(sample_data) %>%
  ggboxplot(x = "tissue", y = "value", color = "tissue", fill = "tissue") +
  scale_color_manual(values = alpha_colors[c(1:group_n)]) +
  scale_fill_manual(values = paste0(alpha_colors[c(1:group_n)], "50")) +
  stat_compare_means() +
  theme_classic() +
  labs(y = "Functional Hill Numbers") +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(color = guide_legend(title = "Tissue"), fill = "none")



### Relationship between alpha diversity and sequencing effort
ggplot(alpha_diversity, aes(x = mapped, y = neutral, label = sample_id)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#e08dde", fill = "#e08dde") +
  geom_point(alpha = 0.5, color = "#6c9ebc") +
  ggrepel::geom_label_repel(max.overlaps = 100, cex = 0.7) +
  labs(
    x = "GBs mapped to MAGs",
    y = "Neutral diversity (effective number of MAGs)"
  ) +
  theme_classic() +
  theme(legend.position = "none")


## 7.3 Beta diversity
beta_colors <- c(
  "#e5bd5b", "#6b7398", "#76b183", "#d57d2c", "#2a2d26", "#f9d4cc", "#3c634e",
  "#ea68c3"
)

beta_q1n <-
  counts_filtered %>%
  column_to_rownames("mag_id") %>%
  hillpair(q = 1, metric = "S")

sample_table_adonis <-
  sample_data %>%
  filter(sample_id %in% labels(beta_q1n)) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

adonis2(
  formula = beta_q1n ~ .,
  data = sample_data %>% column_to_rownames("sample_id"),
  permutations = 999
)


### Beta diversity plot
beta_q1n_nmds <-
  beta_q1n %>%
  metaMDS(., trymax = 500, k = 2, verbosity = FALSE) %>%
  vegan::scores() %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(sample_data, by = "sample_id")

group_n <- 2
beta_q1n_nmds %>%
  mutate(
    x_cen = mean(NMDS1, na.rm = TRUE),
    y_cen = mean(NMDS2, na.rm = TRUE),
    .by = tissue
  ) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = tissue)) +
  scale_color_manual(values = beta_colors) +
  geom_point(size = 2) +
  geom_segment(aes(
    x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2, alpha = 0.2
  )) +
  theme_classic() +
  theme(legend.position = "right", legend.box = "vertical") +
  guides(color = guide_legend(title = "Tissue"))
