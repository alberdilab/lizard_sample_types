source("src/00_load_data.R")
source("src/00_compose_ehi_metadata_file.R")
# source("src/00_parse_reports.R")

# 4.1 Quality mappings ----
reads <-
  ehi_metadata %>%
  mutate(
    low_quality = raw_reads - trimmed_reads,
    unmapped_reads = trimmed_reads - mapped_SceUnd - mapped_mags
  ) %>%
  select(sample_id, low_quality, mapped_SceUnd, mapped_mags, unmapped_reads) %>%
  pivot_longer(-sample_id)

reads %>%
  separate(
    col = "sample_id", into = c("sample", "tissue"), sep = "\\.", remove = FALSE
  ) %>%
  ggplot(aes(x = sample_id, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#CCCCCC", "#178a94", "#ee8080", "#d03161")) +
  facet_wrap(~tissue, scales = "free")
# geom_bar(stat = "identity")


# 4.2 Estimated vs prokaryotic

singlem_data <-
  ehi_metadata %>%
  mutate(
    unmapped_reads = trimmed_reads - mapped_SceUnd - mapped_mags,
    mag_proportion = mapped_mags / (mapped_mags + unmapped_reads),
    singlem_read_fraction = singlem_read_fraction
  ) %>%
  select(sample_id, mag_proportion, singlem_read_fraction) %>%
  mutate(
    mag_proportion = if_else(singlem_read_fraction == 0, 0, mag_proportion),
    singlem_read_fraction = if_else(singlem_read_fraction == 0, NA, singlem_read_fraction),
    singlem_read_fraction = if_else(singlem_read_fraction < mag_proportion, NA, singlem_read_fraction),
    singlem_read_fraction = if_else(singlem_read_fraction > 1, 1, singlem_read_fraction)
  )


singlem_data %>%
  pivot_longer(-sample_id, names_to = "proportion", values_to = "value") %>%
  mutate(
    proportion = factor(
      proportion,
      levels = c("mag_proportion", "singlem_read_fraction")
    )
  ) %>%
  ggplot(aes(x = value, y = sample_id, color = proportion)) +
  geom_line(aes(group = sample_id), color = "#f8a538") +
  geom_point() +
  scale_color_manual(values = c("#52e1e8", "#876b53")) +
  theme_classic() +
  labs(y = "Samples", x = "Prokaryotic fraction") +
  scale_x_continuous(limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 6
    ),
    legend.position = "right"
  )


## Additional sequencing needed

mapped_mags_aim <- 2 * 1e9 / 150
mapped_host_aim <- 5 * 1e9 / 150

# sequence_fractions_required <-
ehi_metadata %>%
  select(
    -singlem_metagenome_size, -singlem_prokaryotic_bases, -singlem_read_fraction
  ) %>%
  mutate(
    mapped_mags = mapped_mags,
    lowqual_reads = raw_reads - trimmed_reads,
    unmapped_bases = trimmed_reads - mapped_SceUnd - mapped_mags,
    mapped_mags_fraction = mapped_mags / raw_reads,
    mapped_mags_difference = mapped_mags_aim - mapped_mags,
    meta_required = mapped_mags_difference / mapped_mags_fraction,
    meta_required = if_else(meta_required < 0, 0, meta_required),
    host_fraction = mapped_SceUnd / raw_reads,
    host_difference = mapped_host_aim / mapped_SceUnd,
    host_required = host_difference / host_fraction,
    host_required = if_else(host_required < 0, 0, host_required)
  ) %>%
  select(sample_id, mapped_mags, unmapped_bases, mapped_SceUnd, lowqual_reads)
# select(sample, mapped_mags, unmapped_bases, host_bases, lowqual_bases, meta_required, host_required)
