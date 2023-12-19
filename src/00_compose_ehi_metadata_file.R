library(tidyverse)

raw_reads <-
  "resources/report/by_step/reads_data/multiqc_general_stats.txt" %>%
  read_tsv() %>%
  select(
    sample_id = Sample,
    raw_reads = `FastQC_mqc-generalstats-fastqc-total_sequences`
  ) %>%
  mutate(
    sample_id = sample_id %>% str_remove_all("_1$") %>% str_remove_all("_2$")
  ) %>%
  summarise(raw_reads = sum(raw_reads), .by = sample_id)

fastp_reads <-
  "resources/report/by_step/preprocessing_data/multiqc_general_stats.txt" %>%
  read_tsv() %>%
  filter(str_detect(Sample, "fastp")) %>%
  select(
    sample_id = Sample,
    trimmed_reads = `FastQC_mqc-generalstats-fastqc-total_sequences`
  ) %>%
  mutate(
    sample_id =
      sample_id %>%
      str_remove_all("_[u12]+$") %>%
      str_remove_all("^fastp \\| ")
  ) %>%
  summarise(trimmed_reads = sum(trimmed_reads), .by = sample_id)

host_mapped <-
  "resources/report/by_step/preprocessing_data/multiqc_general_stats.txt" %>%
  read_tsv() %>%
  filter(!str_detect(Sample, "fastp")) %>%
  select(
    sample_id = Sample,
    host_mapped = `Samtools_mqc-generalstats-samtools-reads_mapped`,
    mapping_total = `Samtools_mqc-generalstats-samtools-raw_total_sequences`
  ) %>%
  mutate(
    host_unmapped = mapping_total - host_mapped
  ) %>%
  filter(!is.na(host_mapped)) %>%
  separate(col = sample_id, into = c("host_name", "sample_id"), sep = " \\| ") %>%
  rename(mapped = host_mapped, unmapped = host_unmapped) %>%
  select(-mapping_total) %>%
  pivot_longer(-host_name:-sample_id) %>%
  mutate(
    name = str_glue("{name}_{host_name}")
  ) %>%
  select(-host_name) %>%
  pivot_wider()

singlem <-
  "resources/singlem/microbial_fraction.tsv" %>%
  read_tsv() %>%
  distinct() %>%
  mutate(
    sample_id = sample %>% str_remove_all("_1$"),
    read_fraction = read_fraction %>% str_remove("%") %>% as.numeric(),
    read_fraction = read_fraction / 100
  ) %>%
  select(
    sample_id,
    singlem_prokaryotic_bases = bacterial_archaeal_bases,
    singlem_metagenome_size = metagenome_size,
    singlem_read_fraction = read_fraction,
  )

nonpareil <-
  "resources/nonpareil/nonpareil.tsv" %>%
  read_tsv() %>%
  select(sample_id, nonpareil_c = C, nonpareil_diversity = diversity)

mag_mapping <-
  "resources/coverm/contig.count.tsv.gz" %>%
  read_tsv() %>%
  pivot_longer(-sequence_id) %>%
  summarise(value = sum(value), .by = "name") %>%
  rename(sample_id = name, mapped_mags = value)

 ehi_metadata <-
   raw_reads %>%
   left_join(fastp_reads) %>%
   left_join(host_mapped) %>%
   left_join(singlem) %>%
   left_join(nonpareil) %>%
   left_join(mag_mapping)

rm(raw_reads, fastp_reads, host_mapped, singlem, nonpareil, mag_mapping)
