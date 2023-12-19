# library(tidyverse)
# library(janitor)
#
#
# # reads ----
#
# reads_fastqc <-
#   read_tsv("resources/report/by_step/reads_data/multiqc_fastqc.txt.xz") %>%
#   clean_names() %>%
#   select(
#     -per_tile_sequence_quality, -basic_statistics, -sequence_duplication_levels,
#     -overrepresented_sequences, -per_base_sequence_quality,
#     -per_sequence_quality_scores, -adapter_content, -file_type,
#     -per_sequence_gc_content, -sequence_length_distribution, -encoding,
#     -per_base_sequence_content, -per_base_n_content
#   )
#
# reads_general <-
#   read_tsv("resources/report/by_step/reads_data/multiqc_general_stats.txt.xz") %>%
#   clean_names() %>%
#   rename_with(~str_remove(., "fast_qc_mqc_generalstats_fastqc_"))
#
# reads <-
#   left_join(
#     reads_general,
#     reads_fastqc %>%
#       select(-total_sequences, -avg_sequence_length, -percent_gc, -median_sequence_length),
#     by = "sample"
#   ) %>%
#   select(-filename) %>%
#   mutate(
#     sample = sample %>% str_remove("_[12]$")
#   ) %>%
#   group_by(sample) %>%
#   summarise(
#     percent_duplicates = mean(percent_duplicates),
#     percent_gc = mean(percent_gc),
#     avg_sequence_length = mean(avg_sequence_length),
#     total_sequences = sum(total_sequences),
#     percent_fails = mean(percent_fails)
#   ) %>%
#   mutate(
#     total_bases = total_sequences * avg_sequence_length
#   ) %>%
#   rename(
#     reads_percent_duplicates = percent_duplicates,
#     reads_percent_gc = percent_gc,
#     reads_avg_sequence_length = avg_sequence_length,
#     reads_total_sequences = total_sequences,
#     reads_percent_fails = percent_fails,
#     reads_total_bases = total_bases,
#   )
#
# # adapters <-
# #   "resources/report/by_step/reads_data/mqc_fa"
# # read_tsv("resources/reports/by_step/reads_data/mqc_fastqc_adapter_content_plot_1.txt") %>%
# #   pivot_longer(-Sample, names_to = "position", values_to = "percent" ) %>%
# #   clean_names() %>%
# #   mutate(
# #     sample = sample %>% str_remove(" - illumina_universal_adapter"),
# #     position = as.integer(position)
# #   ) %>%
# #   arrange(sample, position) %>%
# #   slice(n(), .by = sample) %>%
# #   mutate(
# #     sample = sample %>% str_remove("_[12]$")
# #   ) %>%
# #   summarise(percent_adapters = mean(percent), .by = sample) %>% view()
# #
# # rm(reads_fastqc, reads_general)
#
#
# # fastp ----
#
# fastp_fastqc <-
#   "resources/report/by_step/preprocessing_data/multiqc_fastqc.txt.xz" %>%
#   read_tsv() %>%
#   clean_names() %>%
#   select(
#     sample, percent_gc, avg_sequence_length, total_deduplicated_percentage,
#     total_sequences, median_sequence_length
#   )
#
# fastp_general <-
#   "resources/report/by_step/preprocessing_data/multiqc_general_stats.txt.xz" %>%
#   read_tsv() %>%
#   clean_names() %>%
#   rename_with(~str_remove(., "fast_qc_mqc_generalstats_fastqc_"))
#
# fastp <-
#   left_join(fastp_fastqc, fastp_general) %>%
#   mutate(
#     sample = sample %>% str_remove("_[12]$")
#   ) %>%
#   summarise(
#     fastp_percent_gc = mean(percent_gc),
#     fastp_total_bases = sum(total_sequences * avg_sequence_length),
#     fastp_avg_sequence_length = mean(avg_sequence_length),
#     fastp_total_deduplicated_percentage = mean(total_deduplicated_percentage),
#     fastp_total_sequences = sum(total_sequences),
#     fastp_percent_fails = mean(percent_fails),
#     fastp_percent_duplicates = mean(percent_duplicates),
#     .by = sample
#   )
# rm(fastp_fastqc, fastp_general)
#
#
# # bowtie2 ----
#
#
# human <-
#
#   "resources/reportsby_step/bowtie2_host_Homo_sapiens.GRCh38_data/multiqc_general_stats.txt"
#   read_tsv() %>%
#   clean_names() %>%
#   rename_with(~str_remove(., "samtools_mqc_generalstats_samtools_")) %>%
#   mutate(
#     sample = sample %>% str_remove("Homo_sapiens.GRCh38\\ \\|\\ ")
#   ) %>%
#   select(sample, reads_to_map = raw_total_sequences, human_reads_mapped = reads_mapped)
#
# chicken <-
#   read_tsv("resources/reportsby_step/bowtie2_host_Gallus_gallus.bGalGal1_data/multiqc_general_stats.txt") %>%
#   clean_names() %>%
#   rename_with(~str_remove(., "samtools_mqc_generalstats_samtools_")) %>%
#   mutate(
#     sample = sample %>% str_remove("Gallus_gallus.bGalGal1\\ \\|\\ ")
#   ) %>%
#   select(sample, chicken_reads_mapped = reads_mapped)
#
# pig <-
#   read_tsv("resources/reportsby_step/bowtie2_host_Sus_scrofa.Sscrofa11.1_data/multiqc_general_stats.txt") %>%
#   clean_names() %>%
#   rename_with(~str_remove(., "samtools_mqc_generalstats_samtools_")) %>%
#   mutate(
#     sample = sample %>% str_remove("Sus_scrofa.Sscrofa11.1\\ \\|\\ ")
#   ) %>%
#   select(sample, pig_reads_mapped = reads_mapped)
#
# arabidopsis <-
#   read_tsv("resources/reportsby_step/bowtie2_host_Arabidopsis_thaliana.TAIR10_data/multiqc_general_stats.txt") %>%
#   clean_names() %>%
#   rename_with(~str_remove(., "samtools_mqc_generalstats_samtools_")) %>%
#   mutate(
#     sample = sample %>% str_remove("Arabidopsis_thaliana.TAIR10\\ \\|\\ ")
#   ) %>%
#   select(sample, arabidopsis_reads_mapped = reads_mapped)
#
# # nonpareil
#
# nonpareil <-
#   read_tsv("nonpareil.tsv") %>%
#   clean_names() %>%
#   mutate(estimated_coverage = lr / l_rstar) %>%
#   rename(sample=sample_id)
#
# # singlem
#
# singlem <-
#   read_tsv("singlem.tsv") %>%
#   filter(str_detect(taxonomy, "Root; d__Bacteria")) %>%
#   group_by(sample) %>%
#   summarise(sum_coverage = sum(coverage))
#
#
#
# # mags - coverm
#
# coverm <-
#   read_tsv("coverm_genome_big.count.tsv") %>%
#   pivot_longer(-Genome) %>%
#   mutate(value = as.integer(value)) %>%
#   summarise(counts = sum(value), .by = name) %>%
#   mutate(sample = name %>% str_remove(" Read Count")) %>%
#   select(sample, mags_reads_mapped = counts)
#
#
#
# # big fucking table
#
# bft <-
#   reads %>%
#   left_join(fastp) %>%
#   left_join(nonpareil) %>%
#   left_join(singlem) %>%
#   left_join(human) %>%
#   left_join(chicken) %>%
#   left_join(pig) %>%
#   left_join(arabidopsis) %>%
#   left_join(coverm) %>%
#   mutate(
#     percent_human = human_reads_mapped / reads_to_map * 100,
#     percent_chicken = chicken_reads_mapped / reads_to_map * 100,
#     percent_pig = pig_reads_mapped / reads_to_map * 100,
#     percent_arabidopsis = arabidopsis_reads_mapped / reads_to_map * 100,
#     percent_mags = mags_reads_mapped / reads_to_map * 100,
#     percent_unmapped = 100 - percent_human - percent_chicken - percent_pig - percent_arabidopsis - percent_mags
#   )
#
# write_tsv(bft, "bft.tsv")
