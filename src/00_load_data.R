library(tidyverse)

counts <- read_tsv("resources/coverm/genome.count.tsv.xz")
coverage <- read_tsv("resources/coverm/genome.covered_bases.tsv.xz")
taxonomy <-
  read_tsv("resources/gtdbtk/gtdbtk.summary.tsv.xz") %>%
  select(mag_id = user_genome, classification) %>%
  separate(
    classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";"
  )

checkm2 <-
  read_tsv("resources/checkm2/quality_report.tsv.xz") %>%
  select(
    mag_id = Name, completeness = Completeness, contamination = Contamination,
    genome_size = Genome_Size, gc_content = GC_Content
  )

good_mag_ids <-
  checkm2 %>%
  filter(completeness >= 70, contamination <= 10) %>%
  mutate(mag_id = str_glue("\'{mag_id}\'")) %>%
  pull(mag_id)

tree <- ape::read.tree("resources/gtdbtk/gtdbtk.backbone.bac120.classify.tree")
# tips_to_keep <-
#   tree$tip.label %>%
#   tibble(node_id = .) %>%
#   filter(
#     !str_detect(node_id, "^GB_"),
#     !str_detect(node_id, "^RS_")
#   ) %>%
#   pull(node_id)

tree <-
  tree %>%
  # ape::keep.tip(tips_to_keep) %>%
  ape::keep.tip(good_mag_ids)

tree$tip.label <- tree$tip.label %>% str_remove_all("\\'")

plot(tree)

rm(tips_to_keep)

# colors
ehi_colors <-
  read_tsv("https://raw.githubusercontent.com/earthhologenome/EHI_taxonomy_colour/main/ehi_phylum_colors.tsv")


colors_alphabetic <-
  ehi_colors %>%
  right_join(taxonomy) %>%
  select(phylum, colors) %>%
  unique() %>%
  arrange(phylum) %>%
  pull(colors)

# dram
dram <-
  read_tsv("resources/dram/annotations.tsv.xz") %>%
  rename(gene_id = `...1`, mag_id = fasta, contig_id = scaffold) %>%
  select(mag_id, contig_id, gene_id, everything())


kegg <-
  read_tsv("resources/dram/product.tsv.xz") %>%
  rename(mag_id = genome)

mag_data <-
  taxonomy %>%
  left_join(checkm2)

rm(checkm2)

# Interesting tables
phyla <-
  mag_data %>%
  left_join(ehi_colors) %>%
  arrange(match(mag_id, tree$tip.label)) %>%
  select(phylum, colors) %>%
  unique()

phylum_to_color <-
  mag_data %>%
  left_join(ehi_colors) %>%
  select(phylum, colors) %>%
  unique()

mag_to_phylum <-
  mag_data %>%
  select(mag_id, phylum) %>%
  column_to_rownames("mag_id")

sample_data <-
  counts %>%
  pivot_longer(-sequence_id) %>%
  select(name) %>%
  separate(
    col = name, into = c("sample_name", "tissue"), sep = "\\.", remove = FALSE
  ) %>%
  rename(sample_id = name) %>%
  distinct(sample_id, .keep_all = T)
