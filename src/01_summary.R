source("src/00_load_data.R")


# General data ----
## samples
ncol(counts) - 1

## amount of discarded data
## Amount of host data
## Amount of metagenomic data
## Amount of estimated prokaryotic (singlem)

# Mag statistics ----
## number of mags
nrow(taxonomy)
## mags without species-level annotation
taxonomy %>%
  filter(species == "s__") %>%
  nrow()

## Number of phylums
taxonomy %>%
  pull(phylum) %>%
  unique() %>%
  length()
