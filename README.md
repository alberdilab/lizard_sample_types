# Sample type comparison for lizard gut metagenomics
Repository of the data and analysis procedures for the manuscript:

**Contrasting recovery of metagenome‑assembled genomes and derived microbial communities from lizard fecal and cloacal samples**
Mauricio Hernández, Jorge Langa, Ostaizka Aizpurua, Yendi E. Navarro-Noya, Antton Alberdi

## Bioinformatic procedures

Data processing to generate annotated metagenome-assembled genomes and genome count tables was conducted using the following Snakemake pipeline: [mg_assembly](https://github.com/3d-omics/mg_assembly). Data analysis procedures source from the outputs of this pipeline.

## Analysis procedures

The raw code used for data analysis is in the **Rmd** files stored in the root directory of this repository, while the bookdown-rendered webbook is available at:

[alberdilab.github.io/lizard_sample_types](https://alberdilab.github.io/lizard_sample_types)

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the Rmd documents. Note that the code chunks that require heavy computation have been tuned off using 'eval=FALSE'. To re-render the webbook, you can use the following code:

```r
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```