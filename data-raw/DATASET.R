library(devtools)
library(data.table)
# library(tidyverse)
library(usethis)

## MAQC2 DESeq2 results

maqc_deseq2 <- fread("data-raw/MAQC2_DESeq2_result.tsv", quote = "")

usethis::use_data(maqc_deseq2, overwrite = TRUE)


## version 1
library(Lazy2)
intg_pathways <- readGMT("data-raw/pathway_gmt/KPath_v1_genesets.gmt")
intg_pathinfo <- fread("data-raw/pathway_gmt/KPath_v1_Reference.txt")

usethis::use_data(intg_pathways, overwrite = TRUE)
usethis::use_data(intg_pathinfo, overwrite = TRUE)










