library(devtools)
library(data.table)
# library(tidyverse)
library(usethis)

## MAQC2 DESeq2 results

maqc_deseq2 <- fread("data-raw/MAQC2_DESeq2_result.tsv", quote = "")

usethis::use_data(maqc_deseq2, overwrite = TRUE)
