# PIS: Pathway Impact Score

PIS is an R package for determining differential expression signals in transcriptome analysis. It identifies differentially expressed gene (DEG) which optimizes pathway-level interpretation and facilitate biological insights, suggesting an alternative approach to uniformly applying conventional significance thresholds.

## Features
- Calculate pathway impact scores from differential expression results
- Flexible support for DESeq2, edgeR, and other result formats
- Binning and ranking of genes by statistical metrics

## Installation

```r
# Install from GitHub
devtools::install_github('LittleHeronCodes/PIS')
library(PIS)
```

## Quick Start Example

```r
# Example workflow
library(PIS)

# Assume you have a DESeq2/edgeR result data.frame 'resultDF' with columns 'stat' and 'entGene'
genes_rank <- createGeneRank(resultDF, value.col = "stat", names.col = "entGene")
gene_bins <- binGenesByCntCutoff(genes_rank, max_deg_count = 2000, bin_size = 10)

# Load or define your pathway gene sets as a list (e.g., from a GMT file)
# pathways <- readGMT("pathways.gmt")

# Define gene universe
gspace <- unique(resultDF$entGene)

# Calculate pathway scores
scoresMat <- calculatePathwayScores(gene_bins, gspace, pathways, ef_cut = 2)

# Summarize and get peak results
pisres <- getPeakResults2(gene_bins, scoresMat, verbose = TRUE)
print(pisres)
```

## Main Functions

- `createGeneRank(resultDF, value.col, names.col)`: Create a named vector of gene ranks from a result data frame.
- `binGenesByCntCutoff(genesrank, max_deg_count, bin_size, reverse)`: Bin genes by rank for thresholding.
- `calculatePathwayScores(glist, gspace, ref_geneset, ef_cut, ...)`: Calculate pathway scores for each gene bin.
- `getPeakResults2(geneCntList, scoresMat, verbose)`: Summarize results and identify the optimal threshold.
- `readGMT(path)`: Read gene sets from a GMT file.

## Documentation
See the `man/` directory or use `?function_name` in R for detailed documentation of each function.

## Example Data
Example datasets and pathway files are provided in the `data/` and `data-raw/` directories.

## License
MIT License. See `LICENSE` file for details.

## Author
Yeogha Yoon (<ygyoon17@ewhain.net>)  
[ORCID: 0000-0003-0955-4344](https://orcid.org/0000-0003-0955-4344)

## Citation
If you use PIS in your research, please cite the package and contact the maintainer for citation details.

