#' MAQC2 deseq2 analysis result
#'
#' Human gene ID dataframe containing Entrez,
#' Used for ent2sym and more gene mapper functions.
#'
#' @docType data
#' @usage data(maqc_deseq2)
#' @format A dataframe of 16851 rows and 10 columns:
#' \describe{
#'   \item{entrez}{NCBI Entrez ID, Used in ent2sym}
#' 	 \item{hgnc_id}{HGNC ID}
#' 	 \item{hgnc_gene}{HGNC gene name}
#' 	 \item{hgnc_symbol}{HGNC approved gene symbol, Used in ent2sym}
#' 	 \item{ensembl}{Ensembl ID}
#' 	 \item{gene_type}{gene type}
#' }
#'
#' @source \href{https://www.gencodegenes.org/human/}{Gencode Metadata}
#' @examples
#' data(maqc_deseq2)
"maqc_deseq2"

