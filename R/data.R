#' MAQC2 deseq2 analysis result
#'
#' DESeq2 analysis result of MAQC2 dataset.
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


#' Integrated pathway set
#'
#' Integrated pathway set from KaiPharm. (sources from Reactome, GO)
#'
#' @docType data
#' @usage data(intg_pathways)
#' @format A list of gene sets (human entrez ID). Pathway descriptions in intg_pathinfo.
#'
#' @examples
#' data(intg_pathways)
"intg_pathways"


#' @describeIn intg_pathways
#' Pathway information metadata for intg_pathways
#' @usage data(intg_pathinfo)
#' @format A data.table describing pathways in 'intg_pathways'
#' @examples
#' data(intg_pathinfo)
"intg_pathinfo"

