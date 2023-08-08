#' MAQC2 DESeq2 analysis result
#'
#' Data table of differential expression analysis results using DESeq2 for the Microarray Quality Control 2 (MAQC2) data.
#' Dataset downloaded from NCBI Sequence Read Archive (accession: SRP001847)
#'
#' @docType data
#' @usage data(maqc_deseq2)
#' @format A data table of 16851 rows and 10 columns:
#' \describe{
#'   \item{baseMean}{average expression of each gene across all samples}
#' 	 \item{logFC}{log2 fold change of each gene}
#' 	 \item{lfcSE}{standard error of logFC, calculated by DESeq2}
#' 	 \item{stat}{Wald statistic calculated by DESeq2}
#' 	 \item{P.Value}{raw p-value}
#' 	 \item{adj.P.Val}{FDR adjusted p-value, or q-value}
#' 	 \item{ensGene}{Ensembl ID}
#' 	 \item{entGene}{NCBI Entrez ID}
#' 	 \item{geneSym}{HGNC approved gene symbol}
#' 	 \item{AveExpr}{average log2-expression}
#' }
#' 
#' @source \href{https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP001847}{SRP001847}
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

