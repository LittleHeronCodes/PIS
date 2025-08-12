## Some gene list related functions from Lazy2

#' Get DEGs by cut-offs
#'
#' Extract the list of DEGs by cut-offs from the differential analysis results.
#' @param resultDF result data frame. Column names should include "entGene", "adj.P.Val", "logFC"
#' @param fcos sequence of fold changes to use as cut-offs
#' @param qcos sequence of q-values (FDR adjusted p-value) to use as cut-offs
#' @param colname.qv name for the column in resultDF containing q-values
#' @param colname.lfc name for the column in resultDF containing log2 Fold changes
#' @param colname.gene name for the column in resultDF containing the gene ID (eg. entrez, ensemble, symbol)
#' @return List of gene IDs selected by cutoffs
#' @examples
#' \dontrun{
#' data(maqc_deseq2)
#' fcos <- seq(1.5,4.0,0.1)
#' qcos <- seq(0.01,0.2, 0.01)
#' geneList.conv <- getGenesByCutoffs(maqc_deseq2, fcos, qcos)
#' # head(lapply(geneList.conv$up, head))
#' }
#' 
#' @export

getGenesByCutoffs <- function(
	resultDF, fcos, qcos, 
	colname.qv = "adj.P.Val", 
	colname.lfc = "logFC", 
	colname.gene = "entGene"
) {
	cutoff_grid <- expand.grid(fc = fcos, qv = qcos)
	cutoff_idx <- apply(cutoff_grid, 1, function(v) {
			paste0("q", sprintf("%.2f", v[2]), "_fc", sprintf("%.1f", v[1]))
		})

	## resultDF cleanup ##
	resultDF <- resultDF |> 
		dplyr::rename_at(c(colname.qv, colname.lfc, colname.gene), ~ c("adj.P.Val", "logFC", "geneID"))

	resultLS <- rep(list(resultDF), length(cutoff_idx))
	names(resultLS) <- cutoff_idx
	fcov <- structure(rep(fcos, length(qcos)), names = cutoff_idx)
	qcov <- structure(rep(qcos, each = length(fcos)), names = cutoff_idx)

	# extract gene list
	geneList <- list(up = list(), dn = list(), to = list())
	for (aid in names(resultLS)) {
		resultDF.f <- resultLS[[aid]]
		resultDF.f <- resultDF.f[which(!is.na(resultDF.f$geneID)),]
		
		up_genes <- with(resultDF.f, unique(geneID[which(adj.P.Val < qcov[aid] & logFC >=  log2(fcov[aid]))]))
		dn_genes <- with(resultDF.f, unique(geneID[which(adj.P.Val < qcov[aid] & logFC <= -log2(fcov[aid]))]))

		geneList$up[[aid]] <- up_genes
		geneList$dn[[aid]] <- dn_genes
		geneList$to[[aid]] <- unique(resultDF.f$geneID)
	}
	return(geneList)
}


#' Count number of genes in geneList
#'
#' Lazy function for gene number for geneList.
#' @param geneList Nested DEG list. See example for gene list structure.
#' @examples
#' set.seed(1234)
#' geneList <- list(
#'     up = list(A = sample(letters, 10), B = sample(letters, 5), C = sample(letters, 4)),
#'     dn = list(A = sample(letters, 6), B = sample(letters, 15), C = sample(letters, 7)),
#'     to = list(A = letters, B = letters, C = letters)
#' )
#' geneCount(geneList)
#' @export

geneCount <- function(geneList) {
	sapply(geneList, lengths)
}


#' readGMT
#'
#' Read/write GMT file. A GMT file format is a tab delimited text file containing gene sets.
#' Each line in the GMT file should contain one gene set or pathway, delimited by tab. Genes should start from the 3rd field.
#'
#' @param gmtfile GMT file path
#' @param as.df Return as data frame?
#' @param genelist List of gene set. Should be un-nested level one named list.
#' @param geneset_desc Description meta information for gene set. Either length one or a named vector of same length as glist.
#' @return Gene set dataframe of 2 column or a list of gene sets.
#' @examples
#' \dontrun{
#' readGMT("file/path.gmt")
#' }
#' @export

readGMT <- function(gmtfile, as.df = FALSE) {
	x <- readLines(gmtfile)
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	out <- lapply(res, "[", -c(1:2))
	if (as.df) {
		ont2gene <- utils::stack(out)
		ont2gene <- ont2gene[, c("ind", "values")]
		colnames(ont2gene) <- c("ont", "gene")
		out <- ont2gene
	}
	return(out)
}


#' @describeIn readGMT write gmt file for geneset list.
#' @export

writeGMT <- function(gmtfile, genelist, geneset_desc = "") {
	if (!(is.list(genelist) & all(sapply(genelist, is.character)))) {
		stop("Check genelist format. genelist should be a one-level list of genesets.")
	}
	if (length(names(geneset_desc)) != 0) {
		geneset_desc <- geneset_desc[names(genelist)]
	}

	concat <- sapply(genelist, function(v) paste(v, collapse = "\t"))
	out <- paste0(names(concat), "\t", geneset_desc, "\t", concat)

	writeLines(out, con = gmtfile)
}
