## Some gene list related functions from Lazy2

#' Get DEGs by cut-offs
#'
#' Extract the list of DEGs by cut-offs from the differential analysis results.
#' @param resultDF result data frame. Column names need to include "entGene", "adj.P.Val", "logFC"
#' @param fcos fold change
#' @param qcos q-value
#' @return List of entrez IDs
#' @export

getGenesByCutoffs <- function(resultDF, fcos, qcos) {
    cutoff_idx <- apply(expand.grid(sprintf("fc%.1f", fcos), sprintf("q%.2f", qcos)), 1, function(v) paste0(v[2], "_", v[1]))

    ## resultDF cleanup ##
    ##

    resultLS <- rep(list(resultDF), length(cutoff_idx))
    names(resultLS) <- cutoff_idx
    fcov <- structure(rep(fcos, length(qcos)), names = cutoff_idx)
    qcov <- structure(rep(qcos, each = length(fcos)), names = cutoff_idx)

    # extract gene list (from Lazy2)
    geneList <- list(up = list(), dn = list(), to = list())
    for (aid in names(resultLS)) {
        resultDF.f <- resultLS[[aid]] %>% filter(!is.na(entGene))

        geneList$up[[aid]] <- with(resultDF.f, unique(entGene[which(adj.P.Val < qcov[aid] & logFC >= log2(fcov[aid]))]))
        geneList$dn[[aid]] <- with(resultDF.f, unique(entGene[which(adj.P.Val < qcov[aid] & logFC <= -log2(fcov[aid]))]))
        geneList$to[[aid]] <- unique(resultDF.f$entGene)
    }
    return(geneList)
}


#' Count number of genes in geneList
#' 
#' Lazy function for gene number for geneList. 
#' @param geneList Nested DEG list. See example for gene list structure. 
#' @export
#' @examples
#' set.seed(1234)
#' geneList <- list(
#' 	up = list(A=sample(letters,10), B=sample(letters, 5), C=sample(letters, 4)),
#' 	dn = list(A=sample(letters, 6), B=sample(letters,15), C=sample(letters, 7)),
#' 	to = list(A=letters, B=letters, C=letters)
#' 	)
#' geneCount(geneList)

geneCount <- function(geneList) { sapply(geneList, lengths) }


#' readGMT
#' 
#' GMT file reader. A GMT file format is a tab delimited text file containing gene sets.
#' Each line in the GMT file should contain one geneset or pathway, delimited by tab. Genes should start from the 3rd field.
#' 
#' @param gmtfile GMT file path
#' @param as.df Return as data frame?
#' @param genelist List of gene set. Should be un-nested level one named list.
#' @param geneset_desc Description meta information for gene set. Either length one or a named vector of same length as glist.
#' @return Gene set dataframe of 2 column or list
#' @export
#' @examples
#' \dontrun{
#' readGMT('file/path.gmt')
#' }

readGMT <- function(gmtfile, as.df=FALSE) {
	x <- readLines(gmtfile)
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	out <- lapply(res, "[", -c(1:2))
	if(as.df) {
		ont2gene <- stack(out)
		ont2gene <- ont2gene[, c("ind", "values")]
		colnames(ont2gene) <- c("ont", "gene")
		out <- ont2gene
	}
	return(out)
}


#' @describeIn readGMT write gmt file for geneset list.
#' @export

writeGMT <- function(gmtfile, genelist, geneset_desc="") {
	if( !(is.list(genelist) & all(sapply(genelist, is.character))) ) {
		stop("Check genelist format. genelist should be a one-level list of genesets.")
	}
	if(length(names(geneset_desc)) != 0) {
		geneset_desc <- geneset_desc[names(genelist)]
	}

	concat <- sapply(genelist, function(v) paste(v, collapse="\t"))
	out <- paste0(names(concat), "\t",geneset_desc,"\t", concat)

	writeLines(out, con=gmtfile)
}
