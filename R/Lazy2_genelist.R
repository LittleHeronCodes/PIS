## Some gene list related functions from Lazy2

#' Count number of genes in geneList
#' 
#' Lazy function for gene number for geneList
#' 
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

geneCount <- function(geneList) { sapply(geneList, function(ls) sapply(ls, length)) }

#' readGMT
#' 
#' gmt file reader. A GMT file format is a tab delimited text file containing gene sets.
#' Each line should contain one geneset, delimited by tab. Genes should start from 3rd field.
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

writeGMT <- function(gmtfile, genelist, geneset_desc='') {
	if( !(is.list(genelist) & all(sapply(genelist, is.character))) ) {
		stop('Check genelist format. genelist should be a one-level list of genesets.')
	}
	if(length(names(geneset_desc)) != 0) {
		geneset_desc <- geneset_desc[names(genelist)]
	}

	concat <- sapply(genelist, function(v) paste(v, collapse='\t'))
	out <- paste0(names(concat), '\t',geneset_desc,'\t', concat)

	writeLines(out, con=gmtfile)
}
