## Pathway Impact Analysis Functions

## =======================================================================================
## These functions are for internal use only

#' Pathway score
#'
#' Calculate pathway scores from enrichment factors (used inside calculatePathwayScores)
#' @param efs vector of enrichment factors
#' @param ef_cut minimum enrichment factor cut-offs. Default 2
#' @return vector of pathway enrichment score
#' @export

pathScores2 <- function(efs, ef_cut = 2) {
    scores <- ifelse(efs >= ef_cut, log2(efs), 0)
    return(scores)
}


#' Hypergeometric Test for Gene set
#'
#' Run hypergeometric test for gene set. (used inside calculatePathwayScores)
#' This function is meant to be used inside calculatePathwayScores, 
#' but can be used separately for simple overrepresentation analysis as well.
#'
#' @param query  gene set to query (eg. Differentially Expressed Genes).
#' @param refGMT list of reference gene set (eg. Pathways).
#' @param gspace background gene space. Should contain all genes in query.
#' @param minGeneSet minimum size of gene set. This is used to filter refGMT. Default 10.
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 1.
#' @param verbose print number of filtered entries in refGMT. Default FALSE.
#' @return dataframe of results
#' @export

hypergeoTestForGeneset.simple <- function(query, refGMT, gspace, minGeneSet = 10, ef.psc = 1, verbose = FALSE) {
    # match gene space
    if (!all(query %in% gspace)) {
        # stop(paste(length(setdiff(query, gspace)),'query items were found outside of background space. Check inputs.'))
        query <- intersect(query, gspace)
    }
    if (length(query) == 0) {
        stop("query should be character vector of at least length 1.")
    }

    if (!all(unlist(refGMT) %in% gspace)) {
        refGMT <- lapply(refGMT, function(g) intersect(g, gspace))
    }

    # filter refGMT with less than minimum gene set
    exc <- which(sapply(refGMT, length) < minGeneSet)
    if (length(exc) != 0) {
        refGMT <- refGMT[which(sapply(refGMT, length) >= minGeneSet)]
    }
    if (length(refGMT) == 0) {
        stop("Length of refGMT after filtering for minGeneSet is zero. Set lower minGeneSet or check gene inputs.")
    }

    # hypergeometric test
    N <- length(gspace) # no of balls in urn
    k <- length(query) # no of balls drawn from urn (DEG no)

    qs <- sapply(refGMT, function(x) length(intersect(x, query))) # no of white balls drawn
    ms <- sapply(refGMT, length) # no of white balls in urn

    pvals <- phyper(qs - 1, ms, N - ms, k, lower.tail = FALSE)
    odds <- (qs + ef.psc) / (ms / N * k + ef.psc)
    jacc <- qs / sapply(refGMT, function(x) length(union(x, query)))
    gs.ratio <- paste0(qs, "/", k)
    bg.ratio <- paste0(ms, "/", N)
    enrRes <- data.table(
        ID = names(refGMT), pVal = pvals, oddsRatio = odds, tan = jacc,
        int = qs, gsRatio = gs.ratio, bgRatio = bg.ratio
    )
	
    # Adjust p-value while filtering out 1 values
    pv <- ifelse(enrRes$int == 0, NA, enrRes$pVal)
    enrRes$qVal <- p.adjust(pv, method = "fdr")
    enrRes$qVal <- ifelse(enrRes$int == 0, 1, enrRes$qVal)

	enrRes <- enrRes[,c('ID','pVal','qVal','oddsRatio','tan','int','gsRatio','bgRatio')]
    return(enrRes)
}


### =======================================================================================
## PIS Functions

#' PIS calculation
#'
#' Calculate pathway score for gene list (uses pathScores2 function)
#' Outputs matrix of pathway scores where each row is gene set and column is bin.
#' @param glist list of gene set (result from binGenesByCntCutoff)
#' @param gspace background gene space. (for hypergeoTestForGeneset.simple)
#' @param ref_geneset list of reference gene set (eg. Pathways)
#' @param ef_cut EF cut-off for scoring.
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 0
#' @param minGeneSet minimum size of gene set. This is used to filter refGMT. Default 10.
#' @param ncore no of threads to use in mclapply. (ncore=1 uses lapply)
#' @param verbose print run time. Default FALSE
#' @param ... futher arguments to be passed to the internal hypergeoTestForGeneset.simple function
#' @return Matrix of pathway scores
#' @examples
#' \dontrun{
#' data(maqc_deseq2)
#' fcos <- seq(1.5,4.0,0.1)
#' qcos <- seq(0.01,0.2, 0.01)
#' gspace <- unique(maqc_deseq2$entGene)
#' pathways <- intg_pathways[1:1000]
#' geneList.conv <- getGenesByCutoffs(maqc_deseq2, fcos, qcos)
#' scoresMat <- calculatePathwayScores(geneList.conv$up, gspace, pathways, 2)
#' pisres.up <- getPeakResults2(geneList.conv$up, scoresMat, verbose=TRUE)
#' pisres.up
#' }
#' 
#' @importFrom parallel mclapply
#' @export
calculatePathwayScores <- function(glist, gspace, ref_geneset, 
								   ef_cut = 2, ef.psc = 1, 
								   minGeneSet = 10, 
								   ncore = 1, 
								   verbose = TRUE, 
								   ...) {
    # if( any(c(ef_cut, min.overlap) <= 0 ) ) stop('ef_cut and min.overlap should both be greater than 0.')

    # match gene space
    if (any(!unlist(ref_geneset) %in% gspace)) {
        ref_geneset <- lapply(ref_geneset, function(g) intersect(g, gspace))
    }
	
    if (any(!unlist(glist) %in% gspace)) {
        glist <- lapply(glist, function(g) intersect(g, gspace))
    }

    # filter pathway sets with less than minimum gene set
    exc <- which(sapply(ref_geneset, length) < minGeneSet)
    if (length(exc) != 0) {
        ref_geneset <- ref_geneset[which(sapply(ref_geneset, length) >= minGeneSet)]
    }
	
    if (length(ref_geneset) == 0) {
        stop("Length of ref_geneset after filtering for minGeneSet is zero. Set lower minGeneSet or check gene inputs.")
    }

    # Calculate Pathway scores
    tcheck <- proc.time()
    if (ncore > 1) {
        scoresLS <- mclapply(glist, function(gset) {
            if (length(gset) == 0) {
                out <- structure(rep(0, length(ref_geneset)), names = names(ref_geneset))
            } else {
                enres <- hypergeoTestForGeneset.simple(gset, ref_geneset, gspace = gspace, ef.psc = ef.psc, ...)
                # enres$oddsRatio[which(enres$int <= min.overlap)] <- 0
                enres$score <- pathScores2(efs = enres$oddsRatio, ef_cut = ef_cut)
                out <- structure(enres$score, names = enres$ID)
                out <- out[names(ref_geneset)]
            }
            return(out)
        }, mc.cores = ncore)
    } else {
        scoresLS <- lapply(glist, function(gset) {
            if (length(gset) == 0) {
                out <- structure(rep(0, length(ref_geneset)), names = names(ref_geneset))
            } else {
                enres <- hypergeoTestForGeneset.simple(gset, ref_geneset, gspace = gspace, ef.psc = ef.psc, ...)
                # enres$oddsRatio[which(enres$int <= min.overlap)] <- 0
                enres$score <- pathScores2(efs = enres$oddsRatio, ef_cut = ef_cut)
                out <- structure(enres$score, names = enres$ID)
                out <- out[names(ref_geneset)]
            }
            return(out)
        })
    }
    scoresMat <- do.call(cbind, scoresLS)
    if (verbose) print((proc.time() - tcheck) / 60) # ~1min for glist of length 1000
    return(scoresMat)
}

# calculatePS_internal <- function(gset, gspace,ref_geneset, ef.psc, ef_cut) {
# 	if(length(gset) == 0) {
# 		out <- structure(rep(0, length(ref_geneset)), names=names(ref_geneset))
# 	} else {
# 		enres <- hypergeoTestForGeneset.simple(gset, ref_geneset, gspace = gspace, ef.psc = ef.psc, ...)
# 		enres$score <- pathScores2(efs=enres$oddsRatio, setsizes=enres$int, ef_cut=ef_cut)
# 		# enres$score <- pathScores2(efs=enres$oddsRatio, setsizes=geneset_size[enres$ID], ef_cut=ef_cut)
# 		out <- structure(enres$score, names=enres$ID)
# 	}
# 	return(out)
# }


#' Get Peak Result from calculatePathwayScores
#'
#' Summarize result from calculatePathwayScores and get peak cut-off.
#' @param geneCntList list of genes binned by count (output from binGenesByCntCutoff)
#' @param scoresMat Pathway score matrix (output from calculatePathwayScores)
#' @param verbose print peak scores. Default FALSE
#' @return PIS object list
#' @examples
#' \dontrun{
#' data(maqc_deseq2)
#' fcos <- seq(1.5,4.0,0.1)
#' qcos <- seq(0.01,0.2, 0.01)
#' gspace <- unique(maqc_deseq2$entGene)
#' geneList.conv <- getGenesByCutoffs(maqc_deseq2, fcos, qcos)
#' scoresMat <- calculatePathwayScores(geneList.conv$up, gspace, intg_pathways, 2)
#' pisres.up <- getPeakResults2(geneList.conv$up, scoresMat, verbose=TRUE)
#' pisres.up
#' }
#' 
#' @export

getPeakResults2 <- function(geneCntList, scoresMat, verbose = FALSE) {
    # PIS for each gene count cut-off
	gcntSum <- apply(scoresMat, 2, sum, na.rm = TRUE)

    # gene count cut off at peak
    peak_cnt <- sapply(geneCntList, length)[which.max(abs(gcntSum))]

    # peak_score
    peak_score <- gcntSum[which.max(abs(gcntSum))]

    # gene set in peak
    peak_gset <- geneCntList[[names(peak_cnt)]]

    # scored pathways in peak cut-off
    vv <- sort(abs(scoresMat[, names(peak_cnt)]), decreasing = TRUE)
    scored_pathways <- vv[which(vv != 0)]

    # gene count for each cut-off
    genecnt_cut <- sapply(geneCntList, length)

    # scores by bin
    bin_scores <- data.table(bin = names(genecnt_cut), genecnt = genecnt_cut, bin_scores = gcntSum)

    # PIS result into PISobject
    peakObj <- list(
        peak_cnt = peak_cnt,
        peak_score = peak_score, peak_pathwayCnt = length(scored_pathways),
        peak_gset = peak_gset, scored_pathways = scored_pathways,
        bin_scores = bin_scores
    )

    class(peakObj) <- "PISobj"

    mesg <- paste("peak score :", round(peak_score), "\npeak pathway counts :", peakObj$peak_pathwayCnt)
    if (verbose) message(mesg)

    return(peakObj)
}
