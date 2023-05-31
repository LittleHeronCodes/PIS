## Pathway Impact Analysis Functions

## =======================================================================================
## These functions are for internal use only

#' Pathway score
#'
#' @param efs vector of enrichment factors
#' @param ef_cut minimum enrichment factor cut-offs. Default 2
#' @return vector of pathway scores
#' @export

pathway_scores <- function(efs, ef_cut = 2) {
    scores <- ifelse((efs >= ef_cut) & (efs != 0), log2(efs), 0)
    return(scores)
}


#' Hypergeometric test for gmt style list
#'
#' Gene list overrepresentation analysis against list of gene set. 
#' Runs both hypergeometric test or calculate enrichment factor for gene set. (used inside calculatePathwayScores)
#' 
#'
#' @param query  gene set to query (eg. Differentially Expressed Genes)
#' @param refGMT list of reference gene set (eg. Pathways)
#' @param gspace background gene space. Should contain all genes in query.
#' @param min_geneset minimum size of gene set used to filter refGMT. Default 10
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 1
#' @param verbose print number of filtered entries in refGMT. Default FALSE
#' @return Data frame of gene set analysis results
#' \describe{
#' 	  \item{ID}{: pathway description or ID}
#' 	  \item{pVal}{: hypergeometric test p values from phyper}
#' 	  \item{qVal}{: FDR adjusted p-values}
#' 	  \item{oddsRatio}{: odds ratio}
#' 	  \item{tan}{: tanimoto coefficient (Jaccard index)}
#' 	  \item{int}{: intersected item count}
#' 	  \item{gsRatio}{: gene set ratio (selected genes in gene set / selected genes)}
#' 	  \item{bgRatio}{: background ratio (total genes in gene set / total gene space)}
#' }
#' @return dataframe of results
#' @export

hypergeo_test_geneset <- function(query, refGMT, gspace, min_geneset = 10, ef.psc = 1, verbose = FALSE) {
    # match gene space
    if (!all(query %in% gspace)) {
        query <- intersect(query, gspace)
    }
 
    if (length(query) == 0) {
        stop(
            "query after filtering for gene space is length zero.
			Make sure that genes provided in gspace and query matches in data type and ID (entrez, symbol, etc.).
			Ideally, gspace should contain all query genes."
        )
    }

    if (!all(unlist(refGMT) %in% gspace)) {
        refGMT <- lapply(refGMT, function(g) intersect(g, gspace))
    }

    # filter refGMT with less than minimum gene set
    exc <- which(sapply(refGMT, length) < min_geneset)
    if (length(exc) != 0) {
        refGMT <- refGMT[which(sapply(refGMT, length) >= min_geneset)]
    }
    if (length(refGMT) == 0) {
        stop(
			"Length of refGMT after filtering for min_geneset and gspace is zero. 
			Set lower min_geneset or make sure genes in refGMT and gspace matches in type."
		)
    }

    # overrepresentation test (hypergeomteric, enrichement factor)
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

	enrRes <- enRes[,c('ID','pVal','qVal','oddsRatio','tan','int','gsRatio','bgRatio')]
    return(enrRes)
}



### =======================================================================================
## PIS Functions

#' PIS calculation
#'
#' Calculate pathway score for gene list (uses pathway_scores function)
#' Outputs matrix of pathway scores where each row is gene set and column is bin.
#' @param genelist list of gene set (result from binGenesByCntCutoff)
#' @param gspace background gene space. (for hypergeo_test_geneset)
#' @param ref_geneset list of reference gene set (eg. Pathways)
#' @param ef_cut EF cut-off for scoring.
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 0
#' @param ncore no of threads to use in mclapply. (ncore=1 uses lapply)
#' @param verbose print run time. Default FALSE
#' @param ... futher arguments to be passed to hypergeo_test_geneset
#' @return Matrix of pathway scores
#' @importFrom parallel mclapply
#' @export

calculatePathwayScores <- function(genelist, gspace, ref_geneset, ef_cut = 2, ef.psc = 1, ncore = 4, verbose = TRUE, ...) {

    # match gene space
    if (any(!unlist(ref_geneset) %in% gspace)) {
        ref_geneset <- lapply(ref_geneset, function(g) intersect(g, gspace))
    }
    if (any(!unlist(genelist) %in% gspace)) {
        genelist <- lapply(genelist, function(g) intersect(g, gspace))
    }

    # Calculate Pathway scores
    tcheck <- proc.time()
    if (ncore > 1) {
        enresLS <- mclapply(genelist, function(gset) {
            if (length(gset) != 0) {
                enres <- hypergeo_test_geneset(gset, ref_geneset, gspace = gspace, ef.psc = ef.psc, ...)
                enres$score <- pathway_scores(efs = enres$oddsRatio, ef_cut = ef_cut)
            } else {
				enrRes <- data.table(
					ID = names(ref_geneset), pVal = NA, qVal = NA, oddsRatio = 0, tan = 0, int = 0, gsRatio = "", bgRatio = "")
                # out <- structure(rep(0, length(ref_geneset)), names = names(ref_geneset))
            }
            return(enrRes)
        }, mc.cores = ncore)
    } else {
        enresLS <- lapply(genelist, function(gset) {
            if (length(gset) != 0) {
                enres <- hypergeo_test_geneset(gset, ref_geneset, gspace = gspace, ef.psc = ef.psc, ...)
                enres$score <- pathway_scores(efs = enres$oddsRatio, ef_cut = ef_cut)
            } else {
				enrRes <- data.table(
					ID = names(ref_geneset), pVal = NA, qVal = NA, oddsRatio = 0, tan = 0, int = 0, gsRatio = "", bgRatio = "")
                # out <- structure(rep(0, length(ref_geneset)), names = names(ref_geneset))
            }
            return(enrRes)
        })
    }
	
	# TODO: return all the results instead of just the pathway scores
	scoresLS <- lapply(enresLS, function(enres) {
	    out <- structure(enres$score, names = enres$ID)
        out <- out[names(ref_geneset)]
	})
    scores_mat <- do.call(cbind, scoresLS)


    # runtime check
    if (verbose) print((proc.time() - tcheck) / 60)
    return(scores_mat)
}

# calculatePS_internal <- function(gset, gspace,ref_geneset, ef.psc, ef_cut) {
# 	if(length(gset) == 0) {
# 		out <- structure(rep(0, length(ref_geneset)), names=names(ref_geneset))
# 	} else {
# 		enres <- hypergeo_test_geneset(gset, ref_geneset, gspace = gspace, ef.psc = ef.psc, ...)
# 		enres$score <- pathway_scores(efs=enres$oddsRatio, setsizes=enres$int, ef_cut=ef_cut)
# 		# enres$score <- pathway_scores(efs=enres$oddsRatio, setsizes=geneset_size[enres$ID], ef_cut=ef_cut)
# 		out <- structure(enres$score, names=enres$ID)
# 	}
# 	return(out)
# }


#' Get Peak Result object from calculatePathwayScores
#'
#' Summarize result from calculatePathwayScores and get peak cut-off.
#' @param genelist list of genes binned by cut-off.
#' @param scores_mat Pathway score matrix (output from calculatePathwayScores)
#' @param verbose print peak scores. Default FALSE
#' @return PISobj with peak results.
#' \describe{
#' 	  \item{peak_cnt}{: No. of DEGs selected}
#' 	  \item{peak_score}{: max PIS score}
#' 	  \item{peak_pathwayCnt}{: No. of pathways selected at peak}
#' 	  \item{peak_gset}{: DEGs selected at peak}
#' 	  \item{scored_pathways}{: Pathways and scores selected by PIS}
#' 	  \item{bin_scores}{: scores by each cut-off)}
#' }
#' @export

getPeakResults2 <- function(genelist, scores_mat, verbose = FALSE) {
    # PIS for each gene count cut-off
    gcntSum <- apply(scores_mat, 2, sum)

    # gene count cut off at peak
    peak_cnt <- sapply(genelist, length)[which.max(abs(gcntSum))]

    # peak_score
    peak_score <- gcntSum[which.max(abs(gcntSum))]

    # gene set in peak
    peak_gset <- genelist[[names(peak_cnt)]]

    # scored pathways in peak cut-off
    vv <- sort(abs(scores_mat[, names(peak_cnt)]), decreasing = TRUE)
    scored_pathways <- vv[which(vv != 0)]

    # gene count for each cut-off
    genecnt_cut <- sapply(genelist, length)

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
