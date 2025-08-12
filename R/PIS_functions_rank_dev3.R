#' Create Gene Rank
#'
#' Create gene rank scores from differential analysis result data frames
#' @param resultDF  result dataframe (DESeq2, edgeR)
#' @param value.col  column as character to use as value in gene rank vector (eg. logFC, stat)
#' @param names.col  column as character to use as names of gene rank vector (eg. entrez, ensembl)
#' @return vector of gene rank
#' @examples
#' \dontrun{
#' createGeneRank(resultDF, "stat", "entGene")
#' }
#' @export

createGeneRank <- function(resultDF, value.col, names.col) {
	resultDF <- resultDF[which(!is.na(resultDF[, value.col])), ]
	genesrank <- with(resultDF, structure(get(value.col), names = as.character(get(names.col))))
	genesrank <- sort(genesrank, decreasing = TRUE)
	return(genesrank)
}


#' Bin Genes by rank
#'
#' Create a list of genes for each binned count
#' @param genesrank named vector of genes rank (output from createGeneRank)
#' @param max_deg_count maximum count for DEG
#' @param bin_size gene rank bin size
#' @param reverse reverse order (for down regulated genes)
#' @return list of genes by rank cut-off
#' @export

binGenesByCntCutoff <- function(
	genesrank, max_deg_count = 2000, bin_size = 10, reverse = FALSE
) {
	if (reverse) genesrank <- rev(genesrank)

	# bin sequence index
	seqidx <- c(seq(1, max_deg_count, bin_size), max_deg_count + 1)
	if (is.null(max_deg_count)) {
		seqidx <- c(seq(1, length(genesrank), bin_size), length(genesrank) + 1)
	}

	# list of genes for each count
	geneCntList <- lapply(1:(length(seqidx) - 1), function(ix) {
		gset <- names(genesrank[1:(seqidx[(ix + 1)] - 1)])
	})
	names(geneCntList) <- paste0("cut_", 1:(length(seqidx) - 1))
	return(geneCntList)
}


#' Smooth scores
#'
#' Fit loess smoothing on bin scores
#' @param peakObj peak result object. (Output from getPeakResults2)
#' @param loess.span span option for loess
#' @return peak result object with smoothed scores attribute
#' @export

scoreSmooth <- function(peakObj, loess.span = 0.1) {
	bin_scores <- peakObj$bin_scores

	lw_colname <- paste0("loess_sp.", loess.span)
	fit <- loess(bin_scores ~ genecnt, data = bin_scores, span = loess.span)
	# bin_scores$loess.pred <- fit$fitted
	bin_scores[, lw_colname] <- fit$fitted
	# # library(mgcv)
	# # model <- gam(PIS ~ s(genecnt_cut, bs='cs'), data = plotdf, sp=0.1)
	# # plotdf$predict <- predict(model)

	# predict peak
	yy <- predict(fit, min(bin_scores$genecnt):max(bin_scores$genecnt))

	if (!"smoothed_peak" %in% names(peakObj)) {
		peakObj$smoothed_peak <- list()
	}
	smoothed_df <- data.frame(
		peak_score = max(yy), 
		peak_point = which.max(yy) + min(bin_scores$genecnt)
	)
	peakObj$smoothed_peak[[lw_colname]] <- smoothed_df
	peakObj$bin_scores <- bin_scores

	return(peakObj)
}


#' Random permutation
#'
#' (add description)
#' @param peakObj peak result object. (Output from getPeakResults2)
#' @param gspace total gene space
#' @param ref_geneset list of reference gene set (eg. Pathways)
#' @param ef_cut EF cut-off
#' @param min.overlap minimum overlap
#' @param iter iteration
#' @param use.smoothed use smoothed?
#' @param ncore number of cores
#' @return peak result object with smoothed scores attribute

peakSignifByRandPerm <- function(
	peakObj, gspace, ref_geneset, 
	ef_cut = 2, 
	min.overlap = 5, 
	iter = 1e4, 
	use.smoothed = FALSE, 
	ncore = 1
) {
	if (use.smoothed & ("smoothed_peak" %in% names(peakObj))) {
		smidx <- names(peakObj$smoothed_peak)[1]
		smObj <- peakObj$smoothed_peak[[smidx]]
		cnt <- smObj$peak_point
		score <- smObj$peak_score
	} else {
		if (use.smoothed) {
			message("use.smoothed is TRUE but there is no smoothed attributes. Using raw peaks.")
		}
		use.smoothed <- FALSE
		cnt <- peakObj$peak_cnt
		score <- peakObj$peak_score
	}

	randomls <- lapply(1:iter, function(i) sample(gspace, cnt))
	scoresran <- calculatePathwayScores(
		randomls, gspace, ref_geneset, ef_cut, min.overlap, ncore = ncore
	)
	ran_score_dist <- apply(scoresran, 2, sum)
	signif <- sum(ran_score_dist > score) / length(ran_score_dist)

	if (use.smoothed) {
		peakObj$smoothed_peak[[smidx]] <- cbind(smObj, signif = signif)
	} else {
		peakObj$signif <- signif
	}

	return(peakObj)
}


#' Draw peak plot
#'
#' Draw peak plot
#' @param peakObj peak result object. (Output from getPeakResults2)
#' @param mtitle main title for plot
#' @param lw_colname loess smoothed column name. (from names(peakObj$smoothed_peak))
#' @export

drawGeneCntCutoffPeak2 <- function(peakObj, mtitle = "", lw_colname = NULL) {
	plotdf <- peakObj$bin_scores
	peak_cnt <- peakObj$peak_cnt

	g1 <- ggplot(plotdf, aes(x = genecnt, y = bin_scores)) +
		geom_point() +
		labs(x = "DEG count", y = "PIS", title = mtitle) +
		theme_bw(base_size = 12)

	if (!is.null(lw_colname)) {
		peak_cnt <- peakObj$smoothed_peak[[lw_colname]]$peak_point
		g1 <- g1 +
			geom_line(
				aes(y = get(lw_colname)), 
				data = plotdf, size = 1.5, linetype = 2, col = "red"
			) +
			geom_vline(xintercept = peak_cnt, col = "blue", lty = 2)
	} else {
		g1 <- g1 + geom_vline(xintercept = peak_cnt, col = "blue", lty = 2)
	}
	g1
}

