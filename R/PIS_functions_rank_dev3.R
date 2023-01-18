#' Create Gene Rank 
#'
#' Create gene rank scores from differential analysis result data frames
#' @param resultDF  result dataframe (DESeq2, edgeR)
#' @param value.col  column as character to use as value in gene rank vector (eg. logFC, stat)
#' @param names.col  column as character to use as names of gene rank vector (eg. entrez, ensembl)
#' @return vector of gene rank
#' @export
#' @examples
#' \dontrun{
#' createGeneRank(resultDF, 'stat', 'entGene')
#' }

createGeneRank <- function(resultDF, value.col, names.col) {
	resultDF <- resultDF[which(!is.na(resultDF[,value.col])),]
	genesrank <- with(resultDF3, structure(get(value.col), names=as.character(get(names.col))))
	genesrank <- sort(genesrank, decreasing=TRUE)
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

binGenesByCntCutoff <- function(genesrank, max_deg_count=2000, bin_size=10, reverse=FALSE) {

	if(reverse) genesrank <- rev(genesrank)

	# bin sequence index
	seqidx <- c(seq(1, max_deg_count, bin_size), max_deg_count+1)
	if(is.null(max_deg_count)) seqidx <- c(seq(1,length(genesrank),bin_size), length(genesrank)+1)

	# list of genes for each count
	geneCntList <- lapply(1:(length(seqidx)-1), function(ix) {
		gset = names(genesrank[1:(seqidx[(ix+1)]-1)])
		})
	names(geneCntList) <- paste0('cut_',1:(length(seqidx)-1))
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

	lw_colname <- paste0('loess_sp.',loess.span)
	fit <- loess(bin_scores ~ genecnt, data=bin_scores, span=loess.span)
	# bin_scores$loess.pred <- fit$fitted
	bin_scores[,lw_colname] <- fit$fitted
	# # library(mgcv)
	# # model <- gam(PIS ~ s(genecnt_cut, bs='cs'), data = plotdf, sp=0.1)
	# # plotdf$predict <- predict(model)

	# predict peak
	yy <- predict(fit, min(bin_scores$genecnt):max(bin_scores$genecnt))

	if(!'smoothed_peak' %in% names(peakObj)) peakObj$smoothed_peak <- list()
	peakObj$smoothed_peak[[lw_colname]] <- data.frame(
		peak_score = max(yy), peak_point = which.max(yy) + min(bin_scores$genecnt))

	peakObj$bin_scores <- bin_scores

	return(peakObj)
}


#' Random permutation
#'
#' (add description)
#' @param peakObj peak result object. (Output from getPeakResults2)
#' @param gspace total gene space
#' @param ef_cut EF cut-off
#' @param min.overlap minimum overlap
#' @param iter iteration
#' @param use.smoothed use smoothed?
#' @param ncore number of cores
#' @return peak result object with smoothed scores attribute

peakSignifByRandPerm <- function(peakObj, gspace, ef_cut=2, min.overlap=5, iter=1e4, use.smoothed=FALSE, ncore=1) {
	if(use.smoothed & ('smoothed_peak' %in% names(peakObj))) {
		smidx <- names(peakObj$smoothed_peak)[1]
		smObj <- peakObj$smoothed_peak[[smidx]]
		cnt <- smObj$peak_point
		score <- smObj$peak_score
	} else {
		if(use.smoothed) message('use.smoothed is TRUE but there is no smoothed attributes. Using raw peaks.')
		use.smoothed <- FALSE
		cnt <- peakObj$peak_cnt
		score <- peakObj$peak_score
	}

	randomls <- lapply(1:iter, function(i) sample(gspace, cnt))
	scoresran <- calculatePathwayScores(randomls, gspace, pathways, ef_cut, min.overlap, ncore=ncore) # ~ 1 min
	ran_score_dist <- apply(scoresran, 2, sum)
	signif <- sum(ran_score_dist > score)/length(ran_score_dist)
	
	if(use.smoothed) {
		peakObj$smoothed_peak[[smidx]] <- cbind(smObj, signif = signif)
	} else {
		peakObj$signif <- signif
	}

	return(peakObj)
}



##


# drawGeneCntCutoffPeak <- function(peakObj, mtitle='', logScaleX=FALSE) {
# 	peak_cnt <- peakObj$peak_cnt
# 	peak_score <- peakObj$peak_score
# 	gcntSum <- peakObj$gcntSum
# 	peak_pathwayCnt <- peakObj$peak_pathwayCnt
# 	genecnt_cut <- peakObj$genecnt_cut

# 	tlabel <- paste('peak score :', round(peak_score),'\npeak pathway counts :', peak_pathwayCnt)
# 	# ylim <- sort(c(0,gcntSum[which.max(abs(gcntSum))]))
# 	ylim <- range(c(0, gcntSum))

# 	if(logScaleX) {
# 		plot(log10(genecnt_cut),gcntSum, pch=20, cex=.7, main=mtitle, xlab='DEG count (log10)', ylab='PIS',ylim=ylim)
# 		lines(log10(genecnt_cut),gcntSum, cex=.7, col='grey25')
# 		abline(v=log10(peak_cnt), col='blue', lty=2)
# 		text(x=log10(peak_cnt),y=10, label=paste('cut',peak_cnt), col='blue', adj=c(-0.2,0.5))
# 		text(x=log10(max(genecnt_cut)),y=round(peak_score), label=tlabel, adj = c(0.80,2))	
# 	} else {
# 		plot(genecnt_cut,gcntSum, pch=20, cex=.7,main=mtitle, xlab='DEG count', ylab='PIS',ylim=ylim)
# 		lines(genecnt_cut,gcntSum, cex=.7, col='grey25')
# 		abline(v=peak_cnt, col='blue', lty=2)
# 		text(x=peak_cnt,y=10, label=paste('cut',peak_cnt), col='blue', adj=c(-0.2,0.5))
# 		text(x=max(genecnt_cut),y=round(peak_score), label=tlabel, adj = c(0.80,2))
# 	}	
# }


#' Draw peak plot
#'
#' Draw peak plot
#' @param peakObj peak result object. (Output from getPeakResults2)
#' @param mtitle main title for plot
#' @param lw_colname loess smoothed column name. (from names(peakObj$smoothed_peak))

drawGeneCntCutoffPeak2 <- function(peakObj, mtitle='', lw_colname=NULL, logScaleX=FALSE) {
	plotdf <- peakObj$bin_scores
	peak_cnt <- peakObj$peak_cnt

	g1 <- ggplot(plotdf, aes(x=genecnt, y=bin_scores)) + 
	  geom_point() + 
	  # geom_abline(xintercept=peak_cnt, col='green', lty=2)
	  # geom_smooth(method='gam', se=FALSE)
	  labs(x='DEG count', y='PIS', title=mtitle) +
	  theme_bw(base_size=12)

	if(!is.null(lw_colname)) {
		# lw_colname <- colnames(plotdf)[grep('loess',colnames(plotdf))]
		peak_cnt <- peakObj$smoothed_peak[[lw_colname]]$peak_point
		g1 <- g1 + 
		  geom_line(aes(y = get(lw_colname)), data = plotdf, size = 1.5, linetype = 2, col = "red") +
		  geom_vline(xintercept=peak_cnt, col='blue', lty=2)
	} else {
		g1 <- g1 + geom_vline(xintercept=peak_cnt, col='blue', lty=2)
	}
	g1
}



# drawPeakMean <- function(peakObj, mtitle='', ef_cut=2.0, logScaleX=FALSE) {

# 	# par(mfrow = c(2,1), mar = c(3.1, 4.1, 2.1, 2.1))
# 	## cutoff peak sum
# 	indivpeaks <- apply(scoresMat, 1, function(v) (v >= ef_cut & v == max(v)) + 0 )
# 	aaa <- apply(indivpeaks, 1, sum)

# 	# check duplicates
# 	# table(apply(indivpeaks, 2, function(v) sum(v == 1)))

# 	plot(peakObj$genecnt_cut, aaa, pch = 20, cex=.7, main=mtitle, xlab='DEG count', ylab='no of pathways in peak EF')
# 	abline(v = peakObj$genecnt_cut[which.max(aaa)], col='red',lty=2)
# 	abline(v = peakObj$peak_cnt2, col='blue', lty=2)


# }

