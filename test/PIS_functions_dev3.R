###==============================================###
###  Pathway Impact Score Functions v.0.1.3.9000 ###
###==============================================###
message('Version 0.1.3')

## Functions : getGenesByCutoffs name format (2023-01-09)

## To Do (version 0.1.4)
## Style fix : Fix variable names to something more sensible
## Functions : Data style handling : DESeq2, EdgeR, limma (array)
## Functions : Sort getPeakResults2
## Functions : PISresult into class????
## Documents : Add examples
## Documents : Prepare for package
## There seems to be some performance issues with this function 
## where gene sets with few overlaps run slower somehow...

##=========================================================================


# ## Extract gene rank scores from list of result data frames (not share)
# extractScores <- function(resultsLS, value.col, names.col, ncore=2) {
# 	require(parallel)
# 	out <- mclapply(resultsLS, function(dff) {
# 		dff <- dff[which(!is.na(dff[,names.col])),]
# 		scores <- structure(dff[,value.col], names=dff[,names.col])
# 		return(scores)
# 		}, mc.cores=ncore)
# 	return(out)
# }


#' Gene List Extraction
#' 
#' Extract gene list by cut-offs from DE result
#' @param resultDF
#' @param fcos
#' @param qcos
#' @return List of entrez IDs 
#' @export

getGenesByCutoffs <- function(resultDF, fcos, qcos) {
	cutoff_idx <- apply(expand.grid(sprintf('fc%0.01f',fcos), sprintf('q%g',qcos)),1, function(v) paste0(v[2],'_',v[1]))
	# cutoff_idx <- apply(expand.grid(paste0('fc',fcos), paste0('q',qcos)),1, function(v) paste0(v[2],'_',v[1]))
	
	########## In package version, use this. ###############
	### cutoff_idx <- apply(expand.grid(sprintf('fc%.1f',fcos), sprintf('q%.2f,qcos)),1, function(v) paste0(v[2],'_',v[1])) 

	## resultDF cleanup ##
	##

	resultLS <- rep(list(resultDF), length(cutoff_idx))
	names(resultLS) <- cutoff_idx
	fcov <- structure(rep(fcos, length(qcos)), names = cutoff_idx)
	qcov <- structure(rep(qcos, each=length(fcos)), names = cutoff_idx)

	# extractGeneList from Lazy2	
	geneList <- list(up = list(), dn = list(), to=list())
	for(aid in names(resultLS)) {
		resultDF.f <- resultLS[[aid]] %>% filter( !is.na(entGene) )

		geneList$up[[aid]] <- with(resultDF.f, unique(entGene[which(adj.P.Val < qcov[aid] & logFC >=  log2(fcov[aid]))]) )
		geneList$dn[[aid]] <- with(resultDF.f, unique(entGene[which(adj.P.Val < qcov[aid] & logFC <= -log2(fcov[aid]))]) )
		geneList$to[[aid]] <- unique(resultDF.f$entGene)
	}
	return(geneList)
}

##=======================================================================================
## These functions are for internal use only

## Pathway score function (used inside calculatePathwayScores)
#' @param efs vector of enrichment factors
#' @param ef_cut minimum enrichment factor cut-offs. Default 2

pathScores2 <- function(efs, ef_cut=2) {
	scores <- ifelse(efs >= ef_cut, log2(efs), 0)
	return(scores)
}


#' Hypergeometric Test for Gene set
#'
#' Run hypergeometric test for gene set (used inside calculatePathwayScores)
#'
#' @param query  gene set to query (eg. Differentially Expressed Genes)
#' @param refGMT list of reference gene set (eg. Pathways)
#' @param gspace background gene space. Should contain all genes in query.
#' @param minGeneSet minimum size of gene set. This is used to filter refGMT. Default 10
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 1
#' @param verbose print number of filtered entries in refGMT. Default FALSE

hypergeoTestForGeneset.simple <- function(query, refGMT, gspace, minGeneSet=10, ef.psc=1, verbose=FALSE) {
	require(data.table)

	# match gene space
	if(!all(query %in% gspace)) {
		# stop(paste(length(setdiff(query, gspace)),'query items were found outside of background space. Check inputs.'))
		query <- intersect(query, gspace)
	}
	if(!all(unlist(refGMT) %in% gspace)) {
		refGMT <- lapply(refGMT, function(g) intersect(g,gspace))		
	}

	if(length(query) == 0) stop('query should be character vector of at least length 1.')

	# filter refGMT with less than minimum gene set
	exc <- which(sapply(refGMT, length) < minGeneSet)
	if(length(exc) != 0) {
		if(verbose) {
			if(length(exc) <= 5) {
				mesg <- paste('Ref set no.', paste(exc, collapse=', '), 'had less than', minGeneSet,'genes and were excluded.')
			} else {
				mesg <- paste(length(exc), ' entries in refGMT had less than', minGeneSet,'genes and were excluded.')
			}
			message(mesg)
		}
		refGMT <- refGMT[which(sapply(refGMT, length) >= minGeneSet)]
	}
	if(length(refGMT) == 0) stop('Length of refGMT after filtering for minGeneSet is zero. Set lower minGeneSet or check gene inputs.')

	# hypergeometric test
	N <- length(gspace)								# no of balls in urn
	k <- length(query)								# no of balls drawn from urn (DEG no)
	enrRes <- lapply(refGMT, function(refgenes) {
		q <- length(intersect(refgenes, query))		# no of white balls drawn
		m <- length(intersect(gspace, refgenes)) 	# no of white balls in urn
		# I <- intersect(refgenes, query)

		pVal <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
		odds <- (q + ef.psc) / (m / N * k + ef.psc)
		jacc <- q / length(union(query, refgenes))
		gs.ratio <- paste0(q,'/',k)
		bg.ratio <- paste0(m,'/',N)

		return(data.table(pVal = pVal, oddsRatio=odds, tan = jacc, int=q, gsRatio=gs.ratio, bgRatio=bg.ratio))
		# return(data.table(pVal=pVal, oddsRatio=odds, int=q))
		})
	names(enrRes) <- names(refGMT)

	enrRes <- rbindlist(enrRes, idcol='ID')
	# enrRes$logP <- -log10(enrRes$pVal)

	# enrRes <- enrRes[,c('ID', 'pVal', 'oddsRatio', 'int')]
	return(enrRes)
}


###=======================================================================================
## PIS Functions

## Calculate pathway score for gene list (uses pathScores2 function)
## Outputs matrix of pathway scores where each row is gene set and column is bin.

#' @param glist list of gene set (result from binGenesByCntCutoff)
#' @param gspace background gene space. (for hypergeoTestForGeneset.simple)
#' @param ref_geneset list of reference gene set (eg. Pathways)
#' @param ef_cut EF cut-off for scoring.
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 0
#' @param ncore no of threads to use in mclapply. (ncore=1 uses lapply)
#' @param verbose print run time. Default FALSE

calculatePathwayScores <- function(glist, gspace, ref_geneset, ef_cut=2, ef.psc=1, ncore=4, verbose=TRUE, ...) {

	# if( any(c(ef_cut, min.overlap) <= 0 ) ) stop('ef_cut and min.overlap should both be greater than 0.')

	# match gene space
	if(any(!unlist(ref_geneset) %in% gspace)) {
		ref_geneset <- lapply(ref_geneset, function(g) intersect(g,gspace))		
	}
	if(any(!unlist(glist) %in% gspace)) {
		glist <- lapply(glist, function(g) intersect(g,gspace))		
	}

	# Calculate Pathway scores
	tcheck <- proc.time()
	if(ncore > 1) {
		require(parallel)
		scoresLS <- mclapply(glist, function(gset) {
			if(length(gset) == 0) {
				out <- structure(rep(0, length(ref_geneset)), names=names(ref_geneset))
			} else {
				enres <- hypergeoTestForGeneset.simple(gset, ref_geneset, gspace=gspace, ef.psc=ef.psc, ...)
				# enres$oddsRatio[which(enres$int <= min.overlap)] <- 0
				enres$score <- pathScores2(efs=enres$oddsRatio, ef_cut=ef_cut)
				out <- structure(enres$score, names=enres$ID)
				out <- out[names(ref_geneset)]
			}
			return(out)
			}, mc.cores = ncore)
	} else {
		scoresLS <- lapply(glist, function(gset) {
			if(length(gset) == 0) {
				out <- structure(rep(0, length(ref_geneset)), names=names(ref_geneset))
			} else {
				enres <- hypergeoTestForGeneset.simple(gset, ref_geneset, gspace=gspace, ef.psc=ef.psc, ...)
				# enres$oddsRatio[which(enres$int <= min.overlap)] <- 0
				enres$score <- pathScores2(efs=enres$oddsRatio, ef_cut=ef_cut)
				out <- structure(enres$score, names=enres$ID)
				out <- out[names(ref_geneset)]
			}
			return(out)
			})
	}
	scoresMat <- do.call(cbind, scoresLS)
	if(verbose) print((proc.time() - tcheck)/60)	# ~1min for glist of length 1000
	return(scoresMat)
}

# .calculatePS_internal <- function(gset,gspace,ref_geneset,ef_cut) {
# 	if(length(gset) == 0) {
# 		out <- structure(rep(0, length(ref_geneset)), names=names(ref_geneset))
# 	} else {
# 		enres <- hypergeoTestForGeneset.simple(gset, ref_geneset, gspace=gspace, ef.psc=1)
# 		enres$score <- pathScores2(efs=enres$oddsRatio, setsizes=enres$int, ef_cut=ef_cut)
# 		# enres$score <- pathScores2(efs=enres$oddsRatio, setsizes=geneset_size[enres$ID], ef_cut=ef_cut)
# 		out <- structure(enres$score, names=enres$ID)
# 	}
# 	return(out)
# }


#' Get Peak Result from calculatePathwayScores
#' 
#' Summarize result from calculatePathwayScores and get peak cut-off
#' 
#' @param geneCntList list of genes binned by count (output from binGenesByCntCutoff)
#' @param scoresMat Pathway score matrix (output from calculatePathwayScores)
#' @param verbose print peak scores. Default FALSE

getPeakResults2 <- function(geneCntList, scoresMat, verbose=FALSE) {

	# PIS for each gene count cut-off
	gcntSum <- apply(scoresMat, 2, sum, na.rm=TRUE)

	# gene count cut off at peak
	peak_cnt <- sapply(geneCntList, length)[which.max(abs(gcntSum))]

	# peak_score
	peak_score <- gcntSum[which.max(abs(gcntSum))]

	# gene set in peak
	peak_gset <- geneCntList[[names(peak_cnt)]]

	# scored pathways in peak cut-off
	vv <- sort(abs(scoresMat[,names(peak_cnt)]), decreasing=TRUE)
	scored_pathways <- vv[which(vv != 0)]

	# gene count for each cut-off
	genecnt_cut <- sapply(geneCntList,length)

	# # optimized peak mean value
	# indivpeaks <- apply(scoresMat, 1, function(v) (v >= ef_cut & v == max(v)) + 0 )
	# aaa <- apply(indivpeaks, 1, sum)

	# peak_cnt2 <- mean(sapply(geneCntList, length)[which(aaa==max(aaa))])
	# peak_cnt2 <- mean(sapply(geneCntList, length)[unlist(apply(indivpeaks,2, function(v) which(v == 1)))])

	# scores by bin
	bin_scores <- data.table(bin = names(genecnt_cut), genecnt = genecnt_cut, bin_scores = gcntSum)#, indiv_peak = aaa)

	##=================================================================================
	## THIS NEEDS TO BE MADE TO A CLASS with proper print and keep all parameters
	peakObj <- list(peak_cnt=peak_cnt, #peak_cnt2=peak_cnt2, 
		peak_score=peak_score, peak_pathwayCnt=length(scored_pathways), 
		peak_gset=peak_gset, scored_pathways=scored_pathways, #indiv_peaks=indivpeaks, #smoothed_peak=list(),
		bin_scores = bin_scores )
		# genecnt_cut=genecnt_cut, gcntSum=gcntSum )
	mesg <- paste('peak score :', round(peak_score),'\npeak pathway counts :', peakObj$peak_pathwayCnt)
	if(verbose) message(mesg)
	##=================================================================================

	return(peakObj)
}


## PIS Result Summarize function (testing)
summarizePISResults <- function(scoresMat, genelist, score_mode=NULL) {
	if(!identical(colnames(scoresMat), names(genelist))) stop('sig ids in scoresMat and genelist does not match!')

	pis_with_usual <- data.table(
		sig_id = names(genelist),
		peak_cnt = sapply(genelist, length),							# genecnt: no of DEGs in peak
		peak_score = apply(scoresMat, 2, sum), 							# pis: PIS Score at peak
		peak_pathwayCnt = apply(scoresMat, 2, function(v) sum(v != 0))	# pathwaycnt: no of pathways in peak
		)
	if(!is.null(score_mode)) pis_with_usual$score_mode <- score_mode
	return(pis_with_usual)
}


##=======================================================================================
## Functions for PIS with binned rank (mostly deprecated)


## Create gene rank scores from differential analysis result data frames
#' @param resultDF : result dataframe (DESeq2, edgeR)
#' @param value.col : column as character to use as value in gene rank vector (eg. logFC, stat)
#' @param names.col : column as character to use as names of gene rank vector (eg. entrez, ensembl)

createGeneRank <- function(resultDF, value.col, names.col) {
	resultDF <- resultDF[which(!is.na(resultDF[,value.col])),]
	genesrank <- with(resultDF3, structure(get(value.col), names=as.character(get(names.col))))
	genesrank <- sort(genesrank, decreasing=TRUE)
	return(genesrank)
}


## Genes binning
# Create a list of genes for each binned count
#' @param genesrank : named vector of genes rank (output from createGeneRank)
#' @param max_deg_count : maximum count for DEG
#' @param bin_size : gene rank bin size
#' @param reverse : reverse order (for down regulated genes)

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


## Fit loess smoothing on bin scores
#' @param peakObj peak result object. (Output from getPeakResults2)
#' @param loess.span span option for loess

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


## Draw peak plot
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


## ggplot transparent theme
theme_transparent3 <- function(base_size = 11, base_family = "", base_line_size = base_size/22, 
	base_rect_size = base_size/22) {
	theme_grey(base_size = base_size, base_family = base_family, 
		base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
	theme(
		# theme_bw
		panel.border = element_rect(fill = NA, colour = "grey20"), 
		panel.grid = element_line(colour = "grey92"), 
		panel.grid.minor = element_line(size = rel(0.5)), 

		# transparent background
		plot.background = element_rect(fill = "transparent",colour = NA),
		panel.background = element_rect(fill = "transparent", colour = NA),
		strip.background = element_rect(fill = "transparent", colour = NA, size = 0.7), 
		legend.background = element_rect(fill = "transparent", colour = NA), 
		legend.key = element_blank(),

		panel.ontop=FALSE, complete = TRUE)
}
