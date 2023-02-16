## PIS object S3 methods

PISobj <- setClass("PISobj")

#' PIS object S3 class
#' 
#' Summarize result from calculatePathwayScores and get peak cut-off. 
#' @export

.newPISobj <- function() {

	obj = list(
		peak_cnt = structure(integer(1), names = character(1)), 
		peak_score = structure(numeric(1), names = character(1)), 
		peak_pathwayCnt = integer(1), peak_gset = character(), 
		scored_pathways = structure(numeric(0), names = character(0)),
		bin_scores = data.table( bin = character(0), genecnt = integer(0), bin_scores = numeric(0) )
	)
	obj = structure(obj, class = "PISobj")
	return(obj)
}

# setOldClass("PISobj")

# ## Methods
# setMethod(show, "PISobj", function(object) {
# 	cat("\n")
# 	cat("  PIS Score:", object$peak_score, "\n")
# 	cat("  PIS threshold:", names(object$peak_cnt), "\n")
# 	cat("  No. of DEGs in peak:", object$peak_cnt,"\n")
# 	cat("  No. of enriched pathways in peak:", object$peak_pathwayCnt,"\n")
# 	cat("\n")
# })

# show.PISobj <- function(object) {
# 	cat("\n")
# 	cat("  PIS Score:", object$peak_score, "\n")
# 	cat("  PIS threshold:", names(object$peak_cnt), "\n")
# 	cat("  No. of DEGs in peak:", object$peak_cnt,"\n")
# 	cat("  No. of enriched pathways in peak:", object$peak_pathwayCnt,"\n")
# 	cat("\n")
# }

print.PISobj <- function(object) { 
	cat("\n")
	cat("  PIS Score:", object$peak_score, "\n")
	cat("  PIS threshold:", names(object$peak_cnt), "\n")
	cat("  No. of DEGs in peak:", object$peak_cnt,"\n")
	cat("  No. of enriched pathways in peak:", object$peak_pathwayCnt,"\n")
	cat("\n")
	
}
# print.PISobj <- function(object) { show(object) }

# head.PISobj <- function(object) {
# 	show(object)
# 	out = lapply(object[c("peak_gset", "scored_pathways","bin_scores")], head)
# 	print(out)
# }




