## PIS object S3 methods

PISobj <- setClass("PISobj")

#' Generate PIS object
#' 
#' Generate an empty PIS object or turn a list into PISobj
#' @param obj input object to make into PISobj. Leaving this empty will generate an empty object.
#' @export

new_PISobj <- function(obj = NULL) {
	if( inherits(obj,"PISobj") ) stop("Input is already a PIS object. What are you trying to accomplish?")
	if( is.null(obj) ) {
		obj = list(
			peak_cnt = structure(integer(1), names = character(1)),
			peak_score = structure(numeric(1), names = character(1)),
			peak_pathwayCnt = integer(1), peak_gset = character(),
			scored_pathways = structure(numeric(0), names = character(0)),
			bin_scores = data.table(bin = character(0), genecnt = integer(0), bin_scores = numeric(0))
		)
	} else {
		valid = validate_PISobj(obj)
	}
	obj = structure(obj, class = "PISobj")
	return(obj)
}

#' Validate PIS object
#' 
#' Validate PIS object or check if a list can become PIS object.
#' @param obj PISobj
#' @export

validate_PISobj <- function(obj) {

	# list type check
	if ( !is.list(obj) ) {
		stop("Input should be in list format.", call. = FALSE)
	}
	
	# attribute check
	required <- c("peak_cnt", "peak_score", "peak_pathwayCnt", "peak_gset", "scored_pathways", "bin_scores")
	obj_names <- names(obj)
	if ( !all(required %in% obj_names) ) {
		missing = paste(setdiff(required, obj_names), collapse=", ")
		stop(paste("The following attributes are needed:", missing), call. = FALSE)
	}	

	obj	
}


#' Printing PIS object
#' 
#' Print PIS object.
#' @param x object of 'PISobj'
#' @param ... further arguments passed to or from other methods.
#' @examples 
#' # S3 method for class 'PISobj'
#' x <- new_PISobj()
#' print(x)
#' @export

print.PISobj <- function(x, ...) { 
	cat("\n")
	cat("  PIS Score:", x$peak_score, "\n")
	cat("  PIS threshold:", names(x$peak_cnt), "\n")
	cat("  No. of DEGs in peak:", x$peak_cnt,"\n")
	cat("  No. of enriched pathways in peak:", x$peak_pathwayCnt,"\n")
	cat("\n")
	
}

setMethod(show, "PISobj", function(object) {
	print(object)
})




# head.PISobj <- function(object) {
# 	show(object)
# 	out = lapply(object[c("peak_gset", "scored_pathways","bin_scores")], head)
# 	print(out)
# }





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





