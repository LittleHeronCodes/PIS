## on package load

.onAttach <- function(libname, pkgname) {
	version <- utils::packageDescription(pkgname, fields = "Version")
	msg <- paste0(pkgname, " version ", version)

	today <- format(Sys.time(), "%m%d")
	if(today == "0101") {
		msg <- paste0(msg, "\n\n*-*-*-*-*-*-*-*-*-*\n* Happy New Year! *\n*-*-*-*-*-*-*-*-*-*\n")
	}
	if(today == "0401") {
		msg <- paste0(msg, "\nError: April fools!")
	}
	if(today %in% c("1224","1225") ) {
		msg <- paste0(msg, "\n\n*-*-*-*-*-*-*-*-*-*\n* Merry Christmas! *\n*-*-*-*-*-*-*-*-*-*\n")
	}

	packageStartupMessage(msg)
}


