## on package load

.onAttach <- function(libname, pkgname) {
	version <- packageDescription(pkgname, fields = "Version")
	msg <- paste0(pkgname, " version ", version)

	today <- format(Sys.time(), "%m%d")
	if(today == "0101") {
		msg <- paste0(msg, "\n\n*-*-*-*-*-*-*-*-*-*\n* Happy New Year! *\n*-*-*-*-*-*-*-*-*-*\n")
	}
	if(today == "0401") {
		msg <- paste0(msg, "\nApril fools!")
		# msg <- paste0(msg, "\nSYSTEM SELF DESTRUCT IN 5 MINUTES")
	}
	if(today == "0504") {
		msg <- paste0(msg, "\nMay the Fourth be with you!")
	}
	if(today %in% c("1224","1225") ) {
		msg <- paste0(msg, "\nMerry Christmas!")
	}

	packageStartupMessage(msg)
}


