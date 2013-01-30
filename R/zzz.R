# ###########
# FILE: zzz.R
#

#' @nord
.onLoad <- function(libname, pkgname) {
    library.dynam("niftir", pkgname, libname);
    cat("\nLoading niftir version < 1.0.\n\n")
}

#.noGenerics <- TRUE           # This was a problem, not used.

#' @nord
.onUnload <- function(libpath) {
    library.dynam.unload("niftir", libpath);
}


# \name{unitTests}
# \alias{unitTest.niftir}
# \title{ Unit tests for the package niftir }
# \description{
#     Performs unit tests defined in this package by running
#     \code{example(unitTests.niftir)}. Tests are in \code{runit*.R} files
#     that are located in the '/unitTests' subdirectory or one of its
#     subdirectories ('/inst/unitTests' and subdirectories of the package source).
# }
# \author{Zarrar Shehzad}
# \examples{
# 
#     library(svUnit)
#     library(niftir)
#     
#     clearLog()
#     (runTest(svSuite("package::niftir"), "niftir"))
#     errorlog()
# }
# \keyword{utilities}
