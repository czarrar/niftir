#----------------------------
# Main Nifti Object Functions
#----------------------------

nifti <- function(x, header=NULL, dim=NULL) {
    if (!is.vector(x) && !is.array(x))
        stop("input x must be a vector, matrix, or array")
    
    if (!is.null(dim))
        x <- array(x, dim)
    else
        x <- as.array(x)
    
    if (is.null(header))
        return(as.nifti(x))
    else
        return(as.nifti(x, header))
}

# big.nifti, nifti, nifti4d, niftiXd, array => as.nifti
setGeneric('as.nifti', 
    function(x, header=NULL) standardGeneric('as.nifti')
)

setMethod('as.nifti',
    signature(x='nifti', header='missing'),
    function(x) return(x)
)

setMethod('as.nifti',
    signature(x='big.nifti4d', header='missing'),
    function(x) {
        # Create the nifti array
        arr <- array(0, x@header$dim)
        y <- aperm(x[,])
        arr[x@mask] <- y
        return(as.nifti(arr, x@header))
    }
)

setMethod('as.nifti',
    signature(x='array', header='list'),
    function(x, header) {
        arr <- new("nifti", header=create.header(header))
        arr@.Data <- x
        arr@header$dim <- dim(x)    # ensure consistency between image and header
        return(arr)
    }
)

setMethod('as.nifti',
    signature(x='array', header='missing'),
    function(x) {
        return(as.nifti(x, list(descrip="niftir")))
    }
)

setGeneric('is.nifti', 
    function(x) standardGeneric('is.nifti')
)

setMethod('is.nifti',
    signature(x='nifti'),
    function(x) return(TRUE)
)

setMethod('is.nifti',
    definition=function(x) return(FALSE)
)

setMethod('show',
    signature(object='nifti'),
    function(object) {
        cat("nifti S4 class:\n")
        cat("---------------\n")
        cat("head(Data)\n")
        print(head(object@.Data))
        cat("\nhead(Header)\n")
        print(head(object@header))
    }
)


#-------------------
# Main IO Functions
#-------------------

#' Read data (header/image) from nifti/analyze file
#' 
#' If you have a both a .hdr and .img, you can input the filename of either
#'
#' @author Zarrar Shehzad
#' @param fname nifti/analyze filename
#' @return \code{\link{nifti}} object
#' @note Some of the C code that this function calls was borrowed
#'          from Rniftilib
#' @seealso \code{\link{read.nifti.header}}, \code{\link{write.nifti}} object
#' @examples
#'  fname <- file.path(system.file("data", package="niftir"), "test.nii.gz") # 2x2x2 size file
#'  # or fname <- file.choose()
#'  nim <- read.nifti(fname)
#'  nim[1,1,2]  # returns voxel value at x=1, y=1, z=2
#'  header(nim)  # returns header attributes as a list
read.nifti <- function(fname) {
    x <- .Call("read_nifti", abspath(fname), 1, PACKAGE="niftir")
    as.nifti(x$image, x$header)
}

#' Read header info from nifti/analyze file
#'
#' @author Zarrar Shehzad
#' @param fname nifti/analyze filename
#' @return list with header attributes
#'  (see output of \code{\link{create.header}} for supported attributes)
#' @note Some of the C code that this function calls was borrowed
#'          from Rniftilib
#' @seealso \code{\link{read.nifti}}
read.nifti.header <- function(fname) {
    x <- .Call("read_nifti", abspath(fname), 0, PACKAGE="niftir")
    x
}

#' Read image data from nifti/analyze file
#'
#' @author Zarrar Shehzad
#' @param fname nifti/analyze filename
#' @return array with nifti data
#' @note Some of the C code that this function calls was borrowed
#'          from Rniftilib
#' @seealso \code{\link{read.nifti}}
read.nifti.image <- function(fname) {
    x <- .Call("read_nifti", abspath(fname), 1, PACKAGE="niftir")
    x$image
}


#' @nord
check.outfile <- function(outfile, overwrite, header) {
    if (!is.character(outfile))
        stop("outfile must be a character")
    
    # Check/set the output filename
    if (is.null(outfile)) {
        if (is.null(header$fname))
            stop('Must give an output filename')
        else
            outfile <- header$fname
    }

    # relative => absolute path
    outfile <- abspath(outfile)
    
    # Check existence of output
    if (file.exists(outfile)) {
        # remove it
        if (overwrite) {
            warning("removing file", outfile)
            file.remove(outfile)
        }
        # whoops
        else {
            warning("output file", outfile, "already exists, NOT OVERWRITING")
        }
    }
    if (!file.exists(dirname(outfile)))
        dir.create(dirname(outfile))
    
    return(outfile)
}

#' @nord
setGeneric('write.nifti', 
    function(x, header, mask, outfile=NULL, overwrite=FALSE, odt=NULL, ...)
        standardGeneric('write.nifti')
)

#' @nord
setMethod('write.nifti',
    signature(x='nifti', header='missing', mask='missing'),
    function(x, outfile=NULL, overwrite=FALSE, odt=NULL, ...) {
        header <- autocal(x, x@header, odt)
        header <- create.header(header, ...)
        outfile <- check.outfile(outfile, overwrite, header)
        .Call("write_nifti", header, x@.Data, outfile, PACKAGE="niftir")
    }
)

#' @nord
setMethod('write.nifti',
    signature(x='big.nifti4d', header='missing', mask='missing'),
    function(x, outfile=NULL, overwrite=FALSE, odt=NULL, ...) {
        library(biganalytics)
        header <- autocal(x, header, odt)
        header <- create.header(x@header, ...)
        outfile <- check.outfile(outfile, overwrite, header)
        .Call("write_bignifti", header, x@address, which(x@mask), outfile,
              PACKAGE="niftir")
    }
)

#' @nord
setMethod('write.nifti',
    signature(x='array', header='list', mask='missing'),
    function(x, header, outfile=NULL, overwrite=FALSE, odt=NULL, ...) {
        header <- autocal(x, header, odt)
        header <- create.header(header)
        outfile <- check.outfile(outfile, overwrite, header)
        .Call("write_nifti", header, x, outfile, PACKAGE="niftir")
    }
)

#' @nord
setMethod('write.nifti',
    signature(x='matrix', header='list', mask='logical'),
    function(x, header, mask, outfile=NULL, overwrite=FALSE, odt=NULL, ...) {
        if (sum(mask) != nrow(x))
            stop("number of TRUE elements in mask is not equal to nrow of input x")
        
        header <- autocal(x, header, odt)
        header <- create.header(header, ...)
        if (prod(header$dim[1:3]) != length(mask))
            stop("dimensions in header attribute do not match the length of x")
        header$dim <- c(header$dim[1:3], ncol(x))
        header$pixdim <- c(header$pixdim[1:3], 1)
        
        out <- matrix(NA, prod(header$dim[1:3]), header$dim[4])
        out[mask,] <- t(x)
        out <- as.double(out)
        dim(out) <- header$dim
        
        outfile <- check.outfile(outfile, overwrite, header)
        .Call("write_nifti", header, out, outfile, PACKAGE="niftir")
    }
)

#' @nord
setMethod('write.nifti',
    signature(x='vector', header='list', mask='logical'),
    function(x, header, mask, outfile=NULL, overwrite=FALSE, odt=NULL, ...) {
        if (sum(mask) != length(x))
            stop("number of TRUE elements in mask is not equal to length of input x")
        
        header <- autocal(x, header, odt)
        header <- create.header(header, ...)
        if (prod(header$dim) != length(mask))
            stop("dimensions in header attribute do not match the length of x")
        
        out <- mask*1
        out[mask] <- x
        out <- as.double(out)
        dim(out) <- header$dim
        
        outfile <- check.outfile(outfile, overwrite, header)
        .Call("write_nifti", header, out, outfile, PACKAGE="niftir")
    }
)

#' @nord
setMethod('write.nifti',
    signature(x='vector', header='list', mask='missing'),
    function(x, header, outfile=NULL, overwrite=FALSE, odt=NULL, ...) {
        header <- autocal(x, header, odt)
        header <- create.header(header, ...)
        if (prod(header$dim) != length(x))
            stop("dimensions in header attribute do not match the length of x")
        
        x <- as.double(x)
        dim(x) <- header$dim
        
        outfile <- check.outfile(outfile, overwrite, header)
        .Call("write_nifti", header, x, outfile, PACKAGE="niftir")
    }
)
