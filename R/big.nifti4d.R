#----------------------
# as.big.nifti4d FUNCTIONS
#----------------------

#' Create a big.nifti4d object
#'
#' @name as.big.nifti4d-methods
#'
#' @aliases as.big.nifti4d,nifti-method
#'  as.big.nifti4d,nifti4d-method
#'  as.big.nifti4d,array,list-method
#'  as.big.nifti4d,matrix,list,logical-method
#'  as.big.nifti4d,big.matrix,list,logical-method
#'
#' @sectionMethods
#'  \describe{
#'      \item{\code{signature(x="nifti")}}{...}
#'      \item{\code{signature(x="nifti4d")}}{...}
#'      \item{\code{signature(x="array", header="list")}}{...}
#'      \item{\code{signature(x="matrix", header="list", mask="logical")}}{...}
#'      \item{\code{signature(x="big.matrix", header="list", mask="logical")}}{...}
#'  }
#'
#' @seealso \code{\link{as.big.nifti4d}}
#'
#' @keywords methods


#' Create a big.nifti4d object
#'
#' @name as.big.nifti4d
#'
#' @usage
#'  as.big.nifti4d(x)               # when x is a nifti object or nifti4d object
#'  as.big.nifti4d(x, header)       # when x is an array
#'  as.big.nifti4d(x, header, mask) # when x is a matrix or big.matrix
#'
#' @author Zarrar Shehzad
#' 
#' @param x 4d \code{nifti}, \code{nifti4d}, 4d \code{array}, \code{matrix}, or
#'  \code{big.matrix}
#' @param header list of header attributes (required when \code{x} is an
#'  \code{array} or \code{matrix} or \code{big.matrix})
#' @param mask logical vector specifying which voxels from 3D image are
#'  specified with \code{x} (required when \code{x} is a \code{matrix} or
#'  \code{big.matrix})
#' @param ... Additional arguments that will be passed when creating a big.matrix
#' 
#' @return \code{big.nifti4d} object
#' 
#' @seealso \code{\link{as.big.nifti4d-methods}}
#'
#' @examples
#'  as.big.nifti4d(array(0, c(10,10,10,10)), create.header())
#'  as.big.nifti4d(nifti(0, dim=c(10,10,10,10)))   # should give same thing as above
#' 
#' @keywords methods


#' @nord
setGeneric('as.big.nifti4d', 
    function(x, header=NULL, mask=NULL, ...) standardGeneric('as.big.nifti4d')
)

#' @nord
setMethod('as.big.nifti4d', 
    signature(x='big.nifti4d'),
    function(x) return(x)
)

#' @nord
setMethod('as.big.nifti4d',
    signature(x='nifti', header='missing', mask='missing'),
    function(x, ...) as.big.nifti4d(x@.Data, x@header, ...)
)
 
#' @nord
setMethod('as.big.nifti4d',
    signature(x='array', header='list', mask='missing'),
    function(x, header, ...) {
        if (length(dim(x)) != 4)
            stop("dimensions not equal to 4")
        
        header$dim <- dim(x)  # ensure header/image consistency
        dim(x) <- c(dim(x)[4], prod(dim(x)[1:3]))  # rows = timepoints & cols = voxels
        bigx <- as.big.matrix(x, ...)
        
        as.big.nifti4d(bigx, header, rep(TRUE, prod(header$dim[1:3])))
    }
)

#' @nord
setMethod('as.big.nifti4d',
    signature(x='matrix', header='list', mask='logical'),
    function(x, header, mask, ...) {
        bigx <- as.big.matrix(x, ...)
        as.big.nifti4d(bigx, header, mask, ...)
    }
)

#' @nord
setMethod('as.big.nifti4d',
    signature(x='big.matrix', header='list', mask='logical'),
    function(x, header, mask, ...) {
        # Check input
        if (length(header$dim) != 4)
            stop("dimensions of header attribute not equal to 4")
        if (prod(header$dim[1:3]) != length(mask))
            stop("mask must have length equal to total number of voxels in
                 nifti image")
        if (ncol(x) != sum(mask))
            stop("mask must have as many TRUE elements as columns in input x")
        
        # Get absolute path for backingpath
        if (is.filebacked(x)) {
            args <- list(...)
            bpath <- abspath(args$backingpath)
            bfile <- basename(args$backingfile)
            dfile <- basename(args$descriptorfile)
            # TODO: if file-backed, then save some additional info?
            # with the header, mask, backingpath, etc?
        } else {
            bpath <- ""
            bfile <- ""
            dfile <- ""
        }
        
        # New object
        bigmat <- new("big.nifti4d", header=header, mask=mask, 
            backingpath=bpath, backingfile=bfile, descriptorfile=dfile)
        bigmat@address <- x@address
        
        return(bigmat)
    }
)


#----------------------
# is.big.nifti4d FUNCTIONS
#----------------------

#' @nord
setGeneric('is.big.nifti4d', 
    function(x) standardGeneric('is.big.nifti4d')
)

#' @nord
setMethod('is.big.nifti4d',
    signature(x='big.nifti4d'),
    function(x) return(TRUE)
)

#' @nord
setMethod('is.big.nifti4d',
    definition=function(x) return(FALSE)
)


#----------------------
# read.nifti4d FUNCTION
#----------------------

# Read in header and pointer to nifti class
read.big.nifti_obj <- function(file) {
    obj <- .Call("read_bignifti_header", abspath(file))
    return(obj)
}

# Generate big matrix for later big nifti object
read.big.nifti_gen.bigmat <- function(nifti, tpts=NULL, voxs=NULL, type=NULL, ...) 
{
    # Get new 2d dimensions for original data
    ## note that the setup here is sort-of the transpose of the
    ## original original data
    orig_dim <- nifti$header$dim
    if (length(orig_dim) != 4)
        stop("needs to be a 4 dimensional nifti image")
    orig_dim <- c(orig_dim[4], prod(orig_dim[1:3]))
    
    # Set time-points
    if (is.null(tpts) && !is.null(voxs)) {
        tpts <- 1:orig_dim[1]
        ntpts <- length(tpts)
    } else {
        tpts <- NULL
        ntpts <- orig_dim[1]
    }
    
    # Set voxels
    if (is.null(voxs) && !is.null(tpts)) {
        voxs <- 1:orig_dim[2]
        nvoxs <- length(voxs)
    } else {
        nvoxs <- orig_dim[2]
    }
    
    # Get the data type
    if (is.null(type)) {
        type <- .Call("niftir_datatype_string", nifti$header$datatype)
        type <- switch(type, 
            BINARY = "char", 
            INT8 = "char", 
            UINT8 = "short", 
            INT16 = "short", 
            UINT16 = "integer", 
            INT32 = "integer", 
            UINT32 = "double", 
            FLOAT32 = "double", 
            FLOAT64 = "double", 
            UNKNOWN = "double", 
            stop("Unsupported datatype ", type)
        )
    } else if (!(type %in% c("char", "short", "integer", "double"))) {
        stop("unsupported datatype ", type)
    }
    
    # Create the matrix
    bigmat <- big.matrix(ntpts, nvoxs, type=type, ...)
    
    # Add needed attributes
    attr(bigmat, "read_partial") <- !(is.null(tpts) && is.null(voxs))
    attr(bigmat, "tpts") <- tpts
    attr(bigmat, "voxs") <- voxs
    
    return(bigmat)
}

# Copy over data from file in nifti_obj to bigmat
read.big.nifti_read <- function(nifti_obj, bigmat) {
    # Which time-points and voxels to read in
    read_partial <- attr(bigmat, "read_partial")
    tot_voxs <- prod(nifti_obj$header$dim[1:3])
    
    # Read!
    if (read_partial) {
        tpts <- as.double(attr(bigmat, "tpts"))
        voxs <- as.double(attr(bigmat, "voxs"))
        .Call("read_partial_bignifti_data", nifti_obj$address, bigmat@address, 
              tpts, voxs, tot_voxs, PACKAGE="niftir")
    } else {
        .Call("read_bignifti_data", nifti_obj$address, bigmat@address, PACKAGE="niftir")
    }
    
    return(bigmat)
}

#' Read in a big.nifti4d object from a file
#'
#' @usage read.big.nifti4d(fname, ...)
#'
#' @author Zarrar Shehzad
#' 
#' @param fname character specifying path to analyze/nifti file
#' @param ... Additional arguments passed to \code{\link{big.matrix}}
#'
#' @return \code{big.nifti4d} object
#' 
#' @seealso \code{\link{read.nifti4d}}, \code{\link{as.big.nifti4d}}
#'
#' @examples
#'  # TODO
#' 
#' @keywords methods
read.big.nifti4d <- function(file, ...) {
    read.big.nifti(file, nifti4d=TRUE, ...)
}

read.big.nifti <- function(file, nifti4d=FALSE, type=NULL, ...) {
    # Read in header and pointer to nifti class
    nifti_obj <- read.big.nifti_obj(file)
    
    # Get output big matrix
    bigmat <- read.big.nifti_gen.bigmat(nifti_obj, type=type, ...)
    
    # Read in data
    bigmat <- read.big.nifti_read(nifti_obj, bigmat)
    
    # Clear nifti address
    header <- nifti_obj$header
    rm(nifti_obj); invisible(gc(F, T))
    
    if (nifti4d) {
        # Create mask
        mask <- vector("logical", prod(header$dim[1:3]))
        voxs <- attr(bigmat, "voxs")
        if (is.null(voxs))
            mask[] <- T
        else
            mask[voxs] <- T
        # Convert to big.nifti4d object
        bigmat <- as.big.nifti4d(bigmat, header, mask, ...)
    }
    
    return(bigmat)
}

#' @nord
setMethod("free.memory",
    signature(x="big.nifti4d", backingpath="missing"),
    function(x) {
        if (!is.filebacked(x))
            stop("input to free.memory cannot be a non-filebacked big.matrix")
        # free up memory
        .Call("CDestroyBigMatrix", x@address, PACKAGE="bigmemory")
        gc()
        # reattach matrix
        tmp <- attach.big.matrix(x@descriptorfile)
        x@address <- tmp@address
        # done!
        return(x)
    }
)

#' @nord
setMethod("free.memory",
    signature(x="list", backingpath="missing"),
    function(x) {
        xs <- x
        # check input
        lapply(xs, function(x) {
            if (!is.big.niftiXd(x))
                stop("input to free.memory must be bigniftiXd object")
            if (!is.filebacked(x))
                stop("input to free.memory cannot be a non-filebacked big.matrix")
        })
        # free memory
        lapply(xs, function(x) 
            .Call("CDestroyBigMatrix", x@address, PACKAGE="bigmemory")
        )
        gc()
        # reattach matrices
        for (i in 1:length(xs)) {
            tmp <- attach.big.matrix(xs[[i]]@descriptorfile)
            xs[[i]]@address <- tmp@address
        }
        # done!
        return(xs)
    }
)
