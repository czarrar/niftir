##-------------------
## Masking your niftiXd
##-------------------

#' @nord
read.mask <- function(fname, thresh=0) {
    x <- .Call("read_nifti", abspath(fname), 1, PACKAGE="niftir")
    ret <- as.vector(x$image)
    if (!is.null(thresh))
        ret>thresh
    else
        ret
}


#' Access/set mask info
#' 
#' @name mask
#' @aliases mask<-
#' 
#' @usage
#'  mask(x)
#'  mask(x) <- y
#' 
#' @param x \code{niftiXd} or \code{big.niftiXd} object
#' 
#' @return logical vector for \code{mask(x)}


#' @nord
setGeneric('mask', function(x) standardGeneric('mask'))

#' @nord
setMethod('mask',
    signature(x='big.nifti4d'),
    function(x) x@mask
)

#' @nord
setGeneric('mask<-', function(x, value) standardGeneric('mask<-'))

#' @nord
setMethod('mask<-',
    signature(x='big.nifti4d'),
    function(x, value) {
        if (!is.vector(value) && !is.logical(value))
            stop("mask must be vector and logical")
        x@mask <- value
    }
)

#' Prepare a mask to do masking
#' 
#' Intended as a function called by \code{do.mask}.
#' Ensures that input \code{new_mask} is a vector of logicals for masking.
#' 
#' If \code{orig_mask} is specified, then this implies that \code{new_mask}
#' will be masking an object later that has already been masked. The 
#' \code{new_mask} can then be either the same length as \code{orig_mask} or the 
#' length of the number of TRUE elements in \code{new_mask}.
#' 
#' @title Prepare mask
#' @author Zarrar Shehzad
#'
#' @usage prepare.mask(new_mask, orig_mask=NULL, thresh=0)
#'
#' @param new_mask vector or object that can be converted to a vector
#' @param orig_mask NULL or vector of logical elements
#' @param thresh if new_mask is numeric than threshold at this level
#'
#' @return vector of logicals
#'
#' @seealso \code{\link{do.mask}}
prepare.mask <- function(new_mask, orig_mask=NULL, thresh=0) {
    if (!is.numeric(new_mask) && !is.logical(new_mask))
        stop("new_mask must be a numeric or logical vector")
    if (!is.null(orig_mask) && !is.logical(orig_mask))
        stop("orig_mask must be a logical vector")
    
    ## enforce new_mask is a vector
    new_mask <- as.vector(new_mask)
    
    # 1. Ensure mask is logical
    if (!is.logical(new_mask))
        new_mask <- new_mask>thresh
    
    # 2. Check orig_mask
    if (is.null(orig_mask))
        return(new_mask)
    else if (!is.logical(orig_mask))
        stop("if orig mask specified, then must be vector of logical elements")
        
    # 3. Mask can have 2 possible lengths
    lenmask <- length(new_mask)
    ## a. = length of original mask
    if (lenmask == length(orig_mask)) {
        new_mask <- new_mask & orig_mask
    }    
    ## b. = length of total # of TRUE elements in orig_mask
    else if (lenmask == sum(orig_mask)) {
        tmp <- orig_mask
        tmp[tmp] <- new_mask
        new_mask <- tmp
    }
    ## error
    else {
        stop("mask must have length equal to length of mask slot in x or number 
            of columns of x")
    }
    
    return(new_mask)
}

#' Mask/unmask the voxels of your niftiXd object
#' 
#' Masking a niftiXd object involves including only those columns (voxels) that
#' are specified as TRUE or above a given threshold in the input \code{mask}.
#' 
#' The benefit of using this function over just \code{x[,mask]} is that it will
#' keep a record of the regions that were masked, allowing you to go back to the
#' original structure with \code{\link{do.unmask}} or easily save your file in 
#' the appropriate dimensions with \code{\link{write.nifti}}.
#' 
#' The \code{mask} input argument can have a length that is:
#' (1) equal to the total number of voxels in the original nifti image
#' (2) equal to the number of columns or elements in your input \code{x}
#' 
#' Note that columns or elements in a \code{niftiXd} correspond to different
#' voxels in 3D space.
#' 
#' The \code{do.unmask} function puts back the columns or vector elements that
#' were previously masked. The values of these new elements will be set to 0.
#' 
#' @name do.mask
#' @aliases do.unmask
#' @title Masking/unmasking your niftiXd object
#' @author Zarrar Shehzad
#' 
#' @usage
#'  do.mask(x, mask, thresh=0, output.prefix=NULL)
#'  do.unmask(x, return.niftiXd=TRUE)
#' 
#' @param x \code{niftiXd} object
#' @param mask A vector
#' @param thresh If the mask isn't logical, what values should we threshold it
#'  at (default: 0)
#' @param output.prefix if you have big.niftiXd object that is file-backed,
#'  you can specify this option
#' @param return.niftiXd For \code{do.unmask}, whether or not to return a 
#'  niftiXd type object (doesn't apply to big.niftiXd objects)
#' @param ... Additional options possible for \code{do.mask}
#' 
#' @return A masked \code{niftiXd} or \code{big.niftiXd} object


#' @nord
setGeneric('do.mask', function(x, mask, thresh=0, ...)
           standardGeneric('do.mask'))

#' @nord
setMethod('do.mask',
    signature(x='big.nifti4d', mask='character'),
    function(x, mask, thresh=0, ...) {
        mask <- read.mask(mask)  ## assume 3D
        return(do.mask(x, mask, thresh, ...))
    }
)

#' @nord
setMethod('do.mask',
    signature(x='big.nifti4d', mask='vector'),
    function(x, mask, thresh=0, ...) {
        # Prepare mask
        mask <- prepare.mask(mask, x@mask, thresh)
        n <- sum(mask)
        
        # Get indices of mask
        inds <- which(mask[x@mask])
        
        # Copy over masked data
        y <- deepcopy(x, cols=inds, ...)
        
        # Return a big.nifti4d object
        y <- as.big.nifti4d(y, header=x@header, mask=mask, ...)
        
        return(y)
    }
)

#' @nord
setGeneric('do.unmask', function(x)
           standardGeneric('do.unmask'))

# TODO: make an unmask for big.nifti4d


#-------------------
# Indices of your niftiXd
#-------------------

#' Get the indices of your mask
#' 
#' @name indices
#' @aliases indices-methods
#' @author Zarrar Shehzad
#' 
#' @param x \code{niftiXd} object or mask
#' @param i mask that might be applied to \code{x]}
#'
#' @return vector of column or vector indices in \code{x}

#' @nord
setGeneric('indices', function(x, i) standardGeneric('indices'))

#' @nord
setMethod('indices',
    signature(x='logical', i="logical"),
    function(x, i) {
        i <- prepare.mask(i, x)
        i <- which(i[x])
        return(i)
})

#' @nord
setMethod('indices',
    signature(x='logical', i="integer"),
    function(x, i) {
        oi <- order(i)
        if (oi[1]<1 || oi[length(i)]>length(x))
            stop("Illegal index usage")
        return(i)
})

#' @nord
setMethod('indices',
    signature(x='big.nifti4d', i="logical"),
    function(x, i) {
        indices(x@mask, i)
})
