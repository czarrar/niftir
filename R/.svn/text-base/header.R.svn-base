#----------------------
# Header Business
#----------------------

#' Generates a list with all supported header attributes
#'
#' If you want to specify some header attributes yourself, you can either supply
#' it as a list with the \code{default_header} argument or as keyword=value.
#' Note that if there are duplicate attributes given in both \code{default_header}
#' and as a key=val, then the key=val will take precedance.
#'
#' @author Zarrar Shehzad
#' 
#' @param default_header A list with some or all supported header attributes
#'  (attributes not specified will be set to \code{NULL})
#' @param ... This can be a keyword and argument to be used in place of or in
#'  addition to the \code{default_header} argument.
#' @return list with header attributes
create.header <- function(default_header=NULL, ...) {
    # Create a list with all NULLs and all supported header values
    header_names <- .Call("getHeaderAttributes")
    header <- vector("list", length(header_names))
    names(header) <- header_names
    
    # Check 1st arg
    if (!is.null(default_header) && !is.list(default_header))
        stop("first argument 'default_header' must be a list");
    
    # Add on additional defaults
    default_header <- defaults(list(...), default_header)
    
    # Check user input
    if (is.null(default_header) || length(default_header)==0)
        return(header)
    
    # Get header names
    default_header_names <- names(default_header)
    if (is.null(default_header_names))
        stop("input list must have name attributes")
    
    # Populate header with user supplied values
    n <- length(default_header)
    for (i in 1:n) {
        if (is.null(default_header[[i]]))
            next
        f <- header_names == default_header_names[[i]]
        if (!any(f))
            stop("input header attribute", default_header_names[[i]], "not valid")
        header[f][[1]] <- default_header[[i]]
    }
    
    if (is.null(header$descrip))
        header$descrip <- "niftir"
    
    return(header)
}

get_datatype <- function(name) {
    datatype <- switch(name, 
        char=2,
        short=4,
        int=8,
        float=16,
        double=64,
        "ERROR"
    )
    if (datatype == "ERROR")
        stop("Incorrect datatype ", name)
    return(datatype)
}

autocal <- function(x, header, dt=NULL) {
    header$cal.min <- min(x)
    header$cal.max <- max(x)
    if (is.null(dt)) {
        header$datatype <- switch(typeof(x), 
            character=2,
            integer=8,
            double=64
        )
    } else {
        header$datatype <- get_datatype(dt) 
    }
    return(header)
}

#' @nord
copygeom <- function(x, y) {
    cpattr <- c("slice.duration", "qform.code", "sform.code", "quatern.b", "quatern.c", "quatern.d", "qoffset.x", "qoffset.y", "qoffset.z", "qfac", "scl.slope", "scl.inter", "qto.xyz", "qto.ijk", "sto.xyz", "sto.ijk")
    header <- x@header[cpattr]
    header <- header[sapply(header, function(i) !is.null(i))]
    y@header <- defaults(header, y@header)
    return(y)
}

#' Access/modify the list of header attributes in \code{nifti} object;
#' or accessing/modifying specific header attributes
#'
#' Note that changing the dim (i.e., \code{dim(x) <- c(10,10,10)}) will change
#' both the header attribute 'dim' and the actual dimensions of your data if
#' your input is a \code{nifti} object
#'
#' @name header
#' @aliases dim dim<- pixdim pixdim<-
#' @title accessing/modifying header attributes
#' 
#' @usage
#'  header(x)
#'  dim(x)
#'  pixdim(x)
#'  header(x) <- list(...)
#'  dim(x) <- c(...)
#'  pixdim(x) <- c(...)
#' 
#' @param x \code{nifti} or \code{niftiXd} object
#' 
#' @return list with header attributes or specific attribute or NULL

#' @nord
setGeneric('header', function(x) standardGeneric('header'))

#' @nord
setMethod('header',
    signature(x='header'),
    function(x) x@header
)

#' @nord
setGeneric('header<-', function(x, value) standardGeneric('header<-'))

#' @nord
setMethod('header<-',
    signature(x='header', value='list'),
    function(x, value) x@header <- create.header(value)
)

#' @nord
setGeneric('dim.header', 
    function(x) standardGeneric('dim.header')
)

#' @nord
setMethod('dim.header',
    signature(x='header'),
    function(x) x@header$dim
)

#' @nord
setMethod('dim<-',
    signature(x='nifti'),
    function(x, value) {
        x@header$dim <- value
        .Primitive("dim<-")(x, value)
    }
)

#' @nord
setGeneric('pixdim', 
    function(x) standardGeneric('pixdim')
)

#' @nord
setMethod('pixdim',
    signature(x='header'),
    function(x) x@header$pixdim
)

#' @nord
setGeneric('pixdim<-', 
    function(x, value) standardGeneric('pixdim<-')
)

#' @nord
setMethod('pixdim<-',
    signature(x='header', value='vector'),
    function(x, value) {
        if (length(value) != length(x@header$dim))
            warning("trying to set pixdim attribute with a number of dimensions
                    that don't match the number in the dim attribute")
        x@header$pixdim <- value
    }
)


#--------------------------------------------
# Utitlity functions using stuff in the header
#--------------------------------------------

#' Transform your array coordinates (ijk) into standard space coordinates (xyz), and vice-versa
#' 
#' @aliases xyz2ijk
#' @title Converting between standard and image coordinates
#' @author Zarrar Shehzad
#' 
#' @usage
#'  ijk2xyz(header, ijk, use.qform=T)
#'  xyz2ijk(header, xyz, use.qform=T)
#' 
#' @param header can be \code{list}, \code{nifti}, or \code{niftiXd} object
#' @param ijk coordinates in your array
#' @param xyz coordinates in standard space
#' @param use.qform if \code{TRUE} will use qform transfomation, otherwise will use sform
#' 
#' @return transformed xyz => ijk coordinates or vice-versa (so a vector of length 3)
ijk2xyz <- function(header, ijk, use.qform=T) {
    try(header <- header@header, silent=T)    # in case it is a nifti type object
    ijk <- c(ijk, 1)
    if (use.qform)
        xyz <- ijk %*% t(header$qto.xyz)
    else
        xyz <- ijk %*% t(header$sto.xyz)
    return(xyz[1:3])
}

#' @nord
xyz2ijk <- function(header, xyz, use.qform=T) {
    try(header <- header@header, silent=T)    # in case it is a nifti type object
    xyz <- c(xyz, 1)
    if (use.qform)
        ijk <- xyz %*% t(header$qto.ijk)
    else
        ijk <- xyz %*% t(header$sto.ijk)
    return(ijk[1:3])
}


#-------------------
# Coordinate your niftiXd
#-------------------

#' Get the ijk coordinates of your niftiXd image
#' 
#' This function looks at the dim attribute of your header to create a data frame
#' that gives the ijk coordinates of your niftiXd image in its original nifti space.
#' You can use the rownames of this data frame to locate the vector index of that
#' coordinate (i.e., the column number in the niftiXd object).
#' 
#' @name coords
#' @aliases coords-methods
#' @title get original ijk coordinates of niftiXd object
#' @author Zarrar Shehzad
#' 
#' @param x \code{niftiXd} object
#' @param usemask trim the coordinates returned based on \code{mask} slot of \code{x}
#'
#' @return data frame with 3 columns (i, j, k coordinates) and rownames corresponding to vector indices

#' @nord
setGeneric('coords', function(x, mask=NULL) standardGeneric('coords'))

#' @nord
setMethod('coords',
    signature(x='vector'),
    function(x, mask=NULL) {
        xdim <- x
        if (is.null(xdim))
            stop("coords method requires a header with a dim attribute")
        
        indices <- lapply(xdim[1:3], seq) # list indices in each dimension
        coords <- expand.grid(indices)  # rows=voxels, cols=x,y,z
        colnames(coords) <- c("i", "j", "k")
        
        if (!is.null(mask)) coords <- coords[mask,]
        
        return(coords)
})

#' @nord
setMethod('coords',
    signature(x='header', mask='missing'),
    function(x) {
        if (any(slotNames(x)=="mask"))
            coords(x@header$dim, x@mask)
        else
            coords(x@header$dim)
})

#' @nord
setMethod('coords',
    signature(x='header', mask="logical"),
    function(x, mask) {
        coords(x@header$dim, mask)
})
