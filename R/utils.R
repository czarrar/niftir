niftir.split.indices <- function(from, to, by=NULL, length.out=NULL, squeeze=FALSE) {
    if (!is.null(by) && !is.null(length.out))
        stop("too many arguments")
    
    if (!is.null(length.out)) {
        n <- to - from + 1
        by <- floor(n/length.out)
        starts <- seq(from, to - n%%length.out, by=by)
    } else if (!is.null(by)) {
        starts <- seq(from, to, by=by)
    } else {
        stop("must specify either by or length.out")
    }
    
    ends <- starts + by - 1
    n <- length(ends)
    
    if (squeeze && (ends[n] == (ends[n-1]+1))) {
        starts <- starts[1:(n-1)]
        ends <- ends[1:(n-1)]
        n <- n - 1
    }
    
    ends[n] <- to
    
    return(list(starts=starts, ends=ends, n=n))
}

#' Remove extension from file
#'
#' @param fname filename
#' @param compressed do you want to remove a compression extension (.gz or .zip),
#'  before removing the real extension?
#'
#' @return character with extension removed
#'
#' @examples
#'  rmext('test.nii.gz')    # should get test
#'  rmext('test.nii.gz', FALSE) # should get test.nii
rmext <- function(fname, compressed=T) {
    if (!is.character(fname))
        stop("input argument must be a character")
    if (compressed==T)
        re = "[.]([a-zA-Z]{2,})[.](gz|bz2)$|[.]([a-zA-Z0-9]{1,})$"
    else
        re = "[.]([a-zA-Z0-9]{1,})$"
    sub(re, "", fname)
}

getext <- function(fname, compressed=T) {
    if (!is.character(fname))
        stop("input argument must be a character")
    if (compressed==T)
        re = "[.]([a-zA-Z]{2,})[.](gz|bz2)$|[.]([a-zA-Z0-9]{1,})$"
    else
        re = "[.]([a-zA-Z0-9]{1,})$"
    f <- regexpr(re, fname)
    if (f == -1)
        return("")
    else
        return(substr(fname, f+1, f+attr(f,"match.length")))
}

#' Convert relative path to absolute path
#'
#' This will also deal with expanding a path as in \code{\link{path.expand}}.
#' You must ensure that the directory to your path exists, otherwise you will
#' get an error.
#'
#' @param path character specifying path in your file system
#'
#' @return character with absolute path
#'
#' @examples
#'  abspath('~/test.nii.gz')
#'  curdir <- getwd()
#'  abspath('../test.nii.gz')
#'  file.path(dirname(curdir), 'test.nii.gz')    # equivalent to above
#'  abspath('test.nii.gz')
#'  file.path(curdir, 'test.nii.gz')    # equivalent to above
abspath <- function(path) {
    if (!is.character(path))
        stop("input must be a character")
    dpath <- dirname(path)
    if (!file.exists(dpath))
        stop("can't find base of path: ", path)
    bpath <- basename(path)
    orig_path <- getwd()
    setwd(dpath)
    new_path <- file.path(getwd(), bpath)
    setwd(orig_path)
    return(new_path)
}

#' Z-transform a correlation value
#'
#' @title Z-Transform
#'
#' @param x Correlation value
#' @return Fischer's Z transform value
#'
#' @seealso \code{\link{atanh}}
#' 
#' @examples
#'  r2z(0.1)
#'  r2z(1)
r2z <- function(x) atanh(x)

#' Convert a Fischer's Z value to a correlation value
#'
#' @title Reverse Z-Transfor
#'
#' @param x Fischer's Z transform value
#' @return correlation value
#'
#' @seealso \code{\link{tanh}}
#' 
#' @examples
#'  z2r(r2z(0.1))
#'  z2r(r2z(1))
z2r <- function(x) tanh(x)

#' Wrapper around txtProgressBar
#'
#' Will start the progress bar at 0 and go up to \code{limit} specified.
#' Uses the 3rd style in \code{\link{txtProgressBar}}.
#; 
#' @title Progress Bar
#'
#' @param limit maximum point of the progress bar (integer)
#' @return \code{txtProgressBar} object
#'
#' @seealso \code{\link{update.txtProgressBar}} \code{\link{end.txtProgressBar}}
#' 
#' @examples
#'  n <- 10
#'  pb <- progressbar(10)
#'  for (i in 1:n) {
#'      Sys.sleep(1)
#'      update(pb, i)
#'  }
#' end(pb)
progressbar <- function(limit, message=NULL) {
    if (!is.null(message))
        cat(message)
    pb <- txtProgressBar(min=0, max=limit, style=3)
    update(pb, 0.01)
    pb
}

#' Wrapper around setTxtProgressBar
#; 
#' @title Update Progress Bar
#'
#' @param pb \code{txtProgressBar} object
#' @param i current point that you are at in your task
#'  (must be less than \code{limit} as specified in \code{progressbar})
#' @return invisible number specifying the previous point in your progress
#'
#' @seealso \code{\link{progressbar}} \code{\link{end.txtProgressBar}}
update.txtProgressBar <- function(pb, i) setTxtProgressBar(pb, i)

#' Close txtProgressBar object
#; 
#' @title Close/End Progress Bar
#'
#' @param pb \code{txtProgressBar} object
#'
#' @seealso \code{\link{progressbar}} \code{\link{update.txtProgressBar}}
end.txtProgressBar <- function(pb) { cat("\n"); close(pb) }
