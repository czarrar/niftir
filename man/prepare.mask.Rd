\name{prepare.mask}
\alias{prepare.mask}
\title{Prepare mask}
\usage{prepare.mask(new_mask, orig_mask=NULL, thresh=0)}
\description{Prepare a mask to do masking}
\details{Intended as a function called by \code{do.mask}.
Ensures that input \code{new_mask} is a vector of logicals for masking.

If \code{orig_mask} is specified, then this implies that \code{new_mask}
will be masking an object later that has already been masked. The \code{new_mask}
can then be either the same length as \code{orig_mask} or the length of the
number of TRUE elements in \code{new_mask}.}
\author{Zarrar Shehzad}
\value{vector of logicals}
\seealso{\code{\link{do.mask}}}
\arguments{\item{new_mask}{vector or object that can be converted to a vector}
\item{orig_mask}{NULL or vector of logical elements}
\item{thresh}{if new_mask is numeric than threshold at this level}
}

