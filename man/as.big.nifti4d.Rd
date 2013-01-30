\name{as.big.nifti4d}
\alias{as.big.nifti4d}
\alias{as.big.nifti4d}
\title{Create a big...}
\usage{as.big.nifti4d(x)               # when x is a nifti object or nifti4d object
as.big.nifti4d(x, header)       # when x is an array
as.big.nifti4d(x, header, mask) # when x is a matrix or big.matrix}
\description{Create a big.nifti4d object}
\author{Zarrar Shehzad}
\value{\code{big.nifti4d} object}
\seealso{\code{\link{as.big.nifti4d-methods}}}
\keyword{methods}
\arguments{\item{x}{4d \code{nifti}, \code{nifti4d}, 4d \code{array}, \code{matrix}, or
\code{big.matrix}}
\item{header}{list of header attributes (required when \code{x} is an
\code{array} or \code{matrix} or \code{big.matrix})}
\item{mask}{logical vector specifying which voxels from 3D image are
specified with \code{x} (required when \code{x} is a \code{matrix} or
\code{big.matrix})}
\item{...}{Additional arguments that will be passed when creating a big.matrix}
}
\examples{as.big.nifti4d(array(0, c(10,10,10,10)), create.header())
as.big.nifti4d(nifti(0, dim=c(10,10,10,10)))   # should give same thing as above}

