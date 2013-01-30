\name{nifti}
\alias{nifti}
\alias{as.nifti}
\alias{is.nifti}
\title{Nifti}
\description{
    \code{nifti}: create a nifti object from input\n
    \code{as.nifti}: converts input to become a nifti object\n
    \code{is.nifti}: checks if input is a nifti object
}
\usage{
    nifti(x, header=NULL, dim=NULL)
    as.nifti(x, header=NULL)
    is.nifti(x)
}
\arguments{
    \item{x}{
        a vector, matrix, or array for \code{nifti}; can also be a nifti, big.nifti4d, niftiXd for \code{as.nifti}
    }
    \item{header}{
        a list with the header attributes
    }
    \item{dim}{
        dimensions of input \code{x}
    }
}
\details{
    A \code{nifti} object is an array with an additional header slot.

    When the \code{header} argument is missing or \code{NULL}, a list is created
    with all supported header attributes set to \code{NULL} except for dim which
    is set to the dim of \code{x}.
}
\value{
    \code{nifti} object for \code{nifti} and \code{as.nifti}\n
    \code{TRUE} or \code{FALSE} for \code{is.nifti}
}
\author{Zarrar Shehzad}
\seealso{\code{\link{as.nifti-methods}}}
\keyword{methods}
\examples{
    nim <- nifti(0, list(dim=c(10,10,10)))  # creates a 10x10x10 nifti object with values 0

    tmp <- read.nifti(file.path(system.file("data", package="niftir"), "test.nii.gz"))
    x <- array(0, c(10,10,10))
    nim <- nifti(x, header(tmp))    # OR even better
    nim <- as.nifti(x, header(tmp)) # both create a 10x10x10 nifti object with header from tmp
    is.nifti(nim)   # should be TRUE}
}
