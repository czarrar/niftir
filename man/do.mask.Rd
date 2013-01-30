\name{do.mask}
\alias{do.mask}
\alias{do.mask}
\title{Masking/unmasking your niftiXd object}
\usage{do.mask(x, mask, thresh=0, output.prefix=NULL)
do.unmask(x, return.niftiXd=TRUE)}
\description{Mask/unmask the voxels of your niftiXd object}
\details{Masking a niftiXd object involves including only those columns (voxels) that
are specified as TRUE or above a given threshold in the input \code{mask}.

The benefit of using this function over just \code{x[,mask]} is that it will
keep a record of the regions that were masked, allowing you to go back to the
original structure with \code{\link{do.unmask}} or easily save your file in the
appropriate dimensions with \code{\link{write.nifti}}.

The \code{mask} input argument can have a length that is:
(1) equal to the total number of voxels in the original nifti image
(2) equal to the number of columns or elements in your input \code{x}

Note that columns or elements in a \code{niftiXd} correspond to different
voxels in 3D space.

The \code{do.unmask} function puts back the columns or vector elements that
were previously masked. The values of these new elements will be set to 0.}
\alias{do.unmask}
\author{Zarrar Shehzad}
\value{A masked \code{niftiXd} object for \code{mask} and \code{unmask};}
\arguments{\item{x}{\code{niftiXd} object}
\item{mask}{A vector}
\item{thresh}{If the mask isn't logical, what values should we threshold it
at (default: 0)}
\item{output.prefix}{if you have big.niftiXd object that is file-backed,
you can specify this option}
\item{return.niftiXd}{Whether or not to return a niftiXd type object
(doesn't apply to big.niftiXd objects) for \code{unmask} function only}
}

