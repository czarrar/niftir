\name{read.nifti}
\alias{read.nifti}
\title{Read data (header/image) from nifti/analyze file...}
\usage{read.nifti(fname)
}
\description{Read data (header/image) from nifti/analyze file}
\details{If you have a both a .hdr and .img, you can input the filename of either}
\author{Zarrar Shehzad}
\value{\code{\link{nifti}} object}
\note{Some of the C code that this function calls was borrowed
from Rniftilib}
\seealso{\code{\link{read.nifti.header}}, \code{\link{write.nifti}} object}
\arguments{\item{fname}{nifti/analyze filename}
}
\examples{fname <- file.path(system.file("data", package="niftir"), "test.nii.gz") # 2x2x2 size file
# or fname <- file.choose()
nim <- read.nifti(fname)
nim[1,1,2]  # returns voxel value at x=1, y=1, z=2
header(nim)  # returns header attributes as a list}

