\name{rmext}
\alias{rmext}
\title{Remove extension from file...}
\usage{rmext(fname, compressed=T)
}
\description{Remove extension from file}
\value{character with extension removed}
\arguments{\item{fname}{filename}
\item{compressed}{do you want to remove a compression extension (.gz or .zip),
before removing the real extension?}
}
\examples{rmext('test.nii.gz')    # should get test
rmext('test.nii.gz', FALSE) # should get test.nii}

