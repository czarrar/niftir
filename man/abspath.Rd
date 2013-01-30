\name{abspath}
\alias{abspath}
\title{Convert relative path to absolute path...}
\usage{abspath(path)
}
\description{Convert relative path to absolute path}
\details{This will also deal with expanding a path as in \code{\link{path.expand}}.
You must ensure that the directory to your path exists, otherwise you will
get an error.}
\value{character with absolute path}
\arguments{\item{path}{character specifying path in your file system}
}
\examples{abspath('~/test.nii.gz')
curdir <- getwd()
abspath('../test.nii.gz')
file.path(dirname(curdir), 'test.nii.gz')    # equivalent to above
abspath('test.nii.gz')
file.path(curdir, 'test.nii.gz')    # equivalent to above}

