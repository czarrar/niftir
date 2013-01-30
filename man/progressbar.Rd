\name{progressbar}
\alias{progressbar}
\title{Progress Bar}
\usage{progressbar(limit)
}
\description{Wrapper around txtProgressBar}
\details{Will start the progress bar at 0 and go up to \code{limit} specified.
Uses the 3rd style in \code{\link{txtProgressBar}}.}
\value{\code{txtProgressBar} object}
\seealso{\code{\link{update.txtProgressBar}} \code{\link{end.txtProgressBar}}}
\arguments{\item{limit}{maximum point of the progress bar (integer)}
}
\examples{n <- 10
pb <- progressbar(10)
for (i in 1:n) {
Sys.sleep(1)
update(pb, i)
}
end(pb)}

