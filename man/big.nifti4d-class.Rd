\name{big.nifti4d-class}
\alias{big.nifti4d-class}
\alias{big.nifti4d}
\title{4d nifti => 2d big matrix}
\description{Represents a 4D nifti image as a 2D big matrix}
\details{Assuming that the 4D image is akin to 3D+time, then the columns represent the voxels in 3D while the rows represent the different timepoints.}
\section{Slots}{\describe{\item{\code{header}:}{}\item{\code{mask}:}{}\item{\code{address}:}{}}}
\section{Extends}{\code{\linkS4class{header}}, \code{\linkS4class{big.niftiXd}}}

