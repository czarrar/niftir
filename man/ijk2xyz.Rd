\name{ijk2xyz}
\alias{ijk2xyz}
\title{Converting between standard and image coordinates}
\usage{ijk2xyz(header, ijk, use.qform=T)
xyz2ijk(header, xyz, use.qform=T)}
\description{Transform your array coordinates (ijk) into standard space coordinates (xyz), and vice-versa}
\alias{xyz2ijk}
\author{Zarrar Shehzad}
\value{transformed xyz => ijk coordinates or vice-versa (so a vector of length 3)}
\arguments{\item{header}{can be \code{list}, \code{nifti}, or \code{niftiXd} object}
\item{ijk}{coordinates in your array}
\item{xyz}{coordinates in standard space}
\item{use.qform}{if \code{TRUE} will use qform transfomation, otherwise will use sform}
}

