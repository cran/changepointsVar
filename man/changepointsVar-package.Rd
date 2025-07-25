\docType{package}
\name{changepointsVar-package}
\alias{changepointsVar-package}
\title{
Change-Points Detections for Changes in Variance
}
\description{
This algorithm allows breakpoint detections for changes in variation assuming that the variance function can be described by a piecewise constant function with segments delimited by unknown change-points. The approach is a generalization of the \code{cumSeg} procedure proposed by Muggeo and Adelfio (2011) assuming that testing for stepwise changes in variance of a sequence of Gaussian random variables may be transformed equivalently to the case of testing for changes in mean of the squared residuals (from an estimated linear model that accounts for the mean behavior of the observed signal) assuming a gamma GLM with a log-link function. A variation of \code{\link[lars]{lars}} procedure adapted to the GLM case is considered to discard the spurious change-points on the basis of a generalized version of the BIC. The proposed approach results in a very efficient algorithm even with n large and many change-points to be estimated. Adelfio, G. (2012), Change-point detection for variance piecewise constant models, \emph{Communications in Statistics, Simulation and Computation}, 41:4, 437-448. Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection for genomic sequences of continuous measurements, \emph{Bioinformatics} 27, 161-166.
}
\details{
\tabular{ll}{
Package: \tab jumpointsVar\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2025-07-25\cr
License: \tab GPL-2\cr
}
The function \code{\link{jumpointsVar}} allows to specify the variable in which to look for change-point in variance. The auxiliary function \code{\link{plot.jumpointsVar}} can be used to plot the brekpoints detected from the main algorithm.
}
\author{
Giada Adelfio

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\references{
Adelfio, G. (2012), Change-point detection for variance piecewise constant models, \emph{Communications in Statistics, Simulation and Computation}, 41:4, 437-448

Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection for genomic sequences of continuous measurements, \emph{Bioinformatics} 27, 161-166.}
\keyword{ package }
\examples{
##---- see jumpointsVar documentation ----
}



