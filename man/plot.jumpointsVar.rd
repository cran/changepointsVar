\name{plot.jumpointsVar}
\alias{plot.jumpointsVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for changes in variance
}
\description{
Plots signal with changes in variance and corresponding changepoints
}
\usage{
\method{plot}{jumpointsVar}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object returned by jumpointsVar
}
  \item{\dots}{
additional arguments.
}
}
\details{
This fuction takes a fitted object returned by jumpointsVar and plots the resulting fit with changepoints.
}
\value{
The function simply plot the fit returned by 'jumpointsVar'
}
\references{
Adelfio, G. (2012), Change-point detection for variance piecewise constant models,\emph{Communications in Statistics, Simulation and Computation}, 41:4, 437-448

Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection for genomic sequences of continuous measurements, \emph{Bioinformatics} 27, 161-166.
}
\author{
Giada Adelfio

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{jumpointsVar}}
}
\examples{
##---- see jumpointsVar documentation ----
}
\keyword{ changepoints }
\keyword{ jumpointsVar }
