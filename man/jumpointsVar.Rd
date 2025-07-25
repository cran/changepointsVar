\name{jumpointsVar}
\alias{jumpointsVar}
\title{
Change-point detection for variance piecewise constant models
}
\description{
Detection of change-points for variance of
heteroscedastic Gaussian variables  with piecewise constant
variance function.
}
\usage{
jumpointsVar(y, x, y.res=FALSE, k=min(30, round(length(y)/10)), print.level=0,
             plot.it=FALSE, psi=NULL, round=TRUE, control=fit.control(),
             selection=sel.control())
}
\arguments{
  \item{y}{
the response variable.
}
  \item{x}{
	the 'segmented' variable; if missing simple indices 1,2,... are assumed.
}
\item{y.res}{
	logical:
if  \code{FALSE}, \code{y}  is the observed sequence with piecewise constant variance function;
if \code{TRUE}, \code{y} is the vector of the squared residuals from a fitted linear
model that accounts for the mean behavior of the observed signal with changes in variation; see Details
}
  \item{k}{
the starting number of changepoints. It should be quite larger than the supposed number of (true) changepoints. This argument is ignored if starting values of the changepoints are specified via psi.
}
  \item{print.level}{
the default value is 0 indicating nothing is printed; 1 allows to print some informations during the algorithm; 2 the same as for 1 plus lars information.
}
  \item{plot.it}{
logical indicating if the curve of the criterion choosen has to be displayed.
}
   \item{psi}{
	numeric vector to indicate the starting values for the changepoints. When \code{psi=NULL} (default), \code{k} quantiles are assumed
}
  \item{round}{
	logical: if the change-point values should be rounded
}
  \item{control}{
a list returned by \code{fit.control}
}
  \item{selection}{
a list returned by \code{sel.control}
}
}
\details{
This algorithm allows breakpoint detections for changes in variation assuming that the variance function can be described by a piecewise constant function with segments delimited by unknown change-points.

The approach is a generalization of the \code{cumSeg} procedure proposed by Muggeo and Adelfio (2011) assuming that testing for stepwise changes in variance of a sequence of Gaussian random variables may be transformed equivalently to the case of testing for changes in mean of the squared residuals (from an estimated linear model that accounts for the mean behavior of the observed signal) assuming a gamma GLM with a log-link function.

A variation of \code{\link[lars]{lars}} procedure adapted to the GLM case is considered to discard the spurious change-points on the basis of a generalized version of the BIC.

The proposed approach results in a very efficient algorithm even with n large and many change-points to be estimated.
}
\value{
An object of class \code{jumpointsVar}. It's a list including several components:
 \item{psi}{the estimated changepoints}
 \item{est.means}{the estimated means}
 \item{n.psi}{the estimated number of changepoints}
 \item{psi0}{the initial estimated changepoints (before applying the selection criterion)}
 \item{est.means0}{the initial estimated means (before applying the selection criterion)}
 \item{criterion}{the curve of the selected criterion}
 \item{fitted.values}{the fitted values}
 \item{input}{the input parameters}
 \item{call}{the call function}
}
\references{
Adelfio, G. (2012), Change-point detection for variance piecewise constant models, \emph{Communications in Statistics, Simulation and Computation}, 41:4, 437-448

Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection for genomic sequences of continuous measurements, \emph{Bioinformatics} 27, 161-166.
}
\author{
Giada Adelfio

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\seealso{
 \code{\link{plot.jumpointsVar}}
}
\examples{
set.seed(2)
n = 500
x = 1:n/n
mu = 10+6*sin(3*pi*x)
# if there are two change-points
sigma = c(rep(0.5,.2*n), rep(8,.4*n), rep(3,.4*n))
y = mu + rnorm(n, 0, sigma)
reg = lm(y ~ mu)
h = influence(reg)$hat
r2 = resid(reg)^2/(1-h) + 1
o = jumpointsVar(y=r2, y.res=TRUE, k=30)
o
plot(o)

# if there are no change-points
sigma = 0.5*x
y = mu + rnorm(n, 0, sigma)
reg = lm(y ~ mu)
h = influence(reg)$hat
r2 = resid(reg)^2/(1-h)+1
o = jumpointsVar(y=r2, y.res=TRUE, k=30)
plot(o)

# if the mean behavior of the observed signal is unkown
sigma = c(rep(0.5,.2*n), rep(8,.4*n), rep(3,.4*n))
y = 10 + rnorm(n, 0, sigma)
o = jumpointsVar(y=y, y.res=FALSE, k=30)
plot(o)

}
\keyword{changepoints}

