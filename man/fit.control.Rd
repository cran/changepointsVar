\name{fit.control}
\alias{fit.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function for controlling model fitting
}
\description{
Auxiliary function as user interface for model fitting. Typically
     only used when calling 'jumpointsVar'
     }
\usage{
fit.control(toll=0.001, it.max=10, last=TRUE, scale.res=FALSE,
            maxit.glm=30, h=1, stop.if.error=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toll}{
  positive convergence tolerance.
}
  \item{it.max}{
integer giving the maximal number of iterations.
}
  \item{last}{
Currently ignored.
}
  \item{scale.res}{
logical indicating if the residuals have to be scaled.
}
  \item{maxit.glm}{
Currently ignored.
}
  \item{h}{
Currently ignored.
}
  \item{stop.if.error}{
logical indicating if the algorithm should stop when one or more estimated changepoints
do not assume admissible values. Default is \code{FALSE} which implies automatic
changepoint selection.

}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A list with the arguments as components to be used by 'jumpointsVar'.

}
\author{
Gianluca Sottile
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{\code{\link{jumpointsVar}}}

\keyword{ regression }
