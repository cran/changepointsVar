\name{sel.control}
\alias{sel.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function for controlling model selection
}
\description{
    Auxiliary function as user interface for model selection. Typically
     only used when calling 'jumpointsVar'
}
\usage{
sel.control(type=c("bic", "rss"), S=1, Cn="2*log(log(n))",
            alg=c("lasso", "stepwise"), edf.psi=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
the criterion to be used to perform model selection.
}
  \item{S}{
if \code{type="rss"} the optimal model is selected when the residual sum of squares decreases
by the threshold \code{S}.
}
  \item{Cn}{
if \code{type="bic"} a character string (as a function of \code{'n'}) to specify to generalized BIC.
If \code{Cn=1} the standard BIC is used.
}
  \item{alg}{
which procedure should be used to perform model selection? The value of \code{alg} is
passed to the argument \code{'type'} of \code{lars}.
}
  \item{edf.psi}{
logical indicating if the number of changepoints should be computed in the model df.
}
}
\details{
This function specifies how to perform model seletion, namely how many change points should be selected.
}
\value{
A list with the arguments as components to be used by 'jumpointsVar' and in turn by 'lars'.
}
\author{
Gianluca Sottile
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{\code{\link{jumpointsVar}}, \code{\link[lars]{lars}}}

\keyword{ regression }
