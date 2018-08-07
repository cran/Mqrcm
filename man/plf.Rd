\name{plf}
\alias{plf}
\title{
Basis of a Piecewise Linear Function
}
\description{
Generates \eqn{b_1(p), b_2(p), \ldots}{b1(p), b2(p), \ldots} such that, for \kbd{0 < p < 1},
\deqn{\theta_1*b_1(p) + \theta_2*b_2(p) + \ldots}{\theta1*b1(p) + \theta2*b2(p) + \ldots}
is a piecewise linear function with slopes 
\eqn{(\theta_1, \theta_2, \ldots)}{\theta1, \theta2, \ldots}.
}
\usage{
plf(p, knots)
}
\arguments{
  \item{p}{
    a numeric vector of values between \kbd{0} and \kbd{1}.
  }
  \item{knots}{
    a set of \emph{internal} knots between \kbd{0} and \kbd{1}.
    It can be \kbd{NULL} for no internal knots. 
}
}
\details{
  This function permits computing a piecewise linear function on the unit interval. 
  A different slope holds between each pair of knots, and the function
  is continuous at the knots.
}
\value{
A matrix with one row for each element of \kbd{p}, and \code{length(knots) + 1} columns. 
The knots are returned as \code{attr(, "knots")}.
Any linear combination of the basis matrix is a piecewise linear function where
each coefficient represents the slope in the corresponding sub-interval (see \sQuote{Examples}). 
}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\note{
This function is typically used within a call to \code{\link{iMqr}}.
A piecewise linear function can be used to describe how M-quantile regression coefficients
depend on the order of the quantile.
}
\seealso{
\code{\link{slp}}, for shifted Legendre polynomials.
}
\examples{
  p <- seq(0,1, 0.1)

  a1 <- plf(p, knots = NULL) # returns p

  a2 <- plf(p, knots = c(0.2,0.7))
  plot(p, 3 + 1*a2[,1] - 1*a2[,2] + 2*a2[,3], type = "l") 
    # intercept = 3; slopes = (1,-1,2)
}
\keyword{array}
