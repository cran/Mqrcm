\name{slp}
\alias{slp}
\title{
Shifted Legendre Polynomials
}
\description{
Computes shifted Legendre polynomials.
}
\usage{
slp(p, k = 3, intercept = FALSE)
}
\arguments{
  \item{p}{
  the variable for which to compute the polynomials. Must be \kbd{0 <= p <= 1}.
}
  \item{k}{
  the degree of the polynomial.
}
  \item{intercept}{
  logical. If \kbd{TRUE}, the polynomials include the constant term.
}
}
\details{
Shifted Legendre polynomials (\acronym{SLP}) are orthogonal polynomial functions in (0,1) that can be used
to build a spline basis, typically within a call to \code{\link{iMqr}}.
The constant term is omitted unless \kbd{intercept = TRUE}: for example,
the first two \acronym{SLP} are \code{(2*p - 1, 6*p^2 - 6*p + 1)}, 
but \code{slp(p, k = 2)} will only return \code{(2*p, 6*p^2 - 6*p)}.
}
\value{
An object of class \dQuote{\code{slp}}, i.e., 
a matrix with the same number of rows as \kbd{p}, and with \kbd{k} columns
named \code{slp1, slp2, \ldots} containing the \acronym{SLP} of the corresponding orders.
The value of \kbd{k} is reported as attribute.
}
\references{
Refaat El Attar (2009), \emph{Legendre Polynomials and Functions}, CreateSpace, ISBN 978-1-4414-9012-4.
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\note{
The default for \command{iMqr} is \code{formula.p = ~ slp(p, k = 3)}.
}


\seealso{
\code{\link{plf}}, for piecewise linear functions in the unit interval.
}
\examples{
  p <- seq(0,1,0.1)
  slp(p, k = 1) # = 2*p
  slp(p, k = 1, intercept = TRUE) # = 2*p - 1 (this is the true SLP of order 1)
  slp(p, k = 3) # a linear combination of (p, p^2, p^3), with slp(0,k) = 0
}
\keyword{smooth}

