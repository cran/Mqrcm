\name{psi}
\alias{psi}
\alias{Huber}
\title{
Generate Various Influence Functions for M-Quantiles
}
\description{
Influence function to be passed to \code{\link{iMqr}}.
}
\usage{
Huber(c = 1.345)
}
\arguments{
  \item{c}{
    tuning parameter for Huber's influence function.
  }
}
\details{
  These functions are only meant to be used used within a call to \code{\link{iMqr}}.
}
\value{
A list with the following items:
\item{psi, psi_tau, psi1_tau, rho_tau}{define the influence function.}
\item{par}{the parameters of the influence function, e.g., the value of \code{c} in Huber's function.}
\item{name}{a character string indicating the name of the influence function.}
}
\references{
Huber, P. J. (1981). "Robust Statistics", John Wiley and Sons, New York.
}
\seealso{
\code{\link{iMqr}}
}
\examples{

  # The following are identical:
  # iMqr(y ~ x, psi = "Huber")
  # iMqr(y ~ x, psi = Huber)
  # iMqr(y ~ x, psi = Huber(c = 1.345))

}
\keyword{models}

