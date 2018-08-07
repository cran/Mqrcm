\name{plot.iMqr}
\alias{plot.iMqr}
\title{
Plot M-Quantile Regression Coefficients
}
\description{
Plots M-quantile regression coefficients 
\eqn{\beta(p)} as a function of \eqn{p},
based on a fitted model of class \dQuote{\code{iMqr}}.
}
\usage{
\method{plot}{iMqr}(x, conf.int = TRUE, polygon = TRUE, which = NULL, ask = TRUE, \ldots)
}
\arguments{
  \item{x}{
an object of class \dQuote{\code{iMqr}}, typically the result of a call to \code{\link{iMqr}}.
}
  \item{conf.int}{
  logical. If \kbd{TRUE}, asymptotic 95\% confidence intervals are added to the plot.
}
  \item{polygon}{
  logical. If \kbd{TRUE}, confidence intervals are represented by shaded areas via \code{polygon}.
  Otherwise, dashed lines are used.
}
  \item{which}{
  an optional numerical vector indicating which coefficient(s) to plot. If \kbd{which = NULL}, 
  all coefficients are plotted.
}
\item{ask}{logical. If \kbd{which = NULL} and \kbd{ask = TRUE} (the default), you will be asked
interactively which coefficients to plot.}
\item{\ldots}{additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd, cex.lab, cex.axis, axes, frame.plot}. 
See \code{\link{par}}.}
}
\details{
Using \command{iMqr}, each M-quantile regression coefficient \eqn{\beta(p)} is described by a linear
combination of known parametric functions of \eqn{p}. With this command, a plot of
\eqn{\beta(p)} versus \eqn{p} is created.
}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\seealso{
\code{\link{iMqr}} for model fitting; \code{\link{summary.iMqr}} and \code{\link{predict.iMqr}}
for model summary and prediction.
}
\examples{
  \donttest{
  # using simulated data
  
  n <- 250
  x <- runif(n)
  qy <- function(p,x){p^2 + x*log(p)}
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
     # beta0(p) = p^2
     # beta1(p) = log(p)
  y <- qy(runif(n), x) # to generate y, plug uniform p in qy(p,x)
  
  par(mfrow = c(1,2))
  plot(iMqr(y ~ x, formula.p = ~ slp(p,3)), ask = FALSE) 
  # flexible fit with shifted Legendre polynomials
  }
}
\keyword{methods}

