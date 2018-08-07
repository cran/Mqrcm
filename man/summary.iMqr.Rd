\name{summary.iMqr}
\alias{summary.iMqr}
\title{
Summary After M-Quantile Regression Coefficients Modeling
}
\description{
Summary of an object of class \dQuote{\code{iMqr}}.
}
\usage{
\method{summary}{iMqr}(object, p, cov = FALSE, \ldots)
}
\arguments{
  \item{object}{
an object of class \dQuote{\code{iMqr}}, the result of a call to \code{\link{iMqr}}.
}
  \item{p}{
  an optional vector of quantiles.
}
  \item{cov}{
  logical. If \kbd{TRUE}, the covariance matrix of \eqn{\beta(p)} is reported.
  Ignored if \kbd{p} is missing.
}
\item{\ldots}{for future methods.}
}
\details{
If \kbd{p} is missing, a summary of the fitted model is reported. This includes the
estimated coefficients, their standard errors, and other summaries (see \sQuote{Value}).
If \kbd{p} is supplied, the M-quantile regression coefficients of order \kbd{p} are extrapolated
and summarized.
}
\value{

If \kbd{p} is supplied, 
a standard summary of the estimated M-quantile regression coefficients 
is returned for each value of \kbd{p}. If \kbd{cov = TRUE}, the covariance
matrix is also reported.

If \kbd{p} is missing (the default), a list with the following items:
\item{converged}{logical value indicating the convergence status.}
\item{n.it}{the number of iterations.}
\item{n}{the number of observations.}
\item{free.par}{the number of free parameters in the model.}
\item{coefficients}{the matrix of estimated coefficients. Each row corresponds to 
a covariate, while each column corresponds to an element of \eqn{b(p)},
the set of functions that describe how M-quantile regression coefficients vary with the 
order of the quantile. See \sQuote{Examples}.}
\item{se}{the estimated standard errors.}
\item{test.x}{Wald test for the covariates. Each \emph{row} of \code{coefficients} is tested for nullity.}
\item{test.p}{Wald test for the building blocks of the quantile function. Each \emph{column} of \code{coefficients} is tested for nullity.}
\item{obj.function}{the minimized loss function.}
\item{call}{the matched call.}
}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\seealso{
\code{\link{iMqr}}, for model fitting; \code{\link{predict.iMqr}} and \code{\link{plot.iMqr}},
for predicting and plotting objects of class \dQuote{\code{iMqr}}. 
}

\examples{

# using simulated data
\donttest{
set.seed(1234); n <- 250
x1 <- rexp(n)
x2 <- runif(n)
qy <- function(p,x){qnorm(p)*(1 + x)}
# true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
   # beta0(p) = beta1(p) = qnorm(p)

y <- qy(runif(n), x1) # to generate y, plug uniform p in qy(p,x)
                      # note that x2 does not enter

model <- iMqr(y ~ x1 + x2, formula.p = ~ I(qnorm(p)) + p + I(p^2))
# beta(p) is modeled by linear combinations of b(p) = (1, qnorm(p),p,p^2)

summary(model)
# interpretation: 
  # beta0(p) = model$coef[1,]*b(p)
  # beta1(p) = model$coef[2,]*b(p); etc.
# x2 and (p, p^2) are not significant


summary(model, p = c(0.25, 0.75)) # summary of beta(p) at selected quantiles
}
}
\keyword{methods}
