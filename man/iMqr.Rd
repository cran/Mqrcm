\name{iMqr}
\alias{iMqr}
\title{
M-Quantile Regression Coefficients Modeling
}
\description{
This function implements Frumento and Salvati's (2018) method for M-quantile regression
coefficients modeling (Mqrcm). M-quantile regression coefficients are described
by parametric functions of the order of the quantile.
}
\usage{
iMqr(formula, formula.p = ~ slp(p,3), weights, data, s, 
  psi = "Huber", plim = c(0,1), tol = 1e-6, maxit)
}
\arguments{
  \item{formula}{
    a two-sided formula of the form \code{y ~ x1 + x2 + \ldots}: 
    a symbolic description of the M-quantile regression model.
  }
  \item{formula.p}{
    a one-sided formula of the form \code{~ b1(p, \ldots) + b2(p, \ldots) + \ldots}, 
    describing how M-quantile regression coefficients depend on \kbd{p}, 
    the order of the quantile.
  }
  \item{weights}{
    an optional vector of weights to be used in the fitting process. 
    The weights will always be normalized to sum to the sample size. 
    This implies that, for example, using double weights will \emph{not} 
    halve the standard errors.}
  \item{data}{
    an optional data frame, list or environment containing the variables in \code{formula}.
  }
  \item{s}{an optional 0/1 matrix that permits excluding some model coefficients
  (see \sQuote{Examples}).}
  \item{psi}{a character string indicating the \sQuote{psi} function. Currently,
  only \sQuote{\code{Huber}} is implemented.}
  \item{plim}{the extremes of the estimation interval. You may want to model the M-quantile
  regression coefficients in an interval, say, \code{(a,b)} instead of \code{(0,1)}.}
  \item{tol}{
    convergence criterion for numerical optimization.
  }
  \item{maxit}{
    maximum number of iterations.
  }
}
\details{

  A linear model is used to describe the \code{p}-th conditional M-quantile:
  \deqn{M(p | x) = \beta_0(p) + \beta_1(p)x_1 + \beta_2(p)x_2 + \ldots.}{%
  Q(p | x) = \beta0(p) + \beta1(p)*x1 + \beta2(p)*x2 + \ldots.}

  Assume that each M-quantile regression coefficient can be expressed as a parametric function of \eqn{p} of the form:
  \deqn{\beta(p | \theta) = \theta_{0} + \theta_1 b_1(p) + \theta_2 b_2(p) + \ldots}{%
  \beta(p | \theta) = \theta0 + \theta1*b1(p) + \theta2*b2(p) + \ldots}
  where \eqn{b_1(p), b_2(p, \ldots)}{b1(p), b2(p), \ldots} are known functions of \eqn{p}. 
  If \eqn{q} is the dimension of 
  \eqn{x = (1, x_1, x_2, \ldots)}{x = (1, x1, x2, \ldots)} 
  and \eqn{k} is that of  
  \eqn{b(p) = (1, b_1(p), b_2(p), \ldots)}{b(p) = (1, b1(p), b2(p), \ldots)}, 
  the entire M-conditional quantile function is described by a 
  \eqn{q \times k}{q*k} matrix \eqn{\theta} of model parameters. 
  
  
  Users are required to specify two formulas: \code{formula} describes the regression model,
  while \code{formula.p} identifies the 'basis' \eqn{b(p)}. 
  By default, \code{formula.p = ~ slp(p, k = 3)}, a 3rd-degree shifted 
  Legendre polynomial (see \code{\link{slp}}). Any user-defined function \eqn{b(p, \ldots)}
  can be used, see \sQuote{Examples}.
  
  Estimation of \eqn{\theta} is carried out by minimizing an integrated loss function, 
  corresponding to the
  integral, over \eqn{p}, of the loss function of standard M-quantile regression. This
  motivates the acronym \code{iMqr} (integrated M-quantile regression). The scale parameter
  \code{sigma} is estimated as the minimizer of the log-likelihood of a Generalized 
  Asymmetric Least Informative distribution (Bianchi et al 2017), and is \dQuote{modeled}
  as a piecewise-constant function of the order of the quantile.
  
}
\value{
An object of class \dQuote{\code{iMqr}}, a list containing the following items:
\item{coefficients}{a matrix of estimated model parameters describing the fitted M-quantile function.}
\item{plim}{a vector of two elements indicating the range of estimation.}
\item{call}{the matched call.}
\item{converged}{logical. The convergence status.}
\item{n.it}{the number of iterations.}
\item{obj.function}{the value of the minimized integrated loss function, which does not include the terms that are
only associated with the scale parameter \code{sigma}.}
\item{s}{the used \sQuote{\kbd{s}} matrix.}
\item{psi}{the used \sQuote{\code{psi}} function.}
\item{covar}{the estimated covariance matrix.}
\item{mf}{the model frame used.}
\item{PDF, CDF}{the fitted values of the conditional probability density function (\acronym{PDF})
and cumulative distribution function (\acronym{CDF}). The \acronym{CDF} value should be interpreted as the order
of the M-quantile that corresponds to the observed \code{y} variable, while the \acronym{PDF} is just the first
derivative of the \acronym{CDF}.}

Use \code{\link{summary.iMqr}}, \code{\link{plot.iMqr}}, and \code{\link{predict.iMqr}}
for summary information, plotting, and predictions from the fitted model.
The generic accessory functions \code{coefficients}, \code{formula}, \code{terms}, 
\code{model.matrix}, \code{vcov} are available to extract information from the fitted model.
}
\references{
Frumento, P., Salvati, N. (2018). \emph{Parametric modeling of M-quantile regression coefficient functions with application to small area estimation} [forthcoming].

Bianchi, et al. (2018). \emph{Estimation and testing in M-quantile regression with application to small area estimation}, International Statistical Review, 0(0), p. 1-30.
}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}

\seealso{
\code{\link{summary.iMqr}}, \code{\link{plot.iMqr}}, \code{\link{predict.iMqr}},
for summary, plotting, and prediction, and \code{\link{plf}} and \code{\link{slp}} 
that may be used to define \eqn{b(p)}
to be a piecewise linear function and a shifted Legendre polynomial basis, respectively.
}
\examples{

  ##### Using simulated data in all examples
  ##### NOTE 1: the true quantile and M-quantile functions do not generally coincide
  ##### NOTE 2: the true M-quantile function is usually unknown, even with simulated data
  
  
  ##### Example 1
  
  n <- 250
  x <- runif(n)
  y <- rnorm(n, 1 + x, 1 + x)
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with 
    # beta0(p) = beta1(p) = 1 + qnorm(p)
                              
  # fit the 'true' model: b(p) = (1 , qnorm(p))
  m1 <- iMqr(y ~ x, formula.p = ~ I(qnorm(p)))
  # the fitted M-quantile regression coefficient functions are
    # beta0(p) = m1$coef[1,1] + m1$coef[1,2]*qnorm(p)
    # beta1(p) = m1$coef[2,1] + m1$coef[2,2]*qnorm(p)
  \donttest{
  # a basis b(p) = (1, p), i.e., beta(p) is assumed to be a linear function of p
  m2 <- iMqr(y ~ x, formula.p = ~ p)

  # a 'rich' basis b(p) = (1, p, p^2, log(p), log(1 - p))
  m3 <- iMqr(y ~ x, formula.p = ~ p + I(p^2) + I(log(p)) + I(log(1 - p)))

  # 'slp' creates an orthogonal spline basis using shifted Legendre polynomials
  m4 <- iMqr(y ~ x, formula.p = ~ slp(p, k = 3)) # note that this is the default
  
  # 'plf' creates the basis of a piecewise linear function
  m5 <- iMqr(y ~ x, formula.p = ~ plf(p, knots = c(0.1,0.9)))
  }
  
  summary(m1)
  summary(m1, p = c(0.25,0.5,0.75))
  par(mfrow = c(1,2)); plot(m1, ask = FALSE)
  # see the documentation for 'summary.iMqr' and 'plot.iMqr'

  

  \donttest{
  ##### Example 2 ### excluding coefficients
  
  n <- 250
  x <- runif(n)
  qy <- function(p,x){(1 + qnorm(p)) + (1 + log(p))*x}
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
    # beta0(p) = 1 + qnorm(p) 
    # beta1(p) = 1 + log(p)
  
  y <- qy(runif(n), x) # to generate y, plug uniform p in qy(p,x) 
  iMqr(y ~ x, formula.p = ~ I(qnorm(p)) + I(log(p)))

  # I would like to exclude log(p) from beta0(p), and qnorm(p) from beta1(p)
  # I set to 0 the corresponding entries of 's'

  s <- matrix(1,2,3); s[1,3] <- s[2,2] <- 0
  iMqr(y ~ x, formula.p = ~ I(qnorm(p)) + I(log(p)), s = s)
  }
}
\keyword{models}
\keyword{regression}
