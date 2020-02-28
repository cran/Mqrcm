\name{Mqrcm-package}
\alias{Mqrcm-package}
\docType{package}
\title{
M-Quantile Regression Coefficients Modeling
}
\description{
This package implements Frumento and Salvati (2018) method for M-quantile regression
coefficients modeling (Mqrcm), in which M-quantile regression coefficients are described
by (flexible) parametric functions of the order of the quantile. This permits modeling
the entire conditional M-quantile function of a response variable. 
}
\details{
\tabular{ll}{
Package: \tab Mqrcm\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2020-02-26\cr
License: \tab GPL-2\cr
}
The function \code{\link{iMqr}} permits specifying the regression model.
Two special functions, \code{\link{slp}} and \code{\link{plf}}, are provided to facilitate model building.
The auxiliary functions \code{\link{summary.iMqr}}, \code{\link{predict.iMqr}}, and \code{\link{plot.iMqr}}
can be used to extract information from the fitted model.
}
\author{
Paolo Frumento

Maintainer: Paolo Frumento <paolo.frumento@unipi.it>
}
\references{
Frumento, P., and Salvati, N. (2020). \emph{Parametric modelling of M-quantile regression coefficient functions with application to small area estimation}. Journal of the Royal Statistical Society, Series A, 183(1), 229-250.
}

\keyword{ package }
\examples{

# use simulated data

n <- 250
x <- rexp(n)
y <- runif(n, 0, 1 + x)
model <- iMqr(y ~ x, formula.p = ~ p + I(p^2))
summary(model)
summary(model, p = c(0.1,0.2,0.3))
predict(model, type = "beta", p = c(0.1,0.2,0.3))
predict(model, type = "CDF", newdata = data.frame(x = c(1,2,3), y = c(0.5,1,2)))
predict(model, type = "QF", p = c(0.1,0.2,0.3), newdata = data.frame(x = c(1,2,3)))
predict(model, type = "sim", newdata = data.frame(x = c(1,2,3)))
par(mfrow = c(1,2)); plot(model, ask = FALSE)
}



