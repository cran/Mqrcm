\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only, or not yet documented.
}
\usage{
check.in(mf, formula.p, s, plim)
check.out(theta, S, covar)
start.theta(y,x, weights, bfun, df, yy, s)
iMqr.internal(mf,cl, formula.p, tol = 1e-6, maxit, s, psi, plim)
cov.theta(fit, y,X,Xw, weights, bfun, psi, s, sigma)
iobjfun(theta, y,X,Xw,weights, bfun, psi, sigma, u, H = FALSE, i = FALSE)
sortH(H)
isigma(theta, y,X,weights,bfun,psi)
iMqr.newton(theta, y,X,Xw, weights, bfun,psi, s,sigma, 
                        tol, maxit, safeit, eps0)

num.fun(dx,fx, op = c("int", "der"))
make.bfun(p, x)
apply_bfun(bfun, p, fun = c("bfun", "b1fun"))
p.bisec(theta, y, X, bfun, n.it = 20)
slp.basis(k, intercept)
is.slp(f)

iqr.waldtest(obj)
extract.p(model, p, cov = FALSE)
pred.beta(model, p, se = FALSE)

\method{print}{iMqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{print}{summary.iMqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{terms}{iMqr}(x, \ldots)
\method{model.matrix}{iMqr}(object, \ldots)
\method{vcov}{iMqr}(object, \ldots)
\method{nobs}{iMqr}(object, \ldots)
}
\keyword{internal}

