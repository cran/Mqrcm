#' @importFrom stats integrate splinefun model.response model.weights model.matrix terms model.frame delete.response coef pnorm qnorm qchisq
#' @importFrom stats approxfun sd prcomp lm.wfit pchisq printCoefmat .getXlevels pchisq runif vcov nobs predict pbeta qbeta qexp
#' @importFrom graphics plot points abline polygon
#' @importFrom grDevices adjustcolor
#' @importFrom utils menu setTxtProgressBar txtProgressBar tail getFromNamespace
#' @importFrom Hmisc wtd.quantile
#' @import pch

#' @export
iMqr <- function(formula, formula.p = ~ slp(p,3), weights, data, s, psi = "Huber", plim = c(0,1), tol = 1e-6, maxit){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	iMqr.internal(mf = mf,cl = cl, formula.p = formula.p, tol = tol, maxit = maxit, s = s, psi = psi, plim = plim)
}

# IMPORTANT: in the 'psi' functions, tau is vector of length d, and u is a n*d matrix.
# rho_tau = loss, psi_tau = derivative of rho_tau, psi1_tau = derivative of psi_tau.
#' @export
Huber <- function(c = 1.345){
  
  psi <- function(u,c = 1.345){
    i <- (abs(u) <= c)
    u*i + c*sign(u)*(1 - i)
  }
  
  psi_tau <- function(u, tau, c = 1.345){
    omega <- t(u <= 0)
    2*psi(u,c)*t(tau*(1 - omega) + (1 - tau)*omega)
  }
  
  psi1_tau <- function(u, tau, c = 1.345){
    omega <- t(u <= 0)
    2*(abs(u) <= c)*t(tau*(1 - omega) + (1 - tau)*omega)
  }
  
  rho_tau <- function(u,tau, c = 1.345){ 
    U <- abs(u)
    i <- (U <= c)
    omega <- (u <= 0)
    h <- t(abs(tau - t(omega)))
    2*(((u^2)/2*h)*i + (c*U - (c^2)/2)*h*(1 - i))
  }	
  
  list(psi = psi, psi_tau = psi_tau, psi1_tau = psi1_tau, rho_tau = rho_tau, par = c, name = "Huber")
}


check.in <- function(mf, formula.p, s, plim){

	if(!missing(s) && all(s == 0)){stop("'s' cannot be all zero")}
	explore.s <- function(s, dim){
		if(dim == 2){s <- t(s)}
		out <- 1
		if((r <- nrow(s)) > 1){
			for(j in 2:r){
				done <- FALSE; rj <- s[j,]
				for(h in 1:(j - 1)){
					if(all(rj == s[h,])){out[j] <- out[h]; done <- TRUE}
				}
				if(!done){out[j] <- max(out) + 1}
			}
		}
		out
	}

	# weights

	if(any((weights <- model.weights(mf)) < 0)){stop("negative 'weights'")}
	if(is.null(weights)){weights <- rep.int(1, nrow(mf)); alarm <- FALSE}
	else{
	alarm <- (weights == 0)
	  sel <- which(!alarm)
	  mf <- mf[sel,]
	  weights <- weights[sel]
	  weights <- weights/mean(weights)
	}
	if(any(alarm)){warning("observations with null weight will be dropped", call. = FALSE)}
	if((n <- nrow(mf)) == 0){stop("zero non-NA cases", call. = FALSE)}

	# y

	y <- model.response(mf)

	# x and b(p)

	X <- model.matrix(attr(mf, "terms"), mf); q <- ncol(X)
	termlabelsX <- attr(attr(mf, "terms"), "term.labels")
	assignX <- attr(X, "assign")
	coefnamesX <- colnames(X)
	
	# p1 is used to evaluate the splinefuns. A non-evenly spaced grid, with more values on the tails.
	# p2 is for external use (p.bisec). A grid with p reachable by bisection on the p scale.
	# p3 is the grid to evaluate the integrated objective function
	p1 <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = 1000),2,2)
	p2 <- (1:1023)/1024
	
	if(length(plim) != 2 | plim[1] >= plim[2]){stop("wrong specification of 'plim'")}
	d <- 99
	p3 <- seq.int(plim[1], plim[2], length.out = d + 2)[2:(d + 1)]

	if((use.slp <- is.slp(formula.p))){
		k <- attr(use.slp, "k")
		intercept <- attr(use.slp, "intercept") 	# slp(0) = 0?
		intB <- attr(use.slp, "intB")			        # b(p) includes 1?
		assignB <- (1 - intB):k
		termlabelsB <- paste("slp", 1:k, sep = "")
		coefnamesB <- (if(intB) c("(Intercept)", termlabelsB) else termlabelsB)
		k <- k + intB
	}
	else{
	  B <- model.matrix(formula.p, data = data.frame(p = c(p1,p2,p3)))
	  B1 <- B[1:1000,, drop = FALSE]
	  B2 <- B[1001:2023,, drop = FALSE]
	  B3 <- B[2024:(2023 + d),, drop = FALSE]

		k <- ncol(B)
		assignB <- attr(B, "assign")
		termlabelsB <- attr(terms(formula.p), "term.labels")
		coefnamesB <- colnames(B)
	}
	if(missing(s)){s <- matrix(1,q,k)}
	else{
		if(any(dim(s) != c(q,k))){stop("wrong size of 's'")}
		if(any(s != 0 & s != 1)){stop("'s' can only contain 0 and 1")}
	}

	# x singularities (set s = 0 where singularities occur)
	# x is dropped as in a linear model, irrespective of s.

	vx <- qr(X); selx <- vx$pivot[1:vx$rank]
	if(vx$rank < q){s[-selx,] <- 0}

	# b(p) singularities. Dropped row by row, based on s

	if(!use.slp && qr(B3)$rank < k){
		u <- explore.s(s,1)
		for(j in unique(u)){
			sel <- which(s[which(u == j)[1],] == 1)
			if(length(sel) > 1){
				vbj <- qr(B3[,sel, drop = FALSE])
				if((rj <- vbj$rank) < length(sel)){
					s[u == j, sel[-vbj$pivot[1:rj]]] <- 0
				}
			}
		}
	}

	# location-scale statistics for x, b(p), and y

	ry <- range(y); my <- ry[1]; My <- ry[2]

	sX <- apply(X,2,sd); mX <- colMeans(X)
	intX <- (length((constX <- which(sX == 0 & mX != 0))) > 0)
	varsX <- which(sX > 0); zeroX <- which(sX == 0 & mX == 0)
	sX[constX] <- X[1,constX]; mX[constX] <- 0; sX[zeroX] <- 1
	if(length(constX) > 1){zeroX <- c(zeroX, constX[-1]); constX <- constX[1]}

#sX <- sX - sX + 1
#mX <- mX*0	
	
	if(!use.slp){
		sB <- apply(B3,2,sd); mB <- colMeans(B3)
		intB <- (length((constB <- which(sB == 0 & mB != 0))) > 0); varsB <- which(sB > 0)
		if(length(varsB) == 0){stop("the M-quantile function must depend on p")}
		if(length(constB) > 1){stop("remove multiple constant functions from 'formula.p'")}
		if(any(sB == 0 & mB == 0)){stop("remove zero functions from 'formula.p'")}
		sB[constB] <- B3[1,constB]; mB[constB] <- 0
	}
	else{
		sB <- rep.int(1, k); mB <- rep.int(0, k)
		if(intB){constB <- 1; varsB <- 2:k}
		else{constB <- integer(0); varsB <- 1:k}
	}
#sB <- sB - sB + 1
#mB <- mB*0
	if(all(s[, varsB] == 0)){stop("the M-quantile function must depend on p (wrong specification of 's')")}
	if(!(theta00 <- ((intX & intB) && s[constX, constB] == 1)))
		{my <- 0; My <- sd(y)*5; mX <- rep.int(0,q)}
	else{for(j in varsX){if(any(s[j,] > s[constX,])){mX[j] <- 0}}}
	if(!intB | (intB && any(s[,constB] == 0))){mB <- rep.int(0,k)}

	# Create bfun (only used by post-estimation functions)

	if(!use.slp){
		bfun <- list()
		if(intB){bfun[[constB]] <- function(p, deriv = 0){rep.int(1 - deriv, length(p))}}
		for(j in varsB){bfun[[j]] <- make.bfun(p1,B1[,j])}
		names(bfun) <- coefnamesB
		attr(bfun, "k") <- k
	}
	else{
		bfun <- slp.basis(k - intB, intercept)
		if(!intB){bfun$a[1,1] <- bfun$A[1,1] <- bfun$AA[1,1] <- 0}
		attr(bfun, "intB") <- intB
		B2 <- apply_bfun(bfun,p2, "bfun")
	}

	attr(bfun, "bp") <- B2
	attr(bfun, "p") <- p2

	# first scaling of x, b(p), y
#my <- 0
#My <- 10
	U <- list(X = X, y = y)
	X <- scale(X, center = mX, scale = sX)
	y <- (y - my)/(My - my)*10
	if(!use.slp){B3 <- scale(B3, center = mB, scale = sB)}

	# principal component rotations that I can apply to x and b(p); second scaling

	rotX <- (diag(1,q))
	MX <- rep.int(0,q); SX <- rep.int(1,q)
	if(length(varsX) > 0){
		uX <- explore.s(s,1)
		X_in <- rowSums(s)	
		for(j in unique(uX)){
			sel <- which(uX == j)
			if(intX){sel <- sel[sel != constX]}
			if(length(sel) > 1 && X_in[sel[1]] != 0){ # & 1 == 2){ ############ OCCHIOOOOO ####################
				PC <- prcomp(X[,sel], center = FALSE, scale. = FALSE)
				X[,sel] <- PC$x
				rotX[sel,sel] <- PC$rotation
			}
		}
		MX <- colMeans(X); MX[mX == 0] <- 0
		SX <- apply(X,2,sd); SX[constX] <- 1; SX[zeroX] <- 1
#SX <- SX - SX + 1
#MX <- MX*0
		X <- scale(X, center = MX, scale = SX)
	}

	rotB <- (diag(1,k))
	MB <- rep.int(0,k); SB <- rep.int(1,k)
	if(!use.slp){
		uB <- explore.s(s,2)
		B_in <- colSums(s)	
		for(j in unique(uB)){
			sel <- which(uB == j)
			if(intB){sel <- sel[sel != constB]}
			if(length(sel) > 1 && B_in[sel[1]] != 0){ # & 1 == 2){ ############ OCCHIOOOOO
				PC <- prcomp(B3[,sel], center = FALSE, scale. = FALSE)
				B3[,sel] <- PC$x
				rotB[sel,sel] <- PC$rotation
			}
		}
		MB <- colMeans(B3); MB[mB == 0] <- 0
		SB <- apply(B3,2,sd); SB[constB] <- 1
#SB <- SB - SB + 1
#MB <- MB*0
		B3 <- scale(B3, center = MB, scale = SB)
	}

	# Create a pre-evaluated basis (only used internally)

	p <- p3
	
	if(!use.slp){
	  bp <- B3
	  dp <- p[-1] - p[-d]
	  b1p <- num.fun(dp,bp, "der")
	}
	else{
	  k <- attr(bfun, "k")
	  pp <- matrix(, d, k + 1)
	  pp[,1] <- 1; pp[,2] <- p
	  if(k > 1){for(j in 2:k){pp[,j + 1] <- pp[,j]*p}}
	  bp <- tcrossprod(pp, t(bfun$a))
	  b1p <- cbind(0, tcrossprod(pp[,1:k, drop = FALSE], t(bfun$a1[-1,-1, drop = FALSE])))
	    
	  if(!intB){
	    bp <- bp[,-1, drop = FALSE]
	    b1p <- b1p[,-1, drop = FALSE]
	  }
	}

	bpij <- NULL; for(i in 1:ncol(bp)){bpij <- cbind(bpij, bp*bp[,i])}
	
	internal.bfun <- list(p = p, bp = bp, b1p = b1p, bpij = bpij, deltap = p[2] - p[1])
    attr(internal.bfun, "pfun") <- approxfun(c(p[1], 0.5*(p[-d] + p[-1])),p, method = "constant", rule = 2)

	# output. U = the original variables. V = the scaled/rotated variables.
	# stats.B, stats.X, stats.y = lists with the values use to scale/rotate

	stats.B <- list(m = mB, s = sB, M = MB, S = SB, rot = rotB, const = constB, vars = varsB,
		intercept = intB, term.labels = termlabelsB, assign = assignB, coef.names = coefnamesB)
	stats.X <- list(m = mX, s = sX, M = MX, S = SX, rot = rotX, const = constX, vars = varsX,
		intercept = intX, term.labels = termlabelsX, assign = assignX, coef.names = coefnamesX)
	stats.y <- list(m = my, M = My)

	V <- list(X = X, Xw = X*weights, y = y, weights = weights)
	list(mf = mf, U = U, V = V, stats.B = stats.B, stats.X = stats.X, stats.y = stats.y, 
		internal.bfun = internal.bfun, bfun = bfun, s = s)
}


check.out <- function(theta, S, covar){

	blockdiag <- function(A, d, type = 1){
		h <- nrow(A); g <- d/h
		if(type == 1){
			out <- diag(1,d)
			for(j in 1:g){ind <- (j*h - h  + 1):(j*h); out[ind,ind] <- A}
		}
		else{
			out <- matrix(0,d,d)
			for(i1 in 1:h){
				for(i2 in 1:h){
					ind1 <- (i1*g - g  + 1):(i1*g)
					ind2 <- (i2*g - g  + 1):(i2*g)
					out[ind1, ind2] <- diag(A[i1,i2],g)
				}
			}
			out <- t(out)
		}
		out
	}

	mydiag <- function(x){
		if(length(x) > 1){return(diag(x))}
		else{matrix(x,1,1)}
	}

	th <- cbind(c(theta))
	q <- nrow(theta)
	k <- ncol(theta)
	g <- q*k
	aX <- S$X; ay <- S$y; aB <- S$B
	cX <- aX$const; cB <- aB$const

	##########################

	A <- blockdiag(mydiag(1/aX$S), g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aX$intercept){
		A <- diag(1,q); A[cX,] <- -aX$M; A[cX, cX] <- 1
		A <- blockdiag(A,g)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(mydiag(1/aB$S),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aB$intercept){
		A <- diag(1,k); A[,cB] <- -aB$M; A[cB, cB] <- 1
		A <- blockdiag(A,g,2)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(aX$rot,g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	A <- blockdiag(t(aB$rot),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	##########################

	A <- blockdiag(mydiag(1/aX$s),g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aX$intercept){
		A <- diag(1,q); A[cX,] <- -aX$m/aX$s[cX]; A[cX, cX] <- 1
		A <- blockdiag(A,g)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(mydiag(1/aB$s),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aB$intercept){
		A <- diag(1,k); A[,cB] <- -aB$m/aB$s[cB]; A[cB, cB] <- 1
		A <- blockdiag(A,g,2)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	v <- (ay$M - ay$m)/10
	th <- th*v
	covar <- covar*(v^2)
	theta <- matrix(th,q,k)
	theta[cX,cB] <- theta[cX,cB] + ay$m/aB$s[cB]/aX$s[cX]
	
	list(theta = theta, covar = covar)
}


# integrated M-quantile regression

iMqr.internal <- function(mf,cl, formula.p, tol = 1e-6, maxit, s, psi, plim){

	A <- check.in(mf, formula.p, s, plim); V <- A$V; U <- A$U; s <- A$s
	mf <- A$mf; n <- nrow(mf)
	S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
	attributes(A$bfun) <- c(attributes(A$bfun), S$B)
	bfun <- A$internal.bfun
	
	if(is.character(psi)){psi <- get(psi, mode = "function", envir = parent.frame())}
	if(is.function(psi)){psi <- psi()}
	if(is.null(psi$name)){stop("only 'Huber' psi function is currently implemented")}

	if(missing(maxit)){maxit <- 10 + 5*sum(s)}
	else{maxit <- max(1, maxit)}
	

	yy <- (if((q <- length(S$X$vars)) > 0) qexp(rank(V$y)/(n + 1)) else NULL)
	theta <- start.theta(V$y, V$X, V$weights, bfun, 
	   df = max(5, min(15, round(n/30/(q + 1)))), yy, s = s)
	
	# Remark: given sigma, I perform complete optimization over theta; then I update sigma again.
	# The possibility of updating sigma at each Newton step for theta has been considered, and discarded
	# (as it was slower: estimating sigma is time-consuming)
	for(i in 1:maxit){
		sigma <- isigma(theta, V$y,V$X,V$weights, bfun,psi)
		fit <- iMqr.newton(theta, V$y, V$X, V$Xw, V$weights, 
		          bfun, psi, s,sigma, 
		          tol = tol, maxit = maxit, safeit = 2 + sum(s)*(i < 3), eps0 = 0.1)
		if(max(abs(fit$theta - theta)) < tol){break}
		theta <- fit$theta
	}
	covar <- cov.theta(fit, V$y,V$X,V$Xw, V$weights, bfun, psi, s, sigma)

	# output
  
	attr(mf, "assign") <- S$X$assign
	attr(mf, "stats") <- S
	attr(mf, "all.vars") <- V
	attr(mf, "all.vars.unscaled") <- U
	attr(mf, "Q0") <- covar$Q0
	attr(mf, "internal.bfun") <- bfun
	attr(mf, "bfun") <- A$bfun
	attr(mf, "theta") <- fit$coefficients

	out <- check.out(fit$theta, S, covar = covar$Q)
	v <- (S$y$M - S$y$m)/10
	fit <- list(coefficients = out$theta, plim = plim, #sigma = sigma*v, 
		call = cl, converged = (i < maxit), n.it = i,
		obj.function = fit$L*v + length(V$y)*bfun$deltap*sum(log(sigma) + log(v)), 
		s = s, psi = psi$psi_name,
		covar = out$covar, mf = mf)
	jnames <- c(sapply(attr(A$bfun, "coef.names"), 
		function(x,y){paste(x,y, sep = ":")}, y = S$X$coef.names))
	dimnames(fit$covar) <- list(jnames, jnames)
	dimnames(fit$coefficients) <- dimnames(fit$s) <- list(S$X$coef.names, S$B$coef.names)

	# CDF and PDF, precision ~ 1e-6

	fit$CDF <- p.bisec(fit$coef,U$y,U$X,A$bfun)
	b1 <- apply_bfun(A$bfun, fit$CDF, "b1fun")
	fit$PDF <- 1/c(rowSums((U$X%*%fit$coef)*b1))
	if(any(fit$PDF < 0)){warning("crossing M-quantiles detected (PDF < 0 at some y)")}
	# fit$PDF[attr(fit$CDF, "out")] <- 0 # removed in v1.1
	attributes(fit$CDF) <- attributes(fit$PDF) <- list(names = rownames(mf))

	# finish

	class(fit) <- "iMqr"
	if(!fit$converged){warning("the algorithm did not converge")}
	fit
}

#' @export
print.iMqr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	"\n\n", sep = "")

	cat("Coefficients:\n")
	print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	invisible(x)
}




# Compute either bfun or b1fun. Note that I never need them both in predictions!
apply_bfun <- function(bfun, p, fun = c("bfun", "b1fun")){
  
  k <- attr(bfun, "k")
  n <- length(p)
  
  if(inherits(bfun, "slp.basis")){
    pp <- matrix(, n, k + 1)
    pp[,1] <- 1; pp[,2] <- p
    if(k > 1){for(j in 2:k){pp[,j + 1] <- pp[,j]*p}}
    if(fun == "bfun"){out <- tcrossprod(pp, t(bfun$a))}
    else{out <- cbind(0, tcrossprod(pp[,1:k, drop = FALSE], t(bfun$a1[-1,-1, drop = FALSE])))}
    if(!attr(bfun, "intB")){out <- out[,-1, drop = FALSE]}
  }
  else{
    out <- matrix(,n,k)
    if(fun == "bfun"){for(j in 1:k){out[,j] <- bfun[[j]](p)}}
    else{for(j in 1:k){out[,j] <- bfun[[j]](p, deriv = 1)}}
  }
  out
}



# Bisection for external use. Precision about 1e-6
p.bisec <- function(theta, y,X, bfun, n.it = 20){
  
  n <- length(y); k <- ncol(theta)
  bp <- attr(bfun, "bp")
  p <- attr(bfun, "p")
  Xtheta <- tcrossprod(X, t(theta))
  
  eta <- tcrossprod(Xtheta,bp[512,, drop = FALSE])
  m <- as.integer(512 + sign(y - eta)*256)
  
  for(i in 3:10){
    eta <- .rowSums(Xtheta*bp[m,, drop = FALSE], n, k)
    m <- m + as.integer(sign(y - eta)*(2^(10 - i)))
  }
  m <- p[m]
  
  for(i in 11:n.it){
    bp <- apply_bfun(bfun, m, "bfun")
    delta.m <- y - .rowSums(Xtheta*bp, n,k)
    m <- m + sign(delta.m)/2^i
  }
  
  m <- c(m)
  
  out.l <- which(m == 1/2^n.it)
  out.r <- which(m == 1 - 1/2^n.it)
  m[out.l] <- 0; m[out.r] <- 1
  attr(m, "out") <- c(out.l, out.r)
  attr(m, "out.r") <- out.r
  
  m
}

# for external use, only returns b(p)
make.bfun <- function(p,x){
  n <- length(x)
  x1 <- x[1:(n-1)]
  x2 <- x[2:n]
  if(all(x1 < x2) | all(x1 > x2)){method <- "hyman"}
  else{method <- "fmm"}
  splinefun(p,x, method = method)
}

num.fun <- function(dx,fx, op = c("int", "der")){
  n <- length(dx) + 1
  k <- ncol(fx)
  fL <- fx[1:(n-1),, drop = FALSE]
  fR <- fx[2:n,, drop = FALSE]
  
  if(op == "int"){out <- apply(rbind(0, 0.5*dx*(fL + fR)),2,cumsum)}
  else{
    out <- (fR - fL)/dx
    out <- rbind(out[1,],out)
  }
  out
}




#' @export
plf <- function(p, knots){ # basis of piecewise linear function
	if(is.null(knots)){return(cbind(b1 = p))}
	k <- length(knots <- sort(knots))
	ind1 <- 1
	ind2 <- NULL
	for(j in knots){ind1 <- cbind(ind1, (p > j))}
	for(j in k:1){ind2 <- cbind(ind2, ind1[,j] - ind1[,j+1])}
	ind2 <- cbind(ind1[,k+1], ind2)[,(k + 1):1, drop = FALSE]
	ind1 <- cbind(ind1,0)
	knots <- c(0,knots,1)
	a <- NULL
	for(j in 1:(k + 1)){
		a <- cbind(a, (p - knots[j])*ind2[,j] + (knots[j + 1] - knots[j])*ind1[,j + 1])
	}
	colnames(a) <- paste("b", 1:(k + 1), sep = "")
	attr(a, "knots") <- knots[2:(k+1)]
	a
}

slp.basis <- function(k, intercept){ # shifted Legendre polynomials basis

	K <- k + 1

	# matrix a such that P%*%a is an orthogonal polynomial, P = (1, p, p^2, p^3, ...)

	a <- matrix(0, K,K)
	for(i1 in 0:k){
		for(i2 in 0:i1){
			a[i2 + 1, i1 + 1] <- choose(i1,i2)*choose(i1 + i2, i2)
		}
	}
	a[,seq.int(2, K, 2)] <- -a[,seq.int(2, K, 2)]
	a[seq.int(2, K, 2),] <- -a[seq.int(2, K, 2),]

	# a1 = first derivatives to be applied to P' = (0,1, p, p^2, ...)
	# A = first integral to be applied to PP = (p, p^2, p^3, p^4, ...)
	# AA = second integral to be applied to PPP = (p^2, p^3, p^4, p^5, ...)

	a1 <- A <- AA <- matrix(,K,K)
	for(j in 0:k){
		a1[j + 1,] <- a[j + 1,]*j
		A[j + 1,] <- a[j + 1,]/(j + 1)
		AA[j + 1,] <- A[j + 1,]/(j + 2)
	}

	if(!intercept){a[1,-1] <- A[1,-1] <- AA[1, -1] <- 0}
	out <- list(a = a, a1 = a1, A = A, AA = AA)
	attr(out, "k") <- k
	class(out) <- "slp.basis"
	out
}


#' @export
slp <- function(p, k = 3, intercept = FALSE){
	if((k <- round(k)) < 1){stop("k >= 1 is required")}
	P <- cbind(1, outer(p, seq_len(k), "^"))
	B <- P%*%slp.basis(k, intercept)$a
	colnames(B) <- paste("slp", 0:k, sep = "")
	B <- B[,-1, drop = FALSE]
	attr(B, "k") <- k
	class(B) <- "slp"
	B
}


is.slp <- function(f){
	test.p <- seq(0,1,0.1)
	B <- model.matrix(f, data = data.frame(p = test.p))
	if(nrow(B) == 0){return(FALSE)}
	a <- attr(B, "assign")
	if(any(a > 1)){return(FALSE)}
	B <- B[,a == 1, drop = FALSE]
	k <- ncol(B)
	intercept <- FALSE
	if(any(B != slp(test.p, k = k, intercept))){
		intercept <- TRUE
		if(any(B != slp(test.p, k = k, intercept))){
			return(FALSE)
		}
	}
	out <- TRUE
	attr(out, "k") <- k
	attr(out, "intercept") <- intercept
	attr(out, "intB") <- any(a == 0)
	out
}

# in iobjfun, H is computed in the wrong order. This function restores the correct format of H.
sortH <- function(H){
  q <- sqrt(nrow(H))
  k <- sqrt(ncol(H))
  
  kstep <- seq.int(1, 1 + (k - 1)*k, k)
  qstep <- seq.int(1, 1 + (q - 1)*q, q)
  
  HH <- NULL
  for(h in 1:k){
    for(j in 1:q){
      HH <- rbind(HH, c(H[qstep + (j - 1), kstep + (h - 1)]))
    }
  }
  HH
}

# Note: the loss does not include the part with sigma, n*sum(log(sigma))*dtau
iobjfun <- function(theta, y,X,Xw,weights, bfun, psi, sigma, u, H = FALSE, i = FALSE){
  
  tau <- bfun$p
  dtau <- bfun$deltap
  bp <- bfun$bp
  bpij <- bfun$bpij
  par <- psi$par
  n <- length(y)

  if(missing(u)){
    eta <- tcrossprod(tcrossprod(X, t(theta)), bp)
    u <- t(t(y - eta)/sigma)
    
    l <- .rowSums(psi$rho_tau(u, tau, par), nrow(u),ncol(u))*weights
    L <- sum(l)*dtau
    return(list(L = L, u = u))
  }
  else{
    if(!i){
      g <- psi$psi_tau(u, tau, par)
      G <- -tcrossprod(crossprod(Xw,g), t(bp/sigma*dtau))
      G <- c(G)
    }
    else{
      g <- psi$psi_tau(u, tau, par)
      G.i <- 0
      for(i in 1:length(tau)){
        g.i <- NULL
        for(h in 1:ncol(theta)){g.i <- cbind(g.i, X*g[,i]*bp[i,h]/sigma[i])}
        G.i <- G.i + g.i
      }
      G <- -G.i*dtau
    }
    
    if(!H){return(G)}
    
    q <- nrow(theta); k <- ncol(theta)
    h <- psi$psi1_tau(u, tau, par)
    bpij <- t(bpij/sigma^2*dtau)
       
    H <- NULL
    count <- 0
    for(j in 1:q){
      htemp <- tcrossprod(crossprod(Xw*X[,j],h), bpij)
      H <- rbind(H, htemp)
    }
    return(list(G = G, H = sortH(H)))
  }
}


isigma <- function(theta, y,X,weights,bfun,psi){
  
  O <- function(sigma,tau,psi,u0){
    n <- nrow(u0)
    d <- ncol(u0)
    par <- psi$par
    -n/sigma + (1/sigma^2)*.colSums(psi$psi_tau(t(t(u0)/sigma), tau, par)*u0,n,d)
  }
  
  tau <- bfun$p
  bp <- bfun$bp
  eta <- tcrossprod(tcrossprod(X, t(theta)), bp)
  d <- length(tau)
  n <- length(y)
  
  u0 <- y - eta
  sigmaL <- rep.int(0,d); sigmaR <- .colMeans(abs(u0), n,d)
  oL <- rep.int(Inf,d); oR <- O(sigmaR,tau,psi,u0)
  while(any(oR > 0)){sigmaR <- 2*sigmaR; oR <- O(sigmaR,tau,psi,u0)}
  
  nit <- ceiling(log2(max(sigmaR)*1e+15))
  for(i in 1:nit){
    sigmaC <- (sigmaL + sigmaR)/2
    oC <- O(sigmaC,tau,psi,u0)
    w <- (oC < 0); ww <- (!w)
    oR[w] <- oC[w]; sigmaR[w] <- sigmaC[w]
    oL[ww] <- oC[ww]; sigmaL[ww] <- sigmaC[ww]
    if(max(abs(c(oL,oR))/n) < 1e-6){break}
  }

  (sigmaL + sigmaR)/2
}


iMqr.newton <- function(theta, y,X,Xw, weights, bfun,psi, s,sigma, 
                        tol, maxit, safeit, eps0){
 
  q <- nrow(theta)
  k <- ncol(theta)
  s <- c(s == 1)
  
  LL <- iobjfun(theta, y,X,Xw,weights, bfun, psi, sigma)
  G <- iobjfun(theta, y,X,Xw,weights, bfun, psi, sigma, LL$u)
  
  g <- G[s]; L <- LL$L
  conv <- FALSE
  eps <- eps0

  # Preliminary safe iterations, only g is used
  
  for(i in 1:safeit){
    
    if(conv | max(abs(g)) < tol){break}
    u <- rep.int(0, q*k)
    u[s] <- g
    delta <- matrix(u, q,k)	
    delta[is.na(delta)] <- 0
    cond <- FALSE
    
    while(!cond){
      
      new.theta <- theta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      LL1 <- iobjfun(new.theta, y,X,Xw,weights, bfun, psi, sigma)
      cond <- (LL1$L < L)
      eps <- eps*0.5
    }
    
    if(conv){break}
    
    theta <- new.theta
    LL <- LL1
    G <- iobjfun(theta, y,X,Xw,weights, bfun, psi, sigma, LL$u)
    g <- G[s]; L <- LL$L
    eps <- min(eps*2,0.1)
  }
  
  # Newton-Raphson
  
  alg <- "nr"
  conv <- FALSE
  eps <- 0.1
  G <- iobjfun(theta, y,X,Xw,weights, bfun, psi, sigma, LL$u, H = TRUE)
  g <- G$G[s]; h <- G$H[s,s, drop = FALSE]
  h <- h + diag(0.0001, nrow(h)) # added in version 1.1

  for(i in 1:maxit){
   
    if(conv | max(abs(g)) < tol){break}
    
    ####
    
    H1 <- try(chol(h), silent = TRUE)
    err <- (inherits(H1, "try-error"))
    
    if(!err){
      if(alg == "gs"){alg <- "nr"; eps <- 1}
      delta <- chol2inv(H1)%*%g
    }
    else{
      if(alg == "nr"){alg <- "gs"; eps <- 1}
      delta <- g
    }
    
    u <- rep.int(0, q*k)
    u[s] <- delta
    delta <- matrix(u, q,k)	
    delta[is.na(delta)] <- 0
    cond <- FALSE
    while(!cond){
      new.theta <- theta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      LL1 <- iobjfun(new.theta, y,X,Xw,weights, bfun, psi, sigma)
      cond <- (LL1$L < L)
      eps <- eps*0.5
    }
    
    if(conv){break}
    theta <- new.theta
    LL <- LL1
    G <- iobjfun(theta, y,X,Xw,weights, bfun, psi, sigma, LL$u, H = TRUE)
    g <- G$G[s]; h <- G$H[s,s, drop = FALSE]; L <- LL$L
    
    if(i > 1){eps <- min(eps*10,1)}
    else{eps <- min(eps*10,0.1)}
  }

  list(theta = matrix(theta, q, k), LL = LL, L = L, g = g, h = h, 
       n.it = i, converged = (i < maxit), fullrank = (alg == "nr"))
}





start.theta <- function(y,x, weights, bfun, df, yy, s){

	if(is.null(yy)){p.star <- (rank(y) - 0.5)/length(y)}
	else{
	  pch.fit.ct <- getFromNamespace("pch.fit.ct", ns = "pch")
	  predF.pch <- getFromNamespace("predF.pch", ns = "pch")
	  m0 <- suppressWarnings(pch.fit.ct(y = yy,  
		x = cbind(1,x), w = weights, breaks = df))
	  p.star <- 1 - predF.pch(m0)[,3]
	}

	pfun <- attr(bfun, "pfun")
	p.star <- pfun(p.star)  
	b.star <- bfun$bp[match(p.star, bfun$p),]
	X <- model.matrix(~ -1 + x:b.star)
	X <- X[, c(s) == 1, drop = FALSE]
	start.ok <- FALSE
	while(!start.ok){

		m <- lm.wfit(X, y, weights)
		res <- m$residuals
		start.ok <- all(w <- (abs(res)/sd(res) < 5))
		if(start.ok | sum(w) < 0.5*length(y)){break}
		X <- X[w,, drop = FALSE]
		y <- y[w]
		weights <- weights[w]
	}
	out <- rep.int(0, length(s))
	out[c(s) == 1] <- m$coef
	out <- matrix(out, ncol(x))
	out[is.na(out)] <- 0
	out
}



# Note: if s has some zeroes, the covariance matrix Q will contain some zero-columns and rows,
# while the gradient and hessian will just omit the parameters that are not estimated

cov.theta <- function(fit, y,X,Xw, weights, bfun, psi, s, sigma){

  theta <- fit$theta
	s <- c(s == 1)

	# first derivatives for OP
	
	G.i <- iobjfun(theta, y,X,Xw, weights, bfun, psi, sigma, fit$LL$u, H = FALSE, i = TRUE)[,s, drop = FALSE]

	# Hessian
	
  H <- fit$h # note: this is already computed only where s = 1.
	H <- 0.5*(H + t(H))

	# Covariance matrix

	Omega <- chol2inv(chol(t(G.i*weights)%*%G.i))
	Q0 <- t(H)%*%Omega%*%H
	Q <- chol2inv(chol(Q0))#/(dtau^2)
	U <- matrix(0, length(s), length(s))
	U[s,s] <- Q

	list(Q = U, Q0 = Q0)
}




#' @export
summary.iMqr <- function(object, p, cov = FALSE, ...){

	if(missing(p)){
		mf <- object$mf
		theta <- object$coefficients
		w <- attr(mf, "all.vars")$weights
  
		u <- sqrt(diag(object$covar))
		u <- matrix(u, q <- nrow(theta), r <- ncol(theta))
		dimnames(u) <- dimnames(theta)
		test <- (if(q*r == 1) NULL else iqr.waldtest(object))
		out <- list(converged = object$converged, n.it = object$n.it,
			coefficients = theta, se = u, 
			test.x = test$test.x, test.p = test$test.p)

		out$obj.function <- object$obj.function		
		out$n <- nrow(object$mf)
		out$free.par <- sum(theta != 0)
	}
	else{
		out <- list()
		for(i in 1:length(p)){out[[i]] <- extract.p(object, p[i], cov)}
		names(out) <- paste("p =", p)
		attr(out, "nacoef") <- which(apply(object$coefficients,1, function(v){all(v == 0)}))
	}
	out$call <- object$call
	class(out) <- "summary.iMqr"
	out	
}

#' @export
print.summary.iMqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	if(!is.null(x$coef)){

		nacoef <- which(x$coef == 0)
		x$coef[nacoef] <- x$se[nacoef] <- NA

		cat("converged:", x$converged, "\n")
		cat("n. of iterations:", x$n.it, "\n")
		cat("n. of observations:", x$n, "\n")
		cat("n. of free parameters:", x$free.par, "\n\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Coefficients:\n")
		print.default(format(x$coef, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("Standard errors:\n")
		print.default(format(x$se, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		if(!is.null(x$test.x)){
			cat("Wald test for x:\n")

			printCoefmat(x$test.x, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
			cat("\n\n")
		}

		if(!is.null(x$test.p)){
			cat("Wald test for b(p):\n")
			printCoefmat(x$test.p, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}

		if(!is.null(x$obj.function)){
			cat("\n")
			cat("Minimized loss function:", x$obj.function)
		}
	}

	else{
		nacoef <- attr(x, "nacoef")
		for(j in 1:(length(x) - 1)){
			cat(paste(names(x)[j], "\n"))
			cat("\n")
			cat("Coefficients:\n")
			coe <- x[[j]]$coef; coe[nacoef,] <- NA
			printCoefmat(coe, digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3, 
				P.values = TRUE, has.Pvalue = TRUE)
			cat("\n")

			if(!is.null(x[[j]]$cov)){
				cat("Covar:\n")
				print.default(format(x[[j]]$cov, digits = digits), print.gap = 2L, quote = FALSE)
			}
			cat("\n\n")
		}
	}

	invisible(x)
}



extract.p <- function(model,p, cov = FALSE){

	theta <- model$coefficients
	v <- model$covar
	q <- nrow(theta)
	k <- ncol(theta)

	bfun <- attr(model$mf, "bfun")
	pred.p <- apply_bfun(bfun, p, "bfun")
	beta <- c(pred.p%*%t(theta))

	cov.beta <- matrix(NA,q,q)
	for(j1 in 1:q){
		w1 <- seq.int(j1,q*k,q)
		for(j2 in j1:q){
			w2 <- seq.int(j2,q*k,q)
			cc <- v[w1,w2, drop = FALSE]
			cov.beta[j1,j2] <- cov.beta[j2,j1] <- pred.p%*%cc%*%t(pred.p)
		}
	}
	se <- sqrt(diag(cov.beta))
	z <- beta/se
	out <- cbind(beta, se, z, 2*pnorm(-abs(z)))
	colnames(out) <- c("Estimate", "std.err", "z value", "p(>|z|))")
	rownames(out) <- colnames(cov.beta) <- rownames(cov.beta) <- rownames(theta)
	if(cov){list(coef = out, cov = cov.beta)}
	else{list(coef = out)}
}

#' @export
plot.iMqr <- function(x, conf.int = TRUE, polygon = TRUE, which = NULL, ask = TRUE, ...){

	plot.iqr.int <- function(p,u,j,conf.int,L){
		beta <- u[[j]]$beta
		if(is.null(L$ylim)){
			if(conf.int){y1 <- min(u[[j]]$low); y2 <- max(u[[j]]$up)}
			else{y1 <- min(beta); y2 <- max(beta)}
			L$ylim <- c(y1,y2)
		}
		plot(p, u[[j]]$beta, xlab = L$xlab, ylab = L$ylab, main = L$labels[j], 
		  type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col, axes = L$axes, 
		  frame.plot = L$frame.plot, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
		if(conf.int){
		  if(polygon){
		    yy <- c(u[[j]]$low, tail(u[[j]]$up, 1), rev(u[[j]]$up), u[[j]]$low[1])
		    xx <- c(p, tail(p, 1), rev(p), p[1])
		    polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
		  }
		  else{
		    points(p, u[[j]]$low, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		    points(p, u[[j]]$up, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		  }
		}
	}

	L <- list(...)
	if(is.null(L$xlim)){L$xlim = x$plim}
	if(is.null(L$lwd)){L$lwd <- 2}
	if(is.null(L$col)){L$col <- "black"}
	if(is.null(L$xlab)){L$xlab <- "p"}
	if(is.null(L$ylab)){L$ylab <- "beta(p)"}
	if(is.null(L$cex.lab)){L$cex.lab <- 1}
	if(is.null(L$cex.axis)){L$cex.axis <- 1}
	if(is.null(L$axes)){L$axes <- TRUE}
	if(is.null(L$frame.plot)){L$frame.plot <- TRUE}
	L$labels <- rownames(x$coefficients)
	q <- length(L$labels)

	p <- seq.int(max(0.001,L$xlim[1]), min(0.999,L$xlim[2]), length.out = 100)
	u <- predict.iMqr(x, p = p, type = "beta", se = conf.int)

	if(!is.null(which) | !ask){
		if(is.null(which)){which <- 1:q}
		for(j in which){plot.iqr.int(p,u,j,conf.int,L)}
	}
	else{
		pick <- 1
		while(pick > 0 && pick <= q){
			pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
			if(pick > 0 && pick <= q){plot.iqr.int(p,u,pick,conf.int,L)}
		}
	}
}


# predict function.
# p: default to percentiles for type = "beta". No default for "fitted". Ignored for "CDF".
# se: ignored for type = "CDF"
# x: only for type = "CDF" or type = "fitted"
# y: only for type = "CDF"

#' @export
predict.iMqr <- function(object, type = c("beta", "CDF", "QF", "sim"), newdata, p, se = TRUE, ...){

	if(is.na(match(type <- type[1], c("beta", "CDF", "QF", "sim")))){stop("invalid 'type'")}
	if(type == "beta"){
		if(missing(p)){p <- seq.int(0.01,0.99,0.01)}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
		return(pred.beta(object, p, se))
	}

	mf <- object$mf
	mt <- terms(mf)
	miss <- attr(mf, "na.action")
	nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	xlev <- .getXlevels(mt, mf)
	contr <- attr(mf, "contrasts")

	if(!missing(newdata)){

	  if(type == "CDF"){
	    yn <- as.character(mt[[2]])
	    if(is.na(ind <- match(yn, colnames(newdata))))
	    {stop("for 'type = CDF', 'newdata' must contain the y-variable")}
	  }
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
			{stop("'newdata' must contain all x-variables")}

		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){
			nr <- nrow(newdata)
			if(type == "CDF"){
				out <- data.frame(matrix(NA,nr,3))
				colnames(out) <- c("log.f", "log.F", "log.S")
				rownames(out) <- rownames(newdata)
			}
			else if(type == "QF"){
				out <- data.frame(matrix(NA,nr,length(p)))
				colnames(out) <- paste("p",p, sep = "")
				rownames(out) <- rownames(newdata)
			}
			else{out <- rep.int(NA, nr)}
			return(out)
		}
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
	}
	
	x <- model.matrix(mt, mf, contrasts.arg = contr)

	if(type == "CDF"){
		bfun <- attr(object$mf, "bfun")
		y <- cbind(model.response(mf))[,1]
		Fy <- p.bisec(object$coefficients, y,x, bfun)
		b1 <- apply_bfun(bfun, Fy, "b1fun")
		fy <- 1/c(rowSums((x%*%object$coefficients)*b1))
		fy[attr(Fy, "out")] <- 0
		if(any(fy < 0)){warning("some PDF values are negative (quantile crossing)")}
		CDF <- PDF <- NULL
		CDF[nomiss] <- Fy
		PDF[nomiss] <- fy
		CDF[miss] <- PDF[miss] <- NA
		out <- data.frame(CDF = CDF, PDF = PDF)
		rownames(out)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(out)[miss] <- names(miss)}
		return(out)
	}

	else if(type == "QF"){
		if(missing(p)){stop("please indicate the value(s) of 'p' to compute x*beta(p)")}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}

		fit <- se.fit <- matrix(, length(c(miss,nomiss)), length(p))
		for(j in 1:length(p)){
			fit.beta <- extract.p(object,p[j], cov = se)
			fit[nomiss,j] <- x%*%cbind(fit.beta$coef[,1])
			if(se){se.fit[nomiss,j] <- sqrt(diag(x%*%fit.beta$cov%*%t(x)))}
		}
		fit <- data.frame(fit)
		colnames(fit) <- paste("p",p, sep = "")
		rownames(fit)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(fit)[miss] <- names(miss)}
		if(se){
			se.fit <- data.frame(se.fit)
			colnames(se.fit) <- paste("p",p, sep = "")
			rownames(se.fit)[nomiss] <- rownames(mf)
			if(!is.null(miss)){rownames(se.fit)[miss] <- names(miss)}
			return(list(fit = fit, se.fit = se.fit))
		}
		else{return(fit)}
	}	
	else{
		p <- runif(nrow(x))
		beta <- apply_bfun(attr(object$mf, "bfun"), p, "bfun")%*%t(object$coefficients)
		y <- NULL; y[nomiss] <- rowSums(beta*x); y[miss] <- NA
		return(y)
	}
}

#' @export
terms.iMqr <- function(x, ...){attr(x$mf, "terms")}
#' @export
model.matrix.iMqr <- function(object, ...){
  mf <- object$mf
  mt <- terms(mf)
  model.matrix(mt, mf)
}
#' @export
vcov.iMqr <- function(object, ...){object$covar}
#' @export
nobs.iMqr <- function(object, ...){nrow(object$mf)}


pred.beta <- function(model, p, se = FALSE){

	if(se){
		Beta <- NULL
		SE <- NULL
		for(j in p){
			b <- extract.p(model,j)$coef
			Beta <- rbind(Beta, b[,1])
			SE <- rbind(SE, b[,2])
		}
		out <- list()
		for(j in 1:ncol(Beta)){
			low <- Beta[,j] - 1.96*SE[,j]
			up <- Beta[,j] + 1.96*SE[,j]
			out[[j]] <- data.frame(p = p, beta = Beta[,j], se = SE[,j], low = low, up = up)
		}
		names(out) <- rownames(model$coefficients)
		return(out)
	}
	else{
		theta <- model$coefficients
		beta <- apply_bfun(attr(model$mf, "bfun"), p, "bfun")%*%t(theta)
		out <- list()
		for(j in 1:nrow(theta)){out[[j]] <- data.frame(p = p, beta = beta[,j])}
		names(out) <- rownames(theta)
		return(out)
	}
}


iqr.waldtest <- function(obj){
	bfun <- attr(obj$mf, "bfun")
	ax <- attr(obj$mf, "assign")
	ap <- attr(bfun, "assign")
	theta <- obj$coef
	q <- nrow(theta)
	k <- ncol(theta)
	s <- obj$s
	cc <- obj$covar
	ind.x <- rep(ax,k)
	ind.p <- sort(rep.int(ap,q))
	K <- tapply(rowSums(s), ax, sum)
	Q <- tapply(colSums(s), ap, sum)

	testx <- testp <- NULL

	if(q > 1){
		for(i in unique(ax)){
			theta0 <- c(theta[which(ax == i),])
			w <- which(ind.x == i)
			c0 <- cc[w,w, drop = FALSE]

			theta0 <- theta0[theta0 != 0]
			w <- which(rowSums(c0) != 0)
			c0 <- c0[w,w, drop = FALSE]
			
			if(length(theta0) == 0){tx <- NA}
			else{tx <- t(theta0)%*%chol2inv(chol(c0))%*%t(t(theta0))}
			testx <- c(testx, tx)
		}
		testx <- cbind(testx, df = K, pchisq(testx, df = K, lower.tail = FALSE))
		colnames(testx) <- c("chi-square", "df", "P(> chi)")

		nx <- attr(attr(obj$mf, "terms"), "term.labels")
		if(attr(attr(obj$mf, "terms"), "intercept") == 1){nx <- c("(Intercept)", nx)}
		rownames(testx) <- nx
	}
	if(k > 1){
		for(i in unique(ap)){
			theta0 <- c(theta[,which(ap == i)])
			w <- which(ind.p == i)
			c0 <- cc[w,w, drop = FALSE]

			theta0 <- theta0[theta0 != 0]
			w <- which(rowSums(c0) != 0)
			c0 <- c0[w,w, drop = FALSE]

			if(length(theta0) == 0){tp <- NA}
			else{tp <- t(theta0)%*%chol2inv(chol(c0))%*%t(t(theta0))}
			testp <- c(testp, tp)
		}
		testp <- cbind(testp, df = Q, pchisq(testp, df = Q, lower.tail = FALSE))
		colnames(testp) <- c("chi-square", "df", "P(> chi)")
		np <- attr(bfun, "term.labels")
		if(any(ap == 0)){np <- c("(Intercept)", np)}
		rownames(testp) <- np
	}

	list(test.x = testx, test.p = testp)
}















