## #' @include e_conversion.R



#' @title Symmetric Kiener Distribution K1
#'
#' @description
#' Density, distribution function, quantile function, random generation, 
#' value-at-risk, expected shortfall (+ signed left/right tail mean) 
#' and additional formulae for symmetric Kiener distribution K1. 
#' This distribution is similar to the power hyperbola logistic distribution 
#' but with additional parameters for location (\code{m}) and scale (\code{g}).
#'
#' @param    x    vector of quantiles.
#' @param    q	  vector of quantiles.
#' @param    m	  numeric. The median.
#' @param    g	  numeric. The scale parameter, preferably strictly positive.
#' @param    k	  numeric. The tail parameter, preferably strictly positive.
#' @param    p	  vector of probabilities.
#' @param    lp	  vector of logit of probabilities.
#' @param    n	  number of observations. If length(n) > 1, the length is  
#'                taken to be the number required.
#' @param    log           logical. If TRUE, densities are given in log scale.
#' @param    lower.tail    logical. If TRUE, use p. If FALSE, use 1-p.
#' @param    log.p         logical. If TRUE, probabilities p are given as log(p).
#' @param    signedES      logical. FALSE (default) returns positive numbers for 
#'                         left and right tails. TRUE returns negative number 
#'                         (= \code{ltmkiener1}) for left tail and positive number 
#'                         (= \code{rtmkiener1}) for right tail.
#'
#' @details
#' Kiener distributions use the following parameters, some of them being redundant. 
#' See \code{\link{aw2k}} and \code{\link{pk2pk}} for the formulas and 
#' the conversion between parameters:
#' \itemize{
#'   \item{ \code{m} (mu) is the median of the distribution,. }
#'   \item{ \code{g} (gamma) is the scale parameter. }
#'   \item{ \code{a} (alpha) is the left tail parameter. } 
#'   \item{ \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} 
#'          and describes a global tail parameter. }
#'   \item{ \code{w} (omega) is the right tail parameter. } 
#'   \item{ \code{d} (delta) is the distortion parameter. }
#'   \item{ \code{e} (epsilon) is the eccentricity parameter. }
#' }
#' 
#' Kiener distributions \code{K1(m, g, k, ...)} describe distributions  
#' with symmetric left and right fat tails and with a tail parameter \code{k}. 
#' This parameter is the power exponent mentionned in the Pareto formula and 
#' Karamata theorems. 
#' 
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and is linked for any value of \code{k} to the density at the 
#' median through the relation
#' \deqn{ g * f(x=m, g=g) = \frac{\pi}{4\sqrt{3}} \approx 0.453 }{%
#'        g * f(x=m, g=g) = pi/4/sqrt(3) = 0.453 approximatively}
#'
#' When \code{k = Inf}, \code{g} is very close to \code{sd(x)}.
#' NOTE: In order to match this standard deviation, the value of \code{g} has 
#' been updated from versions < 1.9.0 by a factor 
#' \eqn{ \frac{2\pi}{\sqrt{3}}}{ 2 pi/sqrt(3) }.
#' 
#' The functions \code{dkiener1}, \code{pkiener1} and \code{lkiener1} have an
#' explicit form (whereas \code{dkiener2347}, \code{pkiener2347} and 
#' \code{lkiener2347} have no explicit forms).
#'
#' \code{dkiener1} function is defined for x in (-Inf, +Inf) by: 
#'  \deqn{ 
#'   \begin{array}{l}
#'    y = \frac{1}{k}\frac{\pi}{\sqrt{3}}\frac{(x-m)}{g} \\[4pt]
#'    dkiener1(x,m,g,k) = \pi*\left[2\sqrt{3}\,g\,\sqrt{y^2 +1}
#'                         \left(1+\cosh(k*asinh(y))\right)\right]^{-1}
#'   \end{array}
#'  }{%
#'    y = pi/sqrt(3)*(x-m)/g/k, 
#'    dkiener1(x,m,g,k) = pi/2/sqrt(3)/g /sqrt(y*y+1) /(1+cosh(k*asinh(y)))
#'  }                    
#'
#' \code{pkiener1} function is defined for q in (-Inf, +Inf) by: 
#'  \deqn{ 
#'   \begin{array}{l}
#'    y = \frac{1}{k}\frac{\pi}{\sqrt{3}}\frac{(x-m)}{g} \\[4pt]
#'    pkiener1(q,m,g,k) = 1/(1+exp(-k*asinh(y)))
#'   \end{array}
#'  }{%
#'    y = pi/sqrt(3)*(q-m)/g/k, 
#'    pkiener1(x,m,g,k) = 1/(1+exp(-k*asinh(y)))
#'  }         
#'
#' \code{qkiener1} function is defined for p in (0, 1) by: 
#'  \deqn{ 
#'    qkiener1(p,m,g,k) = m + \frac{\sqrt{3}}{\pi}*g*k*       
#'                            \sinh\left(\frac{logit(p)}{k}\right)
#'  }{%
#'    qkiener1(p,m,g,k) = m + sqrt(3)/pi*g*k*sinh(logit(p)/k) 
#'  }
#'
#' \code{rkiener1} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener1} is the density function calculated from the probability p. 
#' It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dpkiener1(p,m,g,k) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}
#'                          sech\left(\frac{logit(p)}{k}\right) 
#'  }{%
#'    dpkiener1(p,m,g,k) = p*(1-p)*pi/sqrt(3)/g/cosh( logit(p)/k ) 
#'  }
#'
#' \code{dqkiener1} is the derivate of the quantile function calculated from 
#' the probability p. It is defined for p in (0, 1) by: 
#'  \deqn{ 
#'    dqkiener1(p,m,g,k) = \frac{\sqrt{3}}{\pi}\frac{g}{p(1-p)}
#'                         \cosh\left(\frac{logit(p)}{k}\right) 
#'  }{%
#'    dqkiener1(p,m,g,k) = sqrt(3)/pi*g/p/(1-p)*cosh( logit(p)/k ) 
#'  }
#'
#' \code{lkiener1} function is equivalent to kashp function but with additional 
#' parameters \code{m} and \code{g}. Being computed from the x (or q) vector, 
#' it can be compared to the logit of the empirical probability logit(p) 
#' through a nonlinear regression with ordinary or weighted least squares 
#' to estimate the distribution parameters. 
#' It is defined for x in (-Inf, +Inf) by:
#'  \deqn{ 
#'   \begin{array}{l}
#'    y = \frac{1}{k}\frac{\pi}{\sqrt{3}}\frac{(x-m)}{g} \\[4pt]
#'    lkiener1(q,m,g,k) = k*asinh(y)
#'   \end{array}
#'  }{%
#'    lkiener1(x,m,g,k) = k*asinh(pi/sqrt(3)*(x-m)/g/k) 
#'  }
#'
#' \code{dlkiener1} is the density function calculated from the logit of the 
#' probability lp = logit(p). It is defined for lp in (-Inf, +Inf) by: 
#'  \deqn{
#'    dlkiener1(lp,m,g,k) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}
#'                          sech\left(\frac{lp}{k}\right) 
#'  }{%
#'    dlkiener1(lp,m,g,k) = p*(1-p)*pi/sqrt(3)/g/cosh(lp/k) 
#'  }
#'
#' \code{qlkiener1} is the quantile function calculated from the logit of the 
#' probability lp = logit(p). It is defined for lp in (-Inf, +Inf) by: 
#'  \deqn{ 
#'    qlkiener1(p,m,g,k) = m + \frac{\sqrt{3}}{\pi}*g*k*2* \sinh\left(\frac{lp}{k}\right)
#'  }{%
#'    qlkiener1(lp,m,g,k) = m + sqrt(3)/pi*g*k*2*sinh(lp/k)
#'  }
#' 
#' \code{varkiener1} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'  \deqn{ 
#'    varkiener1 <- if\;(p <= 0.5)\;\; (- qkiener1)\;\; else\;\; (qkiener1) 
#'  }{%
#'    varkiener1 <- if (p <= 0.5) (- qkiener1) else (qkiener1) 
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#'
#' \code{ltmkiener1}, \code{rtmkiener1} and \code{eskiener1} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener1} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener1} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'  \deqn{ 
#'    eskiener1 <- if\;(p <= 0.5)\;\; (- ltmkiener1)\;\; else\;\; (rtmkiener1) 
#'  }{%
#'    eskiener1 <- if(p <= 0.5) (- ltmkiener1) else (rtmkiener1)
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#'
#' \code{dtmqkiener1} is the difference between the left tail mean and the quantile 
#' when (p <= 0.5) and the difference between the right tail mean and the quantile 
#' when (p > 0.5). It is in quantile unit and is an indirect measure of the tail curvature.
#' 
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014.  Download it from: 
#' \url{https://www.inmodelia.com/exemples/2014-0627-Rmetrics-Kiener-en.pdf}
#'
#' P. Kiener, Fat tail analysis and package FatTailsR, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' Download it from: 
#' \url{https://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#' 
#' C. Acerbi, D. Tasche, Expected shortfall: a natural coherent alternative to 
#' Value at Risk, 9 May 2001. Download it from: 
#' \url{https://www.bis.org/bcbs/ca/acertasc.pdf}
#' 
#' @seealso 
#' Standardized logistic distribution \code{\link{logisst}}, 
#' asymmetric Kiener distributions K2, K3, K4 and K7  
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}, 
#' \code{\link{kiener7}},  
#' regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' require(graphics)
#' 
#' ### EXAMPLE 1
#' x <- seq(-5, 5, by = 0.1) ; x
#' pkiener1(x, m=0, g=1, k=4)
#' dkiener1(x, m=0, g=1, k=4)
#' lkiener1(x, k=4)
#' plot( x, pkiener1(x, m=0, g=1, k=4), las=1)
#' lines(x, pkiener1(x, m=0, g=1, k=9999))
#' 
#' plot( x, lkiener1(x, m=0, g=1, k=4), las=1)
#' lines(x, lkiener1(x, m=0, g=1, k=9999))
#' 
#' 
#' p <- c(ppoints(11, a = 1), NA, NaN) ; p
#' qkiener1(p, k = 4)
#' dpkiener1(p, k = 4)
#' dqkiener1(p, k=4)
#' 
#' varkiener1(p=0.01, k=4)
#' ltmkiener1(p=0.01, k=4) 
#'  eskiener1(p=0.01, k=4) # VaR and ES should be positive
#' ### END EXAMPLE 1
#' 
#' 
#' ### PREPARE THE GRAPHICS FOR EXAMPLES 2 AND 3
#' xx  <- c(-4,-2, 0, 2, 4)
#' lty <- c( 1, 2, 3, 4, 5, 1)
#' lwd <- c( 2, 1, 1, 1, 1, 1)
#' col <- c("black","green3","cyan3","dodgerblue2","purple2","brown3")
#' lat <- c(-6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9)
#' lgt <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", "logit(0.95)   = 2.9", 
#'          "logit(0.50)   = 0", "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'          "logit(0.001) = -6.9  ")
#' funleg <- function(xy, k) legend(xy, title = expression(kappa), legend = names(k),
#'                   lty = lty, col = col, lwd = lwd, inset = 0.02, cex = 0.8)
#' funlgt <- function(xy) legend(xy, title = "logit(p)", legend = lgt,
#'                               inset = 0.02, cex = 0.6)
#' 
#' ### EXAMPLE 2
#' ### PROBA, DENSITY, LOGIT-PROBA, LOG-DENSITY FROM x
#' x <- seq(-5, 5, by = 0.1) ; head(x, 10)
#' k <- c(9999, 9, 5, 3, 2, 1) ; names(k) <- k
#' 
#' mat11 <- outer(x, k, \(x,k) pkiener1(x, k=k)) ; head(mat11, 10)
#' mat12 <- outer(x, k, \(x,k) dkiener1(x, k=k)) ; mat12
#' mat13 <- outer(x, k, \(x,k) lkiener1(x, k=k)) ; mat13
#' mat14 <- outer(x, k, \(x,k) dkiener1(x, k=k, log=TRUE)) ; mat14
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(x, mat11, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="pkiener1(x, m=0, g=1, k=k)", xlab="", ylab="")
#' 	funleg("topleft", k)
#' 	matplot(x, mat12, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="dkiener1", xlab="", ylab="")
#' 	funleg("topleft", k)
#' 	matplot(x, mat13, type="l", lwd=lwd, lty=lty, col=col, yaxt="n", 
#' 			main="lkiener1", xlab="", ylab="")
#' 	   axis(2, at=lat, las=1)
#' 	funleg("bottomright", k)
#' 	funlgt("topleft")
#' 	matplot(x, mat14, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="log(dkiener1)", xlab="", ylab="")
#' 	funleg("bottom", k)
#' par(op)
#' ### END EXAMPLE 2
#' 
#' 
#' ### EXAMPLE 3
#' ### QUANTILE, DIFF-QUANTILE, DENSITY, LOG-DENSITY FROM p
#' p <- ppoints(1999, a=0) ; head(p, n=10)
#' k <- c(9999, 9, 5, 3, 2, 1) ; names(k) <- k
#' 
#' mat15 <- outer(p, k, \(p,k)  qkiener1(p, k=k)) ; head(mat15, 10)
#' mat16 <- outer(p, k, \(p,k) dqkiener1(p, k=k)) ; head(mat16, 10)
#' mat17 <- outer(p, k, \(p,k) dpkiener1(p, k=k)) ; head(mat17, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(p, mat15, type="l", xlim=c(0,1), ylim=c(-5,5), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="qkiener1(p, m=0, g=1, k=k)", xlab="", ylab="")
#' 	funleg("topleft", k)
#' 	matplot(p, mat16, type="l", xlim=c(0,1), ylim=c(0,40), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="dqkiener1", xlab="", ylab="")
#' 	funleg("top", k)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(0, 0.5), las=1,
#' 		 main="qkiener1, dpkiener1", xlab="", ylab="")
#' 	mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(mat17), 
#' 		   lwd=lwd, lty=1, col=col)
#' 	funleg("topright", k)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(-7, -0.5), las=1,
#' 		 main="qkiener1, log(dpkiener1)", xlab="", ylab="")
#' 	mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(log(mat17)), 
#' 		   lwd=lwd, lty=lty, col=col)
#' 	funleg("bottom", k)
#' par(op)
#' ### END EXAMPLE 3
#' 
#' 
#' ### EXAMPLE 4: PROCESSUS: which processus look credible?
#' ### PARAMETER k VARIES
#' ### RUN SEED ii <- 1 THEN THE cairo_pdf CODE WITH THE 6 SEEDS
#' # cairo_pdf("K1-6x6-stocks-k.pdf")
#' # for (ii in c(1,2016,2018,2022,2023,2024)) {
#' 	ii <- 1
#' 	set.seed(ii)
#' 	p <- sample(ppoints(299, a=0), 299)
#' 	k <- c(9999, 6, 4, 3, 2, 1) ; names(k) <- k
#' 	mat18 <- outer(p, k, \(p,k)  qkiener1(p=p, g=0.85, k=k)) 
#' 	mat19 <- apply(mat18, 2, cumsum)
#' 	title <- paste0(
#' 		"stock_", ii,    
#' 	     ":  k_left = c(", paste(k[1:3], collapse = ", "), ")",
#' 	    ",  k_right = c(", paste(k[4:6], collapse = ", "), ")")
#' 	plot.ts(mat19, ann=FALSE, las=1, 
#' 			mar.multi=c(0,3,0,1), oma.multi=c(3,0,3,0.5))
#' 	mtext(title, outer = TRUE, line=-1.5, font=2)
#' 	plot.ts(mat18, ann=FALSE, las=1, 
#' 			mar.multi=c(0,3,0,1), oma.multi=c(3,0,3,0.5))
#' 	mtext(title, outer=TRUE, line=-1.5, font=2)
#' # }
#' # dev.off()
#' ### END EXAMPLE 4
#' 
#' 
#' 
#' @name kiener1
NULL

#' @export
#' @rdname kiener1
dkiener1 <- function(x, m = 0, g = 1, k = 3.2, log = FALSE) {
	xx  <- (x-m)*pi/sqrt(3)/g
	y   <- xx/k
	kay <- k *asinh(y)
	dx  <- 1/sqrt(y*y + 1) 
	# since v1.9.0
	v     <- pi/2/sqrt(3)/g*dx/(1 + cosh(kay))
	if(log) log(v) else v
}

#' @export
#' @rdname kiener1
pkiener1 <- function(q, m = 0, g = 1, k = 3.2, lower.tail = TRUE, log.p = FALSE) {
	xx <- (q-m)*pi/sqrt(3)/g
	v  <- 1/(1 + exp(-k *asinh(xx/k)))
	if(lower.tail) v <- v else v <- 1-v
	if(log.p) log(v) else v
}

#' @export
#' @rdname kiener1
qkiener1 <- function(p, m = 0, g = 1, k = 3.2, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	# since v1.9.0
	v  <- m + sqrt(3)/pi*g* k*sinh(logit(p)/k)
	v
}

#' @export
#' @rdname kiener1
rkiener1 <- function(n, m = 0, g = 1, k = 3.2) {
	p <- runif(n)
	v <- qkiener1(p, m, g, k, lower.tail = TRUE, log.p = FALSE)
	v
}

#' @export
#' @rdname kiener1
dpkiener1 <- function(p, m = 0, g = 1, k = 3.2, log = FALSE) {
	# since v1.9.0
	v   <- p*(1-p)*pi/sqrt(3)/g/cosh(logit(p)/k)
	if(log) log(v) else v
}

#' @export
#' @rdname kiener1
dqkiener1 <- function(p, m = 0, g = 1, k = 3.2, log = FALSE) {
	# since v1.9.0
	v   <- 1/p/(1-p)*sqrt(3)/pi*g*cosh(logit(p)/k)
	if(log) log(v) else v
}

#' @export
#' @rdname kiener1
lkiener1 <- function(x, m = 0, g = 1, k = 3.2) { 
	# since v1.9.0
	xx <- (x-m)*pi/sqrt(3)/g
	v  <- k * asinh(xx/k) 
	v
}

#' @export
#' @rdname kiener1
dlkiener1 <- function(lp, m = 0, g = 1, k = 3.2, log = FALSE) {
	p <- invlogit(lp)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g /cosh(lp/k)
	if(log) log(v) else v
}

#' @export
#' @rdname kiener1
qlkiener1 <- function(lp, m = 0, g = 1, k = 3.2, lower.tail = TRUE ) {
	if(!lower.tail) lp <- -lp
	# since v1.9.0
	v  <- m + sqrt(3)/pi*g *k *sinh(lp/k)
	v
}

#' @export
#' @rdname kiener1
varkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener1(p[i], m, g, k),
					  qkiener1(p[i], m, g, k))
	}
	va
}

#' @export
#' @rdname kiener1
ltmkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p  <- exp(p)
	ltm <- if (lower.tail) {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k)
			+pbeta(p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k))	
	} else {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k)
			+pbeta(p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k))
	}
	ltm
}

#' @export
#' @rdname kiener1
rtmkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p  <- exp(p)
	rtm <- if (!lower.tail) {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k)
			+pbeta(1-p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k))	
	} else {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k)
			+pbeta(1-p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k))
	}
	rtm
}

#' @export
#' @rdname kiener1
dtmqkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                        lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener1(p[i], m, g, k, lower.tail, log.p) 
			- qkiener1(p[i], m, g, k, lower.tail, log.p),
			rtmkiener1(p[i], m, g, k, lower.tail, log.p) 
			- qkiener1(p[i], m, g, k, lower.tail, log.p))	
	}
	dtmq
}

#' @export
#' @rdname kiener1
eskiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                      lower.tail = TRUE, log.p = FALSE, signedES = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	es  <- p
	for (i in seq_along(p)) {
		if (signedES) {
			es[i] <- ifelse(p[i] <= 0.5, 
					  ltmkiener1(p[i], m, g, k),
					  rtmkiener1(p[i], m, g, k))
		} else {
			es[i] <- ifelse(p[i] <= 0.5, 
				      abs(ltmkiener1(p[i], m, g, k)),
					  abs(rtmkiener1(p[i], m, g, k)))		
		}
	}
	es
}



