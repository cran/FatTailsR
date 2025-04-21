## #' @include f_kiener1.R



#' @title Asymmetric Kiener Distribution K2
#'
#' @description
#' Density, distribution function, quantile function, random generation,
#' value-at-risk, expected shortfall (+ signed left/right tail mean) 
#' and additional formulae for asymmetric Kiener distribution K2.
#'
#' @param    x    vector of quantiles.
#' @param    q    vector of quantiles.
#' @param    m    numeric. The median.
#' @param    g    numeric. The scale parameter, preferably strictly positive.
#' @param    a    numeric. The left tail parameter, preferably strictly positive.
#' @param    w    numeric. The right tail parameter, preferably strictly positive.
#' @param    p    vector of probabilities.
#' @param    lp   vector of logit of probabilities.
#' @param    n    number of observations. If length(n) > 1, the length is  
#'                                  taken to be the number required.
#' @param    log           logical. If TRUE, densities are given in log scale.
#' @param    lower.tail    logical. If TRUE, use p. If FALSE, use 1-p.
#' @param    log.p         logical. If TRUE, probabilities p are given as log(p).
#' @param    signedES      logical. FALSE (default) returns positive numbers for 
#'                         left and right tails. TRUE returns negative number 
#'                         (= \code{ltmkiener4}) for left tail and positive number 
#'                         (= \code{rtmkiener4}) for right tail.
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
#' Kiener distributions \code{K2(m, g, a, w)} are distributions 
#' with asymmetrical left 
#' and right fat tails described by the parameters \code{a} (alpha) for 
#' the left tail and \code{w} (omega) for the right tail. These parameters 
#' correspond to the power exponent that appear in Pareto formula and 
#' Karamata theorems. 
#'
#' As \code{a} and \code{w} are highly correlated, the use of Kiener distributions
#' (\code{K3(..., k, d)} K4 (\code{K4(..., k, e)} is an alternate solution.
#' 
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and is linked for any value of \code{a} and \code{w} to the density at the 
#' median through the relation
#' \deqn{ g * f(x=m, g=g) = \frac{\pi}{4\sqrt{3}} \approx 0.453 }{%
#'        g * f(x=m, g=g) = pi/4/sqrt(3) = 0.453 approximatively}
#'
#' When \code{a = Inf} and \code{w = Inf}, \code{g} is very close to \code{sd(x)}. 
#' NOTE: In order to match this standard deviation, the value of \code{g} has 
#' been updated from versions < 1.9.0 by a factor 
#' \eqn{ \frac{2\pi}{\sqrt{3}}}{ 2 pi/sqrt(3) }.
#' 
#' The functions \code{dkiener2347}, \code{pkiener2347} and \code{lkiener2347} 
#' have no explicit forms. Due to a poor optimization algorithm, their
#' calculations in versions < 1.9 were unreliable. In versions > 1.9, a much better 
#' algorithm was found and the optimization is conducted in a fast way to avoid 
#' a lengthy optimization. The two extreme elements (minimum, maximum) of the 
#' given \code{x} or \code{q} arguments are sent to a second order optimizer that 
#' minimize the residual error of the \code{lkiener2347} function and return the 
#' estimated lower and upper logit values. Then a sequence of logit values of 
#' length 51 times the length of \code{x} or \code{q} is generated between these 
#' lower and upper values and the corresponding quantiles are calculated with 
#' the function \code{qlkiener2347}. These 51 times more numerous quantiles are
#' then compared to the original \code{x} or \code{q} arguments and the closest 
#' values with their associated logit values are selected. The probabilities are then 
#' calculated with the function \code{invlogit} and the densities are calculated 
#' with the function \code{dlkiener2347}. The accuracy of this approach depends  
#' on the sparsity of the initial \code{x} or \code{q} sequences. A 4 digits 
#' accuracy can be expected, enough for most usages. 
#'
#' \code{qkiener2} function is defined for p in (0, 1) by: 
#'  \deqn{ 
#'    qkiener2(p,m,g,a,w) = m + \frac{\sqrt{3}}{\pi}*g*k* 
#'     \left(-exp\left(-\frac{logit(p)}{a} +\frac{logit(p)}{w}\right)\right)
#'  }{%
#'    qkiener2(p,m,g,a,w) = m + sqrt(3)/pi*g*k*(-exp(-logit(p)/a) +exp(logit(p)/w))  
#'  }
#' where k is the harmonic mean of the tail parameters \code{a} and \code{w} 
#' calculated by \eqn{k = aw2k(a, w)}.
#'
#' \code{rkiener2} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener2} is the density function calculated from the probability p. 
#' It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dpkiener2(p,m,g,a,w) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]^{-1} 
#'  }{%
#'    dpkiener2(p,m,g,a,w) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-logit(p)/a)/a +exp(logit(p)/w)/w) 
#'  }
#'
#' \code{dqkiener2} is the derivate of the quantile function calculated from 
#' the probability p. It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dqkiener2(p,m,g,a,w) = \frac{\sqrt{3}}{\pi}\frac{g}{p(1-p)}\frac{k}{2}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]
#'  }{%
#'    dqkiener2(p,m,g,a,w) = sqrt(3)/pi*g/p/(1-p)*k/2*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)
#'  }
#'
#' \code{dlkiener2} is the density function calculated from the logit of the 
#' probability lp = logit(p) defined in (-Inf, +Inf) by: 
#'  \deqn{
#'    dlkiener2(lp,m,g,a,w) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{lp}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{lp}{w}\right) \right]^{-1} 
#'  }{%
#'    dlkiener2(lp,m,g,a,w) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-lp/a)/a +exp(lp/w)/w) 
#'  }
#'
#' \code{qlkiener2} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'  \deqn{ 
#'    qlkiener2(lp,m,g,a,w) = m + \frac{\sqrt{3}}{\pi}*g*k* 
#'     \left(-exp\left(-\frac{lp}{a} +\frac{lp}{w}\right)\right)
#'  }{%
#'    qlkiener2(lp,m,g,a,w) = m + sqrt(3)/pi*g*k*(-exp(-lp/a) +exp(lp/w))  
#'  }
#' 
#' \code{varkiener2} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'  \deqn{ 
#'    varkiener2 <- if\;(p <= 0.5)\;\; (- qkiener2)\;\; else\;\; (qkiener2) 
#'  }{%
#'    varkiener2 <- if (p <= 0.5) (- qkiener2) else (qkiener2) 
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#' 
#' \code{ltmkiener2}, \code{rtmkiener2} and \code{eskiener2} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener2} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener2} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'  \deqn{ 
#'    eskiener2 <- if\;(p <= 0.5)\;\; (- ltmkiener2)\;\; else\;\; (rtmkiener2) 
#'  }{%
#'    eskiener2 <- if(p <= 0.5) (- ltmkiener2) else (rtmkiener2)
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#'
#' \code{dtmqkiener2} is the difference between the left tail mean and the quantile 
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
#' Symmetric Kiener distribution K1 \code{\link{kiener1}}, 
#' asymmetric Kiener distributions K3, K4 and K7
#' \code{\link{kiener3}}, \code{\link{kiener4}}, \code{\link{kiener7}}, 
#' conversion functions \code{\link{aw2k}}, 
#' estimation function \code{\link{fitkienerX}}, 
#' regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' 
#' require(graphics)
#' 
#' ### EXAMPLE 1
#' x <- seq(-5, 5, by = 0.1) ; x
#' pkiener2(x, m=0, g=1, a=2, w=5)
#' dkiener2(x, m=0, g=1, a=2, w=5)
#' lkiener2(x, m=0, g=1, a=2, w=5)
#' plot( x, pkiener2(x, m=0, g=1, a=2, w=5), las=1)
#' lines(x, pkiener1(x, m=0, g=1, k=9999))
#' 
#' plot( x, dkiener2(x, m=0, g=1, a=2, w=5), las=1, type="l", lwd=2)
#' lines(x, dkiener1(x, m=0, g=1, k=9999))
#' 
#' plot( x, lkiener2(x, m=0, g=1, a=2, w=5), las=1)
#' lines(x, lkiener1(x, m=0, g=1, k=9999))
#' 
#' 
#' p <- c(ppoints(11, a = 1), NA, NaN) ; p
#' qkiener2(p, a=2, w=5)
#' dpkiener2(p, a=2, w=5)
#' dqkiener2(p, a=2, w=5)
#' 
#' varkiener2(p=0.01, a=2, w=5)
#' ltmkiener2(p=0.01, a=2, w=5) 
#'  eskiener2(p=0.01, a=2, w=5) # VaR and ES should be positive
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
#' funleg <- function(xy, a) legend(xy, title = expression(alpha), legend = names(a),
#'                   lty = lty, col = col, lwd = lwd, inset = 0.02, cex = 0.8)
#' funlgt <- function(xy) legend(xy, title = "logit(p)", legend = lgt,
#'                               inset = 0.02, cex = 0.6)
#' 
#' ### EXAMPLE 2
#' ### PROBA, DENSITY, LOGIT-PROBA, LOG-DENSITY FROM x
#' x <- seq(-5, 5, by = 0.1) ; x ; length(x)
#' a <- c(9999, 9, 5, 3, 2, 1) ; names(a) <- a
#' 
#' fun1 <- function(a, x) pkiener2(x, a=a, w=5)
#' fun2 <- function(a, x) dkiener2(x, a=a, w=5)
#' fun3 <- function(a, x) lkiener2(x, a=a, w=5)
#' fun4 <- function(a, x) dkiener2(x, a=a, w=5, log=TRUE)
#' 
#' mat11 <- sapply(a, fun1, x) ; head(mat11, 10)
#' mat12 <- sapply(a, fun2, x) ; head(mat12, 10)
#' mat13 <- sapply(a, fun3, x) ; head(mat13, 10)
#' mat14 <- sapply(a, fun4, x) ; head(mat14, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(x, mat11, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="pkiener2(x, m=0, g=1, a=a, w=5)", xlab="", ylab="")
#' 	funleg("topleft", a)
#' 	matplot(x, mat12, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="dkiener2", xlab="", ylab="")
#' 	funleg("topleft", a)
#' 	matplot(x, mat13, type="l", lwd=lwd, lty=lty, col=col, yaxt="n", ylim=c(-10,10),
#' 			main="lkiener2", xlab="", ylab="")
#' 	   axis(2, at=lat, las=1)
#' 	funleg("bottomright", a)
#' 	funlgt("topleft")
#' 	matplot(x, mat14, type="l", lwd=lwd, lty=lty, col=col, ylim=c(-8,0),
#' 			main="log(dkiener2)", xlab="", ylab="")
#' 	funleg("bottom", a)
#' par(op)
#' ### END EXAMPLE 2
#' 
#' 
#' ### EXAMPLE 3
#' ### QUANTILE, DIFF-QUANTILE, DENSITY, LOG-DENSITY FROM p
#' p <- ppoints(1999, a=0) ; head(p, n=10)
#' a <- c(9999, 9, 5, 3, 2, 1) ; names(a) <- a
#' 
#' mat15 <- outer(p, a, \(p,a)  qkiener2(p, a=a, w=5)) ; head(mat15, 10)
#' mat16 <- outer(p, a, \(p,a) dqkiener2(p, a=a, w=5)) ; head(mat16, 10)
#' mat17 <- outer(p, a, \(p,a) dpkiener2(p, a=a, w=5)) ; head(mat17, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(p, mat15, type="l", xlim=c(0,1), ylim=c(-5,5), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="qkiener2(p, m=0, g=1, a=a, w=5)", xlab="", ylab="")
#' 	funleg("topleft", a)
#' 	matplot(p, mat16, type="l", xlim=c(0,1), ylim=c(0,40), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="dqkiener2", xlab="", ylab="")
#' 	funleg("top", a)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(0, 0.5), las=1,
#' 		 main="qkiener2, dpkiener2", xlab="", ylab="")
#' 	invisible(mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(mat17), 
#' 		   lwd=lwd, lty=1, col=col))
#' 	funleg("topright", a)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(-7, -0.5), las=1,
#' 		 main="qkiener2, log(dpkiener2)", xlab="", ylab="")
#' 	invisible(mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(log(mat17)), 
#' 		   lwd=lwd, lty=lty, col=col))
#' 	funleg("bottom", a)
#' par(op)
#' ### END EXAMPLE 3
#' 
#' 
#' ### EXAMPLE 4: PROCESSUS: which processus look credible?
#' ### PARAMETER a VARIES, w=4 IS CONSTANT
#' ### RUN SEED ii <- 1 THEN THE cairo_pdf CODE WITH THE 6 SEEDS
#' # cairo_pdf("K2-6x6-stocks-a.pdf")
#' # for (ii in c(1,2016,2018,2022,2023,2024)) {
#' 	ii <- 1
#' 	set.seed(ii)
#' 	p <- sample(ppoints(299, a=0), 299)
#' 	a <- c(9999, 9, 4, 3, 2, 1) ; names(a) <- a
#' 	mat18 <- outer(p, a, \(p,a)  qkiener2(p=p, g=0.85, a=a, w=4)) 
#' 	mat19 <- apply(mat18, 2, cumsum)
#' 	title <- paste0(
#' 		"stock_", ii,    
#' 	     ":  a_left = c(", paste(a[1:3], collapse = ", "), ")",
#' 	    ",  a_right = c(", paste(a[4:6], collapse = ", "), ")",
#' 		",  w = 4")
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
#' @name kiener2
NULL

#' @export
#' @rdname kiener2
dkiener2 <- function(x, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
	lp <-  lkiener2(x,  m, g, a, w)
	v  <- dlkiener2(lp, m, g, a, w)
	if(log) log(v) else v
}

#' @export
#' @rdname kiener2
pkiener2 <- function(q, m = 0, g = 1, a = 3.2, w = 3.2, 
                     lower.tail = TRUE, log.p = FALSE) {
	lp <- lkiener2(x = q, m, g, a, w)
	v  <- if(lower.tail) invlogit(lp) else 1 - invlogit(lp)
	if(log.p) log(v) else v
}

#' @export
#' @rdname kiener2
qkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	k <- aw2k(a, w)
	# since v1.9.0
	v <- m + sqrt(3)/pi*g*k*(- exp(-logit(p)/a) + exp(logit(p)/w))/2
	v
}

#' @export
#' @rdname kiener2
rkiener2 <- function(n, m = 0, g = 1, a = 3.2, w = 3.2) {
	p <- runif(n)
	v <- qkiener2(p, m, g, a, w)
	v
}

#' @export
#' @rdname kiener2
dpkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
	k <- aw2k(a, w)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g/k/( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)*2
	if(log) log(v) else v
}

#' @export
#' @rdname kiener2
dqkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
# Compute dX/dp
	k <- aw2k(a, w)
	# since v1.9.0
	v <- 1/p/(1-p)*sqrt(3)/pi*g*k*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)/2
	if(log) log(v) else v
}


#' @export
#' @rdname kiener2
lkiener2 <- function(x, m = 0, g = 1, a = 3.2, w = 3.2) {
	funnlslm2 <- function(x, m, g, a, w, lpi) { 
		fn  <- function(lp) x - qlkiener2(lp, m, g, a, w)
		opt <- minpack.lm::nls.lm(par=lpi, fn=fn)
		opt$par
	}
    fuv <- function(u, v) which.min(abs(u-v))[1]
	lk  <- lkiener1(range(x), m, g, k=aw2k(a, w))
	lp2 <- lk*exp(-lk*aw2d(a, w)*g^0.5)
	lr2 <- funnlslm2(range(x), m, g, a, w, range(lp2))
	l5  <- seq(range(lr2)[1], range(lr2)[2],
	           length.out=max(50001, length(x)*51))
	q5  <- qlkiener2(l5, m, g, a, w)
	id5 <- sapply(x, fuv, q5)
	lpi <- l5[id5]
	lpi
}

#' @export
#' @rdname kiener2
dlkiener2 <- function(lp, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
	p <- invlogit(lp)
	k <- aw2k(a, w)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g /k/( exp(-lp/a)/a + exp(lp/w)/w)*2
	if(log) log(v) else v
}

#' @export
#' @rdname kiener2
qlkiener2 <- function(lp, m = 0, g = 1, a = 3.2, w = 3.2, lower.tail = TRUE ) {
	if(!lower.tail) lp <- -lp
	k <- aw2k(a, w)
	# since v1.9.0
	v <- m + sqrt(3)/pi*g *k *(- exp(-lp/a) + exp(lp/w))/2
	v
}

#' @export
#' @rdname kiener2
varkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener2(p[i], m, g, a, w),
					  qkiener2(p[i], m, g, a, w))
	}	
	va
}

#' @export
#' @rdname kiener2
ltmkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	if (log.p) p <- exp(p)
	k  <- aw2k(a, w)
	ltm <- if (lower.tail) {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p, 1-1/a, 1+1/a)*beta(1-1/a, 1+1/a)
			+pbeta(p, 1+1/w, 1-1/w)*beta(1+1/w, 1-1/w))	
	} else {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p, 1+1/a, 1-1/a)*beta(1+1/a, 1-1/a)
			+pbeta(p, 1-1/w, 1+1/w)*beta(1-1/w, 1+1/w))
	}
	ltm
}

#' @export
#' @rdname kiener2
rtmkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	if (log.p) p <- exp(p)
	k  <- aw2k(a, w)
	rtm <- if (!lower.tail) {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1-1/a, 1+1/a)*beta(1-1/a, 1+1/a)
			+pbeta(1-p, 1+1/w, 1-1/w)*beta(1+1/w, 1-1/w))	
	} else {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1+1/a, 1-1/a)*beta(1+1/a, 1-1/a)
			+pbeta(1-p, 1-1/w, 1+1/w)*beta(1-1/w, 1+1/w))
	}
	rtm
}

#' @export
#' @rdname kiener2
dtmqkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener2(p[i], m, g, a, w, lower.tail, log.p) 
			- qkiener2(p[i], m, g, a, w, lower.tail, log.p),
			rtmkiener2(p[i], m, g, a, w, lower.tail, log.p) 
			- qkiener2(p[i], m, g, a, w, lower.tail, log.p))	
	}
	dtmq
}

#' @export
#' @rdname kiener2
eskiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                      lower.tail = TRUE, log.p = FALSE, signedES = FALSE) {
	if (log.p)      p <- exp(p)
	if(!lower.tail) p <- 1-p
	es  <- p
	for (i in seq_along(p)) {
		if (signedES) {
			es[i] <- ifelse(p[i] <= 0.5, 
					  ltmkiener2(p[i], m, g, a, w),
					  rtmkiener2(p[i], m, g, a, w))
		} else {
			es[i] <- ifelse(p[i] <= 0.5, 
				      abs(ltmkiener2(p[i], m, g, a, w)),
					  abs(rtmkiener2(p[i], m, g, a, w)))		
		}
	}
	es
}



