## #' @include g_kiener2.R



#' @title Asymmetric Kiener Distribution K3
#'
#' @description
#' Density, distribution function, quantile function, random generation,
#' value-at-risk, expected shortfall (+ signed left/right tail mean) 
#' and additional formulae for asymmetric Kiener distribution K3.
#'
#' @param    x    vector of quantiles.
#' @param    q    vector of quantiles.
#' @param    m    numeric. The median.
#' @param    g    numeric. The scale parameter, preferably strictly positive.
#' @param    k    numeric. The tail parameter, preferably strictly positive.
#' @param    d    numeric. The distortion parameter between left and right tails.
#' @param    p    vector of probabilities.
#' @param    lp   vector of logit of probabilities.
#' @param    n    number of observations. If length(n) > 1, the length is  
#'                taken to be the number required.
#' @param    log           logical. If TRUE, densities are given in log scale.
#' @param    lower.tail    logical. If TRUE, use p. If FALSE, use 1-p.
#' @param    log.p         logical. If TRUE, probabilities p are given as log(p).
#' @param    signedES      logical. FALSE (default) returns positive numbers for 
#'                         left and right tails. TRUE returns negative number 
#'                         (= \code{ltmkiener3}) for left tail and positive number 
#'                         (= \code{rtmkiener3}) for right tail.
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
#' Kiener distributions \code{K3(m, g, k, d, ...)} are distributions 
#' with asymmetrical left and right fat tails described by a global tail 
#' parameter \code{k} and a distortion parameter \code{d}. 
#' 
#' Distributions K3 (\code{\link{kiener3}}) 
#' with parameters \code{k} (kappa) and \code{d} (delta) and
#' distributions K4 (\code{\link{kiener4}})
#' with parameters \code{k} (kappa) and \code{e} (epsilon))
#' have been created to disantangle the parameters 
#' \code{a} (alpha) and \code{w} (omega) of distributions of 
#' distribution K2 (\code{\link{kiener2}}). 
#' The tiny difference between distributions K3 and K4 (\eqn{d = e/k}) 
#' has not yet been fully evaluated. Both should be tested at that moment.
#' 
#' \code{k} is the harmonic mean of \code{a} and \code{w} and represents a 
#' global tail parameter.
#'
#' \code{d} is a distortion parameter between the left tail parameter
#' \code{a} and the right tail parameter \code{w}.
#' It verifies the inequality: \eqn{-k < d < k} 
#' (whereas \code{e} of distribution K4 verifies \eqn{-1 < e < 1}).
#' The conversion functions (see \code{\link{aw2k}}) are:
#'
#' \deqn{\frac{1}{k} = \frac{1}{2}\left( \frac{1}{a} + \frac{1}{w}\right) }{%
#'                      1/k = (1/a + 1/w)/2}
#' \deqn{ d = \frac{1}{2}\left(-\frac{1}{a} + \frac{1}{w}\right) }{%
#'                      d = (-1/a + 1/w)/2} 
#' \deqn{\frac{1}{a} = \frac{1}{k} - d}{1/a = 1/k - d} 
#' \deqn{\frac{1}{w} = \frac{1}{k} + d}{1/w = 1/k + d}
#'
#' \code{d} (and \code{e}) should be of the same sign than the skewness. 
#' A negative value \eqn{ d < 0 } implies \eqn{ a < w } and indicates a left  
#' tail heavier than the right tail. A positive value \eqn{ d > 0 } implies 
#' \eqn{ a > w } and a right tail heavier than the left tail.  
#' 
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and is linked for any value of \code{k} and \code{d} to the density at the 
#' median through the relation
#' \deqn{ g * dkiener3(x=m, g=g, d=d) = \frac{\pi}{4\sqrt{3}} \approx 0.453 }{%
#'        g * dkiener3(x=m, g=g, d=d) = pi/4/sqrt(3) = 0.453 approximatively}
#'
#' When \code{k = Inf}, \code{g} is very close to \code{sd(x)}. 
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
#' \code{qkiener3} function is defined for p in (0, 1) by:
#'  \deqn{ 
#'    qkiener3(p,m,g,k,d) = m + \frac{\sqrt{3}}{\pi}*g*k*       
#'        \sinh\left(\frac{logit(p)}{k}\right)*exp\left(d*logit(p)\right)
#'  }{%
#'    qkiener3(p,m,g,k,d) = m + sqrt(3)/pi*g*k*sinh(logit(p)/k)*exp(d*logit(p)) 
#'  }
#'
#' \code{rkiener3} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener3} is the density function calculated from the probability p. 
#' The formula is adapted from distribution K2. It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dpkiener3(p,m,g,k,d) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]^{-1} 
#'  }{%
#'    dpkiener3(p,m,g,k,d) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-logit(p)/a)/a +exp(logit(p)/w)/w) 
#'  }
#' with \code{a} and \code{w} defined from \code{k} and \code{d} 
#' with the formula presented above.
#'
#' \code{dqkiener3} is the derivate of the quantile function calculated from 
#' the probability p. The formula is adapted from distribution K2. 
#' It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dqkiener3(p,m,g,k,d) = \frac{\sqrt{3}}{\pi}\frac{g}{p(1-p)}\frac{k}{2}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]
#'  }{%
#'    dqkiener3(p,m,g,k,d) = sqrt(3)/pi*g/p/(1-p)*k/2*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)
#'  }
#' with \code{a} and \code{w} defined above. 
#'
#' \code{dlkiener3} is the density function calculated from the logit of the 
#' probability lp = logit(p) defined in (-Inf, +Inf). The formula is adapted 
#' from distribution K2: 
#'  \deqn{
#'    dlkiener2(lp,m,g,k,e) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{lp}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{lp}{w}\right) \right]^{-1} 
#'  }{%
#'    dlkiener2(lp,m,g,k,e) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-lp/a)/a +exp(lp/w)/w) 
#'  }
#' with \code{a} and \code{w} defined above. 
#'
#' \code{qlkiener3} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'  \deqn{ 
#'    qlkiener3(lp,m,g,k,d) = m + \frac{\sqrt{3}}{\pi}*g*k*       
#'        \sinh\left(\frac{lp}{k}\right)*exp\left(d*lp\right)
#'  }{%
#'    qlkiener3(lp,m,g,k,d) = m + sqrt(3)/pi*g*k*sinh(lp/k)*exp(d*lp) 
#'  }
#' 
#' \code{varkiener3} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'  \deqn{ 
#'    varkiener3 <- if\;(p <= 0.5)\;\; (- qkiener3)\;\; else\;\; (qkiener3) 
#'  }{%
#'    varkiener3 <- if (p <= 0.5) (- qkiener3) else (qkiener3) 
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#' 
#' \code{ltmkiener3}, \code{rtmkiener3} and \code{eskiener3} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener3} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener3} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'  \deqn{ 
#'    eskiener3 <- if\;(p <= 0.5)\;\; (- ltmkiener3)\;\; else\;\; (rtmkiener3) 
#'  }{%
#'    eskiener3 <- if(p <= 0.5) (- ltmkiener3) else (rtmkiener3)
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#'
#' \code{dtmqkiener3} is the difference between the left tail mean and the quantile 
#' when (p <= 0.5) and the difference between the right tail mean and the quantile 
#' when (p > 0.5). It is in quantile unit and is an indirect measure of the tail curvature.
#' 
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014. Download it from:  
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
#' asymmetric Kiener distributions K2, K4 and K7
#' \code{\link{kiener2}}, \code{\link{kiener4}}, \code{\link{kiener7}}, 
#' conversion functions \code{\link{aw2k}},  
#' estimation function \code{\link{fitkienerX}}, 
#' regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' 
#' require(graphics)
#' 
#' ### EXAMPLE 1
#' x <- (-15:15)/3 ; round(x, 2)
#' round(pkiener3(x, m=0, g=1, k=4, d=0.1), 4)
#' round(dkiener3(x, m=0, g=1, k=4, d=0.1), 4)
#' round(lkiener3(x, m=0, g=1, k=4, d=0.1), 4)
#' 
#' plot( x, pkiener3(x, m=0, g=1, k=9999, d=0), las=1, type="l", lwd=2)
#' lines(x, pkiener3(x, m=0, g=1, k=4, d=0.1), col="red")
#' lines(x, pkiener3(x, m=0, g=1, k=4, d=0.25), lwd=1)  # d in [-1:k, 1:k]
#' 
#' plot( x, dkiener3(x, m=0, g=1, k=9999, d=0), las=1, type="l", lwd=2, ylim=c(0,0.6))
#' lines(x, dkiener3(x, m=0, g=1, k=4, d=0.1), col="red")
#' lines(x, dkiener3(x, m=0, g=1, k=4, d=0.25), lwd=1)
#' 
#' plot( x, lkiener3(x, m=0, g=1, k=9999, d=0), las=1, type="l", lwd=2)
#' lines(x, lkiener3(x, m=0, g=1, k=4, d=0.1), col="red")
#' lines(x, lkiener3(x, m=0, g=1, k=4, d=0.25), lwd=1)
#' 
#' 
#' p <- c(ppoints(11, a = 1), NA, NaN) ; p
#' qkiener3(p, k=4, d=0.1)
#' dpkiener3(p, k=4, d=0.1)
#' dqkiener3(p, k=4, d=0.1)
#' 
#' varkiener3(p=0.01, k=4, d=0.1)
#' ltmkiener3(p=0.01, k=4, d=0.1) 
#'  eskiener3(p=0.01, k=4, d=0.1) # VaR and ES should be positive
#' ### END EXAMPLE 1
#' 
#' 
#' ### PREPARE THE GRAPHICS FOR EXAMPLES 2 AND 3
#' xx  <- c(-4,-2, 0, 2, 4)
#' lty <- c( 3, 2, 1, 4, 5, 1)
#' lwd <- c( 1, 1, 2, 1, 1, 1)
#' col <- c("cyan3","green3","black","dodgerblue2","purple2","brown3")
#' lat <- c(-6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9)
#' lgt <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", "logit(0.95)   = 2.9", 
#'          "logit(0.50)   = 0", "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'          "logit(0.001) = -6.9  ")
#' funleg <- function(xy, d) legend(xy, title = expression(delta), legend = names(d),
#'                   lty = lty, col = col, lwd = lwd, inset = 0.02, cex = 0.8)
#' funlgt <- function(xy) legend(xy, title = "logit(p)", legend = lgt,
#'                               inset = 0.02, cex = 0.6)
#' 
#' ### EXAMPLE 2
#' ### PROBA, DENSITY, LOGIT-PROBA, LOG-DENSITY FROM x
#' x <- seq(-5, 5, by = 0.1) ; head(x, 10)
#' d <- c(-0.15, -0.1, 0, 0.1, 0.15, 0.25) ; names(d) <- d
#' 
#' fun1 <- function(d, x) pkiener3(x, k=4, d=d)
#' fun2 <- function(d, x) dkiener3(x, k=4, d=d)
#' fun3 <- function(d, x) lkiener3(x, k=4, d=d)
#' fun4 <- function(d, x) dkiener3(x, k=4, d=d, log=TRUE)
#' 
#' mat11 <- sapply(d, fun1, x) ; head(mat11, 10)
#' mat12 <- sapply(d, fun2, x) ; head(mat12, 10)
#' mat13 <- sapply(d, fun3, x) ; head(mat13, 10)
#' mat14 <- sapply(d, fun4, x) ; head(mat14, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(x, mat11, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="pkiener3(x, m=0, g=1, k=4, d=d)", xlab="", ylab="")
#' 	funleg("topleft", d)
#' 	matplot(x, mat12, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="dkiener3", xlab="", ylab="")
#' 	funleg("topleft", d)
#' 	matplot(x, mat13, type="l", lwd=lwd, lty=lty, col=col, yaxt="n", ylim=c(-9,9),
#' 			main="lkiener3", xlab="", ylab="")
#' 	   axis(2, at=lat, las=1)
#' 	funleg("bottomright", d)
#' 	funlgt("topleft")
#' 	matplot(x, mat14, type="l", lwd=lwd, lty=lty, col=col, ylim=c(-8,0),
#' 			main="log(dkiener3)", xlab="", ylab="")
#' 	funleg("bottom", d)
#' par(op)
#' ### END EXAMPLE 2
#' 
#' 
#' ### EXAMPLE 3
#' ### QUANTILE, DIFF-QUANTILE, DENSITY, LOG-DENSITY FROM p
#' p <- ppoints(1999, a=0) ; head(p, n=10)
#' d <- c(-0.15, -0.1, 0, 0.1, 0.15, 0.25) ; names(d) <- d
#' 
#' mat15 <- outer(p, d, \(p,d)  qkiener3(p, k=4, d=d)) ; head(mat15, 10)
#' mat16 <- outer(p, d, \(p,d) dqkiener3(p, k=4, d=d)) ; head(mat16, 10)
#' mat17 <- outer(p, d, \(p,d) dpkiener3(p, k=4, d=d)) ; head(mat17, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(p, mat15, type="l", xlim=c(0,1), ylim=c(-5,5), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="qkiener3(p, m=0, g=1, k=4, d=d)", xlab="", ylab="")
#' 	funleg("topleft", d)
#' 	matplot(p, mat16, type="l", xlim=c(0,1), ylim=c(0,40), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="dqkiener3", xlab="", ylab="")
#' 	funleg("top", d)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(0, 0.6), las=1,
#' 		 main="qkiener3, dpkiener3", xlab="", ylab="")
#' 	mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(mat17), 
#' 		   lwd=lwd, lty=1, col=col)
#' 	funleg("topright", d)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(-7, -0.5), las=1,
#' 		 main="qkiener3, log(dpkiener3)", xlab="", ylab="")
#' 	mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(log(mat17)), 
#' 		   lwd=lwd, lty=lty, col=col)
#' 	funleg("bottom", d)
#' par(op)
#' ### END EXAMPLE 3
#' 
#' 
#' ### EXAMPLE 4: PROCESSUS: which processus look credible?
#' ### PARAMETER d VARIES, k=4 IS CONSTANT
#' ### RUN SEED ii <- 1 THEN THE cairo_pdf CODE WITH THE 6 SEEDS
#' # cairo_pdf("K3-6x6-stocks-d.pdf")
#' # for (ii in c(1,2016,2018,2022,2023,2024)) {
#' 	ii <- 1
#' 	set.seed(ii)
#' 	p <- sample(ppoints(299, a=0), 299)
#' 	d <- c(-0.1, -0.05, 0, 0.05, 0.1, 0.25) ; names(d) <- d
#' 	mat18 <- outer(p, d, \(p,d)  qkiener3(p=p, g=0.85, k=4, d=d)) 
#' 	mat19 <- apply(mat18, 2, cumsum)
#' 	title <- paste0(
#' 		"stock_", ii,    
#' 	     ":  k = 4",
#' 	     ",  d_left = c(", paste(d[1:3], collapse = ", "), ")",
#' 	    ",  d_right = c(", paste(d[4:6], collapse = ", "), ")")
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
#' @name kiener3
NULL

#' @export 
#' @rdname kiener3
dkiener3 <- function(x, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
	lp <-  lkiener3(x,  m, g, k, d)
	v  <- dlkiener3(lp, m, g, k, d)
	if(log) log(v) else v
} 

#' @export
#' @rdname kiener3
pkiener3 <- function(q, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	lp <- lkiener3(x = q, m, g, k, d)
	v  <- if(lower.tail) invlogit(lp) else 1 - invlogit(lp)
	if(log.p) log(v) else v
} 

#' @export
#' @rdname kiener3
qkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	# since v1.9.0
	v <- m + sqrt(3)/pi*g* k*sinh(logit(p)/k) *exp(d*logit(p))
	v
} 

#' @export
#' @rdname kiener3
rkiener3 <- function(n, m = 0, g = 1, k = 3.2, d = 0) {
	p <- runif(n)
	v <- qkiener3(p, m, g, k, d)
	v
} 

#' @export
#' @rdname kiener3
dpkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
	a <- kd2a(k, d)
	w <- kd2w(k, d)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g/k/( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)*2
	if(log) log(v) else v
} 

#' @export
#' @rdname kiener3
dqkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
# Compute dX/dp
	a <- kd2a(k, d)
	w <- kd2w(k, d)
	# since v1.9.0
	v <- 1/p/(1-p)*sqrt(3)/pi*g*k*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)/2
	if(log) log(v) else v
} 

#' @export
#' @rdname kiener3
lkiener3 <- function(x, m = 0, g = 1, k = 3.2, d = 0) { 
	funnlslm3 <- function(x, m, g, k, d, lpi) { 
		fn  <- function(lp) x - qlkiener3(lp, m, g, k, d)
		opt <- minpack.lm::nls.lm(par=lpi, fn=fn)
		opt$par
	}
    fuv <- function(u, v) which.min(abs(u-v))[1]
	lk  <- lkiener1(range(x), m, g, k)
	lp2 <- lk*exp(-lk*d*g^0.5)
	lr2 <- funnlslm3(range(x), m, g, k, d, range(lp2))
	l5  <- seq(range(lr2)[1], range(lr2)[2],
	           length.out=max(50001, length(x)*51))
	q5  <- qlkiener3(l5, m, g, k, d)
	id5 <- sapply(x, fuv, q5)
	lpi <- l5[id5]
	lpi
}

#' @export
#' @rdname kiener3
dlkiener3 <- function(lp, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
	p <- invlogit(lp)
	a <- kd2a(k, d)
	w <- kd2w(k, d)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g /k/( exp(-lp/a)/a + exp(lp/w)/w)*2
	if(log) log(v) else v
}

#' @export
#' @rdname kiener3
qlkiener3 <- function(lp, m = 0, g = 1, k = 3.2, d = 0, lower.tail = TRUE ) {
	if(!lower.tail) lp <- -lp
	# since v1.9.0
	v  <- m + sqrt(3)/pi*g *k *sinh(lp/k) *exp(d*lp)
	v
} 

#' @export
#' @rdname kiener3
varkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener3(p[i], m, g, k, d),
					  qkiener3(p[i], m, g, k, d))
	}
	va
}

#' @export
#' @rdname kiener3
ltmkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p  <- exp(p)
	ltm <- if (lower.tail) {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p,1+d-1/k, 1-d+1/k)*beta(1+d-1/k, 1-d+1/k)
			+pbeta(p,1+d+1/k, 1-d-1/k)*beta(1+d+1/k, 1-d-1/k))	
	} else {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p,1-d+1/k, 1+d-1/k)*beta(1-d+1/k, 1+d-1/k)
			+pbeta(p,1-d-1/k, 1+d+1/k)*beta(1-d-1/k, 1+d+1/k))
	}
	ltm
} 

#' @export
#' @rdname kiener3
rtmkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                       lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p <- exp(p)
	rtm <- if (!lower.tail) {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1+d-1/k, 1-d+1/k)*beta(1+d-1/k, 1-d+1/k)
			+pbeta(1-p, 1+d+1/k, 1-d-1/k)*beta(1+d+1/k, 1-d-1/k))	
	} else {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1-d+1/k, 1+d-1/k)*beta(1-d+1/k, 1+d-1/k)
			+pbeta(1-p, 1-d-1/k, 1+d+1/k)*beta(1-d-1/k, 1+d+1/k))
	}
	rtm
}

#' @export
#' @rdname kiener3
dtmqkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                        lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener3(p[i], m, g, k, d, lower.tail, log.p) 
			- qkiener3(p[i], m, g, k, d, lower.tail, log.p),
			rtmkiener3(p[i], m, g, k, d, lower.tail, log.p) 
			- qkiener3(p[i], m, g, k, d, lower.tail, log.p)) 
	}
	dtmq
}

#' @export
#' @rdname kiener3
eskiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                      lower.tail = TRUE, log.p = FALSE, signedES = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	es  <- p
	for (i in seq_along(p)) {
		if (signedES) {
			es[i] <- ifelse(p[i] <= 0.5, 
					  ltmkiener3(p[i], m, g, k, d),
					  rtmkiener3(p[i], m, g, k, d))
		} else {
			es[i] <- ifelse(p[i] <= 0.5, 
				      abs(ltmkiener3(p[i], m, g, k, d)),
					  abs(rtmkiener3(p[i], m, g, k, d)))		
		}
	}
	es
}


