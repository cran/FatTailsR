## #' @include h_kiener3.R



#' @title Asymmetric Kiener Distribution K4
#'
#' @description
#' Density, distribution function, quantile function, random generation,
#' value-at-risk, expected shortfall (+ signed left/right tail mean) 
#' and additional formulae for asymmetric Kiener distribution K4.
#'
#' @param    x    vector of quantiles.
#' @param    q    vector of quantiles.
#' @param    m    numeric. The median.
#' @param    g    numeric. The scale parameter, preferably strictly positive.
#' @param    k    numeric. The tail parameter, preferably strictly positive.
#' @param    e    numeric. The eccentricity parameter between left and right tails.
#' @param    p    vector of probabilities.
#' @param    lp   vector of logit of probabilities.
#' @param    n    number of observations. If length(n) > 1, the length is  
#'                taken to be the number required.
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
#' Kiener distributions \code{K4(m, g, k, e, ...)} are distributions 
#' with asymmetrical left and right fat tails described by a global tail 
#' parameter \code{k} and an eccentricity parameter \code{e}. 
#'  
#' Distributions K3 (\code{\link{kiener3}}) 
#' with parameters \code{k} (kappa) and \code{d} (delta) and
#' distributions K4 (\code{\link{kiener4}})
#' with parameters \code{k} (kappa) and \code{e} (epsilon))
#' have been created to disantangle the parameters 
#' \code{a} (alpha) and \code{w} (omega) of distributions K2
#' (\code{\link{kiener2}}). 
#' The tiny difference between distributions K3 and K4 (\eqn{d = e/k}) 
#' has not yet been fully evaluated. Both should be tested at that moment.
#' 
#' \code{k} is the harmonic mean of \code{a} and \code{w} and represents a 
#' global tail parameter.
#'
#' \code{e} is an eccentricity parameter between the left tail parameter
#' \code{a} and the right tail parameter \code{w}.
#' It verifies the inequality: \eqn{-1 < e < 1} 
#' (whereas \code{d} of distribution K3 verifies \eqn{-k < d < k}).
#' The conversion functions (see \code{\link{aw2k}}) are:
#'
#' \deqn{1/k = (1/a + 1/w)/2 }
#' \deqn{  e = (a - w)/(a + w) }
#' \deqn{  a = k/(1 - e) }
#' \deqn{  w = k/(1 + e) }
#'
#' \code{e} (and \code{d}) should be of the same sign than the skewness. 
#' A negative value \eqn{ e < 0 } implies \eqn{ a < w } and indicates a left 
#' tail heavier than the right tail. A positive value \eqn{ e > 0 } implies 
#' \eqn{ a > w } and a right tail heavier than the left tail.    
#' 
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and is linked for any value of \code{k} and \code{e} to the density at the 
#' median through the relation
#' \deqn{ g * dkiener4(x=m, g=g, e=e) = \frac{\pi}{4\sqrt{3}} \approx 0.453 }{%
#'        g * dkiener4(x=m, g=g, e=e) = pi/4/sqrt(3) = 0.453 approximatively}
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
#' \code{qkiener4} function is defined for p in (0, 1) by: 
#'  \deqn{ 
#'    qkiener4(p,m,g,k,e) = m + \frac{\sqrt{3}}{\pi}*g*k*       
#'     \sinh\left(\frac{logit(p)}{k}\right)*exp\left(\frac{e}{k} logit(p)\right)
#'  }{%
#'    qkiener4(p,m,g,k,e) = m + sqrt(3)/pi*g*k*sinh(logit(p)/k)*exp(e/k*logit(p)) 
#'  }
#'
#' \code{rkiener4} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener4} is the density function calculated from the probability p. 
#' The formula is adapted from distribution K2. It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dpkiener4(p,m,g,k,e) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]^{-1} 
#'  }{%
#'    dpkiener4(p,m,g,k,e) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-logit(p)/a)/a +exp(logit(p)/w)/w) 
#'  }
#' with \code{a} and \code{w} defined from \code{k} and \code{e}.
#'
#' \code{dqkiener4} is the derivate of the quantile function calculated from 
#' the probability p. The formula is adapted from distribution K2. 
#' It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dqkiener4(p,m,g,k,e) = \frac{\sqrt{3}}{\pi}\frac{g}{p(1-p)}\frac{k}{2}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]
#'  }{%
#'    dqkiener4(p,m,g,k,e) = sqrt(3)/pi*g/p/(1-p)*k/2*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)
#'  }
#' with \code{a} and \code{w} defined with the formula presented above. 
#'
#' \code{dlkiener4} is the density function calculated from the logit of the
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
#' \code{qlkiener4} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'  \deqn{ 
#'    qlkiener4(lp,m,g,k,e) = m + \frac{\sqrt{3}}{\pi}*g*k*       
#'     \sinh\left(\frac{lp}{k}\right)*exp\left(\frac{e}{k} lp\right)
#'  }{%
#'    qlkiener4(lp,m,g,k,e) = m + sqrt(3)/pi*g*k*sinh(lp/k)*exp(e/k*lp) 
#'  }
#' 
#' \code{varkiener4} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'  \deqn{ 
#'    varkiener4 <- if\;(p <= 0.5)\;\; (- qkiener4)\;\; else\;\; (qkiener4) 
#'  }{%
#'    varkiener4 <- if (p <= 0.5) (- qkiener4) else (qkiener4) 
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#' 
#' \code{ltmkiener4}, \code{rtmkiener4} and \code{eskiener4} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener4} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener4} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'  \deqn{ 
#'    eskiener4 <- if\;(p <= 0.5)\;\; (- ltmkiener4)\;\; else\;\; (rtmkiener4) 
#'  }{%
#'    eskiener4 <- if(p <= 0.5) (- ltmkiener4) else (rtmkiener4)
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#'
#' \code{dtmqkiener4} is the difference between the left tail mean and the quantile 
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
#' asymmetric Kiener distributions K2, K3 and K7
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener7}}, 
#' conversion functions \code{\link{aw2k}}, 
#' estimation function \code{\link{fitkienerX}},
#  regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' require(graphics)
#' 
#' ### EXAMPLE 1
#' x <- seq(-5, 5, by = 0.1) ; round(x, 2)
#' round(pkiener4(x, m=0, g=1, k=4, e=0.1), 4)
#' round(dkiener4(x, m=0, g=1, k=4, e=0.1), 4)
#' round(lkiener4(x, m=0, g=1, k=4, e=0.1), 4)
#' 
#' plot( x, pkiener4(x, m=0, g=1, k=9999, e=0), las=1, type="l", lwd=2)
#' lines(x, pkiener4(x, m=0, g=1, k=4, e=0.5), col="red")
#' lines(x, pkiener4(x, m=0, g=1, k=4, e=1), lwd=1)  # e in [-1, 1]
#' 
#' plot( x, dkiener4(x, m=0, g=1, k=9999, e=0), las=1, type="l", lwd=2, ylim=c(0,0.6))
#' lines(x, dkiener4(x, m=0, g=1, k=4, e=0.5), col="red")
#' lines(x, dkiener4(x, m=0, g=1, k=4, e=1), lwd=1)
#' 
#' plot( x, lkiener4(x, m=0, g=1, k=9999, e=0), las=1, type="l", lwd=2)
#' lines(x, lkiener4(x, m=0, g=1, k=4, e=0.05), col="green")
#' lines(x, lkiener4(x, m=0, g=1, k=4, e=0.5), col="red")
#' lines(x, lkiener4(x, m=0, g=1, k=4, e=1), lwd=1)
#' 
#' 
#' p <- c(ppoints(11, a = 1), NA, NaN) ; p
#' qkiener4(p, k=4, e=0.5)
#' dpkiener4(p, k=4, e=0.5)
#' dqkiener4(p, k=4, e=0.5)
#' 
#' varkiener4(p=0.01, k=4, e=0.5)
#' ltmkiener4(p=0.01, k=4, e=0.5) 
#'  eskiener4(p=0.01, k=4, e=0.5) # VaR and ES should be positive
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
#' funleg <- function(xy, e) legend(xy, title = expression(epsilon), legend = names(e),
#'                   lty = lty, col = col, lwd = lwd, inset = 0.02, cex = 0.8)
#' funlgt <- function(xy) legend(xy, title = "logit(p)", legend = lgt,
#'                               inset = 0.02, cex = 0.6)
#' 
#' ### EXAMPLE 2
#' ### PROBA, DENSITY, LOGIT-PROBA, LOG-DENSITY FROM x
#' x <- seq(-5, 5, by = 0.1) ; head(x, 10)
#' e <- c(-0.5, -0.25, 0, 0.25, 0.50, 1) ; names(e) <- e
#' 
#' fun1 <- function(e, x) pkiener4(x, k=4, e=e)
#' fun2 <- function(e, x) dkiener4(x, k=4, e=e)
#' fun3 <- function(e, x) lkiener4(x, k=4, e=e)
#' fun4 <- function(e, x) dkiener4(x, k=4, e=e, log=TRUE)
#' 
#' mat11 <- sapply(e, fun1, x) ; head(mat11, 10)
#' mat12 <- sapply(e, fun2, x) ; head(mat12, 10)
#' mat13 <- sapply(e, fun3, x) ; head(mat13, 10)
#' mat14 <- sapply(e, fun4, x) ; head(mat14, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(x, mat11, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="pkiener4(x, m=0, g=1, k=4, e=e)", xlab="", ylab="")
#' 	funleg("topleft", e)
#' 	matplot(x, mat12, type="l", lwd=lwd, lty=lty, col=col, 
#' 			main="dkiener4", xlab="", ylab="")
#' 	funleg("topleft", e)
#' 	matplot(x, mat13, type="l", lwd=lwd, lty=lty, col=col, yaxt="n", ylim=c(-9,9),
#' 			main="lkiener4", xlab="", ylab="")
#' 	   axis(2, at=lat, las=1)
#' 	funleg("bottomright", e)
#' 	funlgt("topleft")
#' 	matplot(x, mat14, type="l", lwd=lwd, lty=lty, col=col, ylim=c(-8,0),
#' 			main="log(dkiener4)", xlab="", ylab="")
#' 	funleg("bottom", e)
#' par(op)
#' ### END EXAMPLE 2
#' 
#' 
#' ### EXAMPLE 3
#' ### QUANTILE, DIFF-QUANTILE, DENSITY, LOG-DENSITY FROM p
#' p <- ppoints(1999, a=0) ; head(p, n=10)
#' e <- c(-0.5, -0.25, 0, 0.25, 0.50, 1) ; names(e) <- e
#' 
#' mat15 <- outer(p, e, \(p,e)  qkiener4(p, k=4, e=e)) ; head(mat15, 10)
#' mat16 <- outer(p, e, \(p,e) dqkiener4(p, k=4, e=e)) ; head(mat16, 10)
#' mat17 <- outer(p, e, \(p,e) dpkiener4(p, k=4, e=e)) ; head(mat17, 10)
#' 
#' op <- par(mfcol = c(2,2), mar = c(2.5,3,1.5,1), las=1)
#' 	matplot(p, mat15, type="l", xlim=c(0,1), ylim=c(-5,5), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="qkiener4(p, m=0, g=1, k=4, e=e)", xlab="", ylab="")
#' 	funleg("topleft", e)
#' 	matplot(p, mat16, type="l", xlim=c(0,1), ylim=c(0,40), 
#'             lwd=lwd, lty=lty, col=col, las=1,
#' 			main="dqkiener4", xlab="", ylab="")
#' 	funleg("top", e)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(0, 0.6), las=1,
#' 		 main="qkiener4, dpkiener4", xlab="", ylab="")
#' 	invisible(mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(mat17), 
#' 		   lwd=lwd, lty=1, col=col))
#' 	funleg("topright", e)
#' 	plot(NA, NA, xlim=c(-5, 5), ylim=c(-7, -0.5), las=1,
#' 		 main="qkiener4, log(dpkiener4)", xlab="", ylab="")
#' 	invisible(mapply(matlines, x=as.data.frame(mat15), y=as.data.frame(log(mat17)), 
#' 		   lwd=lwd, lty=lty, col=col))
#' 	funleg("bottom", e)
#' par(op)
#' ### END EXAMPLE 3
#' 
#' 
#' ### EXAMPLE 4: PROCESSUS: which processus look credible?
#' ### PARAMETER e VARIES, k=4 IS CONSTANT
#' ### RUN SEED ii <- 1 THEN THE cairo_pdf CODE WITH THE 6 SEEDS
#' # cairo_pdf("K4-6x6-stocks-e.pdf")
#' # for (ii in c(1,2016,2018,2022,2023,2024)) {
#' 	ii <- 1
#' 	set.seed(ii)
#' 	p <- sample(ppoints(299, a=0), 299)
#' 	e <- c(-0.1, -0.05, 0, 0.05, 0.1, 0.25) ; names(e) <- e
#' 	mat18 <- outer(p, e, \(p,e)  qkiener4(p=p, g=0.85, k=4, e=e)) 
#' 	mat19 <- apply(mat18, 2, cumsum)
#' 	title <- paste0(
#' 		"stock_", ii,
#' 	     ": k = 4", 
#' 		 ",  e_left = c(", paste(e[1:3], collapse = ", "), ")",
#' 	    ",  e_right = c(", paste(e[4:6], collapse = ", "), ")")
#' 	plot.ts(mat19, ann=FALSE, las=1, 
#' 			mar.multi=c(0,3,0,1), oma.multi=c(3,0,3,0.5))
#' 	mtext(title, outer = TRUE, line=-1.5, font=2)
#' 	plot.ts(mat18, ann=FALSE, las=1, 
#' 			mar.multi=c(0,3,0,1), oma.multi=c(3,0,3,0.5))
#' 	mtext(title, outer=TRUE, line=-1.5, font=2)
#' # }
#' # dev.off()
#' 
#' 
#' ### PARAMETER k VARIES, e=0.05 IS CONSTANT
#' # cairo_pdf("K4-6x6-stocks-k.pdf", width=11)
#' # for (ii in c(1,2016,2018,2022,2023,2024)) {
#' 	ii <- 1
#' 	set.seed(ii)
#' 	p <- sample(ppoints(299, a=0), 299)
#' 	k <- c(9999, 6, 4, 3, 2, 1) ; names(k) <- k
#' 	mat20 <- outer(p, k, \(p,k)  qkiener4(p=p, g=0.85, k=k, e=0.05)) 
#' 	mat21 <- apply(mat20, 2, cumsum)
#' 	title <- paste0(
#'  		"stock_", ii,    
#'  	     ":  k_left = c(", paste(k[1:3], collapse = ", "), ")",
#'  	    ",  k_right = c(", paste(k[4:6], collapse = ", "), ")",
#' 		",  e = 0.05")
#' 	plot.ts(mat21, ann=FALSE, las=1, nc=2,
#' 			mar.multi=c(0,3,0,1), oma.multi=c(3,0,3,0.5))
#' 	mtext(title, outer = TRUE, line=-1.5, font=2)
#' 	plot.ts(mat20, ann=FALSE, las=1, nc=2,
#' 			mar.multi=c(0,3,0,1), oma.multi=c(3,0,3,0.5))
#' 	mtext(title, outer=TRUE, line=-1.5, font=2)
#' # }
#' # dev.off()
#' ### END EXAMPLE 4
#' 
#' 
#' @name kiener4
NULL

#' @export 
#' @rdname kiener4
dkiener4 <- function(x, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
	lp <-  lkiener4(x,  m, g, k, e)
	v  <- dlkiener4(lp, m, g, k, e)
	if(log) log(v) else v
}

#' @export
#' @rdname kiener4
pkiener4 <- function(q, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	lp <- lkiener4(x = q, m, g, k, e)
	if(lower.tail) v <- invlogit(lp) else v <- 1 - invlogit(lp)
	if(log.p) log(v) else v
}

#' @export
#' @rdname kiener4
qkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	# since v1.9.0
	v <- m + sqrt(3)/pi*g* k*sinh(logit(p)/k) *exp(e/k*logit(p))
	v
}

#' @export
#' @rdname kiener4
rkiener4 <- function(n, m = 0, g = 1, k = 3.2, e = 0) {
	p <- runif(n)
	v <- qkiener4(p, m, g, k, e)
	v
}

#' @export
#' @rdname kiener4
dpkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
	a <- ke2a(k, e)
	w <- ke2w(k, e)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g/k/( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)*2
	if(log) log(v) else v
}

#' @export
#' @rdname kiener4
dqkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
# Compute dX/dp
	a <- ke2a(k, e)
	w <- ke2w(k, e)
	# since v1.9.0
	v <- 1/p/(1-p)*sqrt(3)/pi*g*k*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)/2
	if(log) log(v) else v
}

## voir fonction NLSstClosestX
#' @export
#' @rdname kiener4
lkiener4 <- function(x, m = 0, g = 1, k = 3.2, e = 0) { 
	funnlslm4 <- function(x, m, g, k, e, lpi) { 
		fn  <- function(lp) x - qlkiener4(lp, m, g, k, e)
		opt <- minpack.lm::nls.lm(par=lpi, fn=fn)
		opt$par
	}
    fuv <- function(u, v) which.min(abs(u-v))[1]
	lk  <- lkiener1(range(x), m, g, k)
	lp2 <- lk*exp(-lk*e/k*g^0.5)
	lr2 <- funnlslm4(range(x), m, g, k, e, range(lp2))
	l5  <- seq(range(lr2)[1], range(lr2)[2],
	           length.out=max(50001, length(x)*51))
	q5  <- qlkiener4(l5, m, g, k, e)
	id5 <- sapply(x, fuv, q5)
	lpi <- l5[id5]
	lpi
}

#' @export
#' @rdname kiener4
dlkiener4 <- function(lp, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
	p <- invlogit(lp)
	a <- ke2a(k, e)
	w <- ke2w(k, e)
	# since v1.9.0
	v <- p*(1-p)*pi/sqrt(3)/g /k/( exp(-lp/a)/a + exp(lp/w)/w)*2
	if(log) log(v) else v
} 

#' @export
#' @rdname kiener4
qlkiener4 <- function(lp, m = 0, g = 1, k = 3.2, e = 0, lower.tail = TRUE ) {
	if(!lower.tail) lp <- -lp
	# since v1.9.0
	v  <- m + sqrt(3)/pi*g *k *sinh(lp/k) *exp(e/k*lp)
	v
} 

#' @export
#' @rdname kiener4
varkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener4(p[i], m, g, k, e),
					  qkiener4(p[i], m, g, k, e))
	}
	va
}

#' @export
#' @rdname kiener4
ltmkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p  <- exp(p)
	ltm <- if (lower.tail) {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
	        -pbeta(p, 1-(1-e)/k, 1+(1-e)/k) * beta(1-(1-e)/k, 1+(1-e)/k) 
	        +pbeta(p, 1+(1+e)/k, 1-(1+e)/k) * beta(1+(1+e)/k, 1-(1+e)/k))
	} else {
		# m+g*k/p*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/p*(
			-pbeta(p, 1+(1-e)/k, 1-(1-e)/k)*beta(1+(1-e)/k, 1-(1-e)/k)
			+pbeta(p, 1-(1+e)/k, 1+(1+e)/k)*beta(1-(1+e)/k, 1+(1+e)/k))
	}
	ltm
}

#' @export
#' @rdname kiener4
rtmkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                       lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p <- exp(p)
	rtm <- if (!lower.tail) {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1-(1-e)/k, 1+(1-e)/k)*beta(1-(1-e)/k, 1+(1-e)/k)
			+pbeta(1-p, 1+(1+e)/k, 1-(1+e)/k)*beta(1+(1+e)/k, 1-(1+e)/k))	
	} else {
		# m+g*k/(1-p)*(
		# since v1.9.2
		m+sqrt(3)/pi/2*g*k/(1-p)*(
			-pbeta(1-p, 1+(1-e)/k, 1-(1-e)/k)*beta(1+(1-e)/k, 1-(1-e)/k)
			+pbeta(1-p, 1-(1+e)/k, 1+(1+e)/k)*beta(1-(1+e)/k, 1+(1+e)/k))
	}
	rtm
}

#' @export
#' @rdname kiener4
dtmqkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener4(p[i], m, g, k, e, lower.tail, log.p) 
			- qkiener4(p[i], m, g, k, e, lower.tail, log.p),
			rtmkiener4(p[i], m, g, k, e, lower.tail, log.p) 
			- qkiener4(p[i], m, g, k, e, lower.tail, log.p)) 
	}
	dtmq
}

#' @export
#' @rdname kiener4
eskiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                      lower.tail = TRUE, log.p = FALSE, signedES = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	es  <- p
	for (i in seq_along(p)) {
		if (signedES) {
			es[i] <- ifelse(p[i] <= 0.5, 
					  ltmkiener4(p[i], m, g, k, e),
					  rtmkiener4(p[i], m, g, k, e))
		} else {
			es[i] <- ifelse(p[i] <= 0.5, 
				  abs(ltmkiener4(p[i], m, g, k, e)),
				  abs(rtmkiener4(p[i], m, g, k, e)))
		}
	}
	es
}

