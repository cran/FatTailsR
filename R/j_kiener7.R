## #' @include i_kiener4.R



#' @title Asymmetric Kiener Distribution K7 (K2)
#'
#' @description
#' Density, distribution function, quantile function, random generation,
#' value-at-risk, expected shortfall (+ signed left/right tail mean) 
#' and additional formulae for asymmetric Kiener distribution K7 = K2.
#' With K7, the vector of parameters is provided as \code{coefk}, usually estimated 
#' with \code{\link{paramkienerX}} (and ~X5,~X7) or \code{\link{regkienerLX}$coefk}. 
#' Main inputs can be supplied as vector (\code{x,q,p}) and matrix (\code{coefk})  
#' and the resulting output is a matrix (useful for simulation).
#' 
#' @param    x      vector of quantiles.
#' @param    q      vector of quantiles.
#' @param    coefk  vector of 7 parameters \code{c(m,g,a,k,w,d,e)} 
#'                  or matrix with 7 columns.
#' @param    p      vector of probabilities.
#' @param    lp     vector of logit of probabilities.
#' @param    n      integer. Number of observations. If length(n) > 1, the length is  
#'                  taken to be the number required.
#' @param    log           logical. If TRUE, densities are given in log scale.
#' @param    lower.tail    logical. If TRUE, use p. If FALSE, use 1-p.
#' @param    log.p         logical. If TRUE, probabilities p are given as log(p).
#' @param    same_p        logical. If FALSE (default), random probabilies are generated 
#'                         on the fly. If TRUE, the same set of random probabilities is 
#'                         used for each line of coefk (if coefk is a matrix). 
#' @param    signedES      logical. FALSE (default) returns positive numbers for 
#'                         left and right tails. TRUE returns negative number 
#'                         (= \code{ltmkiener7}) for left tail and positive number 
#'                         (= \code{rtmkiener7}) for right tail.
#' 
#' @details
#' Kiener distributions use the following parameters, some of them being redundant. 
#' See \code{\link{aw2k}} and \code{\link{pk2pk}} for the formulas and 
#' the conversion between parameters:
#' \itemize{
#'   \item{ \code{m} (mu) is the median of the distribution. }
#'   \item{ \code{g} (gamma) is the scale parameter. }
#'   \item{ \code{a} (alpha) is the left tail parameter. } 
#'   \item{ \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} 
#'          and describes a global tail parameter. }
#'   \item{ \code{w} (omega) is the right tail parameter. } 
#'   \item{ \code{d} (delta) is the distortion parameter. }
#'   \item{ \code{e} (epsilon) is the eccentricity parameter. }
#' }
#' 
#' Kiener distribution \code{K7} is designed after \code{\link{kiener2}} 
#' but uses as input \code{coefk} rather than \code{m}, \code{g}, \code{a} 
#' and \code{w}. 
#'  
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and is linked for any value of \code{a} and \code{w} to the density at the 
#' median through the relation
#' \deqn{ g * dkiener7(x=m, coefk=coefk) = \frac{\pi}{4\sqrt{3}} \approx 0.453 }{%
#'        g * dkiener7(x=m, coefk=coefk) = pi/4/sqrt(3) = 0.453 approximatively}
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
#' \code{qkiener7} function is defined for p in (0, 1) by: 
#'  \deqn{ 
#'    qkiener7(p, coefk) = m + \frac{\sqrt{3}}{\pi}*g*k* 
#'     \left(-exp\left(-\frac{logit(p)}{a} +\frac{logit(p)}{w}\right)\right)
#'  }{%
#'    qkiener7(p, coefk) = m + sqrt(3)/pi*g*k*(-exp(-logit(p)/a) +exp(logit(p)/w))  
#'  }
#' where k is the harmonic mean of the tail parameters \code{a} and \code{w} 
#' calculated by \eqn{k = aw2k(a, w)}.
#'
#' \code{rkiener7} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener7} is the density function calculated from the probability p. 
#' It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dpkiener7(p, coefk) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]^{-1} 
#'  }{%
#'    dpkiener7(p, coefk) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-logit(p)/a)/a +exp(logit(p)/w)/w) 
#'  }
#'
#' \code{dqkiener7} is the derivate of the quantile function calculated from 
#' the probability p. It is defined for p in (0, 1) by: 
#'  \deqn{
#'    dqkiener7(p, coefk) = \frac{\sqrt{3}}{\pi}\frac{g}{p(1-p)}\frac{k}{2}
#'     \left[ +\frac{1}{a}exp\left(-\frac{logit(p)}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{logit(p)}{w}\right) \right]
#'  }{%
#'    dqkiener4(p,m,g,k,e) = sqrt(3)/pi*g/p/(1-p)*k/2*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)
#'  }
#' with \code{a} and \code{w} extracted from \code{coefk}. 
#'
#' \code{dlkiener7} is the density function calculated from the logit of the 
#' probability lp = logit(p) defined in (-Inf, +Inf). The formula is adapted 
#' from distribution K2: 
#'  \deqn{
#'    dlkiener7(lp, coefk) = \frac{\pi}{\sqrt{3}}\frac{p(1-p)}{g}\frac{2}{k}
#'     \left[ +\frac{1}{a}exp\left(-\frac{lp}{a}\right) 
#'            +\frac{1}{w}exp\left( \frac{lp}{w}\right) \right]^{-1} 
#'  }{%
#'    dlkiener7(lp, coefk) = pi/sqrt(3)*p*(1-p)/g*2/k/(+exp(-lp/a)/a +exp(lp/w)/w) 
#'  }
#'
#' \code{qlkiener7} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'  \deqn{ 
#'    qlkiener7(lp, coefk) = m + \frac{\sqrt{3}}{\pi}*g*k* 
#'     \left(-exp\left(-\frac{lp}{a} +\frac{lp}{w}\right)\right)
#'  }{%
#'    qlkiener7(lp, coefk) = m + sqrt(3)/pi*g*k*(-exp(-lp/a) +exp(lp/w))  
#'  }
#' 
#' \code{varkiener7} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'  \deqn{ 
#'    varkiener7 <- if\;(p <= 0.5)\;\; (- qkiener7)\;\; else\;\; (qkiener7) 
#'  }{%
#'    varkiener7 <- if (p <= 0.5) (- qkiener7) else (qkiener7) 
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#' 
#' \code{ltmkiener7}, \code{rtmkiener7} and \code{eskiener7} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener7} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener7} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'  \deqn{ 
#'    eskiener7 <- if\;(p <= 0.5)\;\; (- ltmkiener7)\;\; else\;\; (rtmkiener7) 
#'  }{%
#'    eskiener7 <- if(p <= 0.5) (- ltmkiener7) else (rtmkiener7)
#'  }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than \code{p}.
#'
#' \code{dtmqkiener7} is the difference between the left tail mean and the quantile 
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
#' asymmetric Kiener distributions K2, K3 and K4 
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}, 
#' conversion functions \code{\link{aw2k}}, 
#' estimation function \code{\link{paramkienerX}}, 
#' estimation function \code{\link{fitkienerX}}, 
#' regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' 
#' head(ED <- fatreturns(extractData())) 
#' (coefk  <- paramkienerX(ED, dgts = 3))  
#' x  <- -4
#' xx <- -4:4
#' p  <- 0.1
#' pp <- pprobs2
#' 
#' dkiener7(x)
#' dkiener7(x,  coefk) 
#' dkiener7(xx)
#' dkiener7(xx, coefk)
#' 
#' pkiener7(x)
#' pkiener7(x,  coefk) 
#' pkiener7(xx)
#' pkiener7(xx, coefk)
#' 
#' qkiener7(p)
#' qkiener7(p,  coefk) 
#' qkiener7(pp)
#' qkiener7(pp, coefk)
#' 
#' rkiener7(10)
#' rkiener7(10, coefk)
#' 
#' varkiener7(p)
#' varkiener7(p, coefk)
#' varkiener7(pp)
#' varkiener7(pp, coefk) 
#' 
#' ltmkiener7(p)
#' ltmkiener7(p, coefk)
#' ltmkiener7(pp)
#' ltmkiener7(pp, coefk)
#' 
#' eskiener7(p)
#' eskiener7(p, coefk)
#' eskiener7(pp)
#' eskiener7(pp, coefk) 
#' 
#' 
#' @name kiener7
NULL

#' @export
#' @rdname kiener7
dkiener7 <- function(x, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     log = FALSE) {
    checkcoefk(coefk)
	fundk7 <- function(coefk, x, log) {
		lp <-  lkiener7(x,  coefk)
		v  <- dlkiener7(lp, coefk)
		if (log) {v <- log(v)}
		v
	}
	if (dimdim1(coefk) == 1) { coefk <- t(as.matrix(coefk)) }
	v <- drop(apply(coefk, 1, fundk7, x, log))
	v
}

#' @export
#' @rdname kiener7
pkiener7 <- function(q, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	funpk7 <- function(coefk, q, lower.tail, log.p) {
		lp <- lkiener7(q, coefk)
		v  <- if (lower.tail) invlogit(lp) else 1 - invlogit(lp)
		if (log.p) v <- log(v)
		v
	}
	if (dimdim1(coefk) == 1) { coefk <- t(as.matrix(coefk)) }
	v <- drop(apply(coefk, 1, funpk7, q, lower.tail, log.p))
	v
}

#' @export
#' @rdname kiener7
qkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	funqk7 <- function(p, coefk) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		# since v1.9.0
		v <- m + sqrt(3)/pi*g*k*(- exp(-logit(p)/a) + exp(logit(p)/w))/2
		v
	}
    if (log.p)       p <- exp(p)
    if (!lower.tail) p <- 1-p
	names(p) <- getnamesk(p)$nquantk
	v        <- drop(apply(as.matrix(p), 1, funqk7, coefk))
	v
}

#' @export
#' @rdname kiener7
rkiener7 <- function(n, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     same_p = FALSE) {
    checkcoefk(coefk)
	funqk7 <- function(coefk, p, same_p) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		if (!same_p) p <- runif(length(p))
		# since v1.9.0
		v <- m + sqrt(3)/pi*g*k*(- exp(-logit(p)/a) + exp(logit(p)/w))/2
		v
	}
	if (dimdim1(coefk) == 1) coefk <- t(as.matrix(coefk))
	p <- runif(n)
	v <- drop(apply(coefk, 1, funqk7, p, same_p))
	v
}

#' @export
#' @rdname kiener7
dpkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      log = FALSE) {
    checkcoefk(coefk)
	fundpk7 <- function(p, coefk, log) {
		dck <- dimdimc(coefk)
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		# since v1.9.0
		v <- p*(1-p)*pi/sqrt(3)/g/k/( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)*2
		if(log) log(v) else v
	}
	names(p) <- getnamesk(p)$ndensk
	v        <- drop(apply(as.matrix(p), 1, fundpk7, coefk, log))
	v
}

#' @export
#' @rdname kiener7
dqkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      log = FALSE) {
# Compute dX/dp
    checkcoefk(coefk)
	fundqk7 <- function(p, coefk, log) {
		dck <- dimdimc(coefk)
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		# since v1.9.0
		v <- 1/p/(1-p)*sqrt(3)/pi*g*k*( exp(-logit(p)/a)/a + exp(logit(p)/w)/w)/2
		if(log) log(v) else v
	}
	names(p) <- getnamesk(p)$ndquantk
	v        <- drop(apply(as.matrix(p), 1, fundqk7, coefk, log))
	v
}

#' @export
#' @rdname kiener7
lkiener7 <- function(x, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0)) {
    checkcoefk(coefk)
	funlk7 <- function(coefk, x) {
		funnlslm2 <- function(x, m, g, a, w, lpi) { 
			fn  <- function(lp) x - qlkiener2(lp, m, g, a, w)
			opt <- minpack.lm::nls.lm(par=lpi, fn=fn)
			opt$par
		}
		fuv <- function(u, v) {
			z <- abs(u-v)
			which(z == min(z))[1]
		}
		m   <- coefk[1]
		g   <- coefk[2]
		a   <- coefk[3]
		w   <- coefk[4]
		k   <- coefk[5]
		d   <- coefk[6]
		lk  <- lkiener1(range(x), m, g, k=k)
		lp2 <- lk*exp(-lk*d*g^0.5)
		lr2 <- funnlslm2(range(x), m, g, a, w, range(lp2))
		l5  <- seq(range(lr2)[1], range(lr2)[2],
	               length.out=max(50001, length(x)*51))
		q5  <- qlkiener2(l5, m, g, a, w)
		id5 <- sapply(x, fuv, q5)
		lpi <- l5[id5]
		lpi
	}
	if (dimdim1(coefk) == 1) coefk <- t(as.matrix(coefk))
	drop(apply(coefk, 1, funlk7, x))
}


#' @export
#' @rdname kiener7
dlkiener7 <- function(lp, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      log = FALSE) {
	checkcoefk(coefk)
	fundqk7 <- function(lp, coefk, log) {
		dck <- dimdimc(coefk)
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		p <- invlogit(lp)
		# since v1.9.0
		v <- p*(1-p)*pi/sqrt(3)/g /k/( exp(-lp/a)/a + exp(lp/w)/w)*2
		if(log) log(v) else v
	}
	names(lp) <- paste0("lp", round(lp,1))
	v  <- drop(apply(as.matrix(lp), 1, fundqk7, coefk, log))
	v
}

#' @export
#' @rdname kiener7
qlkiener7 <- function(lp, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      lower.tail = TRUE ) {
    checkcoefk(coefk)
	funqlk7 <- function(lp, coefk, lower.tail) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		if (!lower.tail) {lp <- -lp}
		# since v1.9.0
		v <- m + sqrt(3)/pi*g *k *(- exp(-lp/a) + exp(lp/w))/2
		v
	}
	names(lp) <- paste0("lp", round(lp, 1))
	v  <- drop(apply(as.matrix(lp), 1, funqlk7, coefk, lower.tail))
	v
}

#' @export
#' @rdname kiener7
varkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	funvark7 <- function(p, coefk) {
		if (p <= 0.5) {	- qkiener7(p, coefk)
		}        else {   qkiener7(p, coefk)
		}
	}
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	names(p) <- getnamesk(p)$nvark
	va <- drop(apply(as.matrix(p), 1, funvark7, coefk))
	va
}

#' @export
#' @rdname kiener7
ltmkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                       lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	funltmk7 <- function(p, coefk, lower.tail) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
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
    if (log.p)  {p <- exp(p)}
	names(p) <- getnamesk(p)$nltmk
	ltm      <- drop(apply(as.matrix(p), 1, funltmk7, coefk, lower.tail))
	ltm
}

#' @export
#' @rdname kiener7
rtmkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                       lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	funrtmk7 <- function(p, coefk, lower.tail) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
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
    if (log.p)  {p <- exp(p)}
	names(p) <- getnamesk(p)$nrtmk
	rtm      <- drop(apply(as.matrix(p), 1, funrtmk7, coefk, lower.tail))
	rtm
}

#' @export
#' @rdname kiener7
dtmqkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	fundtmqk7 <- function(p, coefk, lower.tail, log.p) {
		dtmq <- if (p <= 0.5) { 
			ltmkiener7(p, coefk, lower.tail, log.p) 
			- qkiener7(p, coefk, lower.tail, log.p)
		} else {
			rtmkiener7(p, coefk, lower.tail, log.p) 
			- qkiener7(p, coefk, lower.tail, log.p)
		}
		dtmq
	}
	names(p) <- getnamesk(p)$ndtmqk
	dtmq <- drop(apply(as.matrix(p), 1, fundtmqk7, 
	                   coefk, lower.tail, log.p))
	dtmq
}

#' @export
#' @rdname kiener7
eskiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      lower.tail = TRUE, log.p = FALSE, signedES = FALSE) {
    checkcoefk(coefk)
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	funesk7 <- function(p, coefk, signedES) {
		if (signedES) {
			if (p <= 0.5) {	ltmkiener7(p, coefk)
			}        else { rtmkiener7(p, coefk)
			}
		} else {
			if (p <= 0.5) {	abs(ltmkiener7(p, coefk))
			}        else { abs(rtmkiener7(p, coefk))
			}
		}
	}
	names(p) <- getnamesk(p)$nesk
	es <- drop(apply(as.matrix(p), 1, funesk7, coefk, signedES))
	es
}



