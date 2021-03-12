

#' @include i_kiener4.R



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
#' The d, p functions have no explicit forms. They are provided here for 
#' convenience. They are estimated from a reverse optimization on the quantile 
#' function and can be (very) slow, depending the number of points to estimate. 
#' We recommand to use the quantile function as much as possible. 
#' WARNING: Results may become inconsistent when \code{a} or \code{w} are
#' smaller than 1. 
#'
#' \code{qkiener7} function is defined for p in (0, 1) by: 
#'   \deqn{ qkiener7(p, coefk) = 
#'                   m + g * k * (- exp(-logit(p)/a) + exp(logit(p)/w) ) }
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
#'   \deqn{ dpkiener7(p, coefk) = 
#'          p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w }
#'
#' \code{dqkiener7} is the derivate of the quantile function calculated from 
#' the probability p. It is defined for p in (0, 1) by: 
#'   \deqn{ dqkiener7(p, coefk) = 
#'          k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w ) }
#'
#' \code{lkiener7} function is estimated from a reverse optimization and can 
#' be (very) slow depending the number of points to estimate. Initialization 
#' is done by assuming a symmetric distribution \code{\link{lkiener1}} 
#' around the harmonic mean \code{k}, then optimization is performed to 
#' take into account the true values \code{a} and \code{w}. 
#' The result can be then compared to the empirical probability logit(p). 
#' WARNING: Results may become inconsistent when \code{a} or \code{w} are
#' smaller than 1. 
#'
#' \code{dlkiener7} is the density function calculated from the logit of the 
#' probability lp = logit(p).  
#' it is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ dlkiener7(lp, coefk) = 
#'           p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w ) }
#'
#' \code{qlkiener7} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ qlkiener7(lp, coefk) = 
#'           m + g * k * ( - exp(-lp/a) + exp(lp/w) ) }
#' 
#' \code{varkiener7} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'    \deqn{ varkiener7 <- if(p <= 0.5) { - qkiener7 } else { qkiener7 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#' 
#' \code{ltmkiener7}, \code{rtmkiener7} and \code{eskiener7} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener7} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener7} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'    \deqn{ eskiener7 <- if(p <= 0.5) { - ltmkiener7 } else { rtmkiener7 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
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
	return(v)
	}
	if (dimdim1(coefk) == 1) { coefk <- t(as.matrix(coefk)) }
	v <- drop(apply(coefk, 1, fundk7, x, log))
return(v)
}

#' @export
#' @rdname kiener7
pkiener7 <- function(q, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
	funpk7 <- function(coefk, q, lower.tail, log.p) {
		lp <- lkiener7(q, coefk)
		v  <- if (lower.tail) {invlogit(lp)} else {1 - invlogit(lp)}
		if (log.p) {v <- log(v)}
	return(v)
	}
	if (dimdim1(coefk) == 1) { coefk <- t(as.matrix(coefk)) }
	v <- drop(apply(coefk, 1, funpk7, q, lower.tail, log.p))
return(v)
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
		v <- m + g * k * (- exp(-logit(p)/a) + exp(logit(p)/w) )
	return(v)
	}
    if (log.p)       {p <- exp(p)}
    if (!lower.tail) {p <- 1-p}
	names(p) <- getnamesk(p)$nquantk
	v        <- drop(apply(as.matrix(p), 1, funqk7, coefk))
return(v)
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
		if (!same_p) { p <- runif(length(p)) }
		v <- m + g * k * (- exp(-logit(p)/a) + exp(logit(p)/w) )
	return(v)
	}
	if (dimdim1(coefk) == 1) { coefk <- t(as.matrix(coefk)) }
	p <- runif(n)
	v <- drop(apply(coefk, 1, funqk7, p, same_p))
return(v)
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
		v <- p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
		if (log) {v <- log(v)}
	return(v)
	}
	names(p) <- getnamesk(p)$ndensk
	v        <- drop(apply(as.matrix(p), 1, fundpk7, coefk, log))
return(v)
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
		v <- k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
		if (log) {v <- log(v)}
	return(v)
	}
	names(p) <- getnamesk(p)$ndquantk
	v        <- drop(apply(as.matrix(p), 1, fundqk7, coefk, log))
return(v)
}

#' @export
#' @rdname kiener7
lkiener7 <- function(x, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0)) { 
    checkcoefk(coefk)
	funlk7 <- function(coefk, x) {
		## lkiener1 = k * asinh((x - m)/g / 2 / k) 
		lp.ini <- lkiener1(x, coefk[1], coefk[2], coefk[4])
		f      <- function(lp) sum( (x - qlkiener7(lp, coefk))^2 )
		lp.fin <- nlm(f, lp.ini)
		v      <- lp.fin$estimate
	return(v)
	}
	if (dimdim1(coefk) == 1) { coefk <- t(as.matrix(coefk)) }
	v <- drop(apply(coefk, 1, funlk7, x))
return(v)
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
		v <- p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w )
		if (log) {v <- log(v)}
	return(v)
	}
	names(lp) <- paste0("lp",round(lp,1))
	v  <- drop(apply(as.matrix(lp), 1, fundqk7, coefk, log))
return(v)
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
		v <- m + g * k * ( - exp(-lp/a) + exp(lp/w) )
	return(v)
	}
	names(lp) <- paste0("lp", round(lp, 1))
	v  <- drop(apply(as.matrix(lp), 1, funqlk7, coefk, lower.tail))
return(v)
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
    if (log.p)       {p <- exp(p)}
    if (!lower.tail) {p <- 1-p}
	names(p) <- getnamesk(p)$nvark
	va <- drop(apply(as.matrix(p), 1, funvark7, coefk))
return(va)
}

#' @export
#' @rdname kiener7
ltmkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                       lower.tail = TRUE, log.p = FALSE) {
	funltmk7 <- function(p, coefk, lower.tail) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		ltm <- if (lower.tail) {
				m+g*k/p*(
					-pbeta(p, 1-1/a, 1+1/a)*beta(1-1/a, 1+1/a)
					+pbeta(p, 1+1/w, 1-1/w)*beta(1+1/w, 1-1/w))	
			} else {
				m+g*k/p*(
					-pbeta(p, 1+1/a, 1-1/a)*beta(1+1/a, 1-1/a)
					+pbeta(p, 1-1/w, 1+1/w)*beta(1-1/w, 1+1/w))
			}
	return(ltm)
	}
    checkcoefk(coefk)
    if (log.p)  {p <- exp(p)}
	names(p) <- getnamesk(p)$nltmk
	ltm      <- drop(apply(as.matrix(p), 1, funltmk7, coefk, lower.tail))
return(ltm)
}

#' @export
#' @rdname kiener7
rtmkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                       lower.tail = TRUE, log.p = FALSE) {
	funrtmk7 <- function(p, coefk, lower.tail) {
		dck <- dimdimc(coefk)
		m <- switch(dck, "1" = coefk[1], "2" = coefk[,1])
		g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
		a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
		k <- switch(dck, "1" = coefk[4], "2" = coefk[,4])
		w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
		rtm <- if (!lower.tail) {
			m+g*k/(1-p)*(
				-pbeta(1-p, 1-1/a, 1+1/a)*beta(1-1/a, 1+1/a)
				+pbeta(1-p, 1+1/w, 1-1/w)*beta(1+1/w, 1-1/w))	
		} else {
			m+g*k/(1-p)*(
				-pbeta(1-p, 1+1/a, 1-1/a)*beta(1+1/a, 1-1/a)
				+pbeta(1-p, 1-1/w, 1+1/w)*beta(1-1/w, 1+1/w))
		}
	return(rtm)
	}
    checkcoefk(coefk)
    if (log.p)  {p <- exp(p)}
	names(p) <- getnamesk(p)$nrtmk
	rtm      <- drop(apply(as.matrix(p), 1, funrtmk7, coefk, lower.tail))
return(rtm)
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
	return(dtmq)
	}
	names(p) <- getnamesk(p)$ndtmqk
	dtmq <- drop(apply(as.matrix(p), 1, fundtmqk7, 
	                   coefk, lower.tail, log.p))
return(dtmq)
}

#' @export
#' @rdname kiener7
eskiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                      lower.tail = TRUE, log.p = FALSE, signedES = FALSE) {
    checkcoefk(coefk)
    if (log.p)       {p <- exp(p)}
    if (!lower.tail) {p <- 1-p}
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
return(es)
}



