## #' @include j_kiener7.R



#' @title Quantile (VaR) and Expected Shortfall Corrective Functions 
#'
#' @description
#' Quantile functions (or VaR) and Expected Shortfall of Kiener distributions 
#' K1, K2, K3 and K4, usually calculated at pprobs2 = c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99), 
#' can be expressed as: 
#' \enumerate{
#'   \item Quantile of the logit function multiplied by a fat tail 
#'         (c)orrective function \code{ckiener1234};
#'   \item Expected s(h)ortfall of the logistic function multiplied 
#'         by a corrective function \code{hkiener1234}. 
#' }
#' Both functions \code{ckiener1234} and \code{hkiener1234} are independant from 
#' the scale parameter \code{g} and are indirect measures of the tail curvature. 
#' A value close to \code{1} indicates a model similar to the logistic function with  
#' almost no curvature and probably parameter \code{k > 8}. When \code{k} (or \code{a,w}) 
#' decreases, the values of \code{c} and \code{h} increase and indicate some more 
#' pronounced symmetric or asymmetric curvature, depending on values of \code{d,e}. 
#' Note that if \code{(min(a,k,w) <= 1)}, \code{ckiener1234} still exists but 
#' the expected shortfall and \code{hkiener1234} become undefined (\code{NA}).
#' 
#' Some financial applications use threshold values on \code{ckiener1234} or 
#' \code{hkiener1234} to select or discard stocks over time as they become 
#' less or more risky. 
#' 
#' @param    p	         numeric or vector of probabilities. 
#' @param    m           numeric. parameter m used in model K1, K2, K3 and K4.
#' @param    g           numeric. parameter g used in model K1, K2, K3 and K4.
#' @param    k	         numeric. parameter k used in model K1, K3 and K4. 
#' @param    a	         numeric. parameter a used in model K2.
#' @param    w	         numeric. parameter w used in model K2.
#' @param    d           numeric. parameter d used in model K3.
#' @param    e	         numeric. parameter e used in model K4.
#' @param    lower.tail  logical. If TRUE, use p. If FALSE, use 1-p.
#' @param    log.p       logical. If TRUE, probabilities p are given as log(p).
#' @param    coefk	     vector or 7 columns-matrix representing parameters 
#'                       \code{c(m,g,a,k,w,d,e)} obtained from \code{\link{paramkienerX}}.
#' 
#' 
#' @seealso  
#' \code{\link{logit}}, \code{\link{qkiener1}}, \code{\link{qkiener2}}, 
#' \code{\link{qkiener3}}, \code{\link{qkiener4}}, \code{\link{fitkienerX}}.
#' 
#' @name ckiener1234
NULL
#' @export
#' @rdname ckiener1234
hkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	h   <- p
	for (i in seq_along(p)) {
		h[i] <- ifelse(p[i] <= 0.5, 
			(ltmkiener1(p[i], m, g, k) - m) / (ltmlogisst(p[i], m, g) - m),
			(rtmkiener1(p[i], m, g, k) - m) / (rtmlogisst(p[i], m, g) - m))
	}
	h
}
#' @export
#' @rdname ckiener1234
hkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	h   <- p
	for (i in seq_along(p)) {
		h[i] <- ifelse(p[i] <= 0.5, 
			(ltmkiener2(p[i], m, g, a, w) - m) / (ltmlogisst(p[i], m, g) - m),
			(rtmkiener2(p[i], m, g, a, w) - m) / (rtmlogisst(p[i], m, g) - m))
	}
	h
}
#' @export
#' @rdname ckiener1234
hkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	h   <- p
	for (i in seq_along(p)) {
		h[i] <- ifelse(p[i] <= 0.5, 
			(ltmkiener3(p[i], m, g, k, d) - m) / (ltmlogisst(p[i], m, g) - m),
			(rtmkiener3(p[i], m, g, k, d) - m) / (rtmlogisst(p[i], m, g) - m))
	}
	h
}
#' @export
#' @rdname ckiener1234
hkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	h   <- p
	for (i in seq_along(p)) {
		h[i] <- ifelse(p[i] <= 0.5, 
			(ltmkiener4(p[i], m, g, k, e) - m) / (ltmlogisst(p[i], m, g) - m),
			(rtmkiener4(p[i], m, g, k, e) - m) / (rtmlogisst(p[i], m, g) - m))
	}
	h
}
#' @export
#' @rdname ckiener1234
hkiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     lower.tail = TRUE, log.p = FALSE) {
	checkcoefk(coefk)
    dck <- dimdimc(coefk)
    m <- switch(dck, "1" = coefk[1], "2" = coefk[,1]) 
    g <- switch(dck, "1" = coefk[2], "2" = coefk[,2])
    a <- switch(dck, "1" = coefk[3], "2" = coefk[,3]) 
    w <- switch(dck, "1" = coefk[5], "2" = coefk[,5])
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	h <- p
	names(h) <- getnamesk(p)$nhesk 
	for (i in seq_along(p)) {
		h[i] <- ifelse(p[i] <= 0.5, 
			(ltmkiener2(p[i], m, g, a, w) - m) / (ltmlogisst(p[i], m, g) - m),
			(rtmkiener2(p[i], m, g, a, w) - m) / (rtmlogisst(p[i], m, g) - m))
	}
	h
}

#' @export
#' @rdname ckiener1234
ckiener1 <- function(p, k = 3.2, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	lp <- logit(p)
	z  <- k/lp * sinh(lp/k)
	z[which(z == "NaN")] <- 1
	z
}
#' @export
#' @rdname ckiener1234
ckiener2 <- function(p, a = 3.2, w = 3.2, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	lp <- logit(p)
	k <- aw2k(a, w)
	e <- aw2e(a, w)
	z <- k/lp * sinh(lp/k) * exp(lp/k *e)
	z[which(z == "NaN")] <- 1
	z
}
#' @export
#' @rdname ckiener1234
ckiener3 <- function(p, k = 3.2, d = 0, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	lp <- logit(p)
	z  <- k/lp * sinh(lp/k) * exp(lp * d)
	z[which(z == "NaN")] <- 1
	z
}
#' @export
#' @rdname ckiener1234
ckiener4 <- function(p, k = 3.2, e = 0, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	lp <- logit(p)
	z  <- k/lp * sinh(lp/k) * exp(lp/k *e)
	z[which(z == "NaN")] <- 1
	z
}
#' @export
#' @rdname ckiener1234
ckiener7 <- function(p, coefk = c(0, 1, 3.2, 3.2, 3.2, 0, 0), 
                     lower.tail = TRUE, log.p = FALSE) {
    checkcoefk(coefk)
    k  <- switch(dimdimc(coefk), "1" = coefk[4], "2" = coefk[,4]) 
    e  <- switch(dimdimc(coefk), "1" = coefk[7], "2" = coefk[,7])
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
    lp <- logit(p) 
    z  <- k/lp * sinh(lp/k) * exp(lp/k *e)
    z[which(z == "NaN")] <- 1
	z
}



