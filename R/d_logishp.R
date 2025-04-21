## #' @include c_trigohp.R



#' @title Logit and Invlogit Functions
#' 
#' @description
#' The logit and invlogit functions, widely used in this package, are wrappers  
#' of \code{\link{qlogis}} and \code{\link{plogis}} functions. 
#' 
#' 
#' @param    p    numeric. one value or a vector between 0 and 1.
#' @param    x    numeric. one value or a vector of numerics.
#' 
#' @details      \code{logit} function is defined for p in (0, 1) by: 
#'               \deqn{ logit(p) = log( p/(1-p) ) }
#' 
#'               \code{invlogit} function is defined for x in (-Inf, +Inf) by: 
#'               \deqn{ invlogit(x) = exp(x)/(1+exp(x)) = plogis(x) }
#' 
#' @examples     
#' logit( c(ppoints(11, a = 1), NA, NaN) )
#' invlogit( c(-Inf, -10:10, +Inf, NA, NaN) )
#' 
#' @export
#' @name logit
logit    <- function(p) { qlogis(p) } 
#' @export
#' @rdname logit
invlogit <- function(x) { plogis(x) } 


#' @title The Standardized Logistic Distribution
#'
#' @description
#' Density, distribution function, quantile function, random generation, 
#' value-at-risk, left-tail mean, right-tail mean, expected shortfall
#' for the standardized logistic distribution, equivalent to 
#' \code{dpqrlogis(..., scale = g*sqrt(3)/pi)}.
#'
#' @param x	vector of quantiles.
#' @param q	vector of quantiles.
#' @param k	numeric. The tail parameter, preferably strictly positive.  
#'          Can be a vector (see details).
#' @param p	vector of probabilities.
#' @param lp	vector of logit of probabilities.
#' @param n	number of observations. If length(n) > 1, the length is 
#'                                  taken to be the number required.
#' @param log	boolean.
#' @param m     numeric. a central parameter (also used in model K1, K2, K3 and K4).
#' @param g     numeric. a scale parameter (also used in model K1, K2, K3 and K4).
#' @param lower.tail  logical. If TRUE, use p. If FALSE, use 1-p.
#' @param log.p       logical. If TRUE, probabilities p are given as log(p).
#'
#' @details
#' \code{dlogisst} function (log is available) is defined for 
#'                 x in (-Inf, +Inf) by: 
#'   \deqn{ dlogisst(x, m, g) = 
#'         stats::dlogis(x, location = m, scale = g*sqrt(3)/pi) }
#' \code{plogisst} function is defined for q in (-Inf, +Inf) by: 
#'   \deqn{ plogisst(q, m, g) = 
#'          stats::plogis(q, location = m, scale = g*sqrt(3)/pi) }
#' \code{qlogisst} function is defined for p in (0, 1) by: 
#'    \deqn{ qlogisst(p, m, g) = 
#'           stats::qlogis(p, location = m, scale = g*sqrt(3)/pi) }
#' \code{rlogisst} function generates \code{n} random values.
#' 
#' In addition to the classical formats, the prefixes dp, dq, l, dl, ql 
#' are also provided:
#' 
#' \code{dplogisst} function (log is available) is defined for p in (0, 1) by: 
#'      \deqn{ dplogisst(p, m, g) =  p*(1-p)/g*pi/sqrt(3) + m*0 }
#' \code{dqlogisst} function (log is available) is defined for p in (0, 1) by: 
#'      \deqn{ dqlogisst(p, m, g) = 1/p/(1-p)*sqrt(3)/pi*g + m*0 }
#' \code{llogisst} function is defined for x in (-Inf, +Inf) by: 
#'      \deqn{ llogisst(x, m, g) = (x-m)/g*pi/sqrt(3) }
#' \code{dllogisst} function is defined for lp = logit(p) in (-Inf, +Inf) by : 
#'      \deqn{ dllogisst(lp, m, g) = p*(1-p)/g*pi/sqrt(3) }
#' \code{qllogisst} function is defined for lp = logit(p) in (-Inf, +Inf) by : 
#'      \deqn{ qllogisst(lp, m, g) = m + sqrt(3)/pi*g }
#'
#' If k is a vector, then the use of the function \code{\link[base]{outer}} 
#' is recommanded.
#' 
#' Functions \code{eslogis} is the expected shortfall of the logistic function 
#' (times a factor 2). 
#' When \code{p<=0.5}, it is equivalent (times -1) to the left tail mean \code{ltmlogisst}. 
#' When \code{p>0.5}, it is equivalent to the right tail mean \code{rtmlogisst}. 
#' \code{ltmlogisst} and \code{rtmlogisst} are used to calculate the \code{h} parameter 
#' in \code{\link{hkiener1}}, \code{hkiener2}, \code{hkiener3}, \code{hkiener4}.
#' 
#' @seealso Kiener distribution K1 \code{\link{kiener1}} which has 
#' location (\code{m}) and scale (\code{g}) parameters. 
#' 
#' @name logisst
NULL

#' @export
#' @rdname logisst
dlogisst <- function(x, m = 0, g = 1, log = FALSE) { 
	stats::dlogis(x, location = m, scale = g*sqrt(3)/pi, log = log)
}

#' @export
#' @rdname logisst
plogisst <- function(q, m = 0, g = 1, lower.tail = TRUE, log.p = FALSE) {
	stats::plogis(q, location = m, scale = g*sqrt(3)/pi, 
	             lower.tail = lower.tail, log.p = log.p) 
} 
		  
#' @export
#' @rdname logisst
qlogisst <- function(p, m = 0, g = 1, lower.tail = TRUE, log.p = FALSE) {
	stats::qlogis(p, location = m, scale = g*sqrt(3)/pi, 
	              lower.tail = TRUE, log.p = FALSE)
}
	  
#' @export
#' @rdname logisst
rlogisst <- function(n, m = 0, g = 1) { 
	stats::rlogis(runif(n), location = m, scale = g*sqrt(3)/pi) 
}

#' @export
#' @rdname logisst
dplogisst <- function(p, m = 0, g = 1, log = FALSE) {
	v <- p*(1-p)/g*pi/sqrt(3) + m*0
	if(log) log(v) else v
}

#' @export
#' @rdname logisst
dqlogisst <- function(p, m = 0, g = 1, k = 3.2, log = FALSE) {
	v <- 1/p/(1-p)*sqrt(3)/pi*g + m*0
	if(log) log(v) else v
}

#' @export
#' @rdname logisst
llogisst <- function(x, m = 0, g = 1) { 
	(x-m)/g*pi/sqrt(3)
}

#' @export
#' @rdname logisst
dllogisst <- function(lp, m = 0, g = 1, k = 3.2, log = FALSE) {
	p <- invlogit(lp)
	v <- p*(1-p)/g*pi/sqrt(3)
	if(log) log(v) else v
}

#' @export
#' @rdname logisst
qllogisst <- function(lp, m = 0, g = 1, k = 3.2, lower.tail = TRUE ) {
	if(!lower.tail) lp <- -lp
	m + sqrt(3)/pi*g
}

#' @export
#' @rdname logisst
varlogisst <- function(p, m = 0, g = 1, k = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qlogisst(p[i], m, g, k),
					  qlogisst(p[i], m, g, k))
	}
	va
}

#' @export
#' @rdname logisst
ltmlogisst <- function(p, m = 0, g = 1, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	ltm  <- m + sqrt(3)/pi*g/p *( (1-p)*log(1-p) + p*log(p) )
	ltm
}

#' @export
#' @rdname logisst
rtmlogisst <- function(p, m = 0, g = 1, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	rtm  <- m - sqrt(3)/pi*g/(1-p) *( (1-p)*log(1-p) + p*log(p) )
	rtm
}

#' @export
#' @rdname logisst
eslogisst <- function(p, m = 0, g = 1, lower.tail = TRUE, log.p = FALSE) {
	if(log.p)       p <- exp(p)
	if(!lower.tail) p <- 1-p
	es  <- p
	for (i in seq_along(p)) {
		es[i] <- ifelse(p[i] <= 0.5, 
					- ltmlogisst(p[i], m, g),
					  rtmlogisst(p[i], m, g))
	}
	es
}



