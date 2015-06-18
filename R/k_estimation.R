

#' @include j_tailskiener.R



#' @title Five, Seven, Eleven Probabilities
#'
#' @description
#' Extract from a dataset \code{X} a vector of 5, 7 or 11 probabilities 
#' \code{c(1-p1, 0.25, 0.50, 0.75, p1)},   
#' \code{c(1-p1, 1-p2, 0.25, 0.50, 0.75, p2, p1)},   
#' \code{c(1-p1, 1-p2, 1-p3, 0.25, 0.35, 0.50, 0.65, 0.75, p3, p2, p1)}. 
#' p1, p2 and p3 are the most extreme probabilities of the dataset \code{X} 
#' with values finishing either by \code{pp950} or \code{pp975} or \code{pp990}. 
#' 
#' Parameters names take the parameter values if \code{parnames = TRUE}.
#' 
#' @param    X	       numeric. Vector of quantiles.
#' @param    parnames  boolean. Parameter vector with or without names.
#' 
#' @return  
#' A vector of 11 probabilities 
#' \code{c(1-p1, 1-p2, 1-p3, 0.25, 0.35, 0.50, 0.65, 0.75, p3, p2, p1)}.
#' 
#' A vector of 7 probabilities 
#' \code{c(1-p1, 1-p2, 0.25, 0.50, 0.75, p2, p1)}.
#' 
#' A vector of 5 probabilities 
#' \code{c(1-p1, 0.25, 0.50, 0.75, p1)}.
#' 
#' 
#' @seealso      \code{\link{estimkienerX}}.
#' 
#' @examples     
#' 
#' ## DS
#' DS   <- getDSdata()
#' for (j in 1:16) { print(round(elevenprobs(DS[[j]]), 6)) }
#' 
#' ## Choose j in 1:16
#' j    <- 1
#' X    <- sort(DS[[j]])
#' leX  <- logit(eX <- elevenprobs(X))
#' lpX  <- logit(ppoints(length(X), a = 0))
#' plot(X, lpX)
#' abline(h = leX, lty = 3)
#' mtext(eX, side = 4, at = leX, las = 1, line = -3.3)
#' ## end
#' 
#' 
#' @export
#' @name elevenprobs
elevenprobs <- function(X, parnames = FALSE) { 
	listp <- c(as.numeric(c(0.25, 0.5, 1) %o% 10^(-8:-1)), 0.15, 0.20)
	ltX   <- length(X)
	p1    <- listp[findInterval(1/ltX, listp) + 1]
	p2    <- listp[findInterval(1/ltX, listp) + 2]
	p3    <- listp[findInterval(1/ltX, listp) + 3]
	p11   <- c(p1, p2, p3, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p3, 1-p2, 1-p1)
	names(p11) <- as.character(p11)
	if (!parnames) { names(p11) <- NULL }
return(p11)
}
#' 
#' 
#' @export
#' @rdname elevenprobs
sevenprobs <- function(X, parnames = FALSE) { 
	listp <- c(as.numeric(c(0.25, 0.5, 1) %o% 10^(-8:-1)), 0.15, 0.20)
	ltX   <- length(X)
	p1    <- listp[findInterval(1/ltX, listp) + 1]
	p2    <- listp[findInterval(1/ltX, listp) + 2]
	p7    <- c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)
	names(p7) <- NULL
	if (!parnames) { names(p7) <- NULL }
return(p7)
}
#' 
#' 
#' @export
#' @rdname elevenprobs
fiveprobs <- function(X, parnames = FALSE) { 
	listp <- c(as.numeric(c(0.25, 0.5, 1) %o% 10^(-8:-1)), 0.15, 0.20)
	ltX   <- length(X)
	p1    <- listp[findInterval(1/ltX, listp) + 1]
	p5    <- c(p1, 0.25, 0.50, 0.75, 1-p1)
	names(p5) <- NULL
	if (!parnames) { names(p5) <- NULL }
return(p5)
}


#' @title Local Function estimkappa6
#'
#' @description
#' Internal function to estimate parameter k (kappa). 
#' 
#' @param    lg	      logit of probability close to the median. 
#' @param    lp	      logit of extreme probability. 
#' @param    dx1p	  numeric. Delta of quantile m - x(1-p).
#' @param    dx1g	  numeric. Delta of quantile m - x(1-pg).
#' @param    dxg	  numeric. Delta of quantile x(pg) - m.
#' @param    dxp	  numeric. Delta of quantile x(p) - m.
#' @param    maxk	  numeric. Maximum allowed value for k (kappa).
#' 
#' @details
#' Internal function to estimate parameter k (kappa). Called by \code{estimkappa11}.
#' 
#' @name estimkappa6
estimkappa6 <- function(lg, lp, dx1p, dx1g, dxg, dxp, maxk = 10) { 
    h   <- lp/lg
	lgh <- lg * sqrt((-7 +3*h*h) /30) 
	rss <- 1.2 - 1.6 /(-1 + h*h)
	psi <- sqrt(dx1p *dxp /dx1g /dxg) / h
	k   <- if (psi <= 1) { maxk } else {
			lgh /sqrt(-1 + sqrt(1 + rss*(-1 +psi)))
			}
	k   <- if (k < maxk) { k } else { maxk }
	names(k) <- NULL
return(k)
}


#' @title Local Function Checkquantiles
#'
#' @description
#' Internal function to verify that quantiles are all different from
#' each other and correctly ordered.
#' 
#' @param    x	      vector of quantiles
#'  
#' 
#' @name checkquantiles
checkquantiles <- function(x) {
	mat <- outer(x, x, ">=")
	z   <- sum(mat[upper.tri(mat)]) == 0
return(z)
}


#' @title Parameter Estimation of Kiener Distributions From 5, 7 or 11 Quantiles
#'
#' @description
#' Parameter estimation of Kiener Distributions with just 5, 7 or 11 quantiles. 
#' 
#' @param    x5,x7,x11	   vector of 5, 7 or 11 quantiles. 
#' @param    p5,p7,p11	   vector of 5, 7 or 11 probabilities. 
#' @param    ord	       integer. Option for quantile selection and treatment.
#' @param    k	           numeric. parameter k used in model k3 and k4.
#' @param    maxk	       numeric. Maximum allowed value for k (kappa).
#' 
#' @details
#' These functions are called by \code{\link{estimkienerX}} to estimate 
#' the parameters of a Kiener distribution applied to a dataset \code{X}.
#' \code{p11, x11} are obtained with the functions \code{elevenprobs(X)} and \code{quantile(p11)}.
#' \code{p7, x7} are obtained with the functions \code{sevenprobs(X)} and \code{quantile(p7)}.  
#' \code{p5, x5} are obtained with the functions \code{fiveprobs(X)} and \code{quantile(p5)}.  
#' \code{k} is parameter k (kappa) of Kiener distributions.
#' 
#' The treatment of the 11 quantiles is controlled with the option \code{ord} 
#' which can take 12 integer values, \code{ord = 7} being the default: 
#' \itemize{
#'   \item{ 1: c(1-p1, 0.35, 0.50, 0.65, p1)}
#'   \item{ 2: c(1-p2, 0.35, 0.50, 0.65, p2)}
#'   \item{ 3: c(1-p1, 1-p2, 0.35, 0.50, 0.65, p2, p1)}
#'   \item{ 4: c(1-p1, 1-p2, 1-p3, 0.35, 0.50, 0.65, p3, p2, p1)}
#'   \item{ 5: c(1-p1, 0.25, 0.50, 0.75, p1)}
#'   \item{ 6: c(1-p2, 0.25, 0.50, 0.75, p2)}
#'   \item{ 7: c(1-p1, 1-p2, 0.25, 0.50, 0.75, p2, p1)}
#'   \item{ 8: c(1-p1, 1-p2, 1-p3, 0.25, 0.50, 0.75, p3, p2, p1)}
#'   \item{ 9: c(1-p1, 0.25, 0.35, 0.50, 0.65, 0.75, p1)}
#'   \item{10: c(1-p2, 0.25, 0.35, 0.50, 0.65, 0.75, p2)}
#'   \item{11: c(1-p1, 1-p2, 0.25, 0.35, 0.50, 0.65, 0.75, p2, p1)}
#'   \item{12: c(1-p1, 1-p2, 1-p3, 0.25, 0.35, 0.50, 0.65, 0.75, p3, p2, p1)}
#' }
#' 
#' The treatment of the 7 quantiles has only one solution:
#' \itemize{
#'   \item{ 1: c(1-p1, 1-p2, 0.25, 0.50, 0.75, p2, p1)}
#' }
#' 
#' The treatment of the 5 quantiles has only one solution:
#' \itemize{
#'   \item{ 1: c(1-p1, 0.25, 0.50, 0.75, p1)}
#' }
#' 
#' Parameter maxk controls the maximum allowed value for estimated parameter k. 
#' Reasonnable values are \code{maxk = 10, 15, 20}. Default is \code{maxk = 10} 
#' to be consistent with \code{\link{regkienerLX}}. 
#' 
#' @seealso      
#' \code{\link{elevenprobs}}, \code{\link{estimkienerX}}, 
#' \code{\link[stats]{quantile}},
#' \code{\link{qkiener2}}, \code{\link{qkiener3}}, \code{\link{qkiener4}}.
#' 
#' @examples     
#' 
#' require(timeSeries)
#' 
#' ## Choose j in 1:16. Choose ord in 1:12 (5 is default)
#' DS   <- getDSdata()
#' j    <- 5
#' ord  <- 5
#' p11  <- elevenprobs(DS[[j]])
#' x11  <- quantile(DS[[j]], probs = p11, na.rm = TRUE, names = TRUE, type = 6) 
#' estimkiener11(x11, p11, ord)
#' ## end
#' 
#' 
#' @export
#' @name estimkiener11
estimkiener11 <- function(x11, p11, ord = 7, maxk = 10) {
	if (length(x11) != 11) {stop("length(x11) is of wrong size. Must be 11.")}
	if (length(p11) != 11) {stop("length(p11) is of wrong size. Must be 11.")}
	if (!is.element(ord, 1:12)) {stop("ord must be in 1:12")}
	names(x11) <- NULL
	if (checkquantiles(x11)) {
		m  <- x11[6]
		k  <- estimkappa11(x11, p11, ord, maxk)
		d  <- estimdelta11(x11, p11, ord)
		d  <- if (abs(d) > (0.90/k)) {sign(d)*0.90/k} else {d}	
		g  <- estimgamma11(x11, p11, k, ord)
		e  <- kd2e(k, d)
		a  <- ke2a(k, e)
		w  <- ke2w(k, e)
		z  <- c(m, g, a, k, w, d, e) 
	} else {
		z  <- c(NA, NA, NA, NA, NA, NA, NA) 
	}
	names(z) <- c("m", "g", "a", "k", "w", "d", "e")
return(z)
}
#' 
#' 
#' @export
#' @rdname estimkiener11
estimkiener7 <- function(x7, p7, maxk = 10) {
	if (length(x7) != 7) {stop("length(x7) is of wrong size. Must be 7.")}
	if (length(p7) != 7) {stop("length(p7) is of wrong size. Must be 7.")}
	names(x7) <- NULL
	if (checkquantiles(x7)) {
		dx 	<- abs(x7 - x7[4])
		lp7	<- logit(p7)
		m   <- x7[4]
		k   <- ( estimkappa6(lp7[5], lp7[7], dx[1], dx[3], dx[5], dx[7], maxk)
				 +estimkappa6(lp7[5], lp7[6], dx[2], dx[3], dx[5], dx[6], maxk))/2
		d   <- log(dx[7]/dx[1]) /4/lp7[7] + log(dx[6]/dx[2]) /4/lp7[6] 
		d   <- if (abs(d) > (0.90/k)) {sign(d)*0.90/k} else {d}	
		g   <- sqrt(dx[3]*dx[5]) /2/k /sinh(lp7[5] /k)
		e   <- kd2e(k, d)
		a   <- ke2a(k, e)
		w   <- ke2w(k, e)
		z   <- c(m, g, a, k, w, d, e)
	} else {
		z   <- c(NA, NA, NA, NA, NA, NA, NA) 
	} 
	names(z) <- c("m", "g", "a", "k", "w", "d", "e")
return(z)
}
#' 
#' 
#' @export
#' @rdname estimkiener11
estimkiener5 <- function(x5, p5, maxk = 10) {
	if (length(x5) != 5) {stop("length(x5) is of wrong size. Must be 5.")}
	if (length(p5) != 5) {stop("length(p5) is of wrong size. Must be 5.")}
	names(x5) <- NULL
	if (checkquantiles(x5)) {
		dx  <- abs(x5 - x5[3])
		lp5 <- logit(p5)
		m   <- x5[3]
		k   <- estimkappa6(lp5[4], lp5[5], dx[1], dx[2], dx[4], dx[5], maxk)
		d   <- log(dx[5]/dx[1]) /2/lp5[5] 
		d   <- if (abs(d) > (0.90/k)) {sign(d)*0.90/k} else {d}	
		g   <- sqrt(dx[2]*dx[4]) /2/k /sinh(lp5[4] /k)
		e   <- kd2e(k, d)
		a   <- ke2a(k, e)
		w   <- ke2w(k, e)
		z   <- c(m, g, a, k, w, d, e)
	} else {
		z   <- c(NA, NA, NA, NA, NA, NA, NA) 
	}
	names(z) <- c("m", "g", "a", "k", "w", "d", "e")
return(z)
}
#' 
#' 
#' @rdname estimkiener11
estimdelta11 <- function(x11, p11, ord) { 
	if (length(x11) != 11) {stop("length(x11) is of wrong side. Must be 11.")}
	if (length(p11) != 11) {stop("length(p11) is of wrong side. Must be 11.")}
	if (!is.element(ord, 1:12)) {stop("ord must be an integer in 1:12")}
    lp11 <- logit(p11)
	dx   <- abs(x11 - x11[6])
	d1   <- log(dx[11]/dx[1]) /2/lp11[11]
	d2   <- log(dx[10]/dx[2]) /2/lp11[10]
	d3   <- log(dx[9] /dx[3]) /2/lp11[9]
	d12  <- (d1 + d2)/2
	d123 <- (d1 + d2 + d3)/3
	d    <- switch(as.character(ord), 
	       "1"=d1, "2"=d2, "3"=d12, "4"=d123, 
	       "5"=d1, "6"=d2, "7"=d12, "8"=d123, 
	       "9"=d1,"10"=d2,"11"=d12,"12"=d123)
	names(d) <- NULL
return(d)
}
#' 
#' 
#' @rdname estimkiener11
estimkappa11 <- function(x11, p11, ord, maxk = 10) { 
	if (length(x11) != 11) {stop("length(x11) is of wrong size. Must be 11.")}
	if (length(p11) != 11) {stop("length(p11) is of wrong size. Must be 11.")}
	if (!is.element(ord, 1:12)) {stop("ord must be in 1:12")}
	dx   <- abs(x11 - x11[6])
    lp11 <- logit(p11)
	l65  <- lp11[7]
	l75  <- lp11[8]
	l3   <- lp11[9]
	l2   <- lp11[10]
	l1   <- lp11[11]
	k    <- switch(as.character(ord),
	"1"=        estimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
	"2"=        estimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
	"3"= mean(c(estimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    estimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk))), 
	"4"= mean(c(estimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    estimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    estimkappa6(l65, l3, dx[3], dx[5], dx[7], dx[ 9], maxk))), 
	"5"=        estimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
	"6"=        estimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk), 
	"7"= mean(c(estimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    estimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk))), 
	"8"= mean(c(estimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    estimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk), 
			    estimkappa6(l75, l3, dx[3], dx[4], dx[8], dx[ 9], maxk))), 
	"9"= mean(c(estimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    estimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk))), 
   "10"= mean(c(estimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    estimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk))), 
   "11"= mean(c(estimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    estimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    estimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    estimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk))), 
   "12"= mean(c(estimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    estimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    estimkappa6(l65, l3, dx[3], dx[5], dx[7], dx[ 9], maxk), 
			    estimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    estimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk), 
			    estimkappa6(l75, l3, dx[3], dx[4], dx[8], dx[ 9], maxk)))) 
	names(k) <- NULL
return(k)
}
#' 
#' 
#' @rdname estimkiener11
estimgamma11 <- function(x11, p11, k, ord) { 
	dx   <- abs(x11 - x11[6])
    lp11 <- logit(p11)
	g75  <- sqrt(dx[4]*dx[8]) /2/k /sinh(lp11[8] /k)
	g65  <- sqrt(dx[5]*dx[7]) /2/k /sinh(lp11[7] /k)
	gmm  <- (g65 + g75)/2
	g    <- switch(as.character(ord),
	       "1"=g65,  "2"=g65,  "3"=g65,  "4"=g65, 
	       "5"=g75,  "6"=g75,  "7"=g75,  "8"=g75, 
	       "9"=gmm, "10"=gmm, "11"=gmm, "12"=gmm)
	names(g) <- NULL
return(g)
}




#' @title Parameter estimation for Kiener Distributions
#'
#' @description
#' A fast parameter estimation for Kiener distributions of type K2, K3 and K4. 
#' 
#' @param    X	      numeric. Vector of quantiles.
#' @param    ord	  integer. Option for quantile selection and treatment.
#' @param    type	  integer. Option passed to the \code{\link{quantile}} function.
#' @param    maxk	  numeric. Maximum allowed value for k (kappa).
#' @param    parnames boolean. Parameter vector with or without names.
#' @param    dgts     integer. The number of rounded digits. 
#' 
#' @details
#' Asymetric Kiener distributions of type K2, K3, K4 can be estimated with just 
#' 5 quantiles : the median x.50, two quantiles around the median, usually x.35 
#' and x.65 or better x.25 and x.75, two quantiles located just before the extremes 
#' of the dataset, for instance x.01 and x.99 if the dataset \code{X} has more 
#' than 100 points, x.0001 and x.9999 if the dataset \code{X} has more than 
#' 10.000 points.
#'  
#' For averaging purpose, calculation can be done from 5 to 11 probability levels 
#' \code{c(1-p1, 1-p2, 1-p3, 0.25, 0.35, 0.50, 0.65, 0.75, p3, p2, p1)}
#' extracted from \code{X} with the function \code{\link{elevenprobs}}.   
#' p1, p2 and p3 are the most extreme probabilities of the dataset \code{X} 
#' with values finishing either by \code{pp950} or \code{pp975} or \code{pp990}. 
#' The corresponding quantiles, extracted with the standard function 
#' \code{\link{quantile}} are then transfered to \code{\link{estimkiener11}}. 
#' \code{estimkienerX} and \code{estimkienerX5} are wrappers of these functions. 
#' 
#' The probability and quantile selection among these 11 values is controlled 
#' with the option \code{ord} which can take 12 integer values, 
#' \code{ord = 7} being the default: 
#' \itemize{
#'   \item{ 1: c(1-p1, 0.35, 0.50, 0.65, p1)}
#'   \item{ 2: c(1-p2, 0.35, 0.50, 0.65, p2)}
#'   \item{ 3: c(1-p1, 1-p2, 0.35, 0.50, 0.65, p2, p1)}
#'   \item{ 4: c(1-p1, 1-p2, 1-p3, 0.35, 0.50, 0.65, p3, p2, p1)}
#'   \item{ 5: c(1-p1, 0.25, 0.50, 0.75, p1)}
#'   \item{ 6: c(1-p2, 0.25, 0.50, 0.75, p2)}
#'   \item{ 7: c(1-p1, 1-p2, 0.25, 0.50, 0.75, p2, p1)}
#'   \item{ 8: c(1-p1, 1-p2, 1-p3, 0.25, 0.50, 0.75, p3, p2, p1)}
#'   \item{ 9: c(1-p1, 0.25, 0.35, 0.50, 0.65, 0.75, p1)}
#'   \item{10: c(1-p2, 0.25, 0.35, 0.50, 0.65, 0.75, p2)}
#'   \item{11: c(1-p1, 1-p2, 0.25, 0.35, 0.50, 0.65, 0.75, p2, p1)}
#'   \item{12: c(1-p1, 1-p2, 1-p3, 0.25, 0.35, 0.50, 0.65, 0.75, p3, p2, p1)}
#' }
#' 
#' The option \code{type} used by \code{\link{quantile}} can change significantly 
#' the results. Our experience is that \code{type = 6} is appropriate when 
#' \code{k > 1.9} and \code{type = 5} is appropriate when \code{k < 1.9}. 
#' Other types \code{type = 8} and \code{type = 9} can be considered as well. 
#' The other types should be ignored. 
#'  
#' Parameter maxk controls the maximum allowed value for estimated parameter k. 
#' Reasonnable values are \code{maxk = 10, 15, 20}. Default is \code{maxk = 10} 
#' to be consistent with \code{\link{regkienerLX}}. 
#' 
#' \code{estimkienerX5} is a simplified function with predifined values ord = 5, 
#' type = 6, rounded digits \code{digits = c(3, 3, 1, 1, 1, 3, 2) + dgts)} 
#' and parameters names \code{c("m", "g", "a", "k", "w", "d", "e")} that can  
#' be omitted if \code{parnames = FALSE}.
#'  
#' When \code{a, k, w < 1.4}, the estimates become less accurate and the function 
#' \code{\link{regkienerLX}} should rather be preferred.
#'  
#'  
#' @return  
#' A vector of estimated parameters code{c(m, g, a, k, w, d, e)} with  
#' rounded values \code{digits = c(3, 3, 1, 1, 1, 3, 2) + dgts)}.
#' 
#'
#' @references
#' P. Kiener, Fat tail analysis and package FatTailsR - Season 2, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' Download it from: 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#
#' @seealso     \code{\link{regkienerLX}}, \code{\link{fitkienerLX}}.
#' 
#' @examples     
#' 
#' require(timeSeries)
#' 
#' ## Choose j in 1:16
#' DS   <- getDSdata()
#' j    <- 5
#' elevenprobs(DS[[j]])
#' round(estimkienerX(DS[[j]]), 3)
#' round(t(sapply(DS, estimkienerX)), 3)
#' t(sapply(DS, estimkienerX5))
#' ## end
#' 
#' 
#' @export
#' @name estimkienerX
estimkienerX <- function(X, ord = 7, type = 6, maxk = 10, parnames = TRUE) { 
	X   <- sort(as.numeric(X))
	p11 <- elevenprobs(X)
    x11 <- quantile(X, probs = p11, na.rm = TRUE, names = FALSE, type = type) 
	z   <- estimkiener11(x11, p11, ord, maxk)
	if (!parnames) { names(z) <- NULL }
return(z)
}
#' 
#' 
#' @export
#' @rdname estimkienerX
estimkienerX7 <- function(X, dgts = 0, parnames = TRUE) { 
	X   <- sort(as.numeric(X))
	p7  <- sevenprobs(X)
	x7  <- quantile(X, probs = p7, na.rm = TRUE, names = FALSE, type = 6) 
	z   <- estimkiener7(x7, p7)
	z   <- round(z, digits = c(3, 3, 1, 1, 1, 3, 2) + dgts)		
	if (!parnames) { names(z) <- NULL }
return(z)
}
#' 
#' 
#' @export
#' @rdname estimkienerX
estimkienerX5 <- function(X, dgts = 0, parnames = TRUE) { 
	X   <- sort(as.numeric(X))
	p5  <- fiveprobs(X)
	x5  <- quantile(X, probs = p5, na.rm = TRUE, names = FALSE, type = 6) 
	z   <- estimkiener5(x5, p5)
	z   <- round(z, digits = c(3, 3, 1, 1, 1, 3, 2) + dgts)		
	if (!parnames) { names(z) <- NULL }
return(z)
}


