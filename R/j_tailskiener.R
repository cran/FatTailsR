

#' @include i_moments.R



#' @title Fat Tail Corrective Function for Kiener Distributions
#'
#' @description
#' Quantile functions of Kiener distributions K1, K2, K3 and K4 can be expressed 
#' as a combination of logit function and fat tail corrective function. 
#' 
#' @param    p	  vector of probabilities. 
#' @param    k	  numeric. parameter k used in model k1, k3 and k4. 
#' @param    e	  numeric. parameter d used in model k4.
#' @param    a	  numeric. parameter a used in model k2.
#' @param    w	  numeric. parameter w used in model k2.
#' @param    d    numeric. parameter d used in model k3.
#' 
#' @details
#' Quantile functions of Kiener distributions K1, K2, K3 and K4 can be expressed 
#' as a combination of logit function and fat tail corrective function. 
#' 
#' @seealso      \code{\link{logit}}, \code{\link{qlogis}}, \code{\link{qkiener1}}, 
#'               \code{\link{qkiener2}}, \code{\link{qkiener3}}, \code{\link{qkiener4}}.
#' 
#' @export
#' @name ckiener
ckiener <- function(p, k, e) {
	l <- qlogis(p) 
	z <- k/l * sinh(l/k) * exp(l/k *e)
	z[which(z == "NaN")] <- 1
return(z)
}

#' @export
#' @rdname ckiener
ckiener1 <- function(p, k) {
	l <- qlogis(p) 
	z <- k/l * sinh(l/k)
	z[which(z == "NaN")] <- 1
return(z)
}

#' @export
#' @rdname ckiener
ckiener2 <- function(p, a, w) {
	l <- qlogis(p) 
	k <- aw2k(a, w)
	e <- aw2e(a, w)
	z <- k/l * sinh(l/k) * exp(l/k *e)
	z[which(z == "NaN")] <- 1
return(z)
}

#' @export
#' @rdname ckiener
ckiener3 <- function(p, k, d) {
	l <- qlogis(p) 
	z <- k/l * sinh(l/k) * exp(l * d)
	z[which(z == "NaN")] <- 1
return(z)
}

#' @export
#' @rdname ckiener
ckiener4 <- function(p, k, e) {
	l <- qlogis(p) 
	z <- k/l * sinh(l/k) * exp(l/k *e)
	z[which(z == "NaN")] <- 1
return(z)
}

