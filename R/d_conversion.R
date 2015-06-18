

#' @include c_logishp.R



#' @title Local Conversion Functions Between Kiener Distribution Parameters
#'
#' @description
#' Conversion functions between parameters \code{a}, \code{k}, \code{w}, 
#' \code{d}, \code{e} used in Kiener distributions of type II, III and IV.
#' 
#' @param        a a numeric value.
#' @param        k a numeric value.
#' @param        w a numeric value.
#' @param        d a numeric value.
#' @param        e a numeric value.
#' 
#' @details
#' \code{a} (alpha) is the left tail parameter, 
#' \code{w} (omega) is the right tail parameter, 
#' \code{d} (delta) is the distorsion parameter, 
#' \code{e} (epsilon) is the eccentricity parameter. 
#' 
#' \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} and 
#' describes a global tail parameter.
#' 
#' They are defined by:
#' \deqn{ aw2k(a, w) = 2 / (1/a + 1/w) = k }
#' \deqn{ aw2d(a, w) = (-1/a + 1/w) / 2 = d }
#' \deqn{ aw2e(a, w) = (a - w) / (a + w) = e }
#' \deqn{ kd2a(k, d) = 1 / ( 1/k - d) = a }
#' \deqn{ kd2w(k, d) = 1 / ( 1/k + d) = w }
#' \deqn{ ke2a(k, e) = k / (1 - e) = a }
#' \deqn{ ke2w(k, e) = k / (1 + e) = w }
#' \deqn{ ke2d(k, e) = e / k = d }
#' \deqn{ kd2e(k, d) = k * d = e }
#' \deqn{ de2k(k, e) = e / d = k }
#' \code{aw2k} is equivalent to the harmonic mean. 
#'
#' @seealso 
#' The asymmetric Kiener distributions of type II, III and IV:  
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}
#'
#' @examples
#' 
#' aw2k(4, 6); aw2d(4, 6); aw2e(4, 6)
#' outer(1:6, 1:6, aw2k)
#' 
#' @export
#' @name aw2k
           aw2k <- function(a, w) { 2 / (1/a + 1/w) }

#' @export
#' @rdname aw2k
           aw2d <- function(a, w) { (-1/a + 1/w) / 2 }

#' @export
#' @rdname aw2k
           aw2e <- function(a, w) { (a - w) / (a + w) }

#' @export
#' @rdname aw2k
           kd2a <- function(k, d) { 1 / ( 1/k - d) }

#' @export
#' @rdname aw2k
           kd2w <- function(k, d) { 1 / ( 1/k + d) }

#' @export
#' @rdname aw2k
           ke2a <- function(k, e) { k / (1 - e) }

#' @export
#' @rdname aw2k
           ke2w <- function(k, e) { k / (1 + e) }

#' @export
#' @rdname aw2k
           ke2d <- function(k, e) { e / k }

#' @export
#' @rdname aw2k
           kd2e <- function(k, d) { k * d }

#' @export
#' @rdname aw2k
           de2k <- function(d, e) { e / d }



		   

#' @title Global Conversion Function Between Kiener Distribution Parameters
#'
#' @description
#' A conversion function between the parameter vectors of Kiener distributions 
#' of type I: \code{c(m, g, k)}, type II: \code{c(m, g, a, w)},
#' type III: \code{c(m, g, k, d)} and type IV: \code{c(m, g, k, e)} to and from
#' \code{coefk = c(m, g, a, k, w, d, e)} extracted from \code{\link{regkienerLX}}. 
#'  
#' @param     coefk    vectors of numeric of length 3, 4 or 7.
#' @param     model    character. Either "k1", "k2", "k3", "k4", "k7".
#' @param     to       character. Either "k1", "k2", "k3", "k4", "k7".
#' @param     roundk   integer. The rounding applied to the output.	
#' 
#' @details
#' Kiener distributions use the following parameters, some of them being redundant. 
#' See also \code{\link{aw2k}} for the formulas and 
#' the conversion between parameters:
#' \itemize{
#'   \item{ \code{m} (mu) is the median of the distribution,. }
#'   \item{ \code{g} (gamma) is the scale parameter. }
#'   \item{ \code{a} (alpha) is the left tail parameter. } 
#'   \item{ \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} 
#'          and describes a global tail parameter. }
#'   \item{ \code{w} (omega) is the right tail parameter. } 
#'   \item{ \code{d} (delta) is the distorsion parameter. }
#'   \item{ \code{e} (epsilon) is the eccentricity parameter. }
#' }
#' 
#' \code{pk2pk()} performs the conversion between the various representation, from and to:
#' \itemize{
#'   \item{ type I: \code{\link{kiener1}} or code{"k1"}: c(m, g, k) }
#'   \item{ type II: \code{\link{kiener1}} or code{"k2"}: c(m, g, a, w) }
#'   \item{ type III: \code{\link{kiener1}} or code{"k3"}: c(m, g, k, d) }
#'   \item{ type IV: \code{\link{kiener1}} or code{"k4"}: c(m, g, k, e) }
#'   \item{ \code{coefk}, here code{"k7"}, extracted from object of class 
#'          \code{clregk} like \code{\link{regkienerLX}}: c(m, g, a, k, w, d, e) }
#' }
#' 
#' \code{coefk} can take any of the above form. When length(coefk) is 4, 
#' \code{model = "k2", "k3" or "k4"} is required to differentiate the three models.
#' When length(coefk) is 3 or 7, recognition is automatic and  
#' \code{model = "k1" or "k7"} is not required (and is ignored).
#' When a vector of type \code{coefk "k7"} is used as input, the vector is assumed
#' to be correct and there is no check of the consistency between the 
#' parameters \code{a, k, w, d} and \code{e}.
#' 
#' The output may be any of the above forms. Default is \code{"k7" = c(m, g, a, k, w, d, e)} 
#' which is \code{coefk} provided by the regression function \code{\link{regkienerLX}}
#' and is widely used in plotting and subsequent calculations.
#' 
#' An integer rounding parameter is provided trough \code{roundk}. Default is no rounding.
#'
#' @seealso 
#' The local conversion functions \code{\link{aw2k}}, 
#' Kiener distributions of type I, II, III and IV: \code{\link{kiener1}}, 
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}
#'
#' @examples
#'
#' ## Example 1
#' c2 <- c(1, 2, 3, 5)
#' pk2pk(c2, model = "k2", to = "k1") # loose the asymmetry.
#' pk2pk(c2, model = "k2", to = "k2")
#' pk2pk(c2, model = "k2", to = "k3")
#' pk2pk(c2, model = "k2", to = "k4")
#' pk2pk(c2, model = "k2", to = "k4")
#' (c7 <- pk2pk(c2, model = "k2", to = "k7", roundk = 3))
#' pk2pk(c7, model = "k7", to = "k2")
#' 
#' ## Example 2 ("k2" to "k7")
#' b08 <- c(-1, 1, 0.6, 0.8)
#' b10 <- c(-1, 1, 1.0, 1.2)
#' b15 <- c(-1, 1, 1.5, 1.7)
#' b20 <- c(-1, 1, 2.0, 2.2)
#' b25 <- c(-1, 1, 2.5, 2.7)
#' b30 <- c(-1, 1, 3.0, 3.2)
#' b35 <- c(-1, 1, 3.5, 3.7)
#' b40 <- c(-1, 1, 4.0, 4.2)
#' b45 <- c(-1, 1, 4.5, 4.7)
#' round(apply(cbind(b08, b10, b15, b20, b25, b30, b35, b40, b45), 2, pk2pk), 3)
#' 
#' 
#' @export
#' @name pk2pk	   
pk2pk <- function(coefk, model = "k2", to = "k7", roundk = 9) {

coeff    <- if (is.matrix(coefk)) { 
	switch( as.character(ncol(coefk)) , 
	"3"  = cbind(coefk[,1], coefk[,2], coefk[,3], coefk[,3], coefk[,3], rep(0, nrow(coefk)), rep(0, nrow(coefk))),
	"4"  = switch( model , 
			"k2"  = cbind(coefk[,1], coefk[,2], coefk[,3], aw2k(coefk[,3], coefk[,4]), 
					  coefk[,4], aw2d(coefk[,3], coefk[,4]), aw2e(coefk[,3], coefk[,4])) ,
			"k3"  = cbind(coefk[,1], coefk[,2], kd2a(coefk[,3], coefk[,4]), coefk[,3],  
					  kd2w(coefk[,3], coefk[,4]), coefk[,4], kd2e(coefk[,3], coefk[,4])) ,
			"k4"  = cbind(coefk[,1], coefk[,2], ke2a(coefk[,3], coefk[,4]), coefk[,3],  
					  ke2w(coefk[,3], coefk[,4]), ke2d(coefk[,3], coefk[,4]), coefk[,4]) , 
			stop("when ncol(coefk) = 4, model must be either k2, k3 or k4.")
			) ,
	"7"  = coefk ,
	stop("coefk is of wrong size. ncol(coefk) must be a matrix of 3, 4 or 7 columns.")
	)} else { 
	switch( as.character(length(coefk)) , 
	"3"  = c(coefk[1], coefk[2], coefk[3], coefk[3], coefk[3], 0, 0) , 
	"4"  = switch( model , 
			"k2"  = c(coefk[1], coefk[2], coefk[3], aw2k(coefk[3], coefk[4]), 
					  coefk[4], aw2d(coefk[3], coefk[4]), aw2e(coefk[3], coefk[4])) , 
			"k3"  = c(coefk[1], coefk[2], kd2a(coefk[3], coefk[4]), coefk[3],  
					  kd2w(coefk[3], coefk[4]), coefk[4], kd2e(coefk[3], coefk[4])) ,
			"k4"  = c(coefk[1], coefk[2], ke2a(coefk[3], coefk[4]), coefk[3],  
					  ke2w(coefk[3], coefk[4]), ke2d(coefk[3], coefk[4]), coefk[4]) , 
			stop("when length(coefk) = 4, model must be either k2, k3 or k4.")
			) ,
	"7"  = coefk ,
	stop("coefk is of wrong size. length(coefk) must be 3, 4 or 7.")
	)}
zcoefk <- if (is.matrix(coeff)) { 
		switch( to , 
			"k1"  = coeff[,c(1,2,4)], 
			"k2"  = coeff[,c(1,2,3,5)], 
			"k3"  = coeff[,c(1,2,4,6)],
			"k4"  = coeff[,c(1, 2, 4, 7)], 
			"k7"  = coeff, 
			stop("to must be either k1, k2, k3, k4 or k7")
			)} else {
		switch( to , 
			"k1"  = coeff[c(1,2,4)], 
			"k2"  = coeff[c(1,2,3,5)], 
			"k3"  = coeff[c(1,2,4,6)],
			"k4"  = coeff[c(1, 2, 4, 7)], 
			"k7"  = coeff, 
			stop("to must be either k1, k2, k3, k4 or k7")			
			)}
if (is.matrix(zcoefk)) { 
		colnames(zcoefk) <- switch( to , 
			"k1"  = c("m", "g", "k") , 
			"k2"  = c("m", "g", "a", "w") , 
			"k3"  = c("m", "g", "k", "d") ,
			"k4"  = c("m", "g", "k", "e") , 
			"k7"  = c("m", "g", "a", "k", "w", "d", "e") 
			)} else {
		names(zcoefk) <- switch( to , 
			"k1"  = c("m", "g", "k") , 
			"k2"  = c("m", "g", "a", "w") , 
			"k3"  = c("m", "g", "k", "d") ,
			"k4"  = c("m", "g", "k", "e") , 
			"k7"  = c("m", "g", "a", "k", "w", "d", "e") 
			)}
return(round(zcoefk, roundk))
}


