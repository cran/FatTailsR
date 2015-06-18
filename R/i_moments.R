

#' @include h_kiener4.R



#' @title Moments Associated To Kiener Distribution Parameters
#'
#' @description
#' The non-central moments, the central moments, the skewness, the kurtosis 
#' the excess of kurtosis and the cumulants associated to the parameter values 
#' of Kiener distributions of type I, II, III and IV (models k1, k2, k3, k4)
#' or estimated from the dataset. All-in-one vectors are also provided.
#' 
#' @param    n	    	integer. Order of the moment. 
#' @param    coefk    	vector. Parameters of the distribution k1, k2, k3 or k4 (and special "k7"). 
#' @param    model   	character. Model type ("k2", "k3" or "k4") if \code{coefk} is of length 4.
#' @param    x	    	numeric. Vector of quantiles.
#' @param    lengthx	integer. The length of the dataset used to calculate the parameters.
#' 
#' @details 
#' The non-central moments \code{mn: m1, m2, m3, m4,...}, 
#' the central moments \code{un: u1, u2, u3, u4,...} (mu in Greek)
#' and the cumulants \code{kn: k1, k2, k3, k4,...} 
#' (kappa in Greek; not to be confounded with models "k1", "k2", "k3", "k4") 
#' of order n exist for all \eqn{n < min(a, w) < k}. 
#' The standard deviation \code{sd} and the variance \code{u2} exist only
#' if \eqn{2 < min(a, w) < k}. 
#' The skewness \code{sk} exists only if \eqn{3 < min(a, w) < k}. 
#' The kurtosis \code{ku} and the excess kurtosis \code{ke} exist only 
#' if \eqn{4 < min(a, w) < k}. 
#' 
#' \code{coefk} may take five different forms :
#' \itemize{
#'   \item{c(m, g, k) of length 3 for distribution "k1".}
#'   \item{c(m, g, a, w) of length 4 for distribution "k2".}
#'   \item{c(m, g, k, d) of length 4 for distribution "k3".}
#'   \item{c(m, g, k, e) of length 4 for distribution "k4".}
#'   \item{c(m, g, a, k, w, d, e) of length 7 extracted from object of class 
#'         \code{clregk} like \code{regkienerLX} (typically \code{"reg$coefk"}.}
#' }
#' Forms of length 3 and 7 are automatically recognized and do not require 
#' \code{model = "k1"} or \code{"k7"} which is ignored. 
#' Forms of length 4 require an additionnal parameter describing the model used,   
#' either \code{model = "k2", "k3"} or \code{"k4"}. 
#' See \code{\link{pk2pk}} for more details on the parameter conversion function 
#' used within \code{kmoments}.
#' 
#' \code{kmoments} and \code{xmoments} provide all-in-one vectors. 
#' \code{kmoments} calls each of the above functions from order 1 to 4. 
#' \code{xmoments} is the traditional average sum of squares, cubic and power 4 functions 
#' of non-centered and centered values of x, from which NA values have been removed. 
#' Therefore, length of x does not count NA values.
#' 
#' @return  
#' Vectors \code{kmoments} and \code{xmoments} have the following structure:
#' 
#' \item{ku}{Kurtosis.}
#' \item{ke}{Excess kurtosis.}
#' \item{sk}{Skewness.}
#' \item{sd}{Standard deviation.}
#' \item{m1}{Mean.}
#' \item{m2}{Non-central moment of second order.}
#' \item{m3}{Non-central moment of third order.}
#' \item{m4}{Non-central moment of fourth order.}
#' \item{u1}{Central moment of first order. 0}
#' \item{u2}{Central moment of second order. Variance}
#' \item{u3}{Central moment of third order.}
#' \item{u4}{Central moment of fourth order.}
#' \item{k1}{Cumulant of first order. 0}
#' \item{k2}{Cumulant of second order.}
#' \item{k3}{Cumulant of third order.}
#' \item{k4}{Cumulant of fourth order.}
#' \item{lh}{Length of x, from which NA values were removed.}
#' \item{......}{.}
#'  
#'   
#' @seealso    
#' \code{\link{pk2pk}}, \code{\link{kiener1}}, \code{\link{kiener2}},  
#' \code{\link{kiener3}}, \code{\link{kiener4}}, \code{\link{regkienerLX}}.
#' 
#' @examples     
#' 
#' ## Example 1
#' kmoment(5, c(1, 2, 7))
#' kmoment(2, c(-1, 1, 6, 9), model = "k2")
#' kmoment(2, c(-1, 1, 7.2, -0.2/7.2), model = "k3")
#' kmoment(2, c(-1, 1, 7.2, -0.2), model = "k4")
#' kmoment(2, c(-1, 1, 6, 7.2, 9, -0.2/7.2, -0.2))
#' 
#' kcmoment(5, c(1, 2, 7))   # should be 0
#' kcmoment(3, c(-1, 1, 6, 9), model = "k2")
#' kstandev(c(-1, 1, 6, 9), model = "k2")
#' kvariance(c(-1, 1, 6, 9), model = "k2")
#' kskewness(c(-1, 1, 6, 9), model = "k2")
#' kkurtosis(c(-1, 1, 6, 9), model = "k2")
#' kekurtosis(c(-1, 1, 6, 9), model = "k2")
#' 
#' ## Example 2: "k2" and "k7" are preferred input formats for kmoments
#' (m08 <- kmoments(b08 <- c(-1, 1, 0.6, 0.8)))
#' m10 <- kmoments(b10 <- c(-1, 1, 1.0, 1.2))
#' m15 <- kmoments(b15 <- c(-1, 1, 1.5, 1.7))
#' m20 <- kmoments(b20 <- c(-1, 1, 2.0, 2.2))
#' m25 <- kmoments(b25 <- c(-1, 1, 2.5, 2.7))
#' m30 <- kmoments(b30 <- c(-1, 1, 3.0, 3.2))
#' m35 <- kmoments(b35 <- c(-1, 1, 3.5, 3.7))
#' m40 <- kmoments(b40 <- c(-1, 1, 4.0, 4.2))
#' (m45 <- kmoments(b45 <- c(-1, 1, 4.5, 4.7)))
#' (mat <- round(apply(cbind(b08, b10, b15, b20, b25, b30, b35, b40, b45), 2, pk2pk), 3))
#' round(apply(mat, 2, kmoments), 3)
#' round(cbind(m08, m10, m15, m20, m25, m30, m35, m40, m45), 3)
#' 
#' ## Example 3
#' ## Choose a number 1=DAX, 2=SMI, 3=CAC, 4=FTSE
#' ## Print reg$coefk and check that 3 < k and 4 < k
#' j    <- 3 
#' X    <- 100*diff(log(as.numeric(EuStockMarkets[,j])))
#' reg  <- regkienerLX(X)
#' reg$coefk 
#' for (n in 1:4) { print(round(kcmoment(n, reg$coefk), 3)) }
#' round(cbind("k" = kmoments(reg$coefk), "X" = xmoments(X)), 2)
#' 
#' 
#' 
#' @name kmoments
NULL

#' @export
#' @rdname kmoments
kmoment  <- function(n, coefk, model = "k2") {

if (!is.element(n, 1:99)) {stop("n must be an integer 0 < n < 99")}
coeff <- pk2pk(coefk, model = model)
m  <- coeff[1]
g  <- coeff[2]
a  <- coeff[3]
k  <- coeff[4]
w  <- coeff[5] 
if (is.na(min(a, w)) || n >= min(a, w) ) { momk <- NA } else {
	momk <- 0
	for (i in 0:n) {
		for (j in 0:i) { 
			momk <- (momk + choose(n, i) *choose(i, j) 
					 *m^(n-i) *g^(i) *k^(i) *(-1)^(j) 
					 *beta(1 +j/a -(i-j)/w, 1 -j/a +(i-j)/w ) )  
		} 
	} 
}
names(momk) <- paste0("m", n)
return(momk)
}

#' @export 
#' @rdname kmoments
kcmoment  <- function(n, coefk, model = "k2") {

if (!is.element(n, 1:99)) {stop("n must be an integer 0 < n < 99")}
coeff <- pk2pk(coefk, model = model)
m  <- coeff[1]
g  <- coeff[2]
a  <- coeff[3]
k  <- coeff[4]
w  <- coeff[5]
if (is.na(min(a, w)) || n >= min(a, w) ) { redcmomk <- NA } else {
	nu <- -beta(1 -1/a, 1 +1/a) + beta(1 -1/w, 1 +1/w)
	redcmomk  <- 0
	for (i in 0:n) {
		for (j in 0:i) { 
		redcmomk <- (redcmomk + choose(n, i) 
		             *choose(i, j) *(-nu)^(n-i) *(-1)^(j) 
					 *beta(1 +j/a -(i-j)/w, 1 -j/a +(i-j)/w) )
		} 
	}
}  
cmomk  <- redcmomk *g^n *k^n
names(cmomk) <- paste0("u", n)
return(cmomk)
}

#' @export 
#' @rdname kmoments
kstandev <- function(coefk, model = "k2") {
	stdev <- sqrt(kcmoment(2, coefk, model))
	names(stdev) <- "sd"
return(stdev)
}

#' @export 
#' @rdname kmoments
kvariance <- function(coefk, model = "k2") {
	variance <- kcmoment(2, coefk, model)
	names(variance) <- "u2"
return(variance)
}

#' @export 
#' @rdname kmoments
kskewness <- function(coefk, model = "k2") {
	skewk <- kcmoment(3, coefk, model) / kcmoment(2, coefk, model)^(1.5)
	names(skewk) <- "sk"
return(skewk)
}


#' @export 
#' @rdname kmoments
kkurtosis <- function(coefk, model = "k2") {
	kurtk <- kcmoment(4, coefk, model) / kcmoment(2, coefk, model)^(2)
	names(kurtk) <- "ku"
return(kurtk)
}

#' @export 
#' @rdname kmoments
kekurtosis <- function(coefk, model = "k2") {
	kurtk <- ( kcmoment(4, coefk, model) / kcmoment(2, coefk, model)^(2) ) - 3
	names(kurtk) <- "ke"
return(kurtk)
}

#' @export
#' @rdname kmoments
xmoments <- function(x) {
x  <- as.numeric(x[!is.na(x)])
y  <- x - mean(x)
z  <- c( 
	"ku" = mean(y^4) / sd(x)^4 , 
	"ke" =(mean(y^4) / sd(x)^4) - 3 , 
	"sk" = mean(y^3) / sd(x)^3 ,
	"sd" = sd(x) ,
	"m1" = mean(x) ,
	"m2" = mean(x^2) ,
	"m3" = mean(x^3) ,
	"m4" = mean(x^4) ,
	"u1" = mean(y) ,
	"u2" = mean(y^2) ,
	"u3" = mean(y^3) ,
	"u4" = mean(y^4) ,
	"k1" = mean(y) ,
	"k2" = mean(y^2) ,
	"k3" = mean(y^3) ,
	"k4" = mean(y^4) - 3*(mean(y^2))^2 ,
	"lh" = length(x) 
	)
return(z)
}

#' @export
#' @rdname kmoments
kmoments <- function(coefk, model = "k2", lengthx = NA) {

coeff  <- pk2pk(coefk, model = model)
z  <- c(
	"ku" = kkurtosis(coeff) , 
	"ke" = kekurtosis(coeff) , 
	"sk" = kskewness(coeff) ,
	"sd" = sqrt(kcmoment(2, coeff)) ,
	"m1" = kmoment(1, coeff) ,
	"m2" = kmoment(2, coeff) ,
	"m3" = kmoment(3, coeff) ,
	"m4" = kmoment(4, coeff) ,
	"u1" = kcmoment(1, coeff) ,
	"u2" = kcmoment(2, coeff) ,
	"u3" = kcmoment(3, coeff) ,
	"u4" = kcmoment(4, coeff) ,
	"k1" = kcmoment(1, coeff) ,
	"k2" = kcmoment(2, coeff) ,
	"k3" = kcmoment(3, coeff) ,
	"k4" = kcmoment(4, coeff) - 3*(kcmoment(2, coeff))^2 ,
	"lh" = lengthx )
names(z) <- c(
	"ku", 
	"ke", 
	"sk",
	"sd",
	"m1",
	"m2",
	"m3",
	"m4",
	"u1",
	"u2",
	"u3",
	"u4",
	"k1",
	"k2",
	"k3",
	"k4",
	"lh" )
return(z)
}


