## #' @include a_data.R


### 03/04/2016 Copie (mais pas transfert) de FatTailsRplot a FatTailsR 

#' @title Elevate
#'
#' @description
#' A transformation to turn negative prices into positive prices
#' and maintain at the same time the hierachy between all prices. 
#'
#' @param     X    The prices.
#' @param     e    numeric. The focal point of the hyperbola.
#' 
#' @details
#' Negative prices in financial markets, like interest rates in Europe, are a 
#' nightmare as the rough calculation of the returns generates non-sense values.
#' \code{elevate} uses an hyperbola and implements the following formula: 
#'      \deqn{ elevate(x, e) = \frac{x + \sqrt{x^2 + e^2}}{2} }{%
#'             elevate(x, e) = (x + sqrt(x*x + e*e)) / 2 }
#' 
#' There is currently no rule of thumb to calculate \code{e}. 
#' When \eqn{e = NULL}, there is no change and the output is identical to the input.
#' When \eqn{e = 0}, all negative values are turned to 0.
#' 
#' @examples    
#' require(graphics)
#' 
#' X <- (-50:100)/5
#' plot( X, elevate(X, e = 5), type = "l", ylim = c(0, 20) )
#' lines(X, elevate(X, e = 2),   col = 2)
#' lines(X, elevate(X, e = 1),   col = 3)
#' lines(X, elevate(X, e = 0.5), col = 4)
#' lines(X, elevate(X, e = 0),   col = 1)
#' 
#' @export
#' @name elevate
elevate <- function(X, e = NULL) { 
	hyper <- function(X, e) { (X + sqrt(X*X + e*e)) / 2 }
	dimv  <- function(AMV) { 
		if  (is(AMV, "array")) { dim(AMV) } else 
			{ if (is.null(dim(AMV))) { c(length(AMV), 1) } else { dim(AMV) } }
	}

if (is(X, "Date") || is(X, "POSIXt") || is(X, "POSIXct") || is(X, "chron") ) {
		stop("X is of a Date/POSIXt/POSIXct/chron format!") }

if (is.null(e)) { out <- X } 
else {
	if ( is.vector(X) ) { 
		if ( (length(e) != 1) ) { stop("e should be of length 1.") }
		out <- hyper(X, e) } 
	else {
		Y  <-	if (is(X[,1], "Date") || is(X[,1], "POSIXt") || is(X[,1], "POSIXct") || is(X[,1], "chron") ) 
					 { X[,-1, drop = FALSE] } 
				else { X }

		if ( (length(e) != 1) && (length(e) != dimv(Y)[2]) ) { 
				stop("e should be of length 1 or the number of columns of X or X[,-Date].") }

		if ( length(e) == 1 ) { e <- rep(e, dimv(Y)[2]) }

		Z <- hyper(Y, e)

		out  <- if (is(X[,1], "Date") || is(X[,1], "POSIXt") || is(X[,1], "POSIXct") || is(X[,1], "chron") ) {
					data.frame(Date = X[,1], Z, stringsAsFactors = FALSE) 
				} 
				else { Z }
	}
}
return(out)
}


#' @title Simple and Elaborated Prices to Returns
#'
#' @description
#' \code{fatreturns} is an elaborated function to compute prices to returns. 
#' It includes a pre-treatment for negative prices. 
#' It computes either log-returns (default) or percentage-returns. 
#' It handles properly NA values in the input vector, replacing them by 0
#' in the output vector. Doing so, it warrants that the sum of the log-returns 
#' (when selected) is equal to the difference of the log-prices. 
#' It works with vector, matrix, data.frame, timeSeries, xts, zoo, list, list of lists 
#' and even list of vector, data.frame, timeSeries, xts, zoo mixed together.
#' The returned object is of same dimension and same class than the input object 
#' with the first line filled with 0.
#' The results may be as per one, per cent (default), per thousand and per ten thousand. 
#' 
#' \code{logreturns} is an improved version of function \code{100*diff(log(x))} to handle 
#' vector, matrix, data.frame and list. It handles properly the first line and the NA values. 
#' It does not control time, rownames and colnames but may return them. 
#' 
#' @param    x       The prices (vector, data.frame, matrix, timeSeries, xts, zoo, list).
#' @param    log     boolean. log returns or percentage returns.
#' @param    per     character. Either "one", "cent, "thousand", "tenthousand" or 
#'                   "o", "c", "th", "te". Multiply the result by 1, 100, 1000, 10000.
#' @param    e       NULL or positive numeric. NULL is for no change \code{f(x)=x}. 
#'                   A positive numeric designates the focal point of the hyperbola 
#'                   to turn negative prices into positive prices, keeping the hierarchy:  
#'                   \code{f(x)=(x+sqrt(x*x+e*e))/2}. There is currently no rules of thumb 
#'                   for the optimal value of \code{e}. See \code{\link{elevate}}.
#' @param    dfrcol  integer. For data.frame only, designates the column that handles the time 
#'                   and must be processed separately. Use \code{dfrcol = 0} if all columns 
#'                   must be processed and there is no time (or turn the data.frame to a matrix).
#' @param    na.rm   boolean. Replace \code{x[t]=NA} with the previous non-NA value available 
#'                   in the price serie such that \code{(x[t-1], x[t]=x[t-1], x[t+1])} and 
#'                   calculate the returns accordingly. Force 0 in the first line of the returns 
#'                   if \code{x[1]=NA}.
#' 
#' @examples 
#' 
#' fatreturns(extractData())
#' logreturns(extractData())
#' 
#' @export
#' @name  fatreturns
fatreturns  <- function(x, log = TRUE, per = "cent", e = NULL, dfrcol = 1, na.rm = TRUE) {
	fatret  <- function(x, log, per, e, dfrcol, na.rm) {
		calret <- function(x, log, per, e, dfrcol, na.rm) { 
			retret <- function(x, log, per, e, na.rm) { 
				# replaceNA <- function(x) {
					# n   <- length(x)
					# idx <- (1:n)[!is.na(x)]
					# y   <- approx(idx, x[idx], 1:n, method = "constant", rule = 2, f = 0)$y
					# names(y) <- names(x)
				# return(y)
				# }
				##
				per    <- match.arg(per, c("cent", "one", "thousand", "tenthousand"))
				valper <- switch(per, 
					             cent = 100, one = 1, thousand = 1000, tenthousand = 10000)
				##
				x   <- if (na.rm) {replaceNA(x)} else {x}
				x   <- if (is.null(e))  {x} else {(x+sqrt(x*x+e*e))/2} 
				ret <- if (log)   {valper*diff(log(x))} else {valper*diff(x) / x[-length(x)]}
			return(ret)
			}
			##
			if (dimdim1(x) == 1) {ret <- retret(x, log, per, e, na.rm)}
			if (dimdim1(x) == 2) {ret <- apply(x, 2, calret, log, per, e, dfrcol, na.rm)}
			if (dimdim1(x)  < 0) {ret <- lapply(x, fatreturns, log, per, e, dfrcol, na.rm)}
		return(ret)
		}
		##
		z <- calret(x, log, per, e, dfrcol, na.rm)
		if (dimdim1(x) == 1) {
			z        <- c(x[1]*0, z) 
			z[1]     <- 0
			names(z) <- names(x)
			if (is(x, "zoo")) {zoo::index(z) <- zoo::index(x)}
		}
		if (dimdim1(x) == 2) {
			z           <- rbind(x[1,]*0, z)
			z[1,]       <- rep(0, ncol(z))
			rownames(z) <- rownames(x)
			colnames(z) <- colnames(x)
			if (is(x, "zoo")) {zoo::index(z) <- zoo::index(x)}
			if (is(x, "timeSeries")) {timeSeries::time(z) <- timeSeries::time(x)}
		}
	return(z)
	}
	## 
	if (is.data.frame(x)) {
		if (dfrcol != 0) {
			xtime <- x[,dfrcol] 
			xdata <- x[,-dfrcol]
			rdata <- fatret(xdata, log, per, e, dfrcol, na.rm)
			z     <- cbind(xtime, rdata)
			rownames(z) <- rownames(x)
			colnames(z) <- colnames(x)
		} else {
			z     <- fatret(x, log, per, e, dfrcol, na.rm)
			rownames(z) <- rownames(x)
			colnames(z) <- colnames(x)
		}
	} else {
		z     <- fatret(x, log, per, e, dfrcol, na.rm)
	}
return(z)
}

#' @export
#' @rdname fatreturns
logreturns <- function(x) { 
	if (is.null(x) || dimdim1(x) > 2) {
		stop("x is NULL or x is an array. 
		     dimdim1(x) must be 1 or 2, not 0 or 3.")
	}
	if (dimdim1(x) == 1) {
		n    <- length(x)
		idx  <- (1:n)[!is.na(x)]
		y    <- approx(idx, x[idx], 1:n, method = "constant", rule = 2, f = 0)$y
		z    <- c(0, 100*diff(log(y)))
		names(z) <- names(x) 
	}
	if (dimdim1(x) == 2) {z <- apply(x, 2, logreturns)}  
	if (dimdim1(x)  < 0) {z <- lapply(x, logreturns)} 
return(z)
}

#' @export
#' @rdname fatreturns
replaceNA <- function(x) { 
	if (is.null(x) || dimdim1(x) > 2) {
		stop("x is NULL or x is an array. 
		     dimdim1(x) must be 1 or 2, not 0 or 3.")
	}
	if (dimdim1(x) == 1) {
		n    <- length(x)
		idx  <- (1:n)[!is.na(x)]
		z    <- approx(idx, x[idx], 1:n, method = "constant", 
		               rule = 2, f = 0)$y
		names(z) <- names(x) 
	}
	if (dimdim1(x) == 2) {z <- apply(x, 2, replaceNA)}  
	if (dimdim1(x)  < 0) {z <- lapply(x, replaceNA)} 
return(z)
}



