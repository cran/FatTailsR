

#' @include m_laplaceroll.R



#' @title Eleven, Seven, Five Probabilities
#'
#' @description
#' Extract from a dataset \code{X} a vector of 11, 7 or 5 probabilities: 
#' \itemize{
#'   \item{ \code{c(p1, p2, p3, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p3, 1-p2, 1-p1)} }
#'   \item{ \code{c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)} }
#'   \item{ \code{c(p1, 0.25, 0.50, 0.75, 1-p1)} }
#' }
#' where p1, p2 and p3 are the most extreme probabilities with values finishing 
#' by \code{..01}, \code{..025} or \code{..05} that can be extracted from the 
#' dataset \code{X}. Parameters names are displayed if \code{parnames = TRUE}.
#' 
#' From version 1.8-0, p1 and 1-p1 can be associated to the i-th and (N-i)-th element.
#' 
#' @param    X	       numeric. Vector of quantiles.
#' @param    parnames  boolean. Output parameter vector with or without names.
#' @param    i	       integer. The i-th and (N-i)-th elements for which the 
#'                     probabilities p1 and 1-p1 are calculated. If (i == 0), the 
#'                     method used before version 1.8-0 : the extreme finishing 
#'                     by \code{..01}, \code{..025} or \code{..05}.
#' 
#' @seealso  \code{\link{fitkienerX}}, \code{\link{estimkiener11}}.
#' 
#' @examples
#' 
#' require(timeSeries)
#' 
#' ## DS
#' DS  <- getDSdata()
#' for (j in 1:16) { print(round(elevenprobs(DS[[j]]), 6)) }
#' z   <- cbind(t(sapply(DS, elevenprobs)), sapply(DS, length))
#' colnames(z) <- c("p1","p2","p3","p.25","p.35","p.50","p.65","p.75","1-p3","1-p2","1-p1","length")
#' z
#' 
#' ## Choose j in 1:16
#' j   <- 1
#' X   <- sort(DS[[j]])
#' leX <- logit(eX <- elevenprobs(X))
#' lpX <- logit(ppoints(length(X), a = 0))
#' plot(X, lpX)
#' abline(h = leX, lty = 3)
#' mtext(eX, side = 4, at = leX, las = 1, line = -3.3)
#' 
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
	np11  <- character(length(p11))
	for (i in 1:length(np11)) { np11[i] <- format(p11[i], nsmall=2, scientific=FALSE) } 
	names(p11) <- if (parnames) { np11 } else { NULL }
return(p11)
}
#' @export
#' @rdname elevenprobs
sevenprobs <- function(X, parnames = FALSE) { 
	listp <- c(as.numeric(c(0.25, 0.5, 1) %o% 10^(-8:-1)), 0.15, 0.20)
	ltX   <- length(X)
	p1    <- listp[findInterval(1/ltX, listp) + 1]
	p2    <- listp[findInterval(1/ltX, listp) + 2]
	p7    <- c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)
	np7   <- character(length(p7))
	for (i in 1:length(np7)) { np7[i] <- format(p7[i], nsmall=2, scientific=FALSE) } 
	names(p7) <- if (parnames) { np7 } else { NULL }
return(p7)
}
#' @export
#' @rdname elevenprobs
fiveprobs <- function (X, i = 4, parnames = FALSE) {
    X  <- sort(as.numeric(X[is.finite(X)]))
    N  <- length(X)
    if (i == 0) {
        listp <- c(as.numeric(c(0.25, 0.5, 1) %o% 10^(-8:-1)), 
                   0.15, 0.20)
        p1 <- listp[findInterval(1/N, listp) + 1]
        p5 <- c(p1, 0.25, 0.5, 0.75, 1-p1)
    } else { 
        p5 <- c(i/(N+1), 0.25, 0.50, 0.75, (N+1-i)/(N+1))
    }
    if (parnames) {
        np5 <- character(length(p5))
        for (i in 1:length(np5)) {
            np5[i] <- format(p5[i], nsmall = 2, scientific = FALSE)
        } 
        names(p5) <- np5 
    }
return(p5)
}


#' @title Check Quantiles and Probabilities
#'
#' @description
#' Check that quantiles (or probabilities) are all 
#' different from each other and correctly ordered. 
#' If \code{proba = TRUE}, check that values are in range (0, 1). 
#' 
#' @param    x	        vector of quantiles.
#' @param    proba	    boolean. If TRUE, check range (0,1).
#' @param    acceptNA   boolean. If FALSE, NA value are not accepted.
#' @param    STOP       boolean. If an error is encountered, TRUE stops
#'                      the function and returns an error message. 
#'                      FALSE just returns FALSE.
#' 
#' @examples     
#' 
#' lst <- list(
#'   0.8,
#'   c(0.1, 0.5, 0.8),
#'   c(0.1, 0.5, 0.8, 0.2),
#'   c(2, 3, 1),
#'   c(2, 3),
#'   -0.01,
#'   NA,
#'   c(NA, NA),
#'   c(0.1, NA),
#'   c(0.1, NA, 0.5, 0.8),
#'   c(0.1, NA, 0.8, NA, 0.5),
#'   c(12, NA)
#' )
#' 
#' ## Evaluate
#' for (i in seq_along(lst)) { 
#'   cat(i, lst[[i]], " : ",
#'       checkquantiles(lst[[i]], proba = FALSE, STOP = FALSE), 
#'       checkquantiles(lst[[i]], proba = TRUE, STOP = FALSE), 
#'       checkquantiles(lst[[i]], proba = FALSE, acceptNA = TRUE, STOP = FALSE), 
#'       checkquantiles(lst[[i]], proba = TRUE,  acceptNA = TRUE, STOP = FALSE), 
#' 	     "\n") 
#' }
#' 
#' sapply(lst, checkquantiles, proba = TRUE, acceptNA = TRUE, STOP = FALSE)
#' 
#' ## Not run: 
#' checkquantiles(matrix((1:12)/16, ncol=3), proba = TRUE, STOP = FALSE)
#' ## End(Not run)
#' 
#' @export
#' @name checkquantiles
checkquantiles <- function(x, proba = FALSE, acceptNA = FALSE, STOP = TRUE) {
	if (acceptNA) { x  <- x[!is.na(x)] } 
	n   <- length(x)
	nt1 <- (n == 0)
	nt2 <- anyNA(x)
	nt3 <- (dimdim1(x) != 1)
	nt4 <- (!is(x, "numeric"))
	nt5 <- (proba && any(x < 0))
	nt6 <- (proba && any(x > 1))
	nt7 <- any(nt1, nt2, nt3, nt4, nt5, nt6)
	nt1;nt2;nt3;nt4;nt5;nt6;nt7
	z   <- if (nt7) {
			FALSE
		 } else {
			if (n == 1) { 
				TRUE 
			} else {
				(sum(x[2:n] <= x[1:(n-1)]) == 0) 
			}
		 }
	if (STOP && !z) { 
		stop("Error in checkquantiles(): x is not correct (dim, NA, proba, etc...) 
		      or not sorted. Please check.") 
	}
return(z)
}


#' @title Estimation Functions with 5, 7 or 11 Quantiles
#'
#' @description
#' Several functions to estimate the parameters of asymmetric Kiener distributions 
#' with just 5, 7 or 11 quantiles. 
#' 
#' @param    x5,x7,x11	   vector of 5, 7 or 11 quantiles. 
#' @param    p5,p7,p11	   vector of 5, 7 or 11 probabilities. 
#' @param    ord	       integer. Option for probability selection and treatment.
#' @param    maxk          numeric. Maximum value for k (kappa).
#' @param    maxe	       numeric. Maximum value for abs(e) (epsilon).
#'                         Maximum is \code{maxe = 1}.
#' 
#' @details
#' These functions, called by \code{paramkienerX5}, \code{paramkienerX7}, 
#' \code{\link{paramkienerX}}, use 5, 7 or 11 probabilites and quantiles 
#' to estimate the parameters of Kiener distributions.   
#' 
#' \code{p5, x5} are obtained with functions \code{fiveprobs(X)} and \code{quantile(p5)}.  
#' 
#' \code{p7, x7} are obtained with functions \code{sevenprobs(X)} and \code{quantile(p7)}.  
#' 
#' \code{p11, x11} are obtained with functions \code{elevenprobs(X)} and \code{quantile(p11)}.
#' 
#' The extraction of the 11 probabilities is controlled with the option \code{ord} 
#' which can take 12 integer values, \code{ord = 7} being the default. 
#' Small dataset should consider \code{ord = 5} and 
#' large dataset can consider \code{ord = 12}: 
#' \enumerate{
#'   \item{ \code{c(p1, 0.35, 0.50, 0.65, 1-p1)}}
#'   \item{ \code{c(p2, 0.35, 0.50, 0.65, 1-p2)}}
#'   \item{ \code{c(p1, p2, 0.35, 0.50, 0.65, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, p2, p3, 0.35, 0.50, 0.65, 1-p3, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, 0.25, 0.50, 0.75, 1-p1)}}
#'   \item{ \code{c(p2, 0.25, 0.50, 0.75, 1-p2)}}
#'   \item{ \code{c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, p2, p3, 0.25, 0.50, 0.75, 1-p3, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p1)}}
#'   \item{ \code{c(p2, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p2)}}
#'   \item{ \code{c(p1, p2, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, p2, p3, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p3, 1-p2, 1-p1)}}
#' }
#' \code{p5 = fiveprobs(X)} corresponds to \code{c(p1, 0.25, 0.50, 0.75, 1-p1)}.
#' 
#' \code{p7 = sevenprobs(X)} corresponds to \code{c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)}.
#' 
#' The above probabilities are then transfered to the \code{\link{quantile}} function 
#' whose parameter \code{type} can change significantly the extracted quantiles. 
#' Our experience is that \code{type = 6} is appropriate when \code{k > 1.9} and 
#' \code{type = 5} is appropriate when \code{k < 1.9}. 
#' Other types \code{type = 8} and \code{type = 9} can be considered as well. 
#' The other types should be ignored. 
#' (Note: when \code{k < 1.5}, algorithm \code{algo = "reg"} returns better  
#' results).
#' 
#' Parameter maxk controls the maximum allowed value for estimated parameter k. 
#' Reasonnable values are \code{maxk = 10, 15, 20}. Default is \code{maxk = 10} 
#' to be consistent with \code{\link{regkienerLX}}. 
#' 
#' @seealso      
#' \code{\link{elevenprobs}}, \code{\link{paramkienerX}}, \code{\link[stats]{quantile}},
#' \code{\link{roundcoefk}}. 
#' 
#' @examples     
#' 
#' require(timeSeries)
#' 
#' ## Choose j in 1:16. Choose ord in 1:12 (7 is default)
#' j    <- 5
#' ord  <- 5
#' DS   <- getDSdata()
#' p11  <- elevenprobs(DS[[j]])
#' x11  <- quantile(DS[[j]], probs = p11, na.rm = TRUE, names = TRUE, type = 6) 
#' round(estimkiener11(x11, p11, ord), 3)
#' 
#' ## Compare the results obtained with the 12 different values of ord on stock j
#' compare <- function(ord, x11, p11) {estimkiener11(x11, p11, ord)}
#' coefk   <- t(sapply(1:12, compare, x11, p11)) 
#' rownames(coefk) <- 1:12
#' mcoefk  <- apply(coefk, 2, mean) # the mean of the 12 results above
#' roundcoefk(rbind(coefk, mcoefk), 13)
#' 
#' 
#' @export
#' @name estimkiener11
estimkiener11 <- function(x11, p11, ord = 7, maxk = 10) {
	if (length(x11) != 11) {stop("length(x11) is of wrong size. Must be 11.")}
	if (length(p11) != 11) {stop("length(p11) is of wrong size. Must be 11.")}
	if (!is.element(ord, 1:12)) {stop("ord must be in 1:12")}
	names(x11) <- NULL
	if (checkquantiles(x11, STOP = FALSE)) {
		m  <- x11[6]
		k  <- .hestimkappa11(x11, p11, ord, maxk)
		d  <- .hestimdelta11(x11, p11, ord)
		d  <- if (abs(d) > (0.90/k)) {sign(d)*0.90/k} else {d}	
		g  <- .hestimgamma11(x11, p11, k, ord)
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
#' @export
#' @rdname estimkiener11
estimkiener7 <- function(x7, p7, maxk = 10) {
	if (length(x7) != 7) {stop("length(x7) is of wrong size. Must be 7.")}
	if (length(p7) != 7) {stop("length(p7) is of wrong size. Must be 7.")}
	names(x7) <- NULL
	if (checkquantiles(x7, STOP = FALSE)) {
		dx 	<- abs(x7 - x7[4])
		lp7	<- logit(p7)
		m   <- x7[4]
		k   <- ( .hestimkappa6(lp7[5], lp7[7], dx[1], dx[3], dx[5], dx[7], maxk)
				+.hestimkappa6(lp7[5], lp7[6], dx[2], dx[3], dx[5], dx[6], maxk))/2
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
#' @export
#' @rdname estimkiener11
estimkiener5 <- function(x5, p5, maxk = 20, maxe = 0.90) {
    names(x5) <- names(p5) <- NULL
    if (length(x5) != 5) {stop("length(x5) is wrong. Must be 5.")}
    if (length(p5) != 5) {stop("length(p5) is wrong. Must be 5.")}
    if (checkquantiles(x5, STOP = FALSE)) {
        dx  <- abs(x5 - x5[3])
        lp  <- logit(p5)
        m   <- x5[3]
        d   <- log((x5[5]-x5[3])/(x5[3]-x5[1]))/2/lp[5]
        Q   <- (x5[5]-x5[1])/(x5[4]-x5[2])*cosh(d*lp[4])/cosh(d*lp[5])
        fOPT<- function (k, q, p, Q) {
            (Q - sinh(logit(p)/k)/sinh(logit(q)/k))^2 
        }
        k   <- if (Q <= sinh(lp[5]/maxk)/sinh(lp[4]/maxk)) {
                    maxk
               } else { 
                    optimize(fOPT, c(0.1, maxk), tol = 0.0001, 
                             q = p5[4], p = p5[5], Q = Q)$minimum
               }
        k   <- min(k, abs(1/d)*maxe)
        e   <- kd2e(k, d)
        a   <- kd2a(k, d)
        w   <- kd2w(k, d)
        g   <- (x5[4]-x5[2])/4/k/sinh(lp[4]/k)/cosh(d*lp[4])
        z   <- c(m, g, a, k, w, d, e)
    } else {
        z   <- c(NA, NA, NA, NA, NA, NA, NA) 
    }
    names(z) <- c("m", "g", "a", "k", "w", "d", "e")
return(z)
}

 
.hestimdelta11 <- function(x11, p11, ord) { 
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


.hestimkappa11 <- function(x11, p11, ord, maxk = 10) { 
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
	"1"=        .hestimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
	"2"=        .hestimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
	"3"= mean(c(.hestimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    .hestimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk))), 
	"4"= mean(c(.hestimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    .hestimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    .hestimkappa6(l65, l3, dx[3], dx[5], dx[7], dx[ 9], maxk))), 
	"5"=        .hestimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
	"6"=        .hestimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk), 
	"7"= mean(c(.hestimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    .hestimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk))), 
	"8"= mean(c(.hestimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    .hestimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk), 
			    .hestimkappa6(l75, l3, dx[3], dx[4], dx[8], dx[ 9], maxk))), 
	"9"= mean(c(.hestimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    .hestimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk))), 
   "10"= mean(c(.hestimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    .hestimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk))), 
   "11"= mean(c(.hestimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    .hestimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    .hestimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    .hestimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk))), 
   "12"= mean(c(.hestimkappa6(l65, l1, dx[1], dx[5], dx[7], dx[11], maxk), 
			    .hestimkappa6(l65, l2, dx[2], dx[5], dx[7], dx[10], maxk), 
			    .hestimkappa6(l65, l3, dx[3], dx[5], dx[7], dx[ 9], maxk), 
			    .hestimkappa6(l75, l1, dx[1], dx[4], dx[8], dx[11], maxk), 
			    .hestimkappa6(l75, l2, dx[2], dx[4], dx[8], dx[10], maxk), 
			    .hestimkappa6(l75, l3, dx[3], dx[4], dx[8], dx[ 9], maxk)))) 
	names(k) <- NULL
return(k)
}


.hestimgamma11 <- function(x11, p11, k, ord) { 
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


## @title Local Function .hestimkappa6
##
## @description
## Internal function to estimate parameter k (kappa). 
## 
## @param    lg	      logit of probability close to the median. 
## @param    lp	      logit of extreme probability. 
## @param    dx1p	  numeric. Delta of quantile m - x(1-p).
## @param    dx1g	  numeric. Delta of quantile m - x(1-pg).
## @param    dxg	  numeric. Delta of quantile x(pg) - m.
## @param    dxp	  numeric. Delta of quantile x(p) - m.
## @param    maxk	  numeric. Maximum allowed value for k (kappa).
## 
## @details
## Internal function to estimate parameter k (kappa). Called by \code{.hestimkappa11}.

.hestimkappa6 <- function(lg, lp, dx1p, dx1g, dxg, dxp, maxk = 10) { 
    h   <- lp/lg
	lgh <- lg * sqrt((-7 +3*h*h) /30) 
	rss <- 1.2 - 1.6 /(-1 + h*h)
	psi <- sqrt(dx1p *dxp /dx1g /dxg) / h
	k   <- if (psi <= 1) { maxk } else {
			lgh /sqrt(-1 + sqrt(1 + rss*(-1 +psi)))
			}
	# is.na(k) added on 27/07/2015 v 1.2-4, 1.2-5
	k   <- if (is.na(k)) { maxk } else { k }
	k   <- if (k < maxk) { k } else { maxk }
	names(k) <- NULL
return(k)
}


