

#' @include m_regression2.R



#' @title Another Regression Function for Kiener Distributions
#'
#' @description
#' A function to estimate the distribution parameters of a given dataset with 
#' Kiener distributions. The output is a flat vector with the estimated parameters,
#' the moments (from the parameter and from the dataset), some quantiles. 
#' 
#' @param    X       vector of quantiles. 
#' @param    maxk    numeric. The maximum value of tail parameter \code{k}. 
#' @param    app     numeric. The parameter "\code{a}" in the function \code{ppoints}.
#' @param    dgts    vector of length 11. Control the rounding of output parameters.
#' 
#' @details      
#' This function is designed to estimate the parameters of Kiener distributions
#' for a given dataset. 
#' 
#' A typical input is a numeric vector that describes the returns of a stock. 
#' Conversion from a (possible) time series format to a sorted numeric vector 
#' is done automatically and without any check of the initial format. 
#' Missing values \code{NA} are removed. There is no check for \code{NaN}, 
#' \code{-Inf}, \code{+Inf}. 
#' 
#' Empirical probabilities of each point in the sorted dataset is calculated 
#' with the function \code{\link[stats]{ppoints}}. The parameter \code{app} 
#' corresponds to the parameter \code{a} in \code{ppoints} but has been  
#' limited to the range (0, 0.5). Default value is 0 as large datasets are 
#' very common in finance. 
#' 
#' A nonlinear regression is performed with \code{\link[minpack.lm]{nlsLM}} 
#' from the logit of the probabilities \code{logit(p)} over the quantiles X 
#' with the function \code{qlkiener4}. 
#'
#' Tail parameter max values are controlled by \code{maxk}.
#' An upper value \eqn{maxk = 10} is appropriate for datasets
#' of low and medium size, less than 50.000 points. For larger datasets, the
#' upper limit can be extended up to \eqn{maxk = 20}. Such a limit returns 
#' results which are very closed to the logistic distribution, an alternate 
#' distribution which could be more appropriate. Remember 
#' that value \eqn{k < 2} describes distribution with no stable variance and 
#' \eqn{k < 1} describes distribution with no stable mean.
#' 
#' The output is a vector with the regression parameters, some selected moments 
#' estimated from the regression parameters and from the dataset, some quantiles.  
#' See \code{\link{aw2k}} and \code{\link{pk2pk}} for the formulas and 
#' the conversion between the parameters:
#' \itemize{
#'   \item{ \code{m} (mu) is the median of the distribution. }
#'   \item{ \code{g} (gamma) is the scale parameter. }
#'   \item{ \code{a} (alpha) is the left tail parameter. } 
#'   \item{ \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} 
#'          and describes a global tail parameter. }
#'   \item{ \code{w} (omega) is the right tail parameter. } 
#'   \item{ \code{d} (delta) is the distorsion parameter. }
#'   \item{ \code{e} (epsilon) is the eccentricity parameter. }
#' }
#' The selected moments are the following:
#' \itemize{
#'   \item{ \code{m1} the first moment of the distribution, which is the median, thus equal to \code{m}. }
#'   \item{ \code{sd} the standard deviation as the square root of centered moment \code{u2} (mu2 in Greek). }
#'   \item{ \code{sk} the skewness. }
#'   \item{ \code{ke} the excess of kurtosis. }
#'   \item{ \code{x.m1} the mean of the data. }
#'   \item{ \code{x.sd} the square root of the squares of centered data. }
#'   \item{ \code{sk} the skewness estimated from the centered data at power 3. }
#'   \item{ \code{ke} the excess of kurtosis estimated from the centered data at power 4. }
#' }
#' 
#' @return  
#' A vector made of several parts:
#' \item{rdt}{The return between last and first values calculated through \code{sum(x)}. Thus, assume log-returns.} 
#' \item{coefk}{The regression parameters: \code{c(m, g, a, k, w, d, e)}.} 
#' \item{momk}{The moments related to the parameters: \code{c(m1, sd, sk, ke)}.} 
#' \item{momx}{The moments estimated from the dataset: \code{c(x.m1, x.sd, x.sk, x.ke)}.} 
#' \item{quantk}{Some quantiles:  
#'         \code{c("q.0001", "q.0005", "q.001", "q.005", "q.01", "q.05", "q.10",
#'         "q.50", "q.90", "q.95", "q.99", "q.995", "q.999", "q.9995", "q.9999")}.} 
#' \item{logisk}{The quantiles corresponding to the logisitic distribution \code{m + 2g logit(p)}}. 
#' \item{dlogisk}{The difference between the logistic quantiles and the true quantiles.} 
#' \item{.........}{.} 
#' \item{.........}{.} 
#' 
#' @seealso   \code{\link{regkienerLX}}, \code{\link{estimkienerX}}.
#'     
#' 
#' 
#' @examples     
#' 
#' require(graphics)
#' require(minpack.lm)
#' require(timeSeries)
#' 
#' ### Load the datasets and select one number (1-16)
#' DS     <- getDSdata()
#' j      <- 5
#' 
#' ### and run this block
#' X      <- DS[[j]]
#' nameX  <- names(DS)[j]
#' (fit   <- fitkienerLX(X))
#' data.frame(X = round(fit, 2))
#' 
#' ### Full result. Subsetting is recommended
#' t(round(sapply(DS, fitkienerLX), 2))
#' ### End block
#' 
#' @export
#' @name fitkienerLX
fitkienerLX <- function(X, maxk = 10, app = 0, dgts = 9) {

if (app < 0 || app > 0.5) { 
	stop("app (the a of ppoints) must be between 0 and 0.5. 
          Recommended values: 0, 0.375, 0.5.")
	}
if (maxk < 10 || maxk > 20) { 
	stop("maxk must be between 10 and 20. Can be increased with the sample size.")
	}
	
## Data
X        <- sort(as.numeric(X[!is.na(X)])) 
P        <- ppoints(length(X), a = app) 
L        <- logit(P) 
names(X) <- "X"
names(P) <- "P"
names(L) <- "L"
dfrXL    <- data.frame(X, L)

## Initialisation and bounds
Xmean  <- mean(X)
Xmed   <- median(X)
Xs     <- sd(X)

parini <- estimkienerX5(X, parnames = FALSE)
if (anyNA(parini)) { 
	qqq    <- quantile(X, c(0.10, 0.50, 0.90), type = 6)
	gini   <- 0.27*Xs
	dini   <- log(abs(qqq[3]-qqq[2])/abs(qqq[1]-qqq[2])) /4.394
	kini   <- 4
	eini   <- dini*kini
	aini   <- kini/(1-eini)
	wini   <- kini/(1+eini)
	} else {
	gini   <- parini[2]
	aini   <- parini[3]
	kini   <- parini[4]
	wini   <- parini[5]
	dini   <- parini[6]
	eini   <- parini[7]
	}

mink   <- 0.2
gmin   <- 0
kmin   <- mink
emin   <- - (maxk - mink) / (maxk + mink)
gmax   <- Inf
kmax   <- maxk
emax   <- (maxk - mink) / (maxk + mink)

## Regression k4
regk0  <- minpack.lm::nlsLM( X ~ qlkiener4(L, Xmed, g, k, e), 
                 data = dfrXL, 
                 start = list(g = gini, k = kini, e = eini), 
                 lower = c(gmin, kmin, emin), 
                 upper = c(gmax, kmax, emax) 
                )
coefk  <-    c(m = Xmed,
               g = coef(regk0)[1],
               a = ke2a(coef(regk0)[2], coef(regk0)[3]),
               k = coef(regk0)[2],
               w = ke2w(coef(regk0)[2], coef(regk0)[3]),
               d = ke2d(coef(regk0)[2], coef(regk0)[3]),
               e = coef(regk0)[3]
              ) 
names(coefk) <- c("m", "g", "a", "k", "w", "d", "e") 

## Moments
momk		 <- kmoments(coefk, lengthx = length(X))[c("m1", "sd", "sk", "ke")]
momx		 <- xmoments(X)[c("m1", "sd", "sk", "ke", "lh")]
names(momx)  <- c("m1x", "sdx", "skx", "kex", "lh")

## Proba and quantiles
probak    <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.10, 0.5,
               0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
quantk    <- qkiener2(p = probak, m = coefk[1], g = coefk[2], 
                                  a = coefk[3], w = coefk[5] )
logisk         <- qlogis(p = probak, location = coefk[1], scale = 2*coefk[2])
dlogisk        <- quantk - logisk
VaR		       <- Xmean - 2.326*Xs
gaussk	       <- c( VaR, quantk[5] -VaR )  # quantk[5] = quantk["q.01"]
names(quantk)  <- c("q.0001", "q.0005", "q.001", "q.005", "q.01", "q.05", "q.10",
            "q.50", "q.90", "q.95", "q.99", "q.995", "q.999", "q.9995", "q.9999")
names(logisk)  <- c("l.0001", "l.0005", "l.001", "l.005", "l.01", "l.05", "l.10",
            "l.50", "l.90", "l.95", "l.99", "l.995", "l.999", "l.9995", "l.9999")
names(dlogisk) <- c("d.0001", "d.0005", "d.001", "d.005", "d.01", "d.05", "d.10",
            "d.50", "d.90", "d.95", "d.99", "d.995", "d.999", "d.9995", "d.9999")
names(gaussk) <- c("VaR", "g.01")	

## CQE - FTC
tailk         <- ckiener4(p = probak, k = coefk[4], e = coefk[7])
names(tailk)  <- c("c.0001", "c.0005", "c.001", "c.005", "c.01", "c.05", "c.10",
           "c.50", "c.90", "c.95", "c.99", "c.995", "c.999", "c.9995", "c.9999")

## Rendement
ret            <- sum(X, na.rm = TRUE)
names(ret)     <- "ret"

## Final object fitk
fitk          <- round(c(ret, coefk, momk, momx, quantk, logisk, dlogisk, gaussk, tailk), digits = dgts)
return(fitk)
}


