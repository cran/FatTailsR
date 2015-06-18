

#' @include l_regression.R



#' @title A Regression Function for Kiener Distributions
#'
#' @description
#' A function to estimate the distribution parameters of a given dataset with 
#' Kiener distributions of type I, II, III and IV. It performs a nonlinear 
#' regression of the logit of the empirical probabilities logit(p) on 
#' quantiles X. 
#' 
#' @param    X       vector of quantiles. 
#' @param    model   the model used for the regression: "k1", "k2", "k3", "k4". 
#' @param    dgts    vector of length 11. Control the rounding of output parameters.
#' @param    maxk    numeric. The maximum value of tail parameter \code{k}. 
#' @param    mink    numeric. The minimum value of tail parameter \code{k}. 
#' @param    app     numeric. The parameter "\code{a}" in the function 
#'                            \code{ppoints}.
#' 
#' @details      
#' This function is designed to estimate the parameters of Kiener distributions
#' for a given dataset. It encapsulates the four distributions described in
#' this package. 
#' "k1" uses model \code{lqkiener1}, "k2" uses model \code{lqkiener2}, 
#' "k3" uses model \code{lqkiener3} and "k4" uses model \code{lqkiener4}. 
#' 
#' A typical input is a numeric vector that describes the returns of a stock. 
#' Conversion from a (possible) time series format to a sorted numeric vector 
#' is done automatically and without any check of the initial format. 
#' There is also no check of missing values, \code{Na}, \code{NaN}, 
#' \code{-Inf}, \code{+Inf}. 
#' Empirical probabilities of each point in the sorted dataset is calculated 
#' with the function \code{\link[stats]{ppoints}}. The parameter \code{app} 
#' corresponds to the parameter \code{a} in \code{ppoints} but has been  
#' limited to the range (0, 0.5). Default value is 0 as large datasets are 
#' very common in finance. 
#' 
#' A nonlinear regression is performed with \code{\link[minpack.lm]{nlsLM}} 
#' from the logit of the probabilities \code{logit(p)} over the quantiles X 
#' with one of the functions \code{lqkiener1234}. 
#' These functions have been selected as they
#' have an explicit form in the four types (this is unfortunately not the case 
#' for \code{dkiener234}) and return satisfactory results with ordinary least 
#' squares. The median is calculated before the regression and is injected 
#' as a mandatory value in the regression function. 
#'
#' Kiener distributions use the following parameters, some of them being redundant. 
#' See \code{\link{aw2k}} and \code{\link{pk2pk}} for the formulas and 
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
#' Where:
#' \itemize{
#'   \item{c(m, g, k) of length 3 for distribution "k1".}
#'   \item{c(m, g, a, w) of length 4 for distribution "k2".}
#'   \item{c(m, g, k, d) of length 4 for distribution "k3".}
#'   \item{c(m, g, k, e) of length 4 for distribution "k4".}
#'   \item{c(m, g, a, k, w, d, e) of length 7 extracted from object of class 
#'         \code{clregk} like \code{regkienerLX} (typically \code{"reg$coefk"}).}
#' }
#' 
#' Model \code{"k1"} return results with 1+2=3 parameters and describes a 
#' (assumed) symmetric distribution. Parameters \code{d} and \code{e} are set 
#' to 0. Models \code{"k2"}, \code{"k3"} and \code{"k4"} describe asymmetric 
#' distributions. They return results with 1+3=4 parameters.
#' Model "k2" has a very clear parameter definition but unfortunately 
#' parameters \code{a} and \code{w} are highly correlated. 
#' Model \code{"k3"} has the least correlated parameters but the meaning of 
#' the distorsion parameter \code{d}, usually of order 1e-3, is not simple. 
#' Multiplying it by 1000 might be a good choice but has not been done here.
#' 
#' Model \code{"k4"} exhibits a reasonable correlation between each parameter
#' and should be the preferred intermediate model between "k1" and "k2" models.
#' The eccentricity parameter \code{e} is well defined and easy to understand:
#' \eqn{e=(a-w)/(a+w)}, \eqn{a=k/(1-e)} and \eqn{w=k/(1+e)}. It is a relative 
#' measure that could be expressed (and understood) as a percentage. Here, the 
#' order 1e-2 has been preserved. 
#'
#' Tail parameter lower and upper values are controlled by \code{maxk} and 
#' \code{mink}. An upper value \eqn{maxk = 10} is appropriate for datasets
#' of low and medium size, less than 50.000 points. For larger datasets, the
#' upper limit can be extended up to \eqn{maxk = 20}. Such a limit returns 
#' results which are very closed to the logistic distribution, an alternate 
#' distribution which could be more appropriate. The lower limit \code{mink} 
#' is intended to avoid the value \eqn{k=0}. Remind 
#' that value \eqn{k < 2} describes distribution with no stable variance and 
#' \eqn{k < 1} describes distribution with no stable mean.
#' 
#' The output is an object in a flat format of class \code{clregk}. It can be 
#' listed with the function \code{\link{attributes}}. First are the data.frames with
#' the initial data and the estimated results. Second is the result of 
#' the regression \code{regk0} given by \code{\link[minpack.lm]{nlsLM}} from 
#' which a few information have been extracted and listed here. 
#' Third are the regression parameters (without the median) in plain format  
#' (no rounding), the variance-covariance matrix, the variance-covariance 
#' matrix times 1e+6 and the correlation matrix in a rounded format.
#' Note that \code{regk0}, \code{coefk0}, \code{coefk0tt}, \code{vcovk0}, 
#' \code{mcork0} have a polymorphic format and changing parameters that 
#' depend from the selected model: "k1", "k2", "k3", "k4". They should be  
#' used with care in subsequent calculations. 
#' Fourth are the distribution parameters tailored to every model "k1", "k2", 
#' "k3", "k4" plus estimated quantiles at levels: 
#' c(0.001, 0.005, 0.01, 0.05, 0.5, 0.95, 0.99, 0.995, 0.999). 
#' They are intended to subsequent calculations. 
#' Fifth are the same parameters presented in a more readable format thanks 
#' to the vector \code{dgts} which controls the rounding of the parameters in
#' the following order: \code{dgts = c("m", "g", 
#' "a", "k", "w", "d", "e", "vcovk0", "vcovk0m", "mcork0", "quantr")}.
#' Sixth are the estimated quantiles and probabilities of interest stored 
#' in a data.frame format.
#' 
#' 
#' @return  
#' \item{dfrXP}{data.frame. X = initial quantiles. P = empirical probabilities.}
#' \item{dfrXL}{data.frame. X = initial quantiles. L = logit of probabilities.}
#' \item{dfrXR}{data.frame. X = initial quantiles. R = residuals after regression.}
#' \item{dfrEP}{data.frame. E = estimated quantiles. P = probabilities.}
#' \item{dfrEL}{data.frame. E = estimated quantiles. L = logit of probabilities.}
#' \item{dfrED}{data.frame. E = estimated quantiles. 
#'                               D = estimated density (from probabilities).}
#' \item{regk0 }{object of class \code{nls} extracted from 
#'               the regression function \code{\link[minpack.lm]{nlsLM}}.}
#' \item{coefk0}{the regression parameters in plain format. 
#'               The median is out of the regression.} 
#' \item{vcovk0}{rounded variance-covariance matrix.} 
#' \item{vcovk0m}{rounded 1e+6 times variance-covariance matrix.} 
#' \item{mcork0}{rounded correlation matrix.} 
#' \item{coefk }{all parameters in plain format.} 
#' \item{coefk1}{parameters for model "k1".} 
#' \item{coefk2}{parameters for model "k2".} 
#' \item{coefk3}{parameters for model "k3".} 
#' \item{coefk4}{parameters for model "k4".} 
#' \item{quantk}{quantiles of interest.} 
#' \item{coefr }{all parameters in a rounded format.} 
#' \item{coefr1}{rounded parameters for model "k1".} 
#' \item{coefr2}{rounded parameters for model "k2".} 
#' \item{coefr3}{rounded parameters for model "k3".} 
#' \item{coefr4}{rounded parameters for model "k4".} 
#' \item{quantr}{quantiles of interest in a rounded format.} 
#' \item{dfrQkPk}{data.frame. Qk = Estimated quantiles of interest. 
#'                Pk = probabilities.} 
#' \item{dfrQkLk}{data.frame. Qk = Estimated quantiles of interest. 
#'                Lk = Logit of probabilities.} 
#' 
#' @seealso    \code{\link[minpack.lm]{nlsLM}}, \code{\link{laplacegaussnorm}}, 
#'     Kiener distributions of type I, II, III and IV: \code{\link{kiener1}}
#'     \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}.
#'     Other functions: \code{\link{fitkienerLX}}, \code{\link{estimkienerX}}.
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
#' 
#' ### and run this block
#' X      <- DS[[j]]
#' nameX  <- names(DS)[j]
#' reg    <- regkienerLX(X)
#' 
#' ## Plotting
#' lleg   <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", 
#'            "logit(0.95)   = 2.9", "logit(0.50)   = 0", 
#'            "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' pleg   <- c( paste("m =",  reg$coefr4[1]), paste("g  =", reg$coefr4[2]), 
#'              paste("k  =", reg$coefr4[3]), paste("e  =", reg$coefr4[4]) )
#' op     <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' plot(X, type="l", main = nameX)
#' plot(reg$dfrXL, main = nameX, yaxt = "n")
#' axis(2, las=1, at=c(-9.2, -6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9, 9.2))
#' abline(h = c(-4.6, 4.6), lty = 4)
#' abline(v = c(reg$quantk[5], reg$quantk[9]), lty = 4)
#' legend("topleft", legend = lleg, cex = 0.7, inset = 0.02, bg = "#FFFFFF")
#' lines(reg$dfrEL, col = 2, lwd = 2)
#' points(reg$dfrQkLk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' plot(reg$dfrXP, main = nameX)
#' legend("topleft", legend = pleg, cex = 0.9, inset = 0.02 )
#' lines(reg$dfrEP, col = 2, lwd = 2)
#' plot(density(X), main = nameX)
#' lines(reg$dfrED, col = 2, lwd = 2)
#' round(cbind("k" = kmoments(reg$coefk, lengthx = nrow(reg$dfrXL)), "X" = xmoments(X)), 2)
#' 
#' ## Attributes
#' attributes(reg)
#' head(reg$dfrXP)
#' head(reg$dfrXL)
#' head(reg$dfrXR)
#' head(reg$dfrEP)
#' head(reg$dfrEL)
#' head(reg$dfrED)
#' reg$regk0
#' reg$coefk0
#' reg$vcovk0
#' reg$vcovk0m
#' reg$mcork0
#' reg$coefk
#' reg$coefk1
#' reg$coefk2
#' reg$coefk3
#' reg$coefk4
#' reg$quantk
#' reg$coefr
#' reg$coefr1
#' reg$coefr2
#' reg$coefr3
#' reg$coefr4
#' reg$quantr
#' reg$dfrQkPk
#' reg$dfrQkLk
#' reg$fitk
#' ### End block
#' 
#' @export
#' @name regkienerLX
regkienerLX <- function(X, model = "k4", 
                        dgts = c(3, 3, 1, 1, 1, 3, 2, 4, 4, 2, 2),
                        maxk = 10, mink = 0.2, app = 0 ) {

if (app < 0 || app > 0.5) { 
	stop("app (the a of ppoints) must be between 0 and 0.5. 
          Recommended values: 0, 0.375, 0.5.")
	}
if (mink < 0.2 || mink > 2) { 
	stop("mink must be between 0.2 and 2. Value lesser than 1 is for strange 
          distributions!")
	}
if (maxk < 5 || maxk > 20) { 
	stop("maxk must be between 5 and 20. Can be increased with the sample size.")
	}
if (FALSE %in% (dgts %in% 0:6)) { 
	stop("each item of dgts must be in c(0, 1, 2, 3, 4, 5, 6)")
	}
if (length(dgts) != 11) { 
	stop("dgts must be of length 11")
	}
if ( ! model %in% list("k1", "k2", "k3", "k4") ) {
	stop("model must be either k1, k2, k3, k4. Default is k4 (m, g, k, e).")
	}

X        <- sort(as.numeric(X[!is.na(X)])) 
P        <- ppoints(length(X), a = app) 
L        <- logit(P) 
names(X) <- "X"
names(P) <- "P"
names(L) <- "L"
dfrXP    <- data.frame(X, P)
dfrXL    <- data.frame(X, L)
Xmed     <- median(X)
Xmean    <- mean(X)
Xsd      <- sd(X)
parini   <- estimkienerX5(X, parnames = FALSE)
if (anyNA(parini)) { 
	qqq    <- quantile(X, c(0.10, 0.50, 0.90), type = 6)
	gini   <- 0.27*Xsd
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

gmin   <- 0
amin   <- mink
kmin   <- mink
wmin   <- mink
dmin   <- - 1 / mink
emin   <- - (maxk - mink) / (maxk + mink)

gmax   <- Inf
amax   <- maxk
kmax   <- maxk
wmax   <- maxk
dmax   <- 1 / mink
emax   <- (maxk - mink) / (maxk + mink)

if (model == "k1") {
regk0  <- nlsLM( X ~ qlkiener1(L, Xmed, g, k), 
                 data = dfrXL, 
                 start = list(g = gini, k = kini), 
                 lower = c(gmin, kmin), 
                 upper = c(gmax, kmax) 
                ) 
coefk     <- c(m = Xmed,
               g = coef(regk0)[1],
               a = coef(regk0)[2],
               k = coef(regk0)[2],
               w = coef(regk0)[2],
               d = 0,
               e = 0
              ) 
names(coefk)   <- c("m", "g", "a", "k", "w", "d", "e")
}

if (model == "k2") {
regk0  <- nlsLM( X ~ qlkiener2(L, Xmed, g, a, w), 
                 data = dfrXL, 
                 start = list(g = gini, a = aini, w = wini), 
                 lower = c(gmin, amin, wmin), 
                 upper = c(gmax, amax, wmax) 
                )
coefk     <- c(m = Xmed, 
               g = coef(regk0)[1], 
               a = coef(regk0)[2], 
               k = aw2k(coef(regk0)[2], coef(regk0)[3]), 
               w = coef(regk0)[3],
               d = aw2d(coef(regk0)[2], coef(regk0)[3]), 
               e = aw2e(coef(regk0)[2], coef(regk0)[3])
              ) 
names(coefk) <- c("m", "g", "a", "k", "w", "d", "e")
}

if (model == "k3") {
regk0  <- nlsLM( X ~ qlkiener3(L, Xmed, g, k, d), 
                 data = dfrXL, 
                 start = list(g = gini, k = kini, d = dini), 
                 lower = c(gmin, kmin, dmin), 
                 upper = c(gmax, kmax, dmax) 
                )
coefk     <- c(m = Xmed,
               g = coef(regk0)[1],
               a = kd2a(coef(regk0)[2], coef(regk0)[3]),
               k = coef(regk0)[2],
               w = kd2w(coef(regk0)[2], coef(regk0)[3]),
               d = coef(regk0)[3],
               e = kd2e(coef(regk0)[2], coef(regk0)[3])
              ) 
names(coefk) <- c("m", "g", "a", "k", "w", "d", "e") 
}

if (model == "k4") {
regk0  <- nlsLM( X ~ qlkiener4(L, Xmed, g, k, e), 
                 data = dfrXL, 
                 start = list(g = gini, k = kini, e = eini), 
                 lower = c(gmin, kmin, emin), 
                 upper = c(gmax, kmax, emax) 
                )
coefk     <- c(m = Xmed,
               g = coef(regk0)[1],
               a = ke2a(coef(regk0)[2], coef(regk0)[3]),
               k = coef(regk0)[2],
               w = ke2w(coef(regk0)[2], coef(regk0)[3]),
               d = ke2d(coef(regk0)[2], coef(regk0)[3]),
               e = coef(regk0)[3]
              ) 
names(coefk) <- c("m", "g", "a", "k", "w", "d", "e") 
}

# Coefficients in plain format
coefk0    <- coef(regk0)
coefk1    <- c(coefk[1], coefk[2], coefk[4]) 
coefk2    <- c(coefk[1], coefk[2], coefk[3], coefk[5])
coefk3    <- c(coefk[1], coefk[2], coefk[4], coefk[6])
coefk4    <- c(coefk[1], coefk[2], coefk[4], coefk[7])
vcovk0    <- round(vcov(regk0),          dgts[8])
vcovk0m   <- round(vcov(regk0)*1e+6,     dgts[9])
mcork0    <- round(cov2cor(vcov(regk0)), dgts[10])
resik0    <- resid(regk0)
D         <- dlkiener2(lp = L, m = coefk2[1], g = coefk2[2], 
                       a = coefk2[3], w = coefk2[4], log = FALSE )
probak    <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.5,
               0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
quantk    <- qkiener2(p = probak, m = coefk2[1], g = coefk2[2], 
                                  a = coefk2[3], w = coefk2[4] )
names(quantk)  <- c("q.0001", "q.0005", "q.001", "q.005", "q.01", "q.05", 
            "q.50", "q.95", "q.99", "q.995", "q.999", "q.9995", "q.9999")

# data.frame
dfrXR     <- data.frame(X, R = resik0)
dfrEP     <- data.frame(E = fitted(regk0), P)
dfrEL     <- data.frame(E = fitted(regk0), L)
dfrED     <- data.frame(E = fitted(regk0), D = D)
dfrQkPk   <- data.frame(Qk = quantk, Pk = probak)
dfrQkLk   <- data.frame(Qk = quantk, Lk = logit(probak))

# Coefficients in rounded format
coefr     <- round(coefk, dgts[1:7])
coefr1    <- c(coefr[1], coefr[2], coefr[4]) 
coefr2    <- c(coefr[1], coefr[2], coefr[3], coefr[5])
coefr3    <- c(coefr[1], coefr[2], coefr[4], coefr[6])
coefr4    <- c(coefr[1], coefr[2], coefr[4], coefr[7])
quantr    <- round(quantk, dgts[11])

# Moments
momk		 <- kmoments(coefk, lengthx = length(X))[c("m1", "sd", "sk", "ke")]
momx		 <- xmoments(X)[c("m1", "sd", "sk", "ke", "lh")]
names(momx)  <- c("m1x", "sdx", "skx", "kex", "lh")

# Proba and quantiles
probak2      <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.10, 0.5,
                  0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
quantk2      <- qkiener2(p = probak2, m = coefk[1], g = coefk[2], 
                                      a = coefk[3], w = coefk[5] )
logisk         <- qlogis(p = probak2, location = coefk[1], scale = 2*coefk[2])
dlogisk        <- quantk2 - logisk
VaR		       <- Xmean - 2.326*Xsd
gaussk	       <- c( VaR, quantk2[5] -VaR )  # quantk2[5] = quantk2["q.01"]
names(quantk2) <- c("q.0001", "q.0005", "q.001", "q.005", "q.01", "q.05", "q.10",
            "q.50", "q.90", "q.95", "q.99", "q.995", "q.999", "q.9995", "q.9999")
names(logisk)  <- c("l.0001", "l.0005", "l.001", "l.005", "l.01", "l.05", "l.10",
            "l.50", "l.90", "l.95", "l.99", "l.995", "l.999", "l.9995", "l.9999")
names(dlogisk) <- c("d.0001", "d.0005", "d.001", "d.005", "d.01", "d.05", "d.10",
            "d.50", "d.90", "d.95", "d.99", "d.995", "d.999", "d.9995", "d.9999")
names(gaussk)  <- c("VaR", "dg.01")

## CQE - FTC - Tail-korrection
tailk          <- ckiener4(p = probak2, k = coefk[4], e = coefk[7])
names(tailk)   <- c("c.0001", "c.0005", "c.001", "c.005", "c.01", "c.05", "c.10",
            "c.50", "c.90", "c.95", "c.99", "c.995", "c.999", "c.9995", "c.9999")

## Rendement
ret            <- sum(X, na.rm = TRUE)
names(ret)     <- "ret"

## Final object fitk
fitk   <- c(ret, coefk, momk, momx, quantk2, logisk, dlogisk, gaussk, tailk)


# Final objet regk
regk      <- list()
regk$dfrXP     <- dfrXP
regk$dfrXL     <- dfrXL
regk$dfrXR     <- dfrXR
regk$dfrEP     <- dfrEP
regk$dfrEL     <- dfrEL
regk$dfrED     <- dfrED
regk$regk0     <- regk0
regk$coefk0    <- coefk0
regk$vcovk0    <- vcovk0
regk$vcovk0m   <- vcovk0m
regk$mcork0    <- mcork0
regk$coefk     <- coefk
regk$coefk1    <- coefk1
regk$coefk2    <- coefk2
regk$coefk3    <- coefk3
regk$coefk4    <- coefk4
regk$quantk    <- quantk
regk$coefr     <- coefr
regk$coefr1    <- coefr1
regk$coefr2    <- coefr2
regk$coefr3    <- coefr3
regk$coefr4    <- coefr4
regk$quantr    <- quantr
regk$dfrQkPk   <- dfrQkPk
regk$dfrQkLk   <- dfrQkLk
regk$fitk      <- fitk
class(regk)    <- "clregk"


return(regk)
}


