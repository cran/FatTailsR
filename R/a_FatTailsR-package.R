



#' @title Package FatTailsR 
#'
#' @description
#' Kiener distributions of type I, II, III and IV tailored to 
#' distributions which exhibit left and right fat tails like those that 
#' occur in financial markets. These distributions can be used to estimate 
#' with a high accuracy market risks and value-at-risk. Also include power 
#' hyperbolas and power hyperbolic functions. 
#' Download the pdf cited in the references to get an overview of  
#' the theoretical part and several examples on stocks and indices.
#' 
#' @details
#' With so many functions, this package could look fat. But it's not! 
#' It's rather agile and easy to use! It addresses two topics: 
#' Power hyperbolic functions and Kiener distributions with their moments. 
#' The various functions can be assigned to the following groups:
#' \enumerate{
#'   \item Miscellaneous functions for general purpose and power hyperbolas:
#'         \itemize{
#'         \item \code{\link{logit}}, invlogit, \code{\link{kashp}}, dkashp_dx, ashp.
#'         }
#'   \item Power hyperbolas and power hyperbolic functions: 
#'         \itemize{
#'         \item \code{\link{exphp}}, coshp, sinhp, tanhp, sechp, cosechp, 
#'               cotanhp.
#'         }
#'   \item Inverse power hyperbolas and inverse power hyperbolic functions:  
#'         \itemize{
#'         \item \code{\link{loghp}}, acoshp, asinhp, atanhp, asechp, 
#'               acosechp, acotanhp.
#'         }
#'   \item The logishp function and its variants which combines the logistic 
#'         function and the power hyperbolas. Also, the kogit and invkogit: 
#'         \itemize{
#'         \item d, p, q, r, dp, dq, l, dl, ql \code{\link{logishp}}.
#'         \item \code{\link{kogit}}, \code{\link{invkogit}} (added in v1.2-0).
#'         }
#'   \item The Kiener distributions of type I, II, III and IV:
#'         \itemize{
#'         \item d, p, q, r, dp, dq, l, dl, ql \code{\link{kiener1}},
#'         \item d, p, q, r, dp, dq, l, dl, ql \code{\link{kiener2}},
#'         \item d, p, q, r, dp, dq, l, dl, ql \code{\link{kiener3}},
#'         \item d, p, q, r, dp, dq, l, dl, ql \code{\link{kiener4}}.
#'         }
#'   \item Fat tail corrective function for Kiener distributions:
#'         \itemize{
#'         \item \code{\link{ckiener1}}, \code{\link{ckiener2}},
#'               \code{\link{ckiener3}}, \code{\link{ckiener4}}
#'               (added in v1.2-0).
#'         }
#'   \item Conversion functions between parameters pertaining to Kiener 
#'         distributions of type II, III and IV:
#'         \itemize{
#'         \item \code{\link{aw2k}}, aw2d, aw2e, kd2a, kd2w, ke2a, ke2w, 
#'                ke2d, kd2e, de2k.
#'         \item \code{\link{pk2pk}} (added in v1.2-0).
#'         }
#'   \item Regression functions to estimate the distribution parameters 
#'         of a given dataset with Kiener distributions and Laplace-Gauss
#'         normal distribution. \code{regkienerLX} and \code{fitkienerLX} 
#'         use the same regression function. \code{estimkienerX} uses a 
#'         different algorithm:
#'         \itemize{
#'         \item \code{\link{regkienerLX}}, \code{\link{laplacegaussnorm}}.
#'         \item \code{\link{fitkienerLX}} (added in v1.2-0).
#'         \item \code{\link{estimkienerX}} (added in v1.2-0).
#'         }
#'   \item Functions related to \code{estimkienerX} (added in v1.2-0):
#'         \itemize{
#'         \item \code{\link{elevenprobs}}, \code{\link{sevenprobs}}, \code{\link{fiveprobs}}
#'         \item \code{\link{estimkiener11}}, \code{\link{estimkiener7}}, \code{\link{estimkiener5}}.
#'         }
#'   \item Moments of the distribution estimated from the dataset or from 
#'         the regression parameters (added in v1.2-0):
#'         \itemize{
#'         \item \code{\link{kmoments}}, xmoments, kmoment, kcmoment,
#'               kskewness, kkurtosis, kekurtosis.
#'         }
#'   \item One function to aggregate external datasets and one dataset 
#'         in different formats: data.frame, timeSeries, xts, zoo (added in v1.2-0):
#'         \itemize{
#'         \item \code{\link{getDSdata}}.
#'         \item \code{\link{dfData}}, \code{\link{tData}}, \code{\link{xData}}, \code{\link{zData}}.
#'         }
#' }
#' For a quick start on Kiener distributions, jump to the functions 
#' \code{\link{regkienerLX}}, \code{\link{fitkienerLX}}, 
#' \code{\link{estimkienerX}} and run the examples. 
#' Then, download and read the document in pdf format cited in the references 
#' which gives a complete overview on all functions. 
#' Finally, explore the other examples. 
#' 
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014. Download it from: 
#' \url{http://www.inmodelia.com/exemples/2014-0627-Rmetrics-Kiener-en.pdf}
#'
#' P. Kiener, Fat tail analysis and package FatTailsR - Season 2, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' Download it from: 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#'
#' @keywords symbolmath distribution models
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
#' reg    <- regkienerLX(X)
#' lgn    <- laplacegaussnorm(X)
#' lleg   <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", 
#'            "logit(0.95)   = 2.9", "logit(0.50)   = 0", 
#'            "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' pleg   <- c( paste("m =",  reg$coefr4[1]), paste("g  =", reg$coefr4[2]), 
#'              paste("k  =", reg$coefr4[3]), paste("e  =", reg$coefr4[4]) )
#' 
#' ## Main plot
#' op     <- par(mfrow = c(1,1), mgp = c(1.5,0.8,0), mar = c(3,3,2,1))
#' plot(reg$dfrXP, main = nameX)
#' legend("top", legend = pleg, cex = 0.9, inset = 0.02 )
#' lines(reg$dfrEP, col = 2, lwd = 2)
#' points(reg$dfrQkPk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' lines(lgn$dfrXPn, col = 7, lwd = 2)
#' 
#' ## Plot F(X) > 0,97
#' front = c(0.06, 0.39, 0.50, 0.95)
#' par(fig = front, new = TRUE, mgp = c(1.5, 0.6, 0), las = 0)
#' plot( reg$dfrXP[which(reg$dfrXP$P > 0.97),] , pch = 1, xlab = "", ylab = "", main = "F(X) > 0,97" )
#' lines(reg$dfrEP[which(reg$dfrEP$P > 0.97),], type="l", col = 2, lwd = 3 )
#' lines(lgn$dfrXPn[which(lgn$dfrXPn$Pn > 0.97),], type = "l", col = 7, lwd= 2 )
#' points(reg$dfrQkPk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' points(lgn$dfrQnPn, pch = 3, col = 7, lwd = 2, cex = 1)
#' 
#' ## Plot F(X) < 0,03
#' front = c(0.58, 0.98, 0.06, 0.61)
#' par(fig = front, new = TRUE, mgp = c(0.5, 0.6, 0), las = 0 )
#' plot( reg$dfrXP[which(reg$dfrXP$P < 0.03),] , pch = 1, xlab = "", ylab = "", main = "F(X) < 0,03")
#' lines(reg$dfrEP[which(reg$dfrEP$P < 0.03),], type = "l", col = 2, lwd = 3 )
#' lines(lgn$dfrXPn[which(lgn$dfrXPn$Pn < 0.03),], type = "l", col= 7, lwd= 2 )
#' points(reg$dfrQkPk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' points(lgn$dfrQnPn, pch = 3, col = 7, lwd = 2, cex = 1)
#' 
#' ## Moments from the parameters (k) and from the Dataset (X)
#' round(cbind("k" = kmoments(reg$coefk, lengthx = nrow(reg$dfrXL)), "X" = xmoments(X)), 2)
#' attributes(reg)
#' ### End block
#' 
#' @import  timeSeries
#' @import  minpack.lm
#' @import  stats
#' @rdname  FatTailsR
#' @aliases FatTailsR
#' @name    FatTailsR-package
#' @docType package
NULL


