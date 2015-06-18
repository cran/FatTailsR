

#' @include b_trigohp.R



#' @title Logit Function
#'
#' @description
#' The logit function, widely used in this package, is a wrapper of 
#' \code{\link{qlogis}}.
#' 
#' @param        p a numeric value or vector.
#' @details      \code{logit} function is defined for p in (0, 1) by: 
#'               \deqn{ logit(p) = log( p/(1-p) ) }
#' @examples     logit( c(ppoints(11, a = 1), NA, NaN) )
#' @seealso      The equivalent function \code{\link{qlogis}}.
#' @export
#' @name logit
         logit <- function(p) { qlogis(p) } 



#' @title InvLogit Function
#'
#' @description
#' The invlogit function is the inverse of the logit function and a 
#' wrapper of \code{\link{plogis}}.
#' 
#' @param        x a numeric value or vector.
#' @details      \code{invlogit} function is defined for x in (-Inf, +Inf) by: 
#'               \deqn{ invlogit(x) = exp(x)/(1+exp(x)) = plogis(x) }
#' @examples     invlogit( c(-Inf, -10:10, +Inf, NA, NaN) )
#' @seealso      The equivalent function \code{\link{plogis}}.
#' @export
#' @name logit
         invlogit <- function(x) { plogis(x) } 



#' @title The Power Hyperbola Logistic Distribution
#'
#' @description
#' Density, distribution function, quantile function and random generation
#' for the power hyperbola logistic distribution.
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
#'
#' @details
#' \code{dlogishp} function (log is available) is defined for 
#'                 x in (-Inf, +Inf) by: 
#'       \deqn{ dlogishp(x, k) = 
#'              dkashp\_dx(x, k) * plogishp(x, k) * plogishp(-x, k) }
#' \code{invkogit=plogishp} functions are defined for q in (-Inf, +Inf) by: 
#'       \deqn{ invkogit(q, k) = plogishp(q, k) = 1/(1 + exp(- kashp(q, k))) }
#' \code{kogit=qlogishp} functions are defined for p in (0, 1) by: 
#'       \deqn{ kogit(p, k) = qlogishp(p, k) = 2 * k * sinh(logit(p) / k) }
#' \code{rlogishp} function generates \code{n} random values.
#' 
#' In addition to the classical formats, the prefixes dp, dq, l, dl, ql 
#' are also provided:
#' 
#' \code{dplogishp} function (log is available) is defined for p in (0, 1) by: 
#'      \deqn{ dplogishp(p, k = 1) = p * (1 - p) / 2 / cosh( logit(p)/k ) }
#' \code{dqlogishp} function (log is available) is defined for p in (0, 1) by: 
#'      \deqn{ dqlogishp(p, k = 1) = 2 / p / (1 - p) * cosh( logit(p)/k ) }
#' \code{llogishp} function is defined for x in (-Inf, +Inf) by: 
#'      \deqn{ llogishp(x, k) = kashp(x, k) }
#' \code{dllogishp} function is defined for lp = logit(p) in (-Inf, +Inf) by : 
#'      \deqn{ dllogishp(lp, k) = p * (1 - p) / 2 / cosh( lp/k ) }
#' \code{qllogishp} function is defined for lp = logit(p) in (-Inf, +Inf) by : 
#'      \deqn{ qllogishp(lp, k) = 2 * k * sinh(lp / k) }
#'
#' If k is a vector, then the use of the function \code{\link[base]{outer}} 
#' is recommanded.
#'
#' @seealso The Kiener distribution of type I \code{\link{kiener1}} which has 
#' location (\code{m}) and scale (\code{g}) parameters. 
#'
#' @examples
#' 
#' require(graphics)
#' 
#' ### Example 1
#' pp <- c(ppoints(11, a = 1), NA, NaN) ; pp
#' plogishp(-5:5, k = 4)
#' dlogishp(-5:5, k = 4)
#' qlogishp(pp, k = 4)
#' outer(-5:5, 1:6, plogishp)
#' outer(-5:5, 1:6, dlogishp)
#' outer(runif(20), 1:6, qlogishp)
#' 
#' ### Example 2
#' x     <- seq(-15, 15, length.out = 101)
#' k     <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(k) <- k ; k
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 4, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' op    <- par(mfrow = c(2,2), mgp = c(1.5,0.8,0), mar = c(3,3,2,1))
#' 
#' plot(x, plogis(x, scale = 2), type = "b", lwd = 2, ylim = c(0, 1),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", main = "plogishp(x, k)")
#' for (i in 1:length(k)) lines(x, plogishp(x, k = k[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(kappa), legend = c(k, "plogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(x, dlogis(x, scale = 2), type = "b", lwd = 2, xaxs = "i", 
#'      yaxs = "i", xlab = "", ylab = "", main = "dlogishp(x, k)")
#' for (i in 1:length(k)) lines(x, dlogishp(x, k = k[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' 
#' plot(x, x/2, type = "b", lwd = 2, ylim = c(-7.5, 7.5), xaxs = "i", 
#'      yaxs = "i", xlab = "", ylab = "", main = "logit(logishp(h, k))")
#' for (i in 1:length(k)) lines(x, llogishp(x, k = k[i]),  
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' 
#' plot(x, log(dlogis(x, scale = 2)), lwd = 2, type = "b", xaxs = "i", 
#'      yaxs = "i", xlab = "", ylab = "", main = "log(dlogishp(x, k))") 
#' for (i in 1:length(k)) lines(x, dlogishp(x, k = k[i], log = TRUE),  
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' ### End example 2
#' 
#' ### Example 3
#' p <- ppoints(199, a=0)
#' plot(p, qlogis(p, scale = 2), type = "o", lwd = 2, ylim = c(-15, 15),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qlogishp(p, k)")
#' for (i in 1:length(k)) lines(p, qlogishp(p, k = k[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(kappa), legend = c(k, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(p, 2/p/(1-p), type = "o", lwd = 2, xlim = c(0, 1), ylim = c(0, 100),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dqlogishp(p, k)")
#' for (i in 1:length(k)) lines(p, dqlogishp(p, k = k[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("top", title = expression(kappa), legend = c(k, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p*(1-p)/2, type = "o", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 0.14), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qlogishp, dplogishp(p, k)")
#' for (i in 1:length(k)) lines(qlogishp(p, k = k[i]), dplogishp(p, k = k[i]),
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(kappa), legend = c(k, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End example 3
#' 
#' @name logishp
NULL

#' @export
#' @rdname logishp
          dlogishp <- function(x, k = 1, log = FALSE) { 
               v   <- dkashp_dx(x, k) * plogishp(x, k) * plogishp(-x, k)
               if(log) return(log(v)) else return(v) 
               }

#' @export
#' @rdname logishp
          plogishp <- function(q, k = 1) { 1/(1 + exp(- kashp(q, k))) } 

#' @export
#' @rdname logishp
          invkogit <- function(q, k = 1) { 1/(1 + exp(- kashp(q, k))) } 
		  
#' @export
#' @rdname logishp
          qlogishp <- function(p, k = 1) { 2 * k * sinh(logit(p) / k) }

#' @export
#' @rdname logishp
          kogit    <- function(p, k = 1) { 2 * k * sinh(logit(p) / k) }
		  
#' @export
#' @rdname logishp
          rlogishp <- function(n, k = 1) { p <- runif(n) ; qlogishp(p, k) }

#' @export
#' @rdname logishp
         dplogishp <- function(p, k = 1, log = FALSE) { 
                          v  <- p * (1 - p) / 2 / cosh( logit(p)/k )
                          if(log) return(log(v)) else return(v)
                          }

#' @export
#' @rdname logishp
         dqlogishp <- function(p, k = 1, log = FALSE) { 
                          v  <- 2 / p / (1 - p) * cosh( logit(p)/k )
                          if(log) return(log(v)) else return(v)
                          }

#' @export
#' @rdname logishp
          llogishp <- function(x, k = 1) { kashp(x, k) }

#' @export
#' @rdname logishp
         dllogishp <- function(lp, k = 1, log = FALSE) {
                          p <- plogis(lp)       # = invlogit
                          v <- p * (1 - p) / 2 / cosh( lp/k )
                          if(log) return(log(v)) else return(v)
                          }

#' @export
#' @rdname logishp
         qllogishp <- function(lp, k = 1) { 2 * k * sinh(lp / k) }
#                        where lp = logit(p)


