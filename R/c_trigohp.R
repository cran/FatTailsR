## #' @include b_fatreturns.R



#' @title Kashp Function
#'
#' @description
#' \code{kashp}, which stands for kappa times arc-sine-hyperbola-power  
#' is the nonlinear transformation of x at the heart 
#' of power hyperbolas, power hyperbolic functions and symmetric Kiener 
#' distributions.
#' \code{dkashp_dx} is its derivative with respect to \code{x}. 
#' \code{ashp} is provided for convenience.
#'
#' @param x a numeric value, vector or matrix.
#' @param k a numeric value or vector, preferably strictly positive.
#'
#' @details
#' \code{ashp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ ashp(x, k) = asinh(x/k) }
#' \code{kashp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ kashp(x, k) = k * asinh(x/k) }
#' \code{dkashp_dx} function is defined for x in (-Inf, +Inf) by: 
#'   \deqn{ \frac{d}{dx}kashp(x, k) = \frac{1}{\sqrt{(x/k)^2 + 1}} 
#'                                  = \frac{1}{\cosh(ashp(x, k))} }{%
#'          dkashp_dx(x, k) = 1/sqrt(x*x/k/k + 1) = 1/cosh(ashp(x, k)) }
#'
#' If k is a vector, then the use of the function \code{\link[base]{outer}} 
#' is recommanded.
#'
#' The undesired case k=0 returns 0 for kashp and dkashp_dx.
#'
#' 
#' @examples
#' require(graphics)
#' 
#' ### FUNCTIONS kashp, dkashp_dx, ashp
#' xx <- (-3:3)*3
#' x  <- (-9:9) ; names(x) <- x
#' k  <- c(9999, 8, 5, 3, 2, 1) ; names(k) <- k
#' mat1 <- outer(x, k, kashp)    ; mat1
#' mat2 <- outer(x, k,dkashp_dx) ; mat2
#' mat3 <- outer(x, k,  ashp)    ; mat3
#' 
#' ### GRAPHICS
#' op <- par(mfcol = c(2,2), mar = c(3,3,2,1))
#' matplot(x, mat1, type="l", lwd=2, xaxt="n", yaxt="n", main="kashp") 
#' axis(1, at = xx) ; axis(2, at = xx, las = 1)
#' legend("topleft", title = expression(kappa), legend = colnames(mat1),
#'        lty = 1:6, col = 1:6, lwd = 2, inset = 0.02, cex = 0.7)
#' 
#' matplot(x, mat2, type="l", lwd=2, xaxt="n", main="dkashp_dx", las=1, ylim=c(0,1)) 
#' axis(1, at = xx)
#' legend("bottom", title = expression(kappa), legend = colnames(mat1),
#'        lty = 1:6, col = 1:6, lwd = 2, inset = 0.02, cex = 0.7)
#' 
#' matplot(x, mat3, type="l", lwd=2, xaxt="n", main="ashp", las=1) 
#' axis(1, at = xx)
#' legend("topleft", title = expression(kappa), legend = colnames(mat1),
#'        lty = 1:6, col = 1:6, lwd = 2, inset = 0.02, cex = 0.7)
#' par(op)
#' 
#' @export
#' @name kashp
   kashp <- function(x, k = 1) { 
                z <- k * asinh(x / k) 
                z[which(z == "NaN")] <- 0
                z
               } 
#' @export
#' @rdname kashp
 dkashp_dx <- function(x, k = 1) {
                z <- abs(k) / sqrt( x*x + k*k ) 
                z[which(z == "NaN")] <- 0
                z
               }  
#' @export
#' @name kashp
   ashp <- function(x, k = 1) { 
                z <- asinh(x / k) 
                z[which(z == "NaN")] <- 0
                z
               } 

