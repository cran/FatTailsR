

#' @include logishp.R



#' @title Conversion Functions
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
#' \deqn{ kd2e(k, d) = d * k = e }
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
           kd2e <- function(k, d) { d * k }

#' @export
#' @rdname aw2k
           de2k <- function(d, e) { e / d }


