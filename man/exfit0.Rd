% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a_data.R
\docType{data}
\name{exfit0}
\alias{exfit0}
\alias{exfit1}
\alias{exfit2}
\alias{exfit3}
\alias{exfit4}
\alias{exfit5}
\alias{exfit6}
\alias{exfit7}
\title{Parameter Subsets}
\format{
An object of class \code{character} of length 2.

An object of class \code{character} of length 7.

An object of class \code{character} of length 8.

An object of class \code{character} of length 6.

An object of class \code{character} of length 6.

An object of class \code{character} of length 6.

An object of class \code{character} of length 33.

An object of class \code{character} of length 35.
}
\usage{
exfit0

exfit1

exfit2

exfit3

exfit4

exfit5

exfit6

exfit7
}
\description{
Some vectors of parameter names to be used with parameter \code{exfitk} in 
functions regkienerLX(.., exfitk = ...) and \code{fitkienerX(.., exfitk = ...)}
or to subset the vector (or matrix) \code{fitk } obtained after regression 
\code{fitk <- regkienerLX(..)$fitk} or estimation \code{fitk <- fitkienerX(..)}. 
Visit \code{\link{fitkienerX}} for details on each parameter.

\code{exfit0 <- c("lh", "ret")}

\code{exfit1 <- c("m", "g", "a", "k", "w", "d", "e")}

\code{exfit2 <- c("m1", "sd", "sk", "ke", "m1x", "sdx", "skx", "kex")}

\code{exfit3 <- c("q.01", "q.05", "q.95", "q.99", "ltm.025", "rtm.975")}

\code{exfit4 <- c("VaR.01", "VaR.05", "VaR.95", "VaR.99", "ES.025", "ES.975")}

\code{exfit5 <- c("c.01", "c.05", "c.95", "c.99", "h.025", "h.975")}

\code{exfit6 <- c(exfit1, exfit2, exfit3, exfit4, exfit5)}

\code{exfit7 <- c(exfit0, exfit1, exfit2, exfit3, exfit4, exfit5)}
}
\examples{
    

require(minpack.lm)
require(timeSeries)

### Load the datasets and select one number j in 1:16
j      <- 5
DS     <- getDSdata()
(fitk  <- regkienerLX(DS[[j]])$fitk)
fitk[exfit3]
fitkienerX(DS[[j]], exfitk = exfit3)


}
\keyword{datasets}
