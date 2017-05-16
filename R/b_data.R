

#' @include a_FatTailsR-package.R



#' @title Get DS Dataset
#'
#' @description 
#' A function to extract the log-returns 
#' of 16 financial series and time series provided by the packages \code{datasets} 
#' (EuStockMarkets, sunspot.year) and \code{timeSeries} (USDCHF, MSFT, LPP2005REC).
#' The 16 datasets are converted to a list of numeric without any reference 
#' to the original dates. This list is usually called \code{DS}, hence the name.
#' 
#' @details
#' The dataset is usually created by the instruction \code{DS <- getDSdata()}.
#' Then, it is used with a call to DS[[j]] with j in 1:16. 
#' \enumerate{
#'   \item{ "USDCHF" (USDCHF, timeSeries) }
#'   \item{ "MSFT" (MSFT, timeSeries) }
#'   \item{ "DAX" (EuStockMarkets, datasets) }
#'   \item{ "SMI" (EuStockMarkets, datasets) }
#'   \item{ "CAC" (EuStockMarkets, datasets) }
#'   \item{ "FTSE" (EuStockMarkets, datasets) }
#'   \item{ "SBI" (LPP2005REC, timeSeries) }
#'   \item{ "SPI" (LPP2005REC, timeSeries) }
#'   \item{ "SII" (LPP2005REC, timeSeries) }
#'   \item{ "LMI" (LPP2005REC, timeSeries) }
#'   \item{ "MPI" (LPP2005REC, timeSeries) }
#'   \item{ "ALT" (LPP2005REC, timeSeries) }
#'   \item{ "LPP25" (LPP2005REC, timeSeries) }
#'   \item{ "LPP40" (LPP2005REC, timeSeries) }
#'   \item{ "LPP60" (LPP2005REC, timeSeries) }
#'   \item{ "sunspot" (sunspot.year, datasets) }
#' }
#' Note that \code{sunspot.year} is regularly updated with each new version of  
#' \code{R}. The generated dataset is \code{logreturn(sunspot.year + 1000)}.
#' 
#' @seealso    
#' \code{\link[datasets]{EuStockMarkets}}, \code{\link[datasets]{sunspot.year}}, 
#' \code{\link[timeSeries]{TimeSeriesData}}, \code{\link{regkienerLX}}, 
#' \code{\link{fitkienerX}}
#' 
#' @examples    
#' 
#' require(timeSeries) 
#' 
#' getDSdata
#' DS  <- getDSdata()
#' attributes(DS)
#' sapply(DS, length)
#' sapply(DS, head)
#' 
#' @export
#' @name getDSdata
getDSdata <- function() {
DSenv <- new.env(parent = baseenv())
utils::data("USDCHF", "MSFT", "LPP2005REC", package = "timeSeries", envir = DSenv)
utils::data("EuStockMarkets", "sunspot.year", package = "datasets", envir = DSenv)
prices2returns <- function(x) { 100*diff(log(as.numeric(x))) }
DS <- list(
	"USDCHF"	= prices2returns(DSenv$USDCHF),
	"MSFT"		= prices2returns(DSenv$MSFT[,4]),
	"DAX"		= prices2returns(DSenv$EuStockMarkets[,1]),
	"SMI"		= prices2returns(DSenv$EuStockMarkets[,2]),
	"CAC"		= prices2returns(DSenv$EuStockMarkets[,3]),
	"FTSE"		= prices2returns(DSenv$EuStockMarkets[,4]),
	"SBI"		= 100*as.numeric(DSenv$LPP2005REC[,1]),
	"SPI"		= 100*as.numeric(DSenv$LPP2005REC[,2]),
	"SII"		= 100*as.numeric(DSenv$LPP2005REC[,3]),
	"LMI"		= 100*as.numeric(DSenv$LPP2005REC[,4]),
	"MPI"		= 100*as.numeric(DSenv$LPP2005REC[,5]),
	"ALT"		= 100*as.numeric(DSenv$LPP2005REC[,6]),
	"LPP25"		= 100*as.numeric(DSenv$LPP2005REC[,7]),
	"LPP40"		= 100*as.numeric(DSenv$LPP2005REC[,8]),
	"LPP60"		= 100*as.numeric(DSenv$LPP2005REC[,9]),
	"sunspot"	= prices2returns(DSenv$sunspot.year+1000) )
return(DS)
}



#' @title Datasets dfData, mData, tData, xData, zData, extractData : dfData
#'
#' @description
#' A list of datasets in data.frame, matrix, timeSeries, xts and zoo formats. 
#' This is the data.frame format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name dfData
NULL

#' @title Datasets dfData, mData, tData, xData, zData, extractData : mData
#'
#' @description
#' A list of datasets in data.frame, matrix, timeSeries, xts and zoo formats. 
#' This is the matrix format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name mData
NULL

#' @title Datasets dfData, mData, tData, xData, zData, extractData : tData
#'
#' @description
#' A list of datasets in data.frame, matrix, timeSeries, xts and zoo formats. 
#' This is the timeSeries format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name tData
NULL

#' @title Datasets dfData, mData, tData, xData, zData, extractData : xData
#'
#' @description
#' A list of datasets in data.frame, matrix, timeSeries, xts and zoo formats. 
#' This is the xts format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name xData
NULL

#' @title Datasets dfData, mData, tData, xData, zData, extractData : zData
#'
#' @description
#' A list of datasets in data.frame, matrix, timeSeries, xts and zoo formats. 
#' This is the zoo format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name zData
NULL



#' @title Datasets dfData, mData, tData, xData, zData, extractData : extractData
#'
#' @description 
#' dfData, mData, tData, xData, zData are datasets made of lists of data.frame, matrix, 
#' timeSeries, xts and zoo components. They describe prices and returns of 10 financial series
#' used in the documents and demos presented at 8th and 9th R/Rmetrics conferences  
#' (2014, 2015). See the references. 
#' The last serie (CHF, interest rates in Switzerland) exhibits negative prices. 
#' All distributions of logreturns exhibit fat tails.
#' Function \code{extractData} converts subsets of mData, tData, xData, zData.
#' 
#' @param    pr	   	 character. Extract prices or returns: \code{c("p","r","prices","returns")}. 
#' @param    ft    	 character. Output format among \code{c("tss","xts","zoo","dfr","bfr","mat")}. 
#' @param    start   character. Start date.
#' @param    end	 character. End date.
#'
#' @details
#' 10 financial series presented in four different formats for convenient use: 
#' dfData is the complete dataset in data.frame format with dates as row.names 
#' and one or several columns. 
#' tData, xData and zData are lists of timeSeries, xts and zoo formats 
#' and display only one column per stock. 
#' \enumerate{
#'   \item{ "GOLD" from 1999-01-04 to 2013-12-31, dim 3694x1 (df, m, t, x, z). }
#'   \item{ "DEXIA" from 2008-10-27 to 20009-10-26, dim 255x1 (t, x, z), 255x5 (df). }
#'   \item{ "SOCGEN" from 1992-07-20 to 2013-12-31, dim 5445x1 (m, t, x, z), 5445x5 (df). }
#'   \item{ "VIVENDI" from 1992-07-20 to 2013-12-31, dim 5444x1 (m, t, x, z) 5444x5 (df). }
#'   \item{ "EURUSD" from 1999-01-03 to 2013-12-31, dim 3843x1 (m, t, x, z), 3843x4 (df). }
#'   \item{ "VIX" from 2004-01-02 to 2013-12-31, dim 2517x1 (m, t, x, z), 2517x4 (df). }
#'   \item{ "CAC40" from 1988-01-04 to 2013-12-31, dim 6574x1 (m, t, x, z), 6574x4 (df). }
#'   \item{ "DJIA" from 1896-05-26 to 2013-12-31, dim 32064x1 (m, df, t, x, z). }
#'   \item{ "SP500" from 1957-01-02 to 2013-12-31, dim 14350x1 (m, df, t, x, z). }
#'   \item{ "CHF" from 1995-01-02 to 2013-09-13, dim 4880x1 (t, x, z), 4880x8 (df). 
#'             Interst rates in Switzerland. Include negative prices at the end 
#'             of the dataset. Care is required to calculate the returns! 
#'             (Use \code{\link{fatreturns}} and \code{\link{elevate}}). 
#'             See the examples if you decide or need to remove it from the list.}
#' }
#' 
#' Function \code{extractData} extracts 8 financial series from matrix mData 
#' (DEXIA and CHF are not included) and converts them into a 2 dimensions object 
#' with any of the following class:
#' \itemize{
#'   \item{ "tss" is for timeSeries with a timeDate index. }
#'   \item{ "xts" is for xts with a POSIXct index. }
#'   \item{ "zoo" is for zoo with a Date index. }
#'   \item{ "dfr" is the usual R data.frame with the Date as index. }
#'   \item{ "bfr" is a data.frame with the Date in the first column. }
#'   \item{ "mat" is a matrix with the Date in the rownames. }
#' }
#' The \code{start} date must be posterior to "2007-01-01" (default) and the 
#' \code{end} date must be anterior to "2013-12-31" (default).
#'  
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014. Download it from: 
#' \url{http://www.inmodelia.com/exemples/2014-0627-Rmetrics-Kiener-en.pdf}
#'
#' P. Kiener, Fat tail analysis and package FatTailsR, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#' 
#' @seealso 
#' \code{\link{tData}}, \code{\link{xData}}, \code{\link{zData}}, \code{\link{dfData}},
#' \code{\link{getDSdata}}.
#' 
#' @examples    
#' 
#' library(zoo) 
#' library(xts) 
#' library(timeSeries) 
#' 
#' ### dfData, tData, xData, zData : prices only
#' attributes(dfData); attributes(tData); attributes(xData); attributes(zData) 
#' lapply(dfData, head, 3)
#' lapply( mData, head, 3)
#' lapply( tData, head, 3)
#' lapply( xData, head, 3)
#' lapply( zData, head, 3)
#' 
#' ### extractData : prices and logreturns
#' head(ptD <- extractData("p", "tss", "2009-01-01", "2012-12-31")) ; tail(ptD)
#' head(rtD <- extractData("r", "tss")) 
#' head(pxD <- extractData("p", "xts")) 
#' head(rxD <- extractData("r", "xts")) 
#' head(pzD <- extractData("p", "zoo")) 
#' head(rzD <- extractData("r", "zoo")) 
#' head(pbD <- extractData("p", "bfr")) 
#' head(rbD <- extractData("r", "bfr")) 
#' head(pmD <- extractData("p", "mat")) 
#' head(rmD <- extractData("r", "mat")) 
#' 
#' ### Remove item CHF (negative prices) from dfData, tData, xData, zData
#' Z <- dfData[names(dfData)[1:9]]; attributes(Z)
#' Z <- tData[names(tData)[1:9]]; attributes(Z)
#' Z <- xData[names(xData)[1:9]]; attributes(Z)
#' Z <- zData[names(zData)[1:9]]; attributes(Z)
#' 
#' @export
#' @name extractData
extractData <- function(pr = "p", ft = "tss",
						start = "2007-01-01", end = "2013-12-31") {
    if (is.na(charmatch(strtrim(pr, 1)[1], 
              c("p", "r")))) { 
              stop("pr must be p or r, prices or returns.") 
    }
    if (is.na(charmatch(ft, 
              c("tss","xts","zoo","dfr","bfr","mat")))) { 
              stop("ft must be tss, xts, zoo, dfr, bfr or mat.") 
    }
    start3   <- as.Date("2006-12-31")
    end3     <- as.Date("2014-01-01")
    start2   <- as.Date(start)
    end2     <- as.Date(end)
    startend <- as.numeric(c(start3, start2, end2, end3))
    if (checkquantiles(startend, STOP = FALSE)) {
        start1   <- start2
        end1     <- end2
    } else {
        start1   <- start3
        end1     <- end3
        warning("start or end dates were ignored as not in range 2007-01-01 - 2013-12-31")
    }
    mD   <- get("mData")
    DmD  <- as.Date(rownames(mD))
    win  <- DmD >= start1 & DmD <= end1
    mat  <- if (strtrim(pr, 1)[1] == "r") { 
                fatreturns(mD[win,])
            } else { 
                mD[win,]
            }
    z <- switch(ft,
        "tss" = timeSeries::timeSeries(mat, 
                                       as.Date(rownames(mat)), 
                                       units = colnames(mat)),
        "xts" = xts::xts(mat, as.POSIXct(rownames(mat))),
        "zoo" = zoo::zoo(mat, as.Date(rownames(mat))),
        "dfr" = as.data.frame(mat, stringsAsFactors = FALSE),
        "bfr" = data.frame("DATE" = rownames(mat), mat, 
                     row.names = NULL, stringsAsFactors = FALSE),
        "mat" = mat
        ) 
return(z)
}






#' @title Several Vectors of Probabilities
#' 
#' @description
#' Several vectors of probabilities used in FatTailsR. 
#' Remark: pprobs5 <- sort(c(pprobs2, pprobs3, pprobs4)).
#' 
#' pprobs0 <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
#' 
#' pprobs1 <- c(0.01, 0.05, 0.95, 0.99)
#' 
#' pprobs2 <- c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99)
#' 
#' pprobs3 <- c(0.001, 0.0025, 0.005, 0.995, 0.9975, 0.999)
#' 
#' pprobs4 <- c(0.0001, 0.00025, 0.0005, 0.9995, 0.99975, 0.9999)
#'
#' pprobs5 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#' 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999)
#'
#' pprobs6 <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.50, 
#' 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999) 
#'
#' pprobs7 <- c(0.01, 0.025, 0.05, 
#' 0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
#' 0.95, 0.975, 0.99) 
#' 
#' pprobs8 <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#' 0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
#' 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999) 
#' 
#' pprobs9 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#' 0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
#' 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999) 
#' 
#' 
#' @seealso 
#' The conversion function \code{\link{getnamesk}}
#'  
#' @export
#' @name pprobs0
pprobs0 <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99) 
#' @export
#' @rdname pprobs0
pprobs1 <- c(0.01, 0.05, 0.95, 0.99) 
#' @export
#' @rdname pprobs0
pprobs2 <- c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99) 
#' @export
#' @rdname pprobs0
pprobs3 <- c(0.001, 0.0025, 0.005, 0.995, 0.9975, 0.999)
#' @export
#' @rdname pprobs0
pprobs4 <- c(0.0001, 0.00025, 0.0005, 0.9995, 0.99975, 0.9999)
#' @export
#' @rdname pprobs0
pprobs5 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
			 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999) 
#' @export
#' @rdname pprobs0
pprobs6 <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.50, 
             0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999) 
#' @export
#' @rdname pprobs0
pprobs7 <- c(0.01, 0.025, 0.05, 
             0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
             0.95, 0.975, 0.99)
#' @export
#' @rdname pprobs0
pprobs8 <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
             0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
			 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999) 
#' @export
#' @rdname pprobs0
pprobs9 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
             0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
			 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999) 



#' @title Parameter Subsets
#' 
#' @description
#' Some vectors of parameter names to be used with parameter \code{exfitk} in 
#' functions regkienerLX(.., exfitk = ...) and \code{fitkienerX(.., exfitk = ...)}
#' or to subset the vector (or matrix) \code{fitk } obtained after regression 
#' \code{fitk <- regkienerLX(..)$fitk} or estimation \code{fitk <- fitkienerX(..)}. 
#' Visit \code{\link{fitkienerX}} for details on each parameter.
#' 
#' \code{exfit0 <- c("lh", "ret")}
#' 
#' \code{exfit1 <- c("m", "g", "a", "k", "w", "d", "e")}
#' 
#' \code{exfit2 <- c("m1", "sd", "sk", "ke", "m1x", "sdx", "skx", "kex")}
#' 
#' \code{exfit3 <- c("q.01", "q.05", "q.95", "q.99", "ltm.025", "rtm.975")}
#' 
#' \code{exfit4 <- c("VaR.01", "VaR.05", "VaR.95", "VaR.99", "ES.025", "ES.975")}
#' 
#' \code{exfit5 <- c("c.01", "c.05", "c.95", "c.99", "h.025", "h.975")}
#' 
#' \code{exfit6 <- c(exfit1, exfit2, exfit3, exfit4, exfit5)}
#' 
#' \code{exfit7 <- c(exfit0, exfit1, exfit2, exfit3, exfit4, exfit5)}
#' 
#' 
#' @examples     
#' 
#' require(minpack.lm)
#' require(timeSeries)
#' 
#' ### Load the datasets and select one number j in 1:16
#' j      <- 5
#' DS     <- getDSdata()
#' (fitk  <- regkienerLX(DS[[j]])$fitk)
#' fitk[exfit3]
#' fitkienerX(DS[[j]], exfitk = exfit3)
#' 
#' 
#' @export
#' @name exfit0
exfit0 <- c("lh", "ret") 
#' @export
#' @rdname exfit0
exfit1 <- c("m", "g", "a", "k", "w", "d", "e")
#' @export
#' @rdname exfit0
exfit2 <- c("m1", "sd", "sk", "ke", "m1x", "sdx", "skx", "kex")
#' @export
#' @rdname exfit0
exfit3 <- c("q.01", "q.05", "q.95", "q.99", "ltm.025", "rtm.975")
#' @export
#' @rdname exfit0
exfit4 <- c("VaR.01", "VaR.05", "VaR.95", "VaR.99", "ES.025", "ES.975")
#' @export
#' @rdname exfit0
exfit5 <- c("c.01", "c.05", "c.95", "c.99", "h.025", "h.975")
#' @export
#' @rdname exfit0
exfit6 <- c(exfit1, exfit2, exfit3, exfit4, exfit5)
#' @export
#' @rdname exfit0
exfit7 <- c(exfit0, exfit1, exfit2, exfit3, exfit4, exfit5)






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
#'      \deqn{ elevate(x, e) = (x + sqrt(x*x + e*e)) / 2 }
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
#'                   for the optimal value of \code{e}.
#' @param    dfrcol  integer. For data.frame only, designates the column that handles the time 
#'                   and must be processed separately. Use \code{dfrcol = 0} if all columns 
#'                   must be processed and there is no time (or turn the data.frame to a matrix).
#' @param    na.rm   boolean. Replace \code{x[t]=NA} with the previous non-NA value available 
#'                   in the price serie such that \code{(x[t-1], x[t]=x[t-1], x[t+1])} and 
#'                   calculate the returns accordingly. Force 0 in the first line of the returns 
#'                   if \code{x[1]=NA}.
#' 
#' @details
#' Negative prices in financial markets, like interest rates in Europe, are a 
#' nightmare as the rough calculation of the returns generates non-sense values.
#' \code{elevate} uses an hyperbola and implements the following formula: 
#'      \deqn{ elevate(x, e) = (x + sqrt(x*x + e*e)) / 2 }
#' 
#' There is currently no rule of thumb to calculate \code{e}. 
#' When \eqn{e = NULL}, there is no change and the output is identical to the input.
#' When \eqn{e = 0}, all negative values are turned to 0.
#' 
#' @examples 
#' 
#' fatreturns(extractData())
#' logreturns(getDSdata())
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

