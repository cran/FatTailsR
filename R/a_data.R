## #' @include FatTailsR-package.R



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





