

#' @include n_regression3.R



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
#' The dataset is usually called by the instruction \code{DS <- getDSdata()}.
#' Then, it is used with a call to DS[[j]] with j in 1:16. 
#' \enumerate{
#'   \item{ 1: "USDCHF" (USDCHF, timeSeries) }
#'   \item{ 2: "MSFT" (MSFT, timeSeries) }
#'   \item{ 3: "DAX" (EuStockMarkets, datasets) }
#'   \item{ 4: "SMI" (EuStockMarkets, datasets) }
#'   \item{ 5: "CAC" (EuStockMarkets, datasets) }
#'   \item{ 6: "FTSE" (EuStockMarkets, datasets) }
#'   \item{ 7: "SBI" (LPP2005REC, timeSeries) }
#'   \item{ 8: "SPI" (LPP2005REC, timeSeries) }
#'   \item{ 9: "SII" (LPP2005REC, timeSeries) }
#'   \item{10: "LMI" (LPP2005REC, timeSeries) }
#'   \item{11: "MPI" (LPP2005REC, timeSeries) }
#'   \item{12: "ALT" (LPP2005REC, timeSeries) }
#'   \item{13: "LPP25" (LPP2005REC, timeSeries) }
#'   \item{14: "LPP40" (LPP2005REC, timeSeries) }
#'   \item{15: "LPP60" (LPP2005REC, timeSeries) }
#'   \item{16: "sunspot" (sunspot.year, datasets) }
#' }
#' Note that \code{sunspot.year} is regularly updated with each new version of  
#' \code{R}. 
#' 
#' @seealso    
#' \code{\link[datasets]{EuStockMarkets}}, \code{\link[datasets]{sunspot.year}}, 
#' \code{\link[timeSeries]{TimeSeriesData}}, \code{\link{regkienerLX}}, 
#' \code{\link{estimkienerX}}
#' 
#' @examples    
#' 
#' require(timeSeries)
#' getDSdata
#' DS  <- getDSdata()
#' attributes(DS)
#' for (j in 1:16) {print(head(DS[[j]]))}
#' 
#' @export
#' @name getDSdata
getDSdata <- function() {
DSenv <- new.env(parent = baseenv())
data("USDCHF", "MSFT", "LPP2005REC", package = "timeSeries", envir = DSenv)
prices2returns <- function(x) { 100*diff(log(as.numeric(x))) }
DS <- list(
	"USDCHF"	= prices2returns(DSenv$USDCHF),
	"MSFT"		= prices2returns(DSenv$MSFT[,4]),
	"DAX"		= prices2returns(EuStockMarkets[,1]),
	"SMI"		= prices2returns(EuStockMarkets[,2]),
	"CAC"		= prices2returns(EuStockMarkets[,3]),
	"FTSE"		= prices2returns(EuStockMarkets[,4]),
	"SBI"		= 100*as.numeric(DSenv$LPP2005REC[,1]),
	"SPI"		= 100*as.numeric(DSenv$LPP2005REC[,2]),
	"SII"		= 100*as.numeric(DSenv$LPP2005REC[,3]),
	"LMI"		= 100*as.numeric(DSenv$LPP2005REC[,4]),
	"MPI"		= 100*as.numeric(DSenv$LPP2005REC[,5]),
	"ALT"		= 100*as.numeric(DSenv$LPP2005REC[,6]),
	"LPP25"		= 100*as.numeric(DSenv$LPP2005REC[,7]),
	"LPP40"		= 100*as.numeric(DSenv$LPP2005REC[,8]),
	"LPP60"		= 100*as.numeric(DSenv$LPP2005REC[,9]),
	"sunspot"	= prices2returns(sunspot.year+1000) )
return(DS)
}




#' @title Datasets dfData, tData, xData, zData.
#'
#' @description
#' A list of 10 financial series in data.frame, timeSeries, xts and zoo formats 
#' used in the document presented at 8th and 9th Rmetrics conferences (2014, 2015). 
#' See references. The datasets include dates and prices. 
#' All distributions of the returns exhibit fat tails. 
#' The last one (CHF, interest rates in Switzerland) exhibits negative prices!
#'
#' @details
#' A list of 10 financial series presented in four different formats for 
#' convenient use: dfData, tData, xData and zData are respectively of 
#' data.frame format (with dates as row.names), timeSeries, xts and zoo formats. 
#' dfData is the complete dataset whereas tData, xData and zData display are subsets 
#' and display one column per stock. 
#' \enumerate{
#'   \item{ 1: "Gold" from 1999-01-04 to 2013-12-31, dim 3694x1 (df, t, x, z). }
#'   \item{ 2: "Dexia" from 2008-10-27 to 20009-10-26, dim 255x5 (df), dim 255x1 (t, x, z). }
#'   \item{ 3: "SocGen" from 1992-07-20 to 2013-12-31, dim 5445x5 (df), dim 5445x1 (t, x, z). }
#'   \item{ 4: "Vivendi" from 1992-07-20 to 2013-12-31, dim 5444x5 (df), dim 5444x1 (t, x, z). }
#'   \item{ 5: "Eurodollar" from 1999-01-03 to 2013-12-31, dim 3843x4 (df), dim 3843x1 (t, x, z). }
#'   \item{ 6: "VIX" from 2004-01-02 to 2013-12-31, dim 2517x4 (df), dim 2517x1 (t, x, z). }
#'   \item{ 7: "CAC40" from 1988-01-04 to 2013-12-31, dim 6574x4 (df), dim 6574x1 (t, x, z). }
#'   \item{ 8: "DJIA" from 1896-05-26 to 2013-12-31, dim 32064x1 (df, t, x, z). }
#'   \item{ 9: "SP500" from 1957-01-02 to 2013-12-31, dim 14350x1 (df, t, x, z). }
#'   \item{10: "CHF" from 1995-01-02 to 2013-09-13, dim 4880x8 (df), dim 4880x1 (t, x, z). 
#'             Include negative prices (interest rates) at the end of the dataset. 
#'             Care is required to calculate the returns! }
#' }
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
#' @seealso    
#' \code{\link{tData}}, \code{\link{xData}}, \code{\link{zData}}.
#'  
#' @examples    
#' 
#' require(timeSeries)
#' require(xts)
#' attributes(dfData); attributes(tData); attributes(xData); attributes(zData) 
#' for (j in 1:10) {print(head(dfData[[j]]))}
#' for (j in 1:10) {print(head(tData[[j]]))}
#' for (j in 1:10) {print(head(xData[[j]]))}
#' for (j in 1:10) {print(head(zData[[j]]))}
#' 
#' @keywords datasets
#' @docType data
#' @name dfData
NULL

#' @title tData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the timeSeries format. 
#' Visit \code{\link{dfData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name tData
NULL

#' @title xData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the xts format. 
#' Visit \code{\link{dfData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name xData
NULL

#' @title zData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the zoo format. 
#' Visit \code{\link{dfData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name zData
NULL


