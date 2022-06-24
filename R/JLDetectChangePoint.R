#' Johnson-Lindendstrauss Dimension Reduction
#'
#' Computes a random projection with Bernoulli or Gaussian entries of each element in a time series, then calls detectChangePoint to determine where the
#' change point occurs. In most cases of real data, the BFIC will be significantly greater than 3; it is an open
#' question to determine good values of the BFIC for various types of problems.
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param reducedDim The dimension you want to project onto. Should be less than the dimension of the time series. Default is 10
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Default is random matrix of +- 1.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param showplot set to TRUE to see plot of 1-d time series and probability plot.
#' @param showall set to TRUE to see all candidate plots [currently only shows highest BFIC value]
#' @export
#' @import stats
#' @examples
#'\dontrun{
#'
#'data(lennon) #Requires EMD package
#'lennon_ts <- matrix(as.vector(lennon) + rnorm(65536*120,0,1), nrow = 120, byrow = TRUE)
#'lennon_ts[80:120,7500:8000] <- lennon_ts[80:120,7500:8000] + 40
#'image(matrix(lennon_ts[1,], nrow = 256), col = gray(0:100/100))
#'image(matrix(lennon_ts[90,], nrow = 256), col = gray(0:100/100))
#'JLDetectChangePoint(lennon_ts, reducedDim = 10,useGaussian = TRUE)
#'
#'}


JLDetectChangePoint <- function(multiSeries, reducedDim = 5, useGaussian = FALSE, useBFIC = TRUE, setdetail, showplot = TRUE,showall=FALSE) {

  fullDim <- ncol(multiSeries)

  #Create random projection
  if(useGaussian) {
    transMatrix <- matrix(rnorm(fullDim * reducedDim), ncol = reducedDim)
  } else transMatrix <- 1/sqrt(fullDim) * matrix(sample(c(1,-1), fullDim*reducedDim, replace = TRUE), ncol = reducedDim)

  #Reduce Dimensionality of Data
  reducedData <- multiSeries %*% transMatrix
  sd.rd <- .2 * sapply(1:reducedDim, function(x) sd(reducedData[,x]))

  reducedData <- t(t(reducedData) + rnorm(nrow(reducedData * reducedDim), 0, sd.rd))

  if(missing(setdetail)) {
    if(showall){


      #showall currently only plots the graph with highest BFIC value. Need to add all plots in one.
      best_val <- -Inf
      for (number in 1:reducedDim){
        grph <- cpbaywave::detectChangePoint(reducedData[,number], useBFIC = useBFIC, showplot = FALSE)
        val <- grph$value
        if (val > best_val){
          best_val <- val
          best_dim <- number
        }

      }
      cpbaywave::detectChangePoint(reducedData[,best_dim], useBFIC = useBFIC, showplot = showplot)



    }
    else{   #BEST GRAPH CHOICE determined by range of index
      best_rng <- Inf

      for (number in 1:reducedDim){
        grph <- cpbaywave::detectChangePoint(reducedData[,number], useBFIC = useBFIC, showplot = FALSE)
        rng <- max(grph$index)-min(grph$index)
        if (rng < best_rng){
          best_rng <- rng
          best_dim <- number

        }

      }
      cpbaywave::detectChangePoint(reducedData[,best_dim], useBFIC = useBFIC, showplot = showplot)

      }

  }
  else { #same thing just with setdetail included (shows only best choice)

    for (number in 1:reducedDim){
      grph <- cpbaywave::detectChangePoint(reducedData[,number],setdetail, useBFIC = useBFIC, showplot = FALSE)
      rng <- max(grph$index)-min(grph$index)
      if (rng < best_rng){
        best_rng <- rng
        best_dim <- number

      }

    }

    cpbaywave::detectChangePoint(reducedData[,best_dim], setdetail, useBFIC = useBFIC, showplot = showplot)
  }
}
