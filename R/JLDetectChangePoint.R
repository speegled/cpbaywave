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
#' @param showall set to TRUE to see the top three candidate plots based on highest BFIC value
#' @export
#' @import stats
#' @import Rfast
#' @import grid
#' @import gridExtra
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
#'


JLDetectChangePoint <- function(multiSeries, reducedDim = 5, useGaussian = FALSE, useBFIC = TRUE, setdetail, showplot = TRUE,showall=FALSE, fast = TRUE) {

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
      #showall plots top 3 BFIC valued graphs

      best_val <- -Inf
      sec_best_val <- -Inf
      thrd_best_val <- -Inf
      vector <- c()
      for (number in 1:reducedDim){
        grph <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,number], useBFIC = useBFIC, showplot = FALSE))
        val <- grph$value
        vector <- c(vector, val)
      }

      best_vector <- c()
      if (reducedDim > 2){
        for (number in 1:3){
          best_dim <- Rfast::nth(vector, number, descending = T,index.return=TRUE)
          best_vector <- c(best_vector,best_dim)

        }

        print('Plot 1')
        plt1 <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = showplot,showall=TRUE))
        print(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = FALSE))
        print('Plot 2')
        plt2 <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,best_vector[2]], useBFIC = useBFIC, showplot = showplot,showall=TRUE))
        print(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = FALSE))
        print('Plot 3')
        plt3 <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,best_vector[3]], useBFIC = useBFIC, showplot = showplot,showall=TRUE))
        print(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = FALSE))

        gridExtra::grid.arrange(plt1, plt2,plt3, nrow = 1,ncol=3)

      }
      else if (reducedDim == 2){
        for (number in 1:2){
          best_dim <- Rfast::nth(vector, number, descending = T,index.return=TRUE)
          best_vector <- c(best_vector,best_dim)
        }
        print('Plot 1')
        plt1 <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = showplot,showall=TRUE))
        print(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = FALSE))
        print('Plot 2')
        plt2 <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,best_vector[2]], useBFIC = useBFIC, showplot = showplot,showall=TRUE))
        print(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = FALSE))
        gridExtra::grid.arrange(plt1, plt2, nrow = 1,ncol=2)
      }
      else{
        best_dim <- Rfast::nth(vector, 1, descending = T,index.return=TRUE)
        best_vector <- c(best_vector,best_dim)
        print('Plot 1')
        plt1 <- suppressMessages(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = showplot,showall=TRUE))
        print(cpbaywave::detectChangePoint(reducedData[,best_vector[1]], useBFIC = useBFIC, showplot = FALSE))
        gridExtra::grid.arrange(plt1, nrow = 1,ncol=1)
      }





    }

    else{   #BEST GRAPH CHOICE by BFIC value

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

  }
  else { #same thing just with setdetail included (shows only best choice)

    best_val <- -Inf
    for (number in 1:reducedDim){
      grph <- cpbaywave::detectChangePoint(reducedData[,number],setdetail, useBFIC = useBFIC, showplot = FALSE)
      val <- grph$value
      if (val > best_val){
        best_val <- val
        best_dim <- number
      }

    }

    cpbaywave::detectChangePoint(reducedData[,best_dim], setdetail, useBFIC = useBFIC, showplot = showplot)
  }
}
