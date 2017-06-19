#' Bootstrap Change Point Detection with Johnson-Lindenstrauss dimension reduction
#'
#' Bootstraps list of indices at which change point may have occured in a multi-dimensional time series with normal noise.
#' Samples the time component and interpolates between values when there are repeats. This function has not been thoroughly tested, and is not implemented in several key cases.
#'
#'
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param reducedDim The dimension you want to project onto. Should be less than the dimension of the time series. Default is 5
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Set to FALSE for random Bernoulli matrix.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param numRepeat Number of times to repeat random projection. Default is 100.
#' @param numKeep The number of indices that are considered significant at each iteration.
#' @param alpha The significance level of the bootstrap.
#' @param useJL Use Johnson-Lindenstrauss dimension reduction. Not currently implemented for useJL = FALSE.
#' @return List with values
#' \item{prop}{The proportion of times the BFIC of the bootstrapped time series showed significance.}
#' \item{indices}{A table of indices that were in the numKeep most likely change points at least once, together with the number of times that they occured.}
#' \item{sigindices}{The indices which occured more often than would be indicated by chance.}
#' \item{q95}{Estimate for the cutoff of the number of times indices would have to appear to be included in sigindices.}
#' @import stats
#' @import ggplot2
#' @export
#' @examples
#'
#' multiseries <- createTimeSeries(mu1 = rep(0, 100), mu2 = rep(1,100))
#' bootJLDetectChangePoint(multiseries) #True change point is at t = 72
#'
#'
#'
#'


bootJLDetectChangePoint <- function(multiSeries, reducedDim = 5, useGaussian = TRUE, setdetail, useBFIC = TRUE, numRepeat = 100, numKeep = 2, alpha = .05, useJL = TRUE) {

  if(ncol(multiSeries) <= reducedDim && useJL) {
    warning("reduced dimension greater than original, useJL set to FALSE")
    useJL <- FALSE
  }

  if(ncol(multiSeries) == 1) {
    warning("Algorithm Not Yet Implemented For This Case")
    return(1)
  }

  if(numRepeat > 200) {
    warning("numRepeat reset to 200.")
    numRepeat <- 200
  }
  if(missing(setdetail)) setdetail <- 0:(floor(log2(nrow(multiSeries) - 1) - 1))

  if(numKeep > 5) {
    numKeep <- 5
    warning("numKeep reset to 5.")
  }



  if(useJL) {
  N <- nrow(multiSeries)
  M <- ncol(multiSeries)
  results <- data.frame(value = 0, index = rep(0, numKeep))
  for(j in 1:numRepeat) {
    nTrue <- sort(sample(N,N,TRUE))
    a <- multiSeries
    a[1,] <- multiSeries[nTrue[1],]
    for(i in 2:N) {
      if(nTrue[i] == nTrue[i-1]) {
        minIndex <- max(1,nTrue[i]-5)
        maxIndex <- min(N,nTrue[i]+5)
        planck <- sd(multiSeries[minIndex:maxIndex,])
        a[i,] <- multiSeries[nTrue[i],] + rnorm(M,0,planck)
      } else {
        a[i,] <- multiSeries[nTrue[i],]
      }
    }
    scoreIndex1 <- cpbaywave::JLDetectChangePoint(as.matrix(a, ncol = ncol(multiSeries)), reducedDim = reducedDim, useGaussian= useGaussian, setdetail = setdetail, useBFIC = useBFIC, showplot = FALSE)
    scoreIndex <- data.frame(value = rep(scoreIndex1$value, numKeep), index = nTrue[scoreIndex1$index[1:numKeep]])
    results <- rbind(results, scoreIndex)
  }
  results <- results[-numKeep:-1,]
  }

  if(!useJL) {
    warning("Bootstrapping not yet implemented for this case")
    return(1)
  }
  percent <- sum(results$value > 3)/numRepeat/numKeep
  indices <- unlist(results$index)
  alphaIndex <- min(round(1000*(1 - alpha)) ,1000)
  q95 <- sort(replicate(1000, max(table(replicate(numRepeat, sample(sample(N,N/2),numKeep))))))[alphaIndex]
  values <- as.data.frame(table(results$index))
  names(values)[1] <- "Index"
  grid::grid.newpage()
  plot1 <-  ggplot(values, aes_string(x = "Index", y = "Freq")) + geom_bar(stat = "identity") +
    geom_abline(slope = 0, intercept = q95, color = "red", linetype = 2)
  grid::grid.draw(plot1)
  sigindices <- table(indices)[table(indices) > q95]
  return(list(prop = percent, indices = table(indices), sigindices = sigindices, q95 = q95))

}