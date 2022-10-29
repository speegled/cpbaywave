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
#' @param rotate_xaxis set to TRUE if you wish to rotate the values on the x-axis; will not print out all values
#' @param returnPlot set to TRUE if you wish for the plot to be returned.
#' @param fast set to FALSE if you wish to check every possible change point.
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


bootJLDetectChangePoint <- function(multiSeries,
                                    reducedDim = 5,
                                    useGaussian = TRUE,
                                    setdetail,
                                    useBFIC = TRUE,
                                    numRepeat = 100,
                                    numKeep = 2,
                                    alpha = .05,
                                    useJL = TRUE,
                                    rotate_xaxis = FALSE,
                                    returnPlot = FALSE,
                                    fast = TRUE) {

  if(ncol(multiSeries) < reducedDim && useJL) {
    warning("reduced dimension greater than original, useJL set to FALSE")
    useJL <- FALSE
  }

  if(is.na(IsPowerOfTwo(nrow(multiSeries))))
    warning("Algorithm designed for time series length a power of 2.")

  if(ncol(multiSeries) == 1) {
    stop("Algorithm Not Yet Implemented For 1-d Case")
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
      nTrue <- sort(sample(N,N,TRUE))   #Bootsrap resample of time coordinates
      a <- multiSeries
      a[1,] <- multiSeries[nTrue[1],]
      for(i in 2:N) {
        #
        # If we repeated a time sample, then we interpolate a new value.
        #
        if(nTrue[i] == nTrue[i-1]) {
          minIndex <- max(1,nTrue[i]-5)
          maxIndex <- min(N,nTrue[i]+5)
          planck <- sd(multiSeries[minIndex:maxIndex,])
          a[i,] <- multiSeries[nTrue[i],] + rnorm(M,0,planck)
        } else {
          a[i,] <- multiSeries[nTrue[i],]
        }
      }
      scoreIndex1 <- cpbaywave::JLDetectChangePoint(as.matrix(a, ncol = ncol(multiSeries)), reducedDim = reducedDim, useGaussian= useGaussian, setdetail = setdetail, useBFIC = useBFIC, showplot = FALSE, fast = fast)
      scoreIndex <- data.frame(value = rep(scoreIndex1$value, numKeep), index = nTrue[scoreIndex1$index[1:numKeep]])
      results <- rbind(results, scoreIndex)
    }
    results <- results[-numKeep:-1,]
  }

  if(!useJL) {
    stop("Bootstrapping only implemented via JL so far.")
  }

  percent <- sum(results$value > 3)/numRepeat/numKeep #Percent of times BFIC > 3
  indices <- unlist(results$index)
  alphaIndex <- min(round(1000*(1 - alpha)) ,1000)
  #Estimate the significance cut-off. TODO: make this a table lookup.
  q95 <- sort(replicate(1000, max(table(replicate(numRepeat, sample(sample(N,N/2),numKeep))))))[alphaIndex]

  values <- as.data.frame(table(results$index))
  names(values)[1] <- "Index"

  if(rotate_xaxis)
    values$Index <- as.integer(as.character(values$Index))
  plot1 <-  ggplot(values, aes_string(x = "Index", y = "Freq")) + geom_bar(stat = "identity") +
    geom_abline(slope = 0, intercept = q95, color = "red", linetype = 2)
  if(rotate_xaxis)
    plot1 <- plot1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  sigindices <- table(indices)[table(indices) > q95]
  if(returnPlot)
    return(list(prop = percent, indices = table(indices), sigindices = sigindices, q95 = q95, plot = plot1))
  if(!returnPlot) {
    grid::grid.newpage()
    grid::grid.draw(plot1)
    return(list(prop = percent, indices = table(indices), sigindices = sigindices, q95 = q95))
  }

}
