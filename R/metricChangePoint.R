#' Change Point Detection in Metric Spaces
#'
#' Determines the most likely locations of a change point in distribution of elements in a metric space. Forms a new time series which is
#' the pairwise distance matrix; the diagonal is replaced by average value of nearby elements. Then, dimension reduction change
#' point detection is called. If useBootstrap is TRUE, then a bootstrapped confidence interval for possible change points is computed.
#'
#'
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param distance A function that returns the distance between two rows of multiSeries.
#' @param reducedDim The dimension you want to project onto. Should be less than the dimension of the time series. Default is 10.
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Set to FALSE for random Bernoulli matrix.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details except one.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param useBootstrap Set to true if you wish to bootstrap a cutoff for significance of indices.
#' @import ggplot2
#' @import testthat
#' @export
#' @examples
#'
#' dim=10
#' sig1=diag(dim)
#' mu1=rep(0,dim)
#' mu2=rep(1.5,dim)
#' n=128
#' tau=30
#' series1=createTimeSeries(mu1, mu2, sig1, n, tau)
#' sn=sin(2.5*pi*1:n/n)
#' series2=sn+series1
#' myDistance <- function(x,y) {
#'   sum(abs(x - y))
#' }
#' metricChangePoint(series2, myDistance) #True change point is at t = 30
#'
#'










metricChangePoint <- function(multiSeries, distance, useGaussian = TRUE, useBFIC = TRUE, setdetail, reducedDim = 10, useBootstrap = FALSE) {

  if(!is.list(multiSeries)) {
    warning("multiSeries must be a list")
    return (1)
  }

  if(test_that(is.numeric(distance(multiSeries[[1]], multiSeries[[2]])), is_true)) {

    if(missing(setdetail)) setdetail <- 0:(floor(log2(length(multiSeries) - 1) - 1))

    N <- length(multiSeries)
    dists <- matrix(1:N, ncol = N)
    for(i in 1:N) {
      dists <- rbind(dists, sapply(1:N, function(j) distance(multiSeries[[i]], multiSeries[[j]])))
    }
    dists <- dists[-1,]
    for(i in 1:N) {
      for(j in 1:N) {
        if(dists[i,j] == 0 && i > 1 && i < N) dists[i,j] <- 1/2*(dists[i,j-1] + dists[i,j+1])
      }
    }
    dists[1,1] <- dists[1,2]
    dists[N,N] <- dists[N,N - 1]
    if(!useBootstrap) {
      JLDetectChangePoint(multiSeries = dists, reducedDim = reducedDim, useGaussian = useGaussian, setdetail = setdetail, useBFIC = useBFIC)
    } else {
      bootJLDetectChangePoint(multiSeries = multiSeries, useGaussian = useGaussian, setdetail = setdetail, useBFIC = useBFIC)
    }


  } else {
    stop("distance between rows of multiseries non-numeric")
  }
}
