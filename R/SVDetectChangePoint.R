#' Singular Value of the Laplacian Change Point Detection
#'
#' Computes a random projection with Bernoulli or Gaussian entries of each element in a time series, then calls detectChangePoint to determine where the
#' change point occurs. In most cases of real data, the BFIC will be significantly greater than 3; it is an open
#' question to determine good values of the BFIC for various types of problems.
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Default is random matrix of +- 1.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param showplot set to TRUE to see plot of 1-d time series and probability plot.
#' @import stats
#' @import Rfast
#' @import grid
#' @import gridExtra
#'
#'
#'




SVDetectChangePoint <- function(multiSeries, useGaussian = FALSE, useBFIC = TRUE, showplot = TRUE) {


  cpbaywave::detectChangePoint(multiseries, useBFIC = useBFIC, showplot = showplot)

}
