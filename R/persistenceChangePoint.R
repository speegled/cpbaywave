#' Persistence Diagram Change Point Detection
#'
#' Computes the persistence diagrams of a time series of multi-dimensional data. Determines where the most likely change point occurs.
#'
#' The algorithm computes the wasserstein (L^2) distance between each pair of persistence diagrams, and finds the most likely
#' change point in the distances. Note that the multiSeries in this case must be a *list* of matrices with N columns, where each row corresponds
#' to a point in R^N.
#'
#'
#' @param multiSeries A list of matrices, where each row of the matrix corresponds to a point in Euclidean space.
#' @param maxDimension Maximum dimensional topological features to detect in persistence diagram.
#' @param maxScale Maximum scale to consider in persistence diagrams.
#' @param useWasserstein set to False to use bottleneck L^\\infty distance.
#' @param useBootstrap set to TRUE if you want to bootstrap the change point locations. Not advisable unless the time series length is equal to a power of 2 or slighlty less than a power of 2.
#' @param ... additional parameters passed to the distance function
#' @return value The BFIC value associated with change point. Bigger than 3 is considered significant.
#' @return index The five most likely change points.
#' @import stats
#' @import ggplot2
#' @import TDA
#' @export
#' @examples
#'
#' Circles <- lapply(1:30, function(x) {
#' Circle1 <- TDA::circleUnif(60);
#' Circle2 <- TDA::circleUnif(60, r = 2) + 3;
#' rbind(Circle1, Circle2);})
#' Circles2 <- lapply(1:34, function(x) {
#'   Circle1 <- TDA::circleUnif(40);
#'   Circle2 <- TDA::circleUnif(40, r = 2) + 5;
#'   Circle3 <- TDA::circleUnif(40,r = 3) - 10;
#'   rbind(Circle1, Circle2, Circle3);})
#' Circles <- c(Circles, Circles2)
#' persistenceChangePoint(Circles, maxDimension = 1, maxScale = 2) #True change point at t = 31
#'


persistenceChangePoint <- function(multiSeries, maxDimension, maxScale, useWasserstein = TRUE, useBootstrap = FALSE, ...) {
  if(maxDimension > 3) {
    warning("maxdimension reset to 3")
    maxDimension <- 3
  } else if(maxDimension == 2) {
    warning("Consider whether maxdimension 2 is necessary, if time is issue")
  }

  N <- length(multiSeries)

  #
  # Compute the persistence diagrams of each list entry
  #
  Diags <- lapply(1:N, function(i) ripsDiag(X = multiSeries[[i]],
                                             maxdimension = maxDimension,
                                             maxscale = maxScale,
                                             library = "GUDHI",
                                             printProgress = FALSE)$diagram)
  metricChangePoint(Diags, wasserstein, useBootstrap = useBootstrap, ...)
}
