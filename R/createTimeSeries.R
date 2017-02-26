#' Create Time Series
#'
#' This function creates a multi-dimensional nomral time series with normal entries and a possible change point in the mean. Running with no
#' arguments creates five dimensional normal time series with mean 0, mean 1 and identity covariance matrix.
#' @param mu1 A vector of length dimension giving mean before change point, or scalar if constant across dimensions.
#' @param mu2 A vector of length dimesnion giving mean after change point, or scalar if constant across dimensions.
#' @param n The length of the time series to create.
#' @param tau The location of the change point.
#' @param sigma The covariance matrix used for generation. Defaults to the identity.
#' @return A vector of length n containing multivariate normal entries.
#' @export
#' @examples
#' createTimeSeries()

createTimeSeries <- function(mu1, mu2, sigma, n = 128, tau = 72) {

    if (missing(mu1) && missing(mu2) && missing(sigma)) {
        mu1 <- rep(0, 5)
        mu2 <- rep(1, 5)
        sigma <- diag(5)
        dimension <- 5
    }

    if (!missing(sigma))
        dimension <- nrow(as.matrix(sigma))

    if (missing(sigma) && !missing(mu1)) {
        dimension <- length(mu1)
        sigma <- diag(dimension)
    }


    if (length(mu1) != length(mu2) || length(mu1) != nrow(as.matrix(sigma)) || length(mu2) != nrow(as.matrix(sigma)) || nrow(as.matrix(sigma)) != ncol(as.matrix(sigma))) {
        stop("Dimensions of means and/or covariance matrix do not match.")
    }



    if (length(mu1) == 1)
        mu1 <- rep(mu1, dimension)
    if (length(mu2) == 1)
        mu2 <- rep(mu2, dimension)

    s1 = MASS::mvrnorm(tau, mu1, sigma, empirical = FALSE)
    s2 = MASS::mvrnorm(n - tau, mu2, sigma, empirical = FALSE)

    s = rbind(s1, s2)
    return(s)

}
