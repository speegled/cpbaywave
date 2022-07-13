#' Compute Probability Vector Elements
#'
#' This is a helper function for changePointDetect.
#' @param DWTmat The matrix of discrete wavelet transform of data values.
#' @param Qvec The discrete wavelet transform of the idealized time series with shift.
#' @param useBFIC set to TRUE to compute Bayes factor information content.
#' @param isDataVector set to TRUE for 1-d time series
#' @export


computeProb <- function(DWTmat, Qvec, useBFIC = FALSE, isDataVector = FALSE) {

    m <- nrow(DWTmat)
    wid <- ncol(DWTmat)

    if(!isDataVector) {
    # Compute A
    A <- t(DWTmat) %*% DWTmat


    # compute B vector
    B <- as.vector(t(DWTmat) %*% Qvec)

    # compute C scalar
    C <- sum(Qvec^2)

    G <- A - (1/C) * ((B %*% t(B)))

    H <- sum(log(diag(chol(G))))

    } else{
      A <- sum(DWTmat^2)
      B <- sum(DWTmat * Qvec)
      C <- sum(Qvec^2)
      G <- A - (1/C) * B^2
      H <- log(G)
    }
    #GOES TO INF
    BFIC <- (-0.5) * log(C) + (-(m - wid - 1)/2) * log((det(A - (1/C) * ((B %*% t(B)))))) + (wid/2) * log(pi) + lgamma(m/2 + 0.5 - wid/2)



    return(ifelse(useBFIC, BFIC, (-0.5) * log(C) + (-(m - wid - 1)/2) * H))
}
