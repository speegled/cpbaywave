#' Detect Change Point Faster
#'
#' Applies Bayesian-wavelet technique to determine whether a multi-dimensional time series has change point in the mean. Returns the
#' 3-5 most likely change point locations.
#' @param a A vector or matrix representing time series. If matrix, each row is the value at a single time. Numerically unstable if data dimension is greater than about 50. Use JLDetectChangePoint in that case.
#' @param setdetail The detail levels of the wavelet transform to use to detect change in mean. Default is all levels.
#' @param useBFIC Set to true to choose the change point with highest BFIC.
#' @param showplot Set to true to see a plot of the probabilities of a change point at each time, together with a scatterplot of the first dimension versus time.
#' @param padding One of mirror, extension or insertion. Default is insertion.
#' @param slow Set to TRUE if you want to check each point separately, rather than fast search for change point.
#' @return value The value of the BFIC or maximium probability if useBFIC = FALSE. BFIC greater than 3 is evidence that there is a change in mean.
#' @return index A vector giving the 3-5 most likely (or highest IC if useBFIC is TRUE) indices where a change point occurred.
#' @export
#' @import stats
#' @import wavethresh
#' @import utils
#' @import ggplot2
#' @import grid
#' @examples
#' a <- createTimeSeries() #True change point at time 72
#' detectChangePoint(a, 0:6, showplot = TRUE)
#'
#' a <- createTimeSeries(mu1 = 1, mu2 = 0, sigma = 1, n = 100, tau = 55)
#' plot(a[,1])
#' detectChangePoint(a) #Hard problem. True change point at t = 55. Not always possible to detect.
#'
#' #Time series with no change in mean.
#' dim=5
#' sig1=diag(dim)
#' mu1=rep(0,dim)
#' mu2=rep(0,dim)
#' n=128
#' tau=70
#' series=createTimeSeries(mu1, mu2, sig1, n, tau)
#' plot(series[,1])
#' detectChangePoint(series)
#' #value less than 3 indicates no change model is favored
#'
#'
#' #Time series with smooth mean function with jump.
#' dim=10
#' sig1=diag(dim)
#' mu1=rep(0,dim)
#' mu2=rep(1,dim)
#' n=90  #Algorithm doesn't work as well when n not power of 2.
#' tau=20
#' series1=createTimeSeries(mu1, mu2, sig1, n, tau)
#' sn=sin(1*pi*1:n/n) #Compare with sn = sin(3*pi*1:n/n), which is too hard for algorithm
#' series2=sn+series1
#' detectChangePoint(series2, useBFIC = TRUE, showplot = TRUE) #True change point at t = 20.
#'
#'
#' #Time series with smooth mean function, no jump. Illustrates that BFIC unreliable
#' #when there is a smooth underlying mean function.
#' dim=3
#' sig1=diag(dim)
#' mu1=rep(0,dim)
#' mu2=rep(0,dim)
#' n=256
#' tau=105
#' series1=createTimeSeries(mu1, mu2, sig1, n, tau)
#' sn=sin(2*pi*1:n/n)
#' series2=sn+series1
#' plot(series2[,2])
#' detectChangePoint(series2, showplot = TRUE)
#' #value less than 3 indicates no change-model is favored





detectChangePointFast <- function(a, setdetail, useBFIC = TRUE, showplot = FALSE, showall=FALSE, padding = "insertion", slow = FALSE) {

  if(is.vector(a)) a <- as.matrix(a,ncol = 1)
  if(ncol(a) > 100) stop("dimension too high. consider JLdetectChangePoint")

  if(missing(setdetail)) setdetail <- 0:floor(log2(nrow(a) - 1))
  F <- 10
  n <- nrow(a)
  if(slow || n <= 128) {
    detectChangePoint(a, setdetail, useBFIC, showplot, showall, padding)
  } else {
    wid <- ncol(a)
    isDataVector <- FALSE

    if(wid > 50) warning("This computation is likely numerically unstable. Consider JLdetectChangePoint instead.")

    # pad with normal data at top of a, centered at first element of time series
    nxt <- 2^(ceiling(log2(n)))
    pad1 <- nxt - n
    data <- a
    isPadded <- (pad1 > 0)
    #if(browse)
    #  browser()


    #
    # BEGIN Lots of different kinds of padding. Insertion works best on simulations.
    #
    #Padding by extending estimated first entry + some noise.
    if(padding == "extend" && pad1 > 0) {
      xvals <- 1:10
      padVals <- matrix(rep(0, pad1 * wid), ncol = wid)
      for(i in 1:wid) {
        lm.model <- lm(a[1:10,i]~xvals)
        padVals[,i] <- predict(lm.model, newdata = data.frame(xvals = (-1 * pad1 + 1):0)) + rnorm(pad1, 0, mean(sd(lm.model$residuals), sd(a[1:10,i])))
      }
      data <- rbind(padVals, a)
    }
    #data <- rbind(matrix(rnorm(pad1 * wid, a[1,], 0.1), ncol = wid), a)
    if(padding == "mirror" && pad1 > 0)
      data <- rbind(as.matrix(a[pad1:1,], ncol = wid), a)
    #Insertion padding. Inserts similar looking data between randomly chosen points.
    if(padding == "insertion" && pad1 > 0) {
      data <- matrix(rep(0, nxt * wid), ncol = wid)
      insertAfter <- sample(n, pad1)
      for(Cols in 1:wid) {
        valueToInsert <- numeric(0)
        for(ind in insertAfter) {
          minIndex <- max(1,ind-5)
          maxIndex <- min(n,ind+5)
          planck <- sd(a[minIndex:maxIndex,Cols])
          meanVal <- weighted.mean(a[minIndex:maxIndex, Cols], w = 5 - abs(ind - minIndex:maxIndex))
          valueToInsert <- c(valueToInsert, meanVal + rnorm(1, 0,planck))
        }
        #Ninja trick to insert new values into vector.
        id <- c(seq_along(a[,Cols]), insertAfter + 0.5)
        vals <- c(a[,Cols], valueToInsert)
        data[,Cols] <- vals[order(id)]
      }
      pad1 <- 0 #Treat the signal as original; use insertAfter to recalculate indices later
    }
    #
    # END different types of padding.
    #

    # re-establish length
    n <- nxt

    # Get information about wd, could be improved
    genDWT <- (wd(data[, 1], filter.number = F, family = "DaubExPhase"))  #generic DWT
    p <- nlevelsWT(genDWT)
    J <- p

    # Compute details of discrete wavelet transform of data by column and store desired details in DWTmat TODO: write
    # wrapper for accessD that pulls out multiple levels and returns a vector of details unlist(sapply(...)) pulls out
    # the desired levels of detail coefficients and combines them in a vector
    DWTmat <- apply(data, 2, function(x) unlist(sapply(J - setdetail - 1, function(y) accessD(wd(x,
                                                                                                 filter.number = F, family = "DaubExPhase"), y))))

    # Creating idealized data set and its discrete wavelet transform Have 0's followed by 1's with change point in each
    # possible position
    probvec <- rep(-Inf, n - 1)
    skip = ceiling(n/64 + 1)
    probvec[seq(3, n - 3, by = skip)] <- sapply(seq(3, n - 3, by = skip), function(x) {
      tauvec <- c(rep(0, x), rep(1, nxt - x))
      Qvec <- unlist(sapply(J - setdetail - 1, function(y) accessD(wd(tauvec, filter.number = F,
                                                                      family = "DaubExPhase"), y)))
      computeProb(DWTmat, Qvec, useBFIC, isDataVector)
    })

    check_further <- sort(probvec, index.return = T, decreasing = T)$ix[1:5]
    check_further <- unique(as.vector(sapply(check_further, function(x) {
      c(x - 1:(skip - 1), x + 1:(skip - 1))
    })))
    check_further <- check_further[check_further > 3 & check_further < (n - 3)]
    probvec[check_further] <- sapply(check_further, function(x) {
      tauvec <- c(rep(0, x), rep(1, nxt - x))
      Qvec <- unlist(sapply(J - setdetail - 1, function(y) accessD(wd(tauvec, filter.number = F,
                                                                      family = "DaubExPhase"), y)))
      computeProb(DWTmat, Qvec, useBFIC, isDataVector)
    })

    probvec <- ifelse(probvec == -Inf, min(probvec), probvec)
    # If M1 > 3, good evidence of change point -
    M2 = max(probvec)
    m <- nrow(DWTmat)
    # t1 = gamma(m/2 + 0.5)
    A <- t(DWTmat) %*% DWTmat
    t2 = det(A)  #Numerically unstable.
    M1 = lgamma(m/2 + 0.5) + (-m/2 + wid/2) * log(t2)


    ifelse(useBFIC, value <- (M2 - M1 - 0.5 * wid * log(m)), value <- max(probvec))
    indices <- match(head(sort(probvec[(pad1 + 1):n], decreasing = TRUE), 5), probvec[(pad1 + 1):n])

    if(padding == "insertion" && isPadded) {
      n <- nrow(a)
      indices <- ifelse(order(id)[indices] <= n, order(id)[indices], order(id)[indices - 1]) #Picks the point in unpadded time series that corresponds to change point
      probvec <- probvec[order(id) <= n]
      data <- a
    }


    if (showplot) {

      if(showall){
        plt <- plotChangePoint(M1,M2,pad1,probvec,wid,m,n,data,indices,showall= TRUE)
        return(plt)
      }
      else{
        plotChangePoint(M1,M2,pad1,probvec,wid,m,n,data,indices,showall= FALSE)
      }
    }


    # Remove nan and return the 3-5 most likely change points There shouldn't be nan, maybe rather exit with error if
    # there are?
    probvec[is.nan(probvec)] <- min(probvec)


    return(list(value = value, index = unique(indices)))

  }
}
