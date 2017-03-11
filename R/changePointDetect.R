#' Detect Change Point
#'
#' Applies Bayesian-wavelet technique to determine whether a multi-dimensional time series has change point in the mean.
#' @param a A vector or matrix representing time series. If matrix, each row is the value at a single time. Numerically unstable if data dimension is greater than about 50. Use JLDetectChangePoint in that case.
#' @param setdetail The detail levels of the wavelet transform to use to detect change in mean. Default is all levels.
#' @param useBFIC Set to true to choose the change point with highest BFIC.
#' @param showplot Set to true to see a plot of the probabilities of a change point at each time, together with a scatterplot of the first dimension versus time.
#' @return value The value of the BFIC or maximium probability if useBFIC = FALSE. BFIC greater than 3 is evidence that there is a change in mean.
#' @return index A vector giving the 5 most likely (or highest IC if useBFIC is TRUE) indices where a change point occurred.
#' @export
#' @examples
#' a <- createTimeSeries() #True change point at time 72
#' detectChangePoint(a, 0:6, showplot = TRUE)
#'
#' a <- createTimeSeries(mu1 = 1, mu2 = 0, sigma = 1, n = 100, tau = 55)
#' plot(a[,1])
#' detectChangePoint(a)
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
#' dtl=c(0,1,2,3,4,5,6)
#' detectChangePoint(series,dtl)
#' #value less than 3 indicates no change model is favored
#'
#'
#' #Time series with smooth mean function with jump.
#' dim=10
#' sig1=diag(dim)
#' mu1=rep(0,dim)
#' mu2=rep(1.5,dim)
#' n=90
#' tau=20
#' series1=createTimeSeries(mu1, mu2, sig1, n, tau)
#' sn=sin(2*pi*1:n/n)
#' series2=sn+series1
#' plot(series2[,2])
#' dtl=c(0,1,2,3,4,5)
#' detectChangePoint(series2,dtl, useBFIC = TRUE)
#' #value greater than 3 indicates change model is favored
#' detectChangePoint(series2,dtl)
#'
#'
#' #Time series with smooth mean function, no jump
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
#' dtl=c(0,1,2,3,4,5)
#' detectChangePoint(series2,dtl)
#' #value less than 3 indicates no change-model is favored



detectChangePoint <- function(a, setdetail, useBFIC = TRUE, showplot = FALSE) {

    if(is.vector(a)) a <- matrix(a,ncol = 1)
    if(ncol(a) > 100) stop("dimension too high. consider JLdetectChangePoint")

    if(missing(setdetail)) setdetail <- 0:floor(log2(nrow(a) - 1))
    F <- 10
    n <- nrow(a)
    wid <- ncol(a)
    isDataVector <- FALSE

    if(wid > 50) warning("This computation is likely numerically unstable. Consider JLdetectChangePoint instead.")

    #If vector, create matrix with both columns the same. This is hack that should be fixed.
    if(wid == 1) {
      isDataVector <- TRUE
    #a <- matrix(c(a,a + runif(n, -.01,.01)),ncol = 2)
    #wid <- 2
    }
    # pad with normal data at top of a, centered at first element of time series
    nxt <- 2^(ceiling(log2(n)))
    pad1 <- nxt - n
    data <- a
    #Padding by extending first entry + some noise.
    #data <- rbind(matrix(stats::rnorm(pad1 * wid, a[1,], 0.1), ncol = wid), a)

    #Trying mirror padding. Seems better based on simulations.
    if(pad1 > 0) data <- rbind(a[(pad1+1):1,], a)

    #Padding by repeating first entry
    #if(pad1 > 0) data <- rbind(t(replicate(pad1, a[1,])),a)

    #Padding by repeating first entry
    #if(pad1 > 0) data <- rbind(t(replicate(pad1, a[1,])),a)

    # re-establish length
    n <- nxt

    # Get information about wd, could be improved
    genDWT <- (wavethresh::wd(data[, 1], filter.number = F, family = "DaubExPhase"))  #generic DWT
    p <- wavethresh::nlevelsWT(genDWT)
    J <- p

    # Compute details of discrete wavelet transform of data by column and store desired details in DWTmat TODO: write
    # wrapper for accessD that pulls out multiple levels and returns a vector of details unlist(sapply(...)) pulls out
    # the desired levels of detail coefficients and combines them in a vector
    DWTmat <- apply(data, 2, function(x) unlist(sapply(J - setdetail - 1, function(y) wavethresh::accessD(wavethresh::wd(x,
        filter.number = F, family = "DaubExPhase"), y))))

    # Creating idealized data set and its discrete wavelet transform Have 0's followed by 1's with change point in each
    # possible position
    probvec <- sapply(1:(n - 1), function(x) {
        tauvec <- c(rep(0, x), rep(1, nxt - x))
        Qvec <- unlist(sapply(J - setdetail - 1, function(y) wavethresh::accessD(wavethresh::wd(tauvec, filter.number = F,
            family = "DaubExPhase"), y)))
        computeProb(DWTmat, Qvec, useBFIC, isDataVector)
    })




    # If M1 > 3, good evidence of change point
    M2 = max(probvec)
    m <- nrow(DWTmat)
    t1 = gamma(m/2 + 0.5)
    A <- t(DWTmat) %*% DWTmat
    t2 = det(A)  #Numerically unstable.
    M1 = log(t1) + (-m/2 + wid/2) * log(t2)


    ifelse(useBFIC, value <- (M2 - M1 - 0.5 * wid * log(m)), value <- max(probvec))
    indices <- match(utils::head(sort(probvec[(pad1 + 1):n], decreasing = TRUE), 5), probvec[(pad1 + 1):n])

    if (showplot) {
      BFIC <-  M2 - M1 - 0.5 * wid * log(m)
      if(isDataVector) {
        plotData <- data.frame(x = 1:(n - pad1 - 1), y = data[(pad1+1):(n-1),1])
        probData <- data.frame(x = 1:(n - pad1 - 1), probvec = probvec[(pad1+1):(n-1)])
        if(BFIC > 3) {
        plot1 <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(data = plotData[1:indices[1],], mapping = ggplot2::aes(x = x, y = y), method = "loess", se = FALSE) +
          ggplot2::geom_smooth(data = plotData[(indices[1] + 1):(n-pad1 - 1),], mapping = ggplot2::aes(x = x, y = y), method = "loess", se = FALSE)
        plot2 <- ggplot2::ggplot(probData, ggplot2::aes(x = x, y = probvec)) +
          ggplot2::geom_line()
        grid::grid.newpage()
        grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2), size = "last"))
        } else {
          plot1 <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_point() +
            ggplot2::geom_smooth(method = "loess", se = FALSE)
          plot2 <- ggplot2::ggplot(probData, ggplot2::aes(x = x, y = probvec)) +
            ggplot2::geom_line()
          grid::grid.newpage()
          grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2), size = "last"))
          }
      } else {
        plotData <- data.frame(x = 1:(n - pad1 - 1), y = data[(pad1+1):(n-1),1])
        probData <- data.frame(x = 1:(n - pad1 - 1), probvec = probvec[(pad1+1):(n-1)])
        if(BFIC > 3) {
          plot1 <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_point(alpha = .5) +
            ggplot2::geom_line(data = plotData[1:indices[1],], mapping = ggplot2::aes(x = x, y = y), stat = "smooth", method = "loess", span = 1.5, se = FALSE, alpha = 1, color = "blue") +
            ggplot2::geom_line(data = plotData[(indices[1] + 1):(n-pad1 - 1),], mapping = ggplot2::aes(x = x, y = y), stat = "smooth", method = "loess", span = 1.5, se = FALSE, alpha = 1, color = "blue")
          plot2 <- ggplot2::ggplot(probData, ggplot2::aes(x = x, y = probvec)) +
            ggplot2::geom_line()
          grid::grid.newpage()
          grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2), size = "last"))
        } else {
          plot1 <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_point(alpha = 0.5) +
            ggplot2::geom_smooth(method = "loess",span = 1.5, se = FALSE)
          plot2 <- ggplot2::ggplot(probData, ggplot2::aes(x = x, y = probvec)) +
            ggplot2::geom_line()
          grid::grid.newpage()
          grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2), size = "last"))
        }
      }
    }

    # Remove nan and return the 5 most likely change points There shouldn't be nan, maybe rather exit with error if
    # there are?
    probvec[is.nan(probvec)] <- min(probvec)

    indices <- match(utils::head(sort(probvec[(pad1 + 1):n], decreasing = TRUE), 5), probvec[(pad1 + 1):n])

    return(list(value = value, index = indices))
}
