#' Repeated JL Dimension Reduction
#'
#' Repeatedly computes a random projection with Bernoulli entries of each element in a time series, then calls detectChangePoint to determine where the
#' change point occurs. Creates histogram of locations of numKeep most likely change points. 95% of the time, all of the bins will be below the red dashed line
#' if the change points are chosen completely randomly. Note that this is not the case, because the likely change points in any random projection of a fixed
#' time series are likely correlated. As such, this plot requires some experience to read correctly. The current behavior is that change points are only recorded
#' if BFIC is greater than 3.
#'
#' Many times when there is no change point, times in the first few or last few possible locations will appear significant. It is recommended to always ignore
#' times that are in the first 8 or last eight spots when using this technique. Eight is a rough guide, and not written in stone, and seems to depend on the
#' filter length.
#'
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param reducedDim The dimension you want to project onto. Should be less than the dimension of the time series. Default is 10
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Default is random matrix of +- 1.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param numRepeat Number of times to repeat random projection. Default is 100.
#' @param numKeep number of indices to keep. Recommend 2-3 percent of time series length. Max is 5.
#' @param alpha alpha level of the test
#' @return percent The percentage of runs for which BFIC was bigger than 3. (Not useful.)
#' @return indices A table of indices which were included as most likely change points and the number of times they were included.
#' @return sigindices A list of indices which occured more often than would be expected if completely random.
#' @return q95 This is the cutoff for being a significant index.
#' @return results The BFIC and indices of the most likely change points.
#' @import stats
#' @import ggplot2
#' @export
#' @examples
#'
#' mu1 <- rep(0,100)
#' mu2 <- rep(1, 100)
#' sigma <- diag(100)
#' series <- createTimeSeries(mu1, mu2, sigma) #Change point at 72
#' repJLDetectChangePoint(series)
#'
#' mu1 <- rep(0,100)
#' mu2 <- rep(0, 100)
#' sigma <- diag(100)
#' series <- createTimeSeries(mu1, mu2, sigma) #No change point
#' repJLDetectChangePoint(series)
#'
#'
#'
#' \dontrun{
#' data(lennon) #Requires EMD package
#' lennon_ts <- matrix(as.vector(lennon) + rnorm(65536*120,0,1), nrow = 120, byrow = TRUE)
#' lennon_ts[80:128,7500:8000] <- lennon_ts[80:128,7500:8000] + 40
#' image(matrix(lennon_ts[1,], nrow = 256), col = gray(0:100/100))
#' image(matrix(lennon_ts[90,], nrow = 256), col = gray(0:100/100))
#' repJLDetectChangePoint(lennon_ts, numKeep = 3)
#' }
#'
#'

repJLDetectChangePoint <-
  function(multiSeries,
           reducedDim = 5,
           useGaussian = TRUE,
           setdetail,
           useBFIC = TRUE,
           numRepeat = 100,
           numKeep = 2,
           alpha = .05) {

   if (numRepeat > 200) {
      warning("numRepeat reset to 200.")
      numRepeat <- 200
   }


    if (missing(setdetail))
      setdetail <- 0:floor(log2(nrow(multiSeries) - 1))

    if (numKeep > 5) {
      numKeep <- 5
      warning("numKeep reset to 5.")
    }

    sd <- numeric(ncol(multiSeries))

    for (i in 1:ncol(multiSeries)) {
      tempdata <- data.frame(x = 1:nrow(multiSeries), y = multiSeries[, i])
      fit.loess <- loess(y ~ x, data = tempdata)
      sd[i] <- sd(fit.loess$residuals)
    }

    results <- data.frame(value = 0, index = rep(0, numKeep))
    for (i in 1:numRepeat) {
      tempdata <- multiSeries
      tempdata <-
        tempdata + matrix(rnorm(nrow(tempdata) * ncol(tempdata), 0, sd),
                          ncol = ncol(tempdata),
                          byrow = TRUE)
      scoreIndex1 <-
        cpbaywave::JLDetectChangePoint(
          multiSeries,
          reducedDim = reducedDim,
          useGaussian = useGaussian,
          setdetail = setdetail,
          useBFIC = useBFIC,
          showplot = FALSE
        )
      scoreIndex <-
        data.frame(value = rep(scoreIndex1$value, numKeep),
                   index = scoreIndex1$index[1:numKeep])
      results <- rbind(results, scoreIndex)
    }

    results <- results[-numKeep:-1, ]
    percent <- sum(results$value > 3) / numRepeat / numKeep
    indices <- unlist(results$index)
    alphaIndex <- min(round(1000 * (1 - alpha)) , 1000)
    q95 <-
      sort(replicate(1000, max(table(
        replicate(numRepeat, sample(nrow(multiSeries), numKeep))
      ))))[alphaIndex]
    values <- as.data.frame(table(results$index))
    names(values)[1] <- "Index"
    grid::grid.newpage()
    plot1 <-
      ggplot(values, aes_string(x = 'Index', y = 'Freq')) + geom_bar(stat = "identity") +
      geom_abline(
        slope = 0,
        intercept = q95,
        color = "red",
        linetype = 2
      )
    grid::grid.draw(plot1)
    sigindices <- table(indices)[table(indices) > q95]
    return(list(
      percent = percent,
      indices = table(indices),
      sigindices = sigindices,
      q95 = q95,
      results = results
    ))
}


