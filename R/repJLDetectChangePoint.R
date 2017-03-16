#' Repeated JL Dimension Reduction
#'
#' Computes a random projection with Bernoulli entries of each element in a time series, then calls detectChangePoint to determine where the
#' change point occurs. In most cases of real data, the BFIC will be significantly greater than 3; it is an open
#' question to determine good values of the BFIC for various types of problems.
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param reducedDim The dimension you want to project onto. Should be less than the dimension of the time series. Default is 10
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Default is random matrix of +- 1.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param numRepeat Number of times to repeat random projection. Default is 100.
#' @export
#' @examples
#'\dontrun{
#'
#'data(lennon) #Requires EMD package
#'lennon_ts <- matrix(as.vector(lennon) + rnorm(65536*120,0,1), nrow = 120, byrow = TRUE)
#'lennon_ts[80:128,7500:8000] <- lennon_ts[80:128,7500:8000] + 40
#'image(matrix(lennon_ts[1,], nrow = 256), col = gray(0:100/100))
#'image(matrix(lennon_ts[90,], nrow = 256), col = gray(0:100/100))
#'repJLDetectChangePoint(lennon_ts, numKeep = 3)
#'
#'}
#'

repJLDetectChangePoint <- function(multiseries, reducedDim = 5, useGaussian = TRUE, setdetail, useBFIC = TRUE, numRepeat = 100, numKeep = 2, alpha = .05) {

  if(numRepeat > 200) {
    warning("numRepeat reset to 200.")
    numRepeat <- 200
  }
  if(missing(setdetail)) setdetail <- 0:floor(log2(nrow(multiseries) - 1))

  results <- data.frame(value = 0, index = rep(0, numKeep))
  for(i in 1:numRepeat) {
    scoreIndex1 <- cpbaywave::JLDetectChangePoint(multiseries, reducedDim = reducedDim, useGaussian= useGaussian, setdetail = setdetail, useBFIC = useBFIC, showplot = FALSE)
    scoreIndex <- data.frame(value = rep(scoreIndex1$value, numKeep), index = scoreIndex1$index[1:numKeep])
    results <- rbind(results, scoreIndex)
 }
  results <- results[-numKeep:-1,]
  percent <- sum(results$value > 3)/numRepeat/numKeep
  indices <- unlist(results$index)
  alphaIndex <- min(round(1000*(1 - alpha)) ,1000)
  q95 <- sort(replicate(1000, max(table(replicate(numRepeat, sample(nrow(multiseries),numKeep))))))[alphaIndex]
  values <- as.data.frame(table(results$index))
  names(values)[1] <- "Index"
  grid::grid.newpage()
  plot1 <-  ggplot(values, aes(x = Index, y = Freq)) + geom_bar(stat = "identity") +
    geom_abline(slope = 0, intercept = q95, color = "red", linetype = 2)
  grid::grid.draw(plot1)
  sigindices <- table(indices)[table(indices) > q95]
  return(list(percent = percent*5, indices = table(indices), sigindices = sigindices, q95 = q95))

}
