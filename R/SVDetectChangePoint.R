#' Singular Value of the Laplacian Change Point Detection
#'
#' Implements change point detection for time series of graphs by examining the singular values of the Laplacian
#' @param multiSeries The high dimensional time series. Should be a matrix, with each row being an observation of the time series.
#' @param useGaussian Set to TRUE if you want to use a random Gaussian projection. Default is random matrix of +- 1.
#' @param setdetail Optional argument to set the detail level you wish to use. Default is all details.
#' @param useBFIC Optional argument to use BFIC to decide change point location.
#' @param showplot set to TRUE to see plot of 1-d time series and probability plot.
#' @param useNorm use normalized Laplacian
#' @import stats
#' @import Rfast
#' @import grid
#' @import gridExtra
#' @import igraph
#' @export
#'
#'#NOTES:
#'Take in list of time series graphs
#'Calculate the singular values of the laplacian
#'5-6 vertices (more than 50 use JL)
#'expects list of dataframes with source, dest, and weight, where first thing in list is graph at time 1
#'each row an observation (5 vertices is five columns)
#'try ... in function declaration for showall
#'change svd to eigen at some point (user parameter?)


SVDetectChangePoint <- function(timeSeriesList, useNorm = TRUE, useBFIC = TRUE, showplot = TRUE, ...) {

  df_list <- list()
  for (series in 1:length(timeSeriesList)){
    svds <- purrr::map_df(timeSeriesList[[series]], function(val) {
      gg <- igraph::graph_from_data_frame(dplyr::select(timeSeriesList[[series]], source, receive), directed = T)
      E(gg)$weight <- timeSeriesList[[series]]$val
      h <- data.frame(t(svd(graph.laplacian(gg))$d)) #graph fourier transform

    })
    df_list[[length(df_list) + 1]] <- svds
  }
  bb <- data.table::rbindlist(df_list)
  bb <- dplyr::bind_rows(df_list)

  if (length(timeSeriesList) > 50){
    cpbaywave::JLDetectChangePoint(as.matrix(bb),showplot=showplot, useBFIC = useBFIC, ... )
  }
  else{
    cpbaywave::detectChangePoint(as.matrix(bb), useBFIC = useBFIC,showplot=showplot)
  }

}






