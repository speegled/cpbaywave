#'Change point plotter
#'
#'
#'Plots change points. Shows best plot or the top three plots based on BFIC value.
#'
#'@param M1 Max of probability vector [calculated in changePointDetect]
#'@param M2 If M1 > 3, good evidence of change point [calculated in changePointDetect]
#'@param pad1 padding  [calculated in changePointDetect]
#'@param probvec idealized data set and its discrete wavelet transform [calculated in changePointDetect]
#'@param wid number of columns of data (width) [calculated in changePointDetect]
#'@param m row length [calculated in changePointDetect]
#'@param n length [calculated in changePointDetect]
#'@param data change point dataset [calculated in changePointDetect]
#'@param indices indices of change point [calculated in changePointDetect]
#'@param showall set to TRUE to see the top three candidate plots based on highest BFIC value
#'
#'@import stats
#'@import ggpubr
#'@import ggplot2
#'@import grid
#'
#'
#'
#'


plotChangePoint <- function(M1,M2,pad1,probvec,wid,m,n,data,indices,showall=FALSE){

  if (showall){
  BFIC <-  M2 - M1 - 0.5 * wid * log(m)
  plotData <- data.frame(x = 1:(n - pad1 - 1), y = data[(pad1+1):(n-1),1])
  probData <- data.frame(x = 1:(n - pad1 - 1), probvec = probvec[(pad1+1):(n-1)])

  if(BFIC > 3) {
    plot1 <- ggplot(plotData, aes_string(x = 'x', y = 'y')) +
      geom_point(alpha = .5) +
      geom_line(data = plotData[1:indices[1],], mapping = aes_string(x = 'x', y = 'y'), stat = "smooth", method = "loess",span =0.9 ,se = FALSE, alpha = 1, color = "blue") +
      geom_line(data = plotData[(indices[1] + 1):(n-pad1 - 1),], mapping = aes_string(x = 'x', y = 'y'), stat = "smooth", method = "loess",span = 0.6, se = FALSE, alpha = 1, color = "blue")
    plot2 <- ggplot(probData, aes_string(x = 'x', y = 'probvec')) +
      geom_line()

    rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last")




  } else {
    plot1 <- ggplot(plotData, aes_string(x = 'x', y = 'y')) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "loess",se = FALSE)
    plot2 <- ggplot(probData, aes_string(x = 'x', y = 'probvec')) +
      geom_line()
    rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last")


  }


  }
  else{
    BFIC <-  M2 - M1 - 0.5 * wid * log(m)
    plotData <- data.frame(x = 1:(n - pad1 - 1), y = data[(pad1+1):(n-1),1])
    probData <- data.frame(x = 1:(n - pad1 - 1), probvec = probvec[(pad1+1):(n-1)])
    if(BFIC > 3) {
      plot1 <- ggplot(plotData, aes_string(x = 'x', y = 'y')) +
        geom_point(alpha = .5) +
        geom_line(data = plotData[1:indices[1],], mapping = aes_string(x = 'x', y = 'y'), stat = "smooth", method = "loess",span =0.9 ,se = FALSE, alpha = 1, color = "blue") +
        geom_line(data = plotData[(indices[1] + 1):(n-pad1 - 1),], mapping = aes_string(x = 'x', y = 'y'), stat = "smooth", method = "loess",span = 0.6, se = FALSE, alpha = 1, color = "blue")
      plot2 <- ggplot(probData, aes_string(x = 'x', y = 'probvec')) +
        geom_line()
      grid.newpage()
      grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))



    } else {
      plot1 <- ggplot(plotData, aes_string(x = 'x', y = 'y')) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "loess",se = FALSE)
      plot2 <- ggplot(probData, aes_string(x = 'x', y = 'probvec')) +
        geom_line()
      grid.newpage()
      grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))


    }


  }
}


