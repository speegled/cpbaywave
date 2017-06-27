# cpbaywave
(Bayesian) Wavelet Change Point Detection

Detect the change point in mean in high-dimensional time series using dimension reduction JLDetectChangePoint followed by a wavelet technique for estimating the change point. For lower dimensional time series, use detectChangePoint, which does not do dimension reduction. Currently, the Bayesian aspect of this algorithm is not implemented, and the only prior allowed is a uniform prior. 

Change points of time series with values in arbitrary metric spaces are also detected via metricChangePoint. Finally, change points in persistence diagrams of spatial data are detected using persistenceChangePoint. Unfortunately, the TDA library that is required for persistenceChangePoint is not working with R 3.4 for some Mac users. It does work on R 3.3.2 (Sincere Pumpkin Patch), which is why this package has only been tested relative to R 3.3.2. I delayed releasing this in hopes that TDA would get fixed.

To install, use devtools::install_github("speegled/cpbaywave"). 
