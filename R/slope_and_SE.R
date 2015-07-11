##' Calculate the slope and standard error of fluorescence vs incubation time
##' 
##' @param x is a data frame containing columns `fl` and `incTime` (both numerics)


# slope and standard error calculation
# Note that this is superceded by corr_stats, but I didn't bother to update the functions that rely on it

slope_and_SE <- function(x) {
  sm <- summary(lm(fl ~ incTime, data=x))
  data.frame(slopePerHr = sm$coefficients[2, 1], slopeSE = sm$coefficients[2, 2])
}