##' calculate a Gaussian distribution
##' 
##' @param x vector of values on which to calculate the distribution
##' @param mu Mean
##' @param sigma standard deviation
##' @details You would think there would be a function for this in base R, but if there is I can't find it
##' @export
gauss <- function(x, mu, sigma) {
  (1/(sigma * sqrt(2*pi))) * exp(-1*(x-mu)^2/(2*sigma^2))
}
