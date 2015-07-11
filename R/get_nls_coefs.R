##' Gets coefficients and standard errors for an nls object
##' @param x an nls object
##' @value a 1-row data frame. Yes, this is dumb: it should be a vector
##' @export

get_nls_coefs <- function(x) {
  # x is a nls fit
  coeffs <- summary(x)$coefficients
  data.frame(k=coeffs[1, 1],
             k.se=coeffs[1, 2],
             recalcitrant=coeffs[2, 1],
             recalcitrant.se=coeffs[2, 2])
}