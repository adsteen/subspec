##' Gets coefficients and standard errors for an nls object
##' @param x an nls object
##' @value a 1-row data frame. Yes, this is dumb: it should be a vector
##' @export

get_nls_coefs <- function(x) {
  # x is a nls fit
  coeffs <- summary(x)$coefficients
  data.frame(Km=coeffs[1, 1],
             Km.se=coeffs[1, 2],
             Vmax=coeffs[2, 1],
             Vmax.se=coeffs[2, 2])
}