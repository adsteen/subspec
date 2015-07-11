##' Corrects fluorescence for pNA quenching
##' 
##' @description Based on the model $fl_{corrected} = fl_{obs} \times e^{k\times[\text{pNA}]}$
##' @param x a data frame containing a column called origFl
##' @param k an exponential coefficient (e^(kx))
##' @return the same data frame with a corrected column called fl
##' @export


correct_for_pNA_quench <- function(x, k) {
  # Correct fluorescence values for quenching by pNA
  x$origFl <- x$fl
  x$fl <- x$fl * exp(k*x$conc.pNA)
  x
}