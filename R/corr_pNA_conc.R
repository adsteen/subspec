##' Corrects pNA concentrations for my sloppy initial calculation
##' 
##' @description Initially I foolishly calculated pNA concentrations as if the total volume were only 1 mL. In reality, the total volume was 1 + 0.1 + 0.016 + 0.032 
##' @param x data frame with incorrectly calculated pNA concentrations in the `conc.pNA` column
##' @return The data frame with teh conc.pNA column corrected
##' @export

corr_pNA_conc <- function(x) {
  x$conc.pNA <- x$conc.pNA / (1 + 0.1 + 0.016 + 0.032)
  x
}