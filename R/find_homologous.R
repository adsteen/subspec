##' Finds pNA substrates homologous to a AMC substrates (i.e., Leu-pNA for Leu-AMC)
##' @description This function is terrible and if you are paying attention to it you should re-evaluate the choices you have made in life.
##' @param x  A data frame with columns pNA.subs and AMC.substrate (you like the antiparallelism?)
##' @return The same data frame with a new column, `is.homologous`, with boolean for whether `pNA.subs` is homologous to `AMC.substrate`
##' @export

find_homologous <- function(x) {
  
  
  # Right now this is pretty brittle: it will fail if the pNA & AMC substrate prefixes are not exactly 3 characters. 
  #     I should improve this at some point.
  x$is.homologous <- FALSE
  pNA.substrate <- tolower(substr(x$pNA.subs, start=1, stop=3))
  AMC.substrate <- tolower(substr(x$AMC.substrate, start=1, stop=3))
  #browser()
  x$is.homologous[pNA.substrate==AMC.substrate] <- TRUE
  x
}