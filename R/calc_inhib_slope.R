##' Calculate the slope of 1/nM.per.hr ~ conc.pNA; return it with pertinent info (e.g p value)
##' 
##' Caculates the slope (and standard error, rsquared, and pvalue) of inverse activity vs conc.pNA - "m" in equation 4.
##' @description This is just one of three functions that (essentially) replicates the functionality of lm_stats
##' @param x a data frame with nM.per.hr and conc.pNA
##' @export
calc_inhib_slope <- function(x) {
  # Caculates the slope (and standard error, rsquared, and pvalue) of inverse activity vs conc.pNA - "m" in equation 4
  m <- lm(1/nM.per.hr ~ conc.pNA, data=x)
  inhib.slope <- coef(summary(m))[2, 1]
  inhib.slope.SE <- coef(summary(m))[2, 2]
  rsq <- summary(m)$r.squared
  pval <- coef(summary(m))[2, 4]
  data.frame(inhib.slope=inhib.slope, inhib.slope.SE=inhib.slope.SE, rsq=rsq, pval=pval)
}
