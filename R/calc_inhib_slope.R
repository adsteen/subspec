# Calculate the slope of 1/nM.per.hr ~ conc.pNA; return it with pertinent info (e.g p value)
calc_inhib_slope <- function(x) {
  # Caculates the slope (and standard error, rsquared, and pvalue) of inverse activity vs conc.pNA - "m" in equation 4
  m <- lm(1/nM.per.hr ~ conc.pNA, data=x)
  inhib.slope <- coef(summary(m))[2, 1]
  inhib.slope.SE <- coef(summary(m))[2, 2]
  rsq <- summary(m)$r.squared
  pval <- coef(summary(m))[2, 4]
  data.frame(inhib.slope=inhib.slope, inhib.slope.SE=inhib.slope.SE, rsq=rsq, pval=pval)
}
