##' Calculates the exponential correction factor for fluorescence quenching by pNA
##' @description Relies on the data frame pNA_calib.Rdata
##' @param printInfPlot prints plot.
##' @return the exponenetial correction factor $k$
##' @export

calc_pNA_k <- function(printInfPlot=FALSE) {
  ### Determine effect of pNA on specific AMC fluorescences
  
  # Note: This file is for data at 
  #calib <- read.csv("2013-07-18 calib curve with pNA.csv")
  pNA_calib$pNA.subs <- as.factor(pNA_calib$pNA.subs)
  
  meanUnquenched <- mean(pNA_calib$fl[pNA_calib$conc.pNA==0])
  
  pNA_calib$relUnquenched <- pNA_calib$fl / meanUnquenched
  
  ggplot(pNA_calib, aes(x=conc.pNA, y=relUnquenched, colour=buffer)) + geom_point() + #geom_line() +
    aes(ymin=0) +
    geom_smooth(method="lm", se=FALSE) #+
  #facet_wrap(~pNA.subs)
  
  # So, I conclude that there is no particular effect of pNA substrate, but the overall effect is quite large.
  #      Also, I don't really want to deal with doing a proper curve for each substrate, and I think it would
  #      just create noise rather than real signal.
  
  #expForm <- formula(I(relUnquenched ~ A*exp(k*conc.pNA)))
  expForm <- formula(I(relUnquenched ~ exp(-1*k*conc.pNA)))
  
  #expMod <- nls(expForm, calib, start=list(A=1, k=0.01), model=TRUE)
  expMod <- nls(expForm, pNA_calib, start=list(k=0.01), model=TRUE)
  summary(expMod)
  # A=0.9853
  k <- coefficients(expMod)[1]
  # Currently 1.33e-3 +/- 3.4e-5
  
  grid <- 0:360
  myPred <- exp(-1* grid*k)
  myPredDF <- data.frame(x=grid, y=myPred)
  
  # Plot the resutls. They're not so bad.
  pNAinfPlot <- ggplot() + 
    geom_point(data=pNA_calib, aes(x=conc.pNA, y=relUnquenched, colour=pNA.subs)) + 
    geom_line(data=myPredDF, aes(x=x, y=y)) +
    aes(ymin=0) 
  if(printInfPlot) {
    print(pNAinfPlot)
  }
  
  
  return(k)
  
}