##' Analyzes and plots all saturation curves
##' 
##' @description Relies on 4 data files in /data/: sat_arg_TN, sat_leu_TN, sat_curves_IMS, sat_curve_calib
##' @param print_plot Whether to print plots
##' @param save_plot Whether to save plots
##' @param fn filename for saved plots
##' @export

make_sat_curves <- function(print_plot=TRUE, save_plot=FALSE, fn=NA) {
  
  ###########
  # arg_AMC, TN river
  ###########
  #sat_arg_TN <- subset(read.csv("2012-11-27ArgvsArgpna.csv"), Conc.Argp==0)
  sat_arg_TN$incTime <- (sat_arg_TN$Time - min(sat_arg_TN$Time))*24
  attr(sat_arg_TN, "unit") <- "hours"
  
  sat_arg_TN$AMC.subs <- "arg-AMC"
  sat_arg_TN$location <- "TN River"
  
  sat_arg_TN <- rename(sat_arg_TN, c("Conc.Arg" = "conc.AMC", "Fl" = "fl"))
  
  ###########
  # leu-AMC, TN river
  ###########
  #sat_leu_TN <- read.csv("2013-02-10 leuamcstandardcurve.csv")
  sat_leu_TN$time <- ymd_hm(paste("2013-02-10", sat_leu_TN$Time))
  sat_leu_TN$incTime <- as.numeric(sat_leu_TN$time - min(sat_leu_TN$time))/3600
  attr(sat_leu_TN$incTime, "unit") <- "hours"
  
  sat_leu_TN$AMC.subs <- "leu-AMC"
  sat_leu_TN$location <- "TN River"
  
  # rename sat_leu_TN columns appropriately
  sat_leu_TN <- rename(sat_leu_TN, c("Conc.leuamc" = "conc.AMC", "Fl" = "fl"))
  
  ###########
  # arg_AMC, IMS
  ###########
  #sat_curves_IMS <- read.csv("../2013-03-25 Tuesday sat curve.csv")
  sat_curves_IMS$R.time <- ymd_hm(paste("2013-03-25", sat_curves_IMS$time))
  sat_curves_IMS$incTime <- as.numeric(sat_curves_IMS$R.time - min(sat_curves_IMS$R.time))/3600
  attr(sat_curves_IMS$incTime, "unit") <- "hours"
  sat_curves_IMS$location <- "Bogue Sound"
  sat_curves_IMS <- rename(sat_curves_IMS, c("conc.pNA" = "conc.AMC")) #The pNA name is an artifact
  #   of the fact that I used the same worksheet to record sat curve data and inhibition data
  sat_curves_IMS$AMC.subs <- as.character(sat_curves_IMS$AMC.subs)
  
  #########
  # Bind the data
  #########
  cols_to_keep <- c("conc.AMC", "fl", "incTime", "AMC.subs", "location")
  sat_list <- list(sat_arg_TN[ , cols_to_keep],
                   sat_leu_TN[sat_leu_TN$Live.dead=="live", cols_to_keep],
                   sat_curves_IMS[ , cols_to_keep])
  
  sat_data <- rbind.fill(sat_list)
  
  # Rename the substrates
  sat_data$AMC.subs[sat_data$AMC.subs=="arg-AMC"] <- "Arg-AMC"
  sat_data$AMC.subs[sat_data$AMC.subs=="leu-AMC"] <- "Leu-AMC"
  sat_data$AMC.subs[sat_data$AMC.subs=="pro-AMC"] <- "Pro-AMC"
  
  # Calc slopes
  sat_slopes <- ddply(sat_data, c("conc.AMC", "AMC.subs", "location"), slope_and_SE)
  
  #### Calibration
  sat_curve_calib$AMC.nM <- sat_curve_calib$AMC.uM*1000
  calib_slope <- coefficients(lm(fl~AMC.nM, data=sat_curve_calib))[2] # slope of the calibration curve
  
  # Perform the calibration
  sat_slopes$nM.per.hr <- sat_slopes$slopePerHr / calib_slope
  sat_slopes$se.nM.per.hr <- sat_slopes$slopeSE / calib_slope

  # Fit MM curves
  mm_form <- formula(I(nM.per.hr ~ (Vmax * conc.AMC)/(Km + conc.AMC)))
  fit_list <- dlply(sat_slopes, c("AMC.subs", "location"), function(x) nls(mm_form, data=x, start=list(Vmax=100, Km=100)))
  
  g <- data.frame(conc.AMC=0:360)
  data.frame(conc.AMC=g$conc.AMC, nM.per.hr = predict(fit_list[[1]], newdata=data.frame(g)))
  mm_fits <- ldply(fit_list, function(x) data.frame(conc.AMC=g$conc.AMC, nM.per.hr = predict(x, newdata=data.frame(g))))
  
  # Create the plot
  p_sat <- ggplot(sat_slopes, aes(x=conc.AMC, y=nM.per.hr)) + geom_point() + 
    geom_errorbar(aes(ymin=nM.per.hr - se.nM.per.hr, ymax=nM.per.hr + se.nM.per.hr), width=5) +
    geom_line(data=mm_fits, aes(x=conc.AMC, y=nM.per.hr)) +
    xlab(expression(paste("[S], ", mu, "M"))) + 
    ylab(expression(paste(v[0], ", nM ", hr^{-1}))) +
    facet_grid(AMC.subs ~ location, scales="free_y") 
  
  if (print_plot) {
    print(p_sat)
  }
  if (save_plot) {
    if (is.na(fn)) {
      fn <- paste(path, "all_sat_curves.png", sep="")
    }
    ggsave(fn, p_sat, height=4, width=3.5, units="in", type="cairo")
  }
  
  p_sat
}