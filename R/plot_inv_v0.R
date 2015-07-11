# Plot the 1/v0 plots for each site & substrate
plot_inv_v0 <- function(x, wd=5, fs = 12, print_plot=FALSE, save_plot=FALSE, fn=NA, sz=1, spacing=0.1, n_cols=2) {
  
  
  p <- ggplot(x, aes(x=conc.pNA, y=1/slopePerHr, ymax=1/(slopePerHr-slopeSE), ymin=1/(slopePerHr+slopeSE))) +
    geom_point(size=sz) +
    geom_errorbar(width=wd) +
    geom_smooth(method="lm", colour="black") +
    xlab(expression(paste("[I], ", mu, "M"))) +
    ylab(expression(paste("1/", v[0]))) +
    #scale_x_continuous(breaks=c(0, 200)) +
    #facet_grid(pNA.subs ~ AMC.substrate+location, scales="free") +
    facet_wrap(~pNA.subs, scales="free", ncol=n_cols) +
    #expand_limits(y=0) +
    theme(text=element_text(size=fs),
          axis.text.x=element_text(angle=-45, hjust=0),
          panel.margin=unit(c(spacing), "inches"))
  
  # Facet the plot as defined by the rows parameters
  
  
  #if (is.na(rows)) {
  #  p <- p + facet_wrap(~pNA.subs) 
  #} else {
  #  p <- p + facet_wrap(~pNA.subs, nrow=rows) 
  #}
  
  ## Add a title, if requested (supp figs need titles, main ms fig does not)
  #if (use_title) {
  #  if (is.na(ti_text)) {
  #    ti_text <- paste(x$AMC.substrate[1], "in", x$location)
  #  }
  #  p <- p + ggtitle(ti_text)
  #}
  
#   # Print the plot to screen
#   if(print) {
#     print(p)
#   }
#   
#   # Save the plot, using either a user-supplied filename or an auto-generated filename
#   if(save) {
#     if(is.na(fn)) {
#       fn <- paste(path, x$.id[1], ".png", sep="")
#     }
#     png(fn, width=6, height=6, res=300, units="in", type="cairo")
#     print(p)
#     dev.off()
#   }
  #browser()
  p
}