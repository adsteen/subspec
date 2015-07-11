##' Makes fig 4
##' 
##' @param d data frame 
##' @export

plot_all_inv_v0 <- function(d, print_plot=TRUE, save_plot=FALSE, fn=NA, height=4, width=7, dpi=300, ...) {
  # Make a plot of all inverse v0 vs inhibitor concentration
  # browser()
  p_all_inv_v0 <- ggplot(d, aes(x=conc.pNA, y=1/nM.per.hr)) + 
    geom_point(size=1) +
    geom_errorbar(aes(ymin=1/(nM.per.hr + se.nM.per.hr), ymax=1/(nM.per.hr - se.nM.per.hr))) +
    geom_smooth(method="lm", colour="black") +
    #scale_colour_manual(values=c("#56A0D3", "#F77F00")) +
    #scale_fill_manual(values=c("#56A0D3", "#F77F00")) +
    scale_x_continuous(breaks=c(0, 100, 200), labels=c(0, 100, 200)) +
    expand_limits(y=0) +
    xlab(expression(paste("[I], ", mu, "M"))) +
    ylab(expression(paste(1 / v[0], ", nM ", hr^{-1}))) +
    facet_grid(location + AMC.substrate ~ pNA.subs, scales="free_y") +
    theme(legend.position="top",
          axis.text.x = element_text(angle=-45, hjust=0)) 
  if (print_plot) {
    print(p_all_inv_v0)
  }
  
  #browser()
  if (save_plot) {
    if (is.na(fn)) {
      fn <- paste(path, "all_inv_v0.png", sep="")
    }
    ggsave(fn, height = height, width=width, units="in", dpi=myDPI, type="cairo")
  }
  
  p_all_inv_v0
}