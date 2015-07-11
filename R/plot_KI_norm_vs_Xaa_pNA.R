##' Make a faceted dotplot of KI/KI,ref for each Xaa-pNA, Yaa-AMC and location
##' 
##' @param x The underlying data (all_rel_Ki)
##' @param vertical Whether to build a vertical or gridded plot (vertical is not up-to-date) 
##' @param print_plot Whether to print the plot
##' @param save_plot Whether to save the plot
##' @param fn Filename for the saved plot
##' @param sz 'Size' (i.e. thickness) of the lines
##' @return The plot object, which is discarded by write_paper()


plot_KI_norm_vs_Xaa_pNA <- function(x, vertical = FALSE, print_plot=TRUE, save_plot=TRUE, fn=NA, sz=0.25, ...) {
  # Make a plot of normalized KI values vs pNA identity at each site & for each 
  
  
  #x$signif <- TRUE
  #x$signif[x$pval < 0.05] <- FALSE

  
  # Points with insignificant slopes and error bars that span 1 will be reset to 1
  x$low.err[x$low.err < 1 & x$rel.Ki > 1 & x$pval > 0.05] <- 1
  
  # Buld the plots: I went back and forth about how I wanted the plot to look
  if(vertical) {
    # Make the plot single-column dotplot
   
    p_KI_vs_Xaa_pNA <- ggplot(x, aes(x=pNA.subs, y=rel.Ki, shape=is.homologous)) + 
      geom_errorbar(aes(ymin=low.err, ymax=hi.err), width=0.2) +
      geom_point(fill="white") +
      scale_shape_manual(values=c(1, 19), guide=FALSE) + 
      geom_hline(yintercept=1) +
      xlab("inhibitor") +
      scale_y_log10() +
      coord_cartesian(ylim=10^c(-1.1, 1.1)) +
      #coord_flip() +
      #facet_wrap( ~ location + AMC.substrate, nrow=5) 
      facet_grid(AMC.substrate ~ location )
    
  } else {
    # Make the plot as a barplot
    p_KI_vs_Xaa_pNA <- ggplot(x, aes(x=pNA.subs, y=rel.Ki, shape=is.homologous)) + 
      #geom_blank() +###
      #geom_ribbon(aes(ymin=0.001, ymax=0), colour="black", alpha=1) +
      geom_errorbar(aes(ymin=low.err, ymax=hi.err), width=0.2) +
      geom_point(fill="white") +
      scale_shape_manual(values=c(1, 19), guide=FALSE) + 
      geom_hline(yintercept=1, size=sz) +
      xlab("inhibitor") +
      scale_y_log10() +
      coord_cartesian(ylim=10^c(-1.1, 1.4)) +
      annotation_logticks(sides="l", 
                          short=unit(0.05, "cm"), mid=unit(0.1, "cm"), long=unit(0.15, "cm"),
                          size=sz) +
      #facet_wrap( AMC.substrate ~ location, nrow=3, scales="free") +
      facet_grid( AMC.substrate ~ location, scales="free") +
      theme(axis.text.x=element_text(angle=-45, hjust=0))
  }
  #browser()
  p_KI_vs_Xaa_pNA <- p_KI_vs_Xaa_pNA +
    #geom_hline(yintercept=1) +
    #scale_x_discrete("inhibitor", breaks=levels(all_rel_Ki$pNA.subs), labels=substr(levels(all_rel_Ki$pNA.subs), start=1, stop=3)) +
    scale_x_discrete("inhibitor", breaks=levels(x$pNA.subs), labels=substr(levels(x$pNA.subs), start=1, stop=3)) +
    ylab(expression(paste(K[I]/K["I, ref"]))) +
    theme(axis.ticks=element_line(size=sz))
  
#   if (print_plot) {
#     print(p_KI_vs_Xaa_pNA)
#   }
#   
#   if (save_plot) {
#     if (is.na(fn)) {
#       fn <- paste(path, "KI_vs_Xaa_pNA.png", sep="")
#     }
#     ggsave(fn, p_KI_vs_Xaa_pNA, ..., units="in", type="cairo")
#   }
  p_KI_vs_Xaa_pNA
}