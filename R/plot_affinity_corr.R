##' Make plots of affinity correlations to DI, MW, and delta G of reaction
##' 
##' @param x should be `inhib_data` from script
##' @param label_size
##' @param print_plots NOT USED
##' @param save_plots NOT USED
##' @param fn for saved plots NOT USED

#plot_affinity_corr <- function(x, label_size=2, print_plots=TRUE, save_plots=FALSE, fn=NA, hjust=0, vjust=0, sz=0.25, inhib=4, letter_size=3, ...) {
plot_affinity_corr <- function(x, label_size=2, hjust=0, vjust=0, sz=0.25, inhib=4, letter_size=3, ...) {
  # x shoulf be 'inhib_data' from script
  
  break_setter = function(n = break_n) {
    function(lims) {pretty(x = as.numeric(lims), n = n)}
  }
  
  x$log.rel.Ki <- log10(x$rel.Ki)
  corr_df_MW <- ddply(x, c("location", "AMC.substrate"), function(x) corr_stats(x, "MW", "log.rel.Ki"))
  corr_df_MW$p.text <- p_val_labeller(corr_df_MW$pval)
  corr_df_MW$r.text <- paste("r^2==", signif(corr_df_MW$rsq, 2))
  corr_df_MW$lab <- paste("atop(", corr_df_MW$p.text, ",", corr_df_MW$r.text, ")", sep="")
  
  # Plot Ki/KI,ref vs MW
  p_MW_corr <- ggplot(x, aes(x=MW, y=rel.Ki, label=AA.abbrev)) + 
    geom_smooth(method="lm", colour="black") +
    geom_text(size=letter_size) + 
    xlab("molecular weight") +
    ylab(expression(K[I]/K["I,ref"])) +
    scale_y_log10(limits=c(10^floor(log10(min(x$rel.Ki,na.rm=TRUE))), 10^ceiling(log10(max(x$rel.Ki,na.rm=TRUE))))) +
    #scale_y_log10(breaks = break_setter(n=break_n)) +
    scale_x_continuous(limits=c(60, 185), breaks=c(80, 120, 160)) +
    annotation_logticks(sides="l", 
                        short=unit(0.05, "cm"), mid=unit(0.1, "cm"), long=unit(0.15, "cm"),
                        size=sz) +
    geom_text(data=corr_df_MW, aes(x=Inf, y=Inf, label=lab), parse=T, size=label_size, vjust=vjust, hjust=hjust) +
    facet_wrap(~location + AMC.substrate, scales="free", nrow=1) +
    theme(panel.margin=unit(1, "mm"),
          axis.ticks=element_line(size=sz))
  
  corr_df_DI <- ddply(x, c("location", "AMC.substrate"), function(x) corr_stats(x, "Dauwe.score.2", "log.rel.Ki"))
  corr_df_DI$p.text <- p_val_labeller(corr_df_DI$pval)
  corr_df_DI$r.text <- paste("r^2==", signif(corr_df_DI$rsq, 2))
  corr_df_DI$lab <- paste("atop(", corr_df_DI$p.text, ",", corr_df_DI$r.text, ")", sep="")

  # Plot Ki/Ki,ref vs DI
  p_DI_corr <- ggplot(x, aes(x=Dauwe.score.2, y=rel.Ki, label=AA.abbrev)) + 
    geom_smooth(method="lm", colour="black") + 
    geom_text(size=letter_size) +
    xlab("Dauwe DI loading") +
    ylab(expression(K[I]/K["I,ref"])) +
    annotation_logticks(sides="l",
                        short=unit(0.05, "cm"), mid=unit(0.1, "cm"), long=unit(0.15, "cm"),
                        size=sz) +
    scale_x_continuous(limits=c(-0.17, 0.22)) +
    scale_y_log10(limits=c(10^floor(log10(min(x$rel.Ki,na.rm=TRUE))), 10^ceiling(log10(max(x$rel.Ki,na.rm=TRUE))))) +
    #scale_y_log10(breaks=break_setter(n=break_n)) +
    geom_text(data=corr_df_DI, aes(x=Inf, y=Inf, label=lab), parse=T, size=label_size, vjust=vjust, hjust=hjust) +
    facet_wrap(~location + AMC.substrate, scales="free", nrow=1) +
    theme(panel.margin=unit(1, "mm"),
          axis.ticks=element_line(size=sz))
  
  # Plot Ki/KI, ref vs delta GoR
  p_thermo <- ggplot(x, aes(x=deltaGr, y=rel.Ki, label=AA.abbrev)) + 
    geom_smooth(method="lm", colour="black") + 
    geom_text(size=letter_size) +
    xlab(expression(paste(Delta, G[r], ", kJ ", mol^{-1}))) +
    ylab(expression(K[I]/K["I,ref"])) +
    annotation_logticks(sides="l",
                        short=unit(0.05, "cm"), mid=unit(0.1, "cm"), long=unit(0.15, "cm"),
                        size=sz) +
    #scale_x_continuous(limits=c(-0.17, 0.22)) +
    scale_y_log10(limits=c(10^floor(log10(min(x$rel.Ki,na.rm=TRUE))), 10^ceiling(log10(max(x$rel.Ki,na.rm=TRUE))))) +
    #scale_y_log10(breaks=break_setter(n=break_n)) +
    geom_text(data=corr_df_DI, aes(x=Inf, y=Inf, label=lab), parse=T, size=label_size, vjust=vjust, hjust=hjust) +
    facet_wrap(~location + AMC.substrate, scales="free", nrow=1) +
    theme(panel.margin=unit(1, "mm"),
          axis.ticks=element_line(size=sz))
  #browser()
  return(list(p_MW_corr=p_MW_corr, p_DI_corr=p_DI_corr, p_thermo=p_thermo))
}