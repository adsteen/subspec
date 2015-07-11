##' Makes Fig 1, a conceptual figure about substrate specificity
##' 
##' @description Returns a conceptual figure
##' @return a ggplot2 object
##' @export


conceptual_fig <- function() {
  
  # Define x values and positions for peaks; these are universal
  x=seq(0,20,length=2000)
  mu1 <- 5
  mu2 <- 10
  mu3 <- 15
  
  # Define sigmas (SD) for the 'specific' plot
  sig1_spec <- 1.5
  sig2_spec <- 1.2
  sig3_spec <- 1.5
  
  # Create upper values for the specific enzyme ribbons
  spec1 <- gauss(x, mu=mu1, sigma=sig1_spec)
  spec2 <- gauss(x, mu=mu2, sigma=sig2_spec)
  spec3 <- gauss(x, mu=mu3, sigma=sig3_spec)
  
  # Combine these upper values into data frame, melt to long format
  df_spec <- data.frame(x=x, y1=spec1, y2=spec2, y3=spec3)
  dfm_spec <- melt(df_spec, id.vars="x")
  
  # Define parameters for 'substrate' plot (lines)
  x_spec <- c(3, 3, 14, 14, 11, 11)
  y_spec <- c(0,gauss(3, mu1, sig1_spec), 0, gauss(14, mu3, sig3_spec), 0, gauss(11, mu2, sig2_spec))
  l <- gl(3, 2, 6, labels=levels(dfm_spec$variable))
  subs_spec <- data.frame(xP=x_spec, yP=y_spec, lP=l)
  
  ####
  ## Make plot with non-specific enzymes
  ####
  
  sig1_nons <- 2
  sig2_nons <- 8
  sig3_nons <- 3
  
  nons1 <- 0.1*gauss(x, mu=mu1, sigma=sig1_nons)
  nons2 <- gauss(x, mu=mu2, sigma=sig2_nons)
  nons3 <- 0.1*gauss(x, mu=mu3, sigma=sig3_nons)
  
  df_nons <- data.frame(x=x, y1=nons1, y2=nons2, y3=nons3)
  dfm_nons <- melt(df_nons, id.vars="x")
  
  x_nons <- c(3, 3, 14, 14, 11, 11)
  y_nons <- c(0,gauss(3, mu2, sig2_nons), 0, gauss(14, mu2, sig2_nons), 0, gauss(11, mu2, sig2_nons))
  l <- gl(3, 2, 6, labels=levels(dfm_nons$variable))
  subs_nons <- data.frame(xP=x_nons, yP=y_nons, lP=l)
  
  ############
  # Build a single data frame, from which to make the whole plot and add a single label
  ############
  
  dfm_all <- ldply(list(spec=dfm_spec, nons=dfm_nons), identity)
  subs_all <- ldply(list(spec=subs_spec, nons=subs_nons), identity)
  
  dfm_all$.id <- factor(dfm_all$.id, levels=c("spec", "nons"), 
                        labels=c("specific", "unspecific"), ordered=TRUE)
  subs_all$.id <- factor(subs_all$.id, levels=c("spec", "nons"), 
                        labels=c("specific", "unspecific"), ordered=TRUE)
  
  grayscale <- c("gray10", "gray50", "gray75")
  colorscale <- c("black", "black", "black")
  
  ###########
  # build the plot
  ###########
  p <- ggplot() +
    geom_ribbon(data=dfm_all, aes(x=x, ymax=value, ymin=0, fill=variable), alpha=0.4) +
    scale_fill_discrete(name="", breaks=c("y1", "y2", "y3"), labels=c("enzyme 1", "enzyme 2", "enzyme 3")) +
    #scale_fill_manual(name="", values=grayscale, breaks=c("y1", "y2", "y3"), labels=c("enzyme 1", "enzyme 2", "enzyme 3")) +
    #geom_line(data=subs_all, aes(x=xP, y=yP, colour=lP), size=0.75) +
    geom_line(data=subs_all, aes(x=xP, xend=xP, y=yP, colour=lP), size=0.75) +
    geom_point(data=subs_all, aes(x=xP, y=yP, colour=lP), size=1) +
    scale_color_brewer(type="qual", palette="Set1", 
                        name="fluorogenic\nsubstrate proxies",
                        labels=c("proxy A", "proxy B", "proxy C"),
                        guide=FALSE) +
    scale_x_continuous("range of substrates", breaks=c(3, 11, 14), labels=c("", "", "")) + 
    ylab(expression(V[max])) +
    #scale_y_continuous(expression(V[max]), breaks=NULL, limits=c(0,max(dfm_all_norm$value))) +
    facet_wrap(~.id, scales="free_y") +
    theme(text=element_text(size=6),
          legend.position="top",
          panel.border=element_blank(),
          strip.background=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
         # legend.margin=unit(1, "mm"),
          legend.key.size=unit(0.1, "in"),
          legend.key=element_rect(colour="white")) # These seem to have no effect. Don't know why
  
  p
}
