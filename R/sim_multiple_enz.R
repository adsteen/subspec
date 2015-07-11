sim_multiple_enz <- function(print_plot=FALSE, save_plot=FALSE, fn=NA, height=1.5, width=3.33, myDPI=300) {
  dotsize=1
  
  inhib <- function(I, S=100, Km=100, Ki=100, Vmax=1) {
    (Vmax*S) / ((1 + I/Ki) + S)
  }
  
  I <- 0:360
  
  # Just one enzyme & one inhibitor
  #Idf <- data.frame(I=I, invV0 = 1/inhib(I))
  #ggplot(Idf, aes(x=I, y=invV0)) + geom_line()
  
  # Ten enzymes, each with:
  #    The same Vmax (0.1)
  #    A different Km (uniform random from 10 to 1000)
  #    A different Ki (uniform random from 10 to 1000)
  
  set.seed(2)
  Kms <- runif(10, min=100, max=1000)
  Kis <- runif(10, min=0.1, max=10)
  
  
  varInhibs <- data.frame(
    I = I,
    v01 = inhib(I, Km=Kms[1], Ki=Kis[1]),
    v02 = inhib(I, Km=Kms[2], Ki=Kis[2]),
    v03 = inhib(I, Km=Kms[3], Ki=Kis[3]),
    v04 = inhib(I, Km=Kms[4], Ki=Kis[4]),
    v05 = inhib(I, Km=Kms[5], Ki=Kis[5]),
    v06 = inhib(I, Km=Kms[6], Ki=Kis[6]),
    v07 = inhib(I, Km=Kms[7], Ki=Kis[7]),
    v08 = inhib(I, Km=Kms[8], Ki=Kis[8]),
    v09 = inhib(I, Km=Kms[9], Ki=Kis[9]),
    v10 = inhib(I, Km=Kms[10], Ki=Kis[10])
  )
  
  varInhibs$sum <- rowSums(varInhibs[ , -1])
  
  # Make plot of individual enzyme simulated kinetics
  varInhibsM <- melt(varInhibs, id.vars="I", variable.name="enzyme", value.name="v0")
  varInhibsM$invV0 <- varInhibsM 
  fig1a <- ggplot(varInhibsM[-which(varInhibsM$enzyme=="sum"), ], aes(x=I, y=v0, group=enzyme)) + geom_line() +
    ylab(expression(paste(v[0], ", arb. units"))) +
    xlab(expression(paste("[I], ", mu, "M")))

  # Make a new 'grid' of all the concentartions that I actually sampled
  varInhibsSamp <- varInhibs[varInhibs$I %in% seq(from=0, to=360, by=20), ] / (1 + 0.1 + 0.016 + 0.036)
  
  
  fig1b <- ggplot(varInhibsSamp, aes(x=I, y= 1/sum)) + geom_point(size=dotsize) + 
    geom_smooth(method="lm", colour="black", linetype=1, se=FALSE) +
    ylab(expression(paste(1/v[obs], ", arb. units"))) +
    xlab(expression(paste("[I], ", mu, "M")))  
  #ggtitle("Treating multipleisofunctional enzymes as one\nleads to negligable error")
  
  #theme_set(old_theme)
  
  vp1 <- viewport(x=0.25, y=0.5, height=1, width=0.5)
  vp2 <- viewport(x=0.75, y=0.5, height=1, width=0.5)
  
  if(print_plot) {
    print(fig1a, vp=vp1)
    print(fig1b, vp=vp2)
  }
  if(save_plot) {
    if(is.na(fn)) {
      fn <- paste(path, "kinetic_simulation_plot.tiff", sep="")
    }
    tiff(fn, height=height, width=width, units="in", res=myDPI, type="cairo", compression="lzw")
    print(fig1a, vp=vp1)
    print(fig1b, vp=vp2)
    dev.off()
  }
  
  return(list(fig1a, fig1b))
  
}