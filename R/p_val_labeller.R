##' Returns nicely formatted p values
##' 
##' @param pval A vector of pvalues
##' @param sigfigs The number of significant figures to which to round each p value
##' @param cutoff The number below which the pvalue will be reported as "p<cutoff"
##' @param plotmath Whether to format for plotmath or plain text. Plotmath uses "p==", which will display onscreen as "p=" if it is interpreted as plotmath
##' @export

p_val_labeller <- function(pval, sigfigs=2, cutoff=0.001, plotmath=TRUE) {
  ##########
  # Round the numeric p values 
  ##########
  p_signif <- rep(NA, length(pval))
  
  # Use for loop b/c vectorized form will give more sigfigs for some elements than others
  for (i in 1:length(pval)) {
    p_signif[i] <- signif(pval[i], sigfigs)
  }
  
  ########
  # Create text labels
  ########
  
  # Plotmath displays '==' as '='
  if (plotmath) {
    p_text <- paste("p==", p_signif, sep="")
  } else {
    p_text <- paste("p=", p_signif, sep="")
  }
  
  # Replace values less than cutoff with "p<cutoff"
  p_text[pval < cutoff] <- paste("p<", cutoff, sep="")
  p_text
}