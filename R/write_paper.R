##' Master function to replicate all analyses in the manuscript
##' 
##' @param path Where to save plot files
##' @param print_plots Whether to print plots to the default graphics device (usually your screen)
##' Note that this may have no effect in command-line R.
##' @param save_plots Whether to save plots.
##' @import ggplot2
##' @import plyr
##' @import reshape2
##' @import lubridate
##' @import gridExtra
##' @import grid
##' @export
##' @return NA. The point of this function is to create the plots & calculate the statistics (which are printed to the console)

write_paper <- function(path="", print_plots=TRUE, save_plots=FALSE) {
  # This function re-creates all data analysis used in Steen, Vazin & Wilhelm (submitted)
  #  "Substrate specificity of aquatic extracellular peptidases"
  
  # Useful when running from inside the function
  print_plots=TRUE
  save_plots=FALSE
  
  
  # Used in debugging
  path <- ""
  
  # Load required packages
  require(ggplot2)
  require(plyr)
  require(reshape2)
  require(lubridate)
  require(gridExtra)
  require(grid)
  

  source("subspec/R/calc_inhib_slope.R")
  source("subspec/R/calc_pNA_k.R")
  source("subspec/R/conceptual_fig.R")
  source("subspec/R/corr_pNA_conc.R")
  source("subspec/R/corr_stats.R")
  source("subspec/R/correct_for_pNA_quench.R")
  source("subspec/R/find_homologous.R")
  source("subspec/R/gauss.R")
  source("subspec/R/make_sat_curves.R")
  source("subspec/R/plot_affinity_corr.R")
  source("subspec/R/plot_all_inv_v0.R")
  source("subspec/R/plot_inv_v0.R")
  source("subspec/R/plot_KI_norm_vs_Xaa_pNA.R")
  source("subspec/R/read_inhib_data.R")
  source("subspec/R/sim_multiple_enz.R")
  source("subspec/R/slope_and_SE.R")
  source("subspec/R/p_val_labeller.R")
  #source("R/write_paper.R")
  
  
  # Set some constants relevant to plotting (relevant when this was a script)
  #save_plots <- FALSE
  #print_plots <- TRUE #This is implemented spottily, in a few of the plot functions but not in others
  
  # Set graphical theme
  theme_set(theme_bw() + theme(#text = element_text(family="Times"),
                               panel.grid.major=element_blank(),
                               panel.grid.minor=element_blank(),
                               text=element_text(size=8),
                               strip.background=element_rect(fill="white"),
                               panel.border=element_rect(colour="black"),
                               strip.background=element_rect(colour="black"),
                               panel.margin=unit(1, "mm")))
  
  # For posters
  #theme_set(theme_get() + theme(text=element_text(size=16), strip.text=element_text(size=12)))
  
  myDPI <- 900
  singleColumn <- 3.19 #inches
  doubleColumn <- 6.65 
  #column_and_a_half <- 5.25
  
  letter_size <- 3 # For Fig 4, where AAs are represented as letters using geom_text
  
  #############
  # Fig 1: Conceptual figure
  #############
  conceptual_plot <- conceptual_fig() 
  if (print_plots) {
    print(conceptual_plot)
  }
  if (save_plots) {
    ggsave(paste(path, "plots/fig1.tiff", sep=""), height=2.5, width=singleColumn, unit="in", dpi=myDPI, compression="lzw", type="cairo") # Must be shrunk in illustrator
  }

  ############
  # Fig 2: Enzyme kinetic simlation
  ############
  # Returns 2 figs in a list
  fig2_list <- sim_multiple_enz(print_plot=print_plots, save_plot=TRUE, height=1.5, width=singleColumn, fn="subspec/plots/fig2.tiff")
  
  ##########
  # Fig 3: Saturation curves 
  ##########
  fig3 <- make_sat_curves(print_plot=FALSE, save_plot=FALSE)
  if(print_plots) {
    print(fig3)
  }
  if(save_plots) {
    ggsave("subspec/plots/fig3.tiff", fig3, height=4, width=singleColumn, dpi=900, compression="lzw", type="cairo")
  }
  
  ##########
  # Data processing for inhibition experiment
  ##########
  ###  Data setup: 
  
  # Read the data by hand, put it into a list
  raw_data_list <- read_inhib_data()
  
  # Correct the pNA conc
  raw_data_list <- llply(raw_data_list, corr_pNA_conc)
  
  ### Correct raw fluorescence values for pNA quenching
  # Determine exponential correction factor for pNA quenching
  
  k <- calc_pNA_k()
  
  # Perform the correction & store results in correctedRawList
  
  corrected_raw_list <- llply(raw_data_list, correct_for_pNA_quench, k=k)
  
  ### Calculate rates of fluorescence production & calibrate 
  
  # Actually calculate the slope and standard error
  all_slopes <- ldply(corrected_raw_list, function(x) ddply(x, c("pNA.subs", "conc.pNA"), slope_and_SE))
  
  # Put a place and substrate label on each row of all_slopes
  all_slopes$AMC.substrate <- paste(substr(all_slopes$.id, start=1, stop=3), "-AMC", sep="")
  all_slopes$location <- substr(all_slopes$.id, start=4, stop=5)
  all_slopes$location[all_slopes$location=="TN"] <- "TN River"
  all_slopes$location[all_slopes$location=="NC"] <- "Bogue Sound"
  
  # Read in all calibration data
  #TN_calib <- read.csv("ProAMCcalcurve.csv")
  TN_slope <- coefficients(lm(Fl ~ Conc.Pro, data=TN_calib))[2] * 1000 #Confirm that it is correct to multiply by 1000
  #IMS_calib <- read.csv("2013-03-26 calib curve.csv")
  IMS_slope <- coefficients(lm(fl~conc.AMC.nM, data=IMS_calib))[2]
  
  # Perform the calibration
  all_slopes$nM.per.hr <- NA
  all_slopes$se.nM.per.hr <- NA
  all_slopes$nM.per.hr[all_slopes$location == "TN River"] <- all_slopes$slopePerHr[all_slopes$location == "TN River"] / TN_slope
  all_slopes$se.nM.per.hr[all_slopes$location == "TN River"] <- all_slopes$slopeSE[all_slopes$location == "TN River"] / TN_slope
  all_slopes$nM.per.hr[all_slopes$location == "Bogue Sound"] <- all_slopes$slopePerHr[all_slopes$location == "Bogue Sound"] / IMS_slope
  all_slopes$se.nM.per.hr[all_slopes$location == "Bogue Sound"] <- all_slopes$slopeSE[all_slopes$location == "Bogue Sound"] / IMS_slope
  
  # Drop the killed controls
  all_slopes <- all_slopes[-which(all_slopes$pNA.subs=="killed"), ]
  
  # Refactor the pNA substance to contain AA names
  all_slopes$pNA.subs <- factor(all_slopes$pNA.subs, levels=1:12, 
                                labels=paste(c("Ala", "Arg", "Glu", "Gly", "Ile", "Leu", "Lys", 
                                               "Met", "Phe", "Pro", "Tyr", "Val"), "-pNA", sep=""), ordered=TRUE)
  
  # Change AMC substrate names
  all_slopes$AMC.substrate <- as.character(factor(all_slopes$AMC.substrate, 
                                           levels=unique(all_slopes$AMC.substrate),
                                           labels=c("Arg-AMC", "Leu-AMC", "Pro-AMC")))
  
  # Drop some obviously bad data points
  # (note: some of these are probably transcription errors, some are typos and some are shot noise from the instrument whcih sometimes reports unreasonably high values)
  all_slopes <- all_slopes[-which(all_slopes$AMC.substrate=="Pro-AMC" &
                        all_slopes$pNA.subs=="Arg-pNA" &
                        all_slopes$conc.pNA==unique(all_slopes$conc.pNA)[9]), ]
  all_slopes <- all_slopes[-which(all_slopes$AMC.substrate=="Pro-AMC" &
                        all_slopes$pNA.subs=="Gly-pNA" &
                        all_slopes$conc.pNA==unique(all_slopes$conc.pNA)[7]), ]
  all_slopes <- all_slopes[-which(all_slopes$AMC.substrate=="Leu-AMC" &
                        all_slopes$location=="Bogue Sound" &
                        all_slopes$pNA.subs=="Val-pNA" &
                        all_slopes$conc.pNA==unique(all_slopes$conc.pNA)[5]), ]
  
  # Looks like I don't use the function plot_all_inv_v0
  p_all_inv_v0 <- ggplot(all_slopes, aes(x=conc.pNA, y=1/nM.per.hr)) + 
    geom_point(size=0.75) +
    geom_errorbar(aes(ymax=1/(nM.per.hr-se.nM.per.hr), ymin=1/(nM.per.hr+se.nM.per.hr))) +
    geom_smooth(method="lm", colour="black") +
    #xlab(expression(paste("[I], ", , mu, "M ", hr^{-1}))) +
    xlab(expression(paste("[I], ", , mu, "M "))) +
    ylab(expression(paste(1/v[0], ", h ", nM^{-1}))) +
    facet_grid(AMC.substrate + location ~ pNA.subs, scales="free") + 
    theme(text=element_text(size=9),
          axis.text.x=element_text(angle=-45, hjust=0)) 
  
  if (print_plots) {
    print(p_all_inv_v0)
  }
  if(save_plots) {
    ggsave("plots/fig4.tiff", p_all_inv_v0, height=5, width=doubleColumn, units="in", dpi=myDPI, type="cairo")
  }
  
#   #######
#   # Example plot of 1/v0 vs [I] (arg at Bogue Sound)
#   #######
#   
#   #p_inv_v0_arg_Bogue <- plot_inv_v0(subset(all_slopes, location=="Bogue Sound" & AMC.substrate=="arg-AMC"))
#   #p_inv_v0 <- plot_inv_v0(all_slopes, spacing=0.1)
#   
#   
# #   if (print_plots) {
# #     print(p_inv_v0_arg_Bogue)
# #   }
# #   if(save_plots) {
# #     ggsave(paste(path, "inv_v0_arg_Bogue.png", sep=""), p_inv_v0_arg_Bogue, height=7, width=singleColumn, units="in", dpi=myDPI, type="cairo")
# #   }
#   
#   ##########
#   # Calculate relative kI for each combination of Yaa-AMC and Xaa-pNA]
#   ###########
#   
#   # Calculate the "inhibition slopes" (i.e., slopes of 1/v0 vs [I])
#   all_inhib_slopes <- ddply(all_slopes, c("location", "AMC.substrate", "pNA.subs"), calc_inhib_slope)
#   
#   
#   # label the "homologous" substrates
#   all_inhib_slopes <- find_homologous(all_inhib_slopes)
#   
#   # Calculate relative Ki (i.e., Ki / Ki,ref)
#   all_rel_Ki <- ddply(all_inhib_slopes, c("location", "AMC.substrate"), transform, 
#                       rel.Ki=inhib.slope[is.homologous==TRUE] / inhib.slope)
#   all_rel_Ki$rel.Ki.SE <- all_rel_Ki$inhib.slope.SE * all_rel_Ki$rel.Ki / all_rel_Ki$inhib.slope 
#   
#   #######
#   # Determine the quantification limit
#   #   The theory here is that the slope SE is a reasonable estiamte of the minimum slope one could measure
#   #   The quantification limit is therefore the mean of the slope[homologous] / slope.SE
#   #######
#   quant_limit_data <- ddply(all_rel_Ki, c("location", "AMC.substrate"), transform, quant.limit = inhib.slope[is.homologous==TRUE] / inhib.slope.SE)
#   quant_limit <- mean(quant_limit_data$quant.limit)
#   
#   # Filter out "bad" relKi values: where inhibSlope < 0 or insn't significant (for now, that just means negative)
#   #all_rel_Ki$rel.Ki[all_rel_Ki$inhib.slope < 0] <- max(all_rel_Ki$rel.Ki) 
#   all_rel_Ki$rel.Ki[all_rel_Ki$inhib.slope < 0] <- quant_limit # This is arbitrary
#   print(paste("'Bad' KI values are set to ", max(all_rel_Ki$rel.Ki), sep=""))
#   
#   # Set errorbars for Ki plots
#   all_rel_Ki$low.err <- all_rel_Ki$rel.Ki - all_rel_Ki$rel.Ki.SE
#   all_rel_Ki$hi.err <- all_rel_Ki$rel.Ki + all_rel_Ki$rel.Ki.SE
#   all_rel_Ki$low.err[all_rel_Ki$low.err <= 0] <- 0.01
#   
#   
#   #########
#   # Fig 4: Dotplot of all normalized affinities
#   #########
#   
#   #invisible(plot_KI_norm_vs_Xaa_pNA(all_rel_Ki, vertical=FALSE, print_plot = print_plots, save_plot=save_plots, height=5, width=singleColumn))
#   p_KI_vs_Xaa_pNA <- plot_KI_norm_vs_Xaa_pNA(all_rel_Ki, vertical=FALSE, print_plot = print_plots, save_plot=save_plots, height=5, width=singleColumn)
# 
#   
#   if (print_plots) {
#     print(p_KI_vs_Xaa_pNA)
#   }
#   
#   if (save_plots) {
#     #if (is.na(fn)) {
#       fn <- paste(path, "KI_vs_Xaa_pNA.png", sep="")
#     #}
#     ggsave(fn, p_KI_vs_Xaa_pNA, height=5, width=column_and_a_half, units="in", dpi=300, type="cairo")
#   }
#   
#   
#   ### Make a table of the actual relative Ki data for the supplemental
#   rel_KI_table <- all_rel_Ki[ , c("location", "AMC.substrate", "rel.Ki", "rel.Ki.SE")]
#   rel_KI_table$rel.Ki <- signif(rel_KI_table$rel.Ki, digits=2)
#   rel_KI_table$rel.Ki.SE <- signif(rel_KI_table$rel.Ki.SE, digits=2)
#   
#   if(save_plots) {
#     pdf(paste(path, "tableS1_all_relative_kI_page1.pdf", sep=""), width=8.5, height=11, onefile=FALSE)
#     grid.table(rel_KI_table[1:36,])
#     dev.off()
#     
#     pdf(paste(path, "tableS1_all_relative_kI_page2.pdf", sep=""), width=8.5, height=11, onefile=FALSE)
#     grid.table(rel_KI_table[37:nrow(all_rel_Ki), ])
#     dev.off()
#   }
#   
#   
#   
#   ########
#   # Compare to amino acid properties
#   ########
#   
#   # Merge amino acid characteristics into table of relative Ki values
#   AA_char <- read.csv("AA characteristics.csv")
#   all_rel_Ki$AA <- substr(all_rel_Ki$pNA.subs, start=1, stop=3)
#   inhib_data <- merge(all_rel_Ki, AA_char, by="AA")
#   
#   
#   # Make correlation plot for affinity vs MW and affinity vs DI
#   add_corr_plot_list <- plot_affinity_corr(inhib_data, print_plots=TRUE, save_plots=FALSE, vjust=1.3, hjust=1.5, letter_size=letter_size) #Saving the plot from the function does not work for some unknown reason
#   
#   #####
#   # NEed to run delta_G_affinity.R _after_ this entire script
#   
#   #####
#   
#   ####
#   # Figure 6:
#   ####
#   
#   # Boolean for AAs that are in my dataset
#   AA_char$in.my.dataset <- AA_char$AA %in% inhib_data$AA
#   
#   MW_DI_corr_data <- corr_stats(AA_char, "MW", "Dauwe.score.2") # Full data set
#   MW_DI_subs_corr_data <- corr_stats(subset(AA_char, in.my.dataset==TRUE & AA.abbrev != "R"), "MW", "Dauwe.score.2")
#   
#   
#   ###############
#   #
#   # Stats: relKi, DI loadings, and MW
#   #
#   ###############
#   
#   # Test correlation between rel.Ki and hydrophobicity
#   m_Ki_by_hyd <- lm(rel.Ki ~ hydrophobicity, data=inhib_data)
#   #plot(m_Ki_by_hyd) #close-enough to being homoskedastic, I think 
#   print(summary(m_Ki_by_hyd))
#   
#   # Test correlation between rel.Ki and hydrophobicity, rank-transformed
#   print(cor.test(x=inhib_data$hydrophobicity, y=inhib_data$rel.Ki, data=inhib_data, method="spearman"))
#   
#   # Test whether polarity affects rel.Ki
#   m_pol <- lm(rel.Ki ~ polarity, data=inhib_data)
#   #plot(m_pol) #lm assumptions are more-or-less supported, I think 
#   summary(m_pol)
#   
#   # Test Spearman correlation between relKi and Dauwe score (full dataset)
#   cor.test(x=inhib_data$Dauwe.score.2, y=inhib_data$rel.Ki, method="spearman")
#   nrow(inhib_data[!is.na(inhib_data$Dauwe.score.2), ])
#   
#   # Test linear correlation between log10(relKi) and Dauwe score
#   m1 <- lm(log10(rel.Ki) ~ Dauwe.score.2, data=inhib_data)
#   print(summary(m1))
#   
#   # Test Spearman correlation between relKi and MW
#   cor.test(x=inhib_data$MW, y=inhib_data$rel.Ki, method="spearman")
#   nrow(inhib_data[!is.na(inhib_data$MW), ])
#   
#   # Test linear correlation between log10(relKi) and MW
#   m2 <- lm(rel.Ki ~ MW, data=inhib_data)
#   print(summary(m2))
#   
#   # Test correlation between DI and MW  -have to use AAchar, since all_rel_Ki has the same datapoints repeated 5x
#   m3 <- lm(Dauwe.score.2 ~ MW, data=subset(AA_char, in.my.dataset==TRUE))
#   print(summary(m3)) # p = 0.16
#   
#   m3v2 <- lm(Dauwe.score.2 ~ MW, data=AA_char)
#   print(summary(m3v2))
#   
#   #######
#   # Create figure 4, saturation curves
#   #######
# make_sat_curves(print_plot=print_plots, save_plot=save_plots)
# 
# 
# ###########
# # FROM SCRIPT delta_g_affinity.R
# ###########
# AA_thermo <- read.csv("data/AA_thermo.csv")
# 
# inhib_data <- merge(inhib_data, AA_thermo)
# 
# ddply(inhib_data, c("location", "AMC.substrate"), corr_stats, xvar="rel.Ki", yvar="deltaGr")
# 
# summary(lm(rel.Ki~deltaGr, data=inhib_data))
# summary(lm(log10(rel.Ki)~deltaGr, data=inhib_data))
# cor.test(x=inhib_data$deltaGr, y=inhib_data$rel.Ki, method="spearman")
# nrow(inhib_data)
# 
# d_ply(inhib_data, c("AMC.substrate", "location"), function(x) print(summary(lm(rel.Ki~deltaGr, data=x))))
# # They're generally not significant on their own, but the slope is always positive
# 
# 
# cor.test(x=inhib_data$Dauwe.score.2, y=inhib_data$rel.Ki, method="spearman")
# 
# # Check correlations of DI and MW, deltaG and MW
# cor.test(x=inhib_data$MW, y=inhib_data$deltaGr, method="spearman")
# cor.test(x=inhib_data$MW, y=inhib_data$deltaGr, method="spearman")
# cor.test(x=inhib_data$deltaGr, y=inhib_data$Dauwe.score.2, method="spearman")
# p_D_MW <- ggplot(inhib_data, aes(x=Dauwe.score.2, y=MW, label=AA.abbrev)) + 
#   geom_text() + 
#   geom_smooth(method="lm", colour="black") +
#   xlab("Dauwe loading") +
#   ylab("molecular weight")
# if(print_plots) {
#   print(p_D_MW)
# }
# 
# p_MW_deltaG <- ggplot(inhib_data, aes(x=MW, y=deltaGr, label=AA.abbrev)) + 
#   geom_text() + 
#   geom_smooth(method="lm", colour="black") +
#   xlab("Molecular weight") +
#   ylab(expression(paste(Delta, italic(G[r])))) 
# 
# p_deltaG_DI <- ggplot(inhib_data, aes(x=deltaGr, y=Dauwe.score.2, label=AA.abbrev)) + 
#   geom_text() + 
#   geom_smooth(method="lm", colour="black") +
#   xlab(expression(paste(Delta, italic(G[r])))) +
#   ylab("DI loading")
# 
# # #tiff("plots/probably_wont_use_DI_MW_deltaG_correlations.tiff", height=2.5, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# # tiff("../L_and_O_submission/figures/DI_MW_deltaG_correlations.tiff", height=2.5, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# # grid.arrange(p_D_MW, p_MW_deltaG, p_deltaG_DI, nrow=1)
# # dev.off()
# 
# 
# 
# # p_dG <- ggplot(inhib_data, aes(x=deltaGr, y=rel.Ki)) + 
# #   geom_text(aes(label=AA.abbrev)) + 
# #   geom_smooth(method="lm", colour="black") +
# #   xlab(expression(paste(Delta, italic(G[r])))) +
# #   ylab(expression(K[I]/K["I,ref"])) 
# # print(p_dG)
# # ggsave("../../../Funding/2013/Steen CDEBI Oct 2013/administration/KI_DI.eps", p_dG + theme(text=element_text(size=12)), height=3, width=4, units="in", dpi=300)
# 
# # p_MW <- ggplot(inhib_data, aes(x=MW, y=rel.Ki)) + 
# #   geom_text(aes(label=AA.abbrev)) + 
# #   geom_smooth(method="lm", colour="black") + 
# #   xlab("MW") +
# #   ylab(expression(K[I]/K["I,ref"])) 
# # print(p_MW)
# 
# # p_DI <- ggplot(inhib_data, aes(x=Dauwe.score.2, y=rel.Ki)) + 
# #   geom_text(aes(label=AA.abbrev)) + 
# #   geom_smooth(method="lm", colour="black") +
# #   xlab("DI") +
# #   ylab(expression(K[I]/K["I,ref"])) 
# # print(p_DI)
# 
# # tiff("../AME_submission/plots/relKI.tiff", height=3, width=6.65, units="in", res=300)
# # grid.arrange(p_MW, p_DI, p_dG, nrow=1)
# # dev.off()
# 
# # ### Save versions for talk
# # p_dG_talk <- p_dG + 
# #   theme(text=element_text(size=20))
# # ggsave("plots/Ki_v_dG_for_talk.png", height=5, width=6, units="in", dpi=300, type="cairo")
# # 
# # p_MW_talk <- p_MW + 
# #   theme(text=element_text(size=20))
# # ggsave("plots/Ki_v_MW_for_talk.png", height=5, width=6, units="in", dpi=300, type="cairo")
# # 
# # p_DI_talk <- p_DI + 
# #   theme(text=element_text(size=20))
# # ggsave("plots/Ki_v_DI_for_talk.png", height=5, width=6, units="in", dpi=300, type="cairo")
# 
# ###########
# # Test significance of correlations between KI and MW, deltaG and DI
# ###########
# 
# summary(lm(rank(rel.Ki)~rank(deltaGr), data=inhib_data))
# cor.test(x=inhib_data$deltaGr, y=inhib_data$rel.Ki, method="spearman")
# 
# ##########
# # Make hte "affinity plots"
# ##########
# affinity_plots <- plot_affinity_corr(inhib_data, n=3)
# print(affinity_plots$p_thermo)
# 
# affinity_plots <- llply(affinity_plots, function(x) x+theme(strip.text=element_text(size=7)))
# #######
# # Fig 7
# #######
# if(save_plots) {
#   png("~/Desktop/new_affinity_plot.png", height=6, width=7.25, units="in", res=300, type="cairo")
#   #png("new_affinity_plot.png", height=6, width=7.25, units="in", res=300, type="cairo")
#   grid.arrange(affinity_plots[[1]], affinity_plots[[2]], affinity_plots[[3]], nrow=3)
#   dev.off()
# }
# #vp1 <- viewport(height=1/3, width=1, x=1/2, y=5/6)
# #tiff("plots/new_affinity_plot.tiff", height=6, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# tiff("new_affinity_plot.tiff", height=6, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# grid.arrange(affinity_plots[[1]], affinity_plots[[2]], affinity_plots[[3]], nrow=3)
# dev.off()
# 
# id_short <- inhib_data[ , c("AA.abbrev", "rel.Ki", "MW", "deltaGr", "Dauwe.score.2")]
# idm <- melt(id_short, id.vars=c("AA.abbrev", "rel.Ki"))
# 
# affinity_plot <- function(d, xvar, sz=2.5, x_lab) {
#   ggplot(d, aes_string(x=xvar, y="rel.Ki", label="AA.abbrev")) + 
#     geom_text(size=sz) + 
#     geom_smooth(method="lm", se=TRUE, colour="black") +
#     xlab(x_lab) + 
#     ylab(expression(K[I]/K["I,ref"])) 
# }
# 
# ap1 <- affinity_plot(inhib_data, xvar="Dauwe.score.2", x_lab="DI loading")
# ap2 <- affinity_plot(inhib_data, xvar="MW", x_lab="MW")
# ap3 <- affinity_plot(inhib_data, xvar="deltaGr", x_lab=expression(paste(Delta,italic(G[r]))))
# 
# tiff("plots/affinity_plots_grouped.tiff", height=2, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# grid.arrange(ap1, ap2, ap3, nrow=1)
# dev.off()
# 
# 
# cor.test(x=inhib_data$rel.Ki, y=inhib_data$Dauwe.score.2, method="spearman")
# cor.test(x=inhib_data$rel.Ki, y=inhib_data$MW, method="spearman")
# cor.test(x=inhib_data$rel.Ki, y=inhib_data$deltaGr, method="spearman")
# 
# summary(lm(Dauwe.score.2 ~ rel.Ki, data=inhib_data))
# summary(lm(MW ~ rel.Ki, data=inhib_data))
# summary(lm(deltaGr ~ rel.Ki, data=inhib_data))
# 
# 
# 
# 
# 
#   
  #return(list(inhib_data=inhib_data))
  NA
 
}