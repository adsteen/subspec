# # CHeck correlation between delta G and Ki/KI,ref
# 
# source("R/write_paper.R")
# prev_paper <- write_paper(print_plots=FALSE, save_plots=FALSE)
# inhib_data <- prev_paper$inhib_data
# 
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
#   xlab(expression(paste(Delta, italic(G[r])), ", kJ ", mol^{-1})) +
#   ylab("DI loading")
# 
# # # Make plot - for L& O version that wasn't published
# # tiff("../L_and_O_submission/figures/DI_MW_deltaG_correlations.tiff", height=2.5, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# # grid.arrange(p_D_MW, p_MW_deltaG, p_deltaG_DI, nrow=1)
# # dev.off()
# 
# 
# 
# # ggplot(inhib_data, aes(x=rank(deltaGr), y=rank(rel.Ki))) + 
# #   geom_point() + 
# #   geom_smooth(method="lm")
# 
# summary(lm(rank(rel.Ki)~rank(deltaGr), data=inhib_data))
# cor.test(x=inhib_data$deltaGr, y=inhib_data$rel.Ki, method="spearman")
# 
# affinity_plots <- plot_affinity_corr(inhib_data, n=3)
# # print(affinity_plots$p_thermo)
# 
# affinity_plots <- llply(affinity_plots, function(x) x+theme(strip.text=element_text(size=7)))
# 
# 
# if(print_plots) {
#   grid.arrange(affinity_plots[[1]], affinity_plots[[2]], affinity_plots[[3]], nrow=3)
# }
# 
# if(save_plots) {
#   tiff("Fig7.tiff", height=6, width=7.25, units="in", res=300, compression="lzw", type="cairo")
#   grid.arrange(affinity_plots[[1]], affinity_plots[[2]], affinity_plots[[3]], nrow=3)
#   dev.off()
# }
# 
# 
# # id_short <- inhib_data[ , c("AA.abbrev", "rel.Ki", "MW", "deltaGr", "Dauwe.score.2")]
# # idm <- melt(id_short, id.vars=c("AA.abbrev", "rel.Ki"))
# 
# # affinity_plot <- function(d, xvar, sz=2.5, x_lab) {
# #   ggplot(d, aes_string(x=xvar, y="rel.Ki", label="AA.abbrev")) + 
# #     geom_text(size=sz) + 
# #     geom_smooth(method="lm", se=TRUE, colour="black") +
# #     xlab(x_lab) + 
# #     ylab(expression(K[I]/K["I,ref"])) 
# # }
# 
# # ap1 <- affinity_plot(inhib_data, xvar="Dauwe.score.2", x_lab="DI loading")
# # ap2 <- affinity_plot(inhib_data, xvar="MW", x_lab="MW")
# # ap3 <- affinity_plot(inhib_data, xvar="deltaGr", x_lab=expression(paste(Delta,italic(G[r]))))
# # 
# # tiff("plots/affinity_plots_grouped.tiff", height=2, width=7.25, units="in", res=300, compression="lzw", type="cairo")
# # grid.arrange(ap1, ap2, ap3, nrow=1)
# # dev.off()
# 
# # Correlation tests of Ki/Ki,ref vs DI loading, MW and free energy of oxidation
# cor.test(x=inhib_data$rel.Ki, y=inhib_data$Dauwe.score.2, method="spearman")
# cor.test(x=inhib_data$rel.Ki, y=inhib_data$MW, method="spearman")
# cor.test(x=inhib_data$rel.Ki, y=inhib_data$deltaGr, method="spearman")
# 
# summary(lm(Dauwe.score.2 ~ rel.Ki, data=inhib_data))
# summary(lm(MW ~ rel.Ki, data=inhib_data))
# summary(lm(deltaGr ~ rel.Ki, data=inhib_data))
# 
