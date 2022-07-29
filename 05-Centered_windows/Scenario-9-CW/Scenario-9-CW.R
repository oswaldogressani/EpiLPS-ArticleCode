#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 9)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")
library("xlsx")

S <- 100 # Total number of epidemics considered
seedval <- 1905

#------------------------------------------------------------------------------#
#               Scenario 9 "Wiggly then flat R(t)" SI-MERS (NegBin DGP)        #
#------------------------------------------------------------------------------#

si_mers <- round(EpiEstim::discr_si(k = seq(1,20), mu = 6.8, sigma=4.1),3)
si_mers <- si_mers/sum(si_mers)


simcheck <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                   dist="negbin", overdisp = 50, verbose = TRUE, plotsim = TRUE)


# LPSMAP, LPSMALA and EpiEstim with 7 days window
pdf(file = "Figures/LPSMAP_S9_7d.pdf", width = 15, height = 4.6)
sim9_LPSMAP90_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 50, ci_level = 0.90,
                           seed = seedval, themetype = "gray", epidays = 60,
                           midwindow = TRUE)
dev.off()

sim9_LPSMALA90_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 50, seed = seedval,
                                themetype = "gray", epidays = 60,
                                midwindow = TRUE)

sim9_LPSMAP95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                                K = 40, method = "LPSMAP", dist = "negbin",
                                overdisp = 50, ci_level = 0.95,
                                seed = seedval, themetype = "gray", 
                                epidays = 60, midwindow = TRUE)

sim9_LPSMALA95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                                 K = 40, method = "LPSMALA", chain_length = 3000,
                                 burn = 1000, ci_level = 0.95, dist = "negbin",
                                 overdisp = 50, seed = seedval,
                                 themetype = "gray", epidays = 60, 
                                 midwindow = TRUE)

# EpiEstim with 3 days window
sim9_3d90_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 50, ci_level = 0.90,
                            seed = seedval, themetype = "gray", epidays = 60,
                            slidewindow = 2, midwindow = TRUE)

sim9_3d95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 50, ci_level = 0.95,
                            seed = seedval, themetype = "gray", epidays = 60,
                            slidewindow = 2, midwindow = TRUE)
  

#---- Extracting plots
pdf(file = "Figures/S9_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim9_LPSMAP90_mers$inciplot+ggplot2::ggtitle("Incidence (Scenario 9)")
dev.off()
pdf(file = "Figures/S9_LPSMAP_mers.pdf", width = 6, height = 4.5) 
sim9_LPSMAP90_mers$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories")
dev.off()
pdf(file = "Figures/S9_LPSMALA_mers.pdf", width = 6, height = 4.5) 
sim9_LPSMALA90_mers$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories")
dev.off()
pdf(file = "Figures/S9_EpiEstim7d_mers.pdf", width = 6, height = 4.5) 
sim9_LPSMAP90_mers$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories")
dev.off()
pdf(file = "Figures/S9_EpiEstim3d_mers.pdf", width = 6, height = 4.5) 
sim9_3d90_mers$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories")
dev.off()


png(file = "Figures/Scenario9_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim9_LPSMAP90_mers$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 9)"),
                        sim9_LPSMAP90_mers$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories"),
                        sim9_LPSMALA90_mers$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories"),
                        sim9_LPSMAP90_mers$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories"),
                        sim9_3d90_mers$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for MERS SI
Scenario9_mers <- matrix(0, nrow = 4, ncol = 6)
rownames(Scenario9_mers) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d")
colnames(Scenario9_mers) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario9_mers[1,1:3] <- c(mean(sim9_LPSMAP90_mers$simul_summary[,1]),
                      mean(sim9_LPSMAP90_mers$simul_summary[,3]),
                      mean(sim9_LPSMAP90_mers$simul_summary[,5]))
Scenario9_mers[1,4] <- mean(sim9_LPSMAP95_mers$simul_summary[,5])
Scenario9_mers[1,5] <- mean(sim9_LPSMAP90_mers$ciwidthepilps)
Scenario9_mers[1,6] <- mean(sim9_LPSMAP95_mers$ciwidthepilps)

# LPSMALA
Scenario9_mers[2,1:3] <- c(mean(sim9_LPSMALA90_mers$simul_summary[,1]),
                      mean(sim9_LPSMALA90_mers$simul_summary[,3]),
                      mean(sim9_LPSMALA90_mers$simul_summary[,5]))
Scenario9_mers[2,4] <- mean(sim9_LPSMALA95_mers$simul_summary[,5])
Scenario9_mers[2,5] <- mean(sim9_LPSMALA90_mers$ciwidthepilps)
Scenario9_mers[2,6] <- mean(sim9_LPSMALA95_mers$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario9_mers[3,1:3] <- c(mean(sim9_LPSMAP90_mers$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim9_LPSMAP90_mers$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim9_LPSMAP90_mers$simul_summary[,6],
                           na.rm = TRUE))
Scenario9_mers[3,4] <- mean(sim9_LPSMAP95_mers$simul_summary[,6],
                            na.rm = TRUE)
Scenario9_mers[3,5] <- mean(sim9_LPSMAP90_mers$ciwidthepiestim,
                            na.rm = TRUE)
Scenario9_mers[3,6] <- mean(sim9_LPSMAP95_mers$ciwidthepiestim,
                            na.rm = TRUE)

# EpiEstim (3 days sliding window)
Scenario9_mers[4,1:3] <- c(mean(sim9_3d90_mers$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim9_3d90_mers$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim9_3d90_mers$simul_summary[,6],
                           na.rm = TRUE))
Scenario9_mers[4,4] <- mean(sim9_3d95_mers$simul_summary[,6],
                            na.rm = TRUE)
Scenario9_mers[4,5] <- mean(sim9_3d90_mers$ciwidthepiestim,
                            na.rm = TRUE)
Scenario9_mers[4,6] <- mean(sim9_3d95_mers$ciwidthepiestim,
                            na.rm = TRUE)


# Extract tables
S9mers <- round(Scenario9_mers, 3)
write.table(S9mers, file="S9mersNegBin.txt", row.names = TRUE, col.names = TRUE)
write.xlsx(S9mers, file="S9mersNegBin.xls")
S9mers

# Save environment
save.image(file="Scenario9-CW.RData")




































