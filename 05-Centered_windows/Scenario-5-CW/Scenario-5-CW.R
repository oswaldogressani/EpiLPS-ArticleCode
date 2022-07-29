#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 5)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")
library("xlsx")

S <- 100 # Total number of epidemics considered
seedval <- 1998

#------------------------------------------------------------------------------#
#               Scenario 5 "Constant R(t)" SI-SARS (NegBin DGP)                #
#------------------------------------------------------------------------------#

si_sars <- c(0.001,0.012,0.043,0.078,0.104,0.117,0.116,0.108,0.094,0.078,0.063,
             0.049,0.038,0.028,0.021,0.015,0.011,0.008,0.005,0.004,0.003,0.002,
             0.001,0.001)

# Above si obtained with discr_si function of EpiEstim package
# si_alter <- round(EpiEstim::discr_si(k = seq(1,24), mu = 8.4, sigma=3.8),3)
# si_alter
# sum(si_alter)
simcheck <- episim(serial_interval = si_sars, endepi = 40, Rpattern = 1,
                   Rconst = 1.3, dist="negbin", overdisp = 5, 
                   verbose = TRUE, plotsim = TRUE)


# LPSMAP, LPSMALA and EpiEstim with 7 days window
pdf(file = "Figures/LPSMAP_S5_7d.pdf", width = 15, height = 4.6)
sim5_LPSMAP90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 5, ci_level = 0.90, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40,
                           midwindow = TRUE)
dev.off()

sim5_LPSMALA90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 5, Rconst = 1.3, seed = seedval,
                                themetype = "gray", epidays = 40,
                                midwindow = TRUE)

sim5_LPSMAP95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                                K = 40, method = "LPSMAP", dist = "negbin",
                                overdisp = 5, ci_level = 0.95, Rconst = 1.3,
                                seed = seedval, themetype = "gray", epidays = 40,
                                midwindow = TRUE)

sim5_LPSMALA95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.95, dist = "negbin",
                                overdisp = 5, Rconst = 1.3, seed = seedval,
                                themetype = "gray", epidays = 40,
                                midwindow = TRUE)

# EpiEstim with 3 days window
sim5_3d90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.90, Rconst = 1.3,
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 2, midwindow = TRUE)

sim5_3d95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.95, Rconst = 1.3,
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 2, midwindow = TRUE)
  

#---- Extracting plots
pdf(file = "Figures/S5_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim5_LPSMAP90_sars$inciplot+ggplot2::ggtitle("Incidence (Scenario 5)")
dev.off()
pdf(file = "Figures/S5_LPSMAP_sars.pdf", width = 6, height = 4.5) 
sim5_LPSMAP90_sars$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories")
dev.off()
pdf(file = "Figures/S5_LPSMALA_sars.pdf", width = 6, height = 4.5) 
sim5_LPSMALA90_sars$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories")
dev.off()
pdf(file = "Figures/S5_EpiEstim7d_sars.pdf", width = 6, height = 4.5) 
sim5_LPSMAP90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories")
dev.off()
pdf(file = "Figures/S5_EpiEstim3d_sars.pdf", width = 6, height = 4.5) 
sim5_3d90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories")
dev.off()


png(file = "Figures/Scenario5_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim5_LPSMAP90_sars$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 5)"),
                        sim5_LPSMAP90_sars$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories"),
                        sim5_LPSMALA90_sars$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories"),
                        sim5_LPSMAP90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories"),
                        sim5_3d90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for SARS SI
Scenario5_sars <- matrix(0, nrow = 4, ncol = 6)
rownames(Scenario5_sars) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d")
colnames(Scenario5_sars) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario5_sars[1,1:3] <- c(mean(sim5_LPSMAP90_sars$simul_summary[,1]),
                      mean(sim5_LPSMAP90_sars$simul_summary[,3]),
                      mean(sim5_LPSMAP90_sars$simul_summary[,5]))
Scenario5_sars[1,4] <- mean(sim5_LPSMAP95_sars$simul_summary[,5])
Scenario5_sars[1,5] <- mean(sim5_LPSMAP90_sars$ciwidthepilps)
Scenario5_sars[1,6] <- mean(sim5_LPSMAP95_sars$ciwidthepilps)

# LPSMALA
Scenario5_sars[2,1:3] <- c(mean(sim5_LPSMALA90_sars$simul_summary[,1]),
                      mean(sim5_LPSMALA90_sars$simul_summary[,3]),
                      mean(sim5_LPSMALA90_sars$simul_summary[,5]))
Scenario5_sars[2,4] <- mean(sim5_LPSMALA95_sars$simul_summary[,5])
Scenario5_sars[2,5] <- mean(sim5_LPSMALA90_sars$ciwidthepilps)
Scenario5_sars[2,6] <- mean(sim5_LPSMALA95_sars$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario5_sars[3,1:3] <- c(mean(sim5_LPSMAP90_sars$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim5_LPSMAP90_sars$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim5_LPSMAP90_sars$simul_summary[,6],
                           na.rm = TRUE))
Scenario5_sars[3,4] <- mean(sim5_LPSMAP95_sars$simul_summary[,6],
                            na.rm = TRUE)
Scenario5_sars[3,5] <- mean(sim5_LPSMAP90_sars$ciwidthepiestim,
                            na.rm = TRUE)
Scenario5_sars[3,6] <- mean(sim5_LPSMAP95_sars$ciwidthepiestim,
                            na.rm = TRUE)

# EpiEstim (3 days sliding window)
Scenario5_sars[4,1:3] <- c(mean(sim5_3d90_sars$simul_summary[,2], na.rm = TRUE),
                      mean(sim5_3d90_sars$simul_summary[,4],  na.rm = TRUE),
                      mean(sim5_3d90_sars$simul_summary[,6], na.rm = TRUE))
Scenario5_sars[4,4] <- mean(sim5_3d95_sars$simul_summary[,6], na.rm = TRUE)
Scenario5_sars[4,5] <- mean(sim5_3d90_sars$ciwidthepiestim, na.rm = TRUE)
Scenario5_sars[4,6] <- mean(sim5_3d95_sars$ciwidthepiestim, na.rm = TRUE)

# Extract tables
S5sars <- round(Scenario5_sars, 3)
write.table(S5sars, file="S5sarsNegBin.txt", row.names = TRUE, col.names = TRUE)
write.xlsx(S5sars, file="S5sarsNegBin.xls")
S5sars

# Save environment
save.image(file="Scenario5-CW.RData")



































