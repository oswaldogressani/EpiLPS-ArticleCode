#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 6)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")
library("xlsx")

S <- 100 # Total number of epidemics considered
seedval <- 1870

#------------------------------------------------------------------------------#
#               Scenario 6 "Step R(t)" SI-SARS (NegBin DGP)                    #
#------------------------------------------------------------------------------#

si_sars <- c(0.001,0.012,0.043,0.078,0.104,0.117,0.116,0.108,0.094,0.078,0.063,
             0.049,0.038,0.028,0.021,0.015,0.011,0.008,0.005,0.004,0.003,0.002,
             0.001,0.001)

# Above si obtained with discr_si function of EpiEstim package
# si_alter <- round(EpiEstim::discr_si(k = seq(1,24), mu = 8.4, sigma=3.8),3)
# si_alter
# sum(si_alter)
simcheck <- episim(serial_interval = si_sars, endepi = 40, Rpattern = 2,
                   dist="negbin", overdisp = 5, verbose = TRUE, plotsim = TRUE)


# LPSMAP, LPSMALA and EpiEstim with 7 days window
sim6_LPSMAP90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 5, ci_level = 0.90,
                           seed = seedval, themetype = "gray", epidays = 40,
                           midwindow = TRUE)

sim6_LPSMALA90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 5, seed = seedval,
                                themetype = "gray", epidays = 40,
                                midwindow = TRUE)

sim6_LPSMAP95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                                K = 40, method = "LPSMAP", dist = "negbin",
                                overdisp = 5, ci_level = 0.95, 
                                seed = seedval, themetype = "gray", 
                                epidays = 40, midwindow = TRUE)

sim6_LPSMALA95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.95, dist = "negbin",
                                overdisp = 5, seed = seedval,
                                themetype = "gray", epidays = 40,
                                midwindow = TRUE)

# EpiEstim with 3 days window
sim6_3d90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.90,
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 2, midwindow = TRUE)

sim6_3d95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.95, 
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 2, midwindow = TRUE)
  

#---- Extracting plots
pdf(file = "Figures/S6_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim6_LPSMAP90_sars$inciplot+ggplot2::ggtitle("Incidence (Scenario 6)")
dev.off()
pdf(file = "Figures/S6_LPSMAP_sars.pdf", width = 6, height = 4.5) 
sim6_LPSMAP90_sars$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories") +
  ggplot2::xlim(c(8,42))
dev.off()
pdf(file = "Figures/S6_LPSMALA_sars.pdf", width = 6, height = 4.5) 
sim6_LPSMALA90_sars$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories") +
  ggplot2::xlim(c(8,42))
dev.off()
pdf(file = "Figures/S6_EpiEstim7d_sars.pdf", width = 6, height = 4.5) 
sim6_LPSMAP90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories") +
  ggplot2::xlim(c(8,42))
dev.off()
pdf(file = "Figures/S6_EpiEstim3d_sars.pdf", width = 6, height = 4.5) 
sim6_3d90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories") +
  ggplot2::xlim(c(8,42))
dev.off()


png(file = "Figures/Scenario6_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim6_LPSMAP90_sars$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 6)"),
                        sim6_LPSMAP90_sars$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories") +
                          ggplot2::xlim(c(8,42)),
                        sim6_LPSMALA90_sars$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories") +
                          ggplot2::xlim(c(8,42)),
                        sim6_LPSMAP90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories") +
                          ggplot2::xlim(c(8,42)),
                        sim6_3d90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories") +
                          ggplot2::xlim(c(8,42)),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for SARS SI
Scenario6_sars <- matrix(0, nrow = 4, ncol = 6)
rownames(Scenario6_sars) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d")
colnames(Scenario6_sars) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario6_sars[1,1:3] <- c(mean(sim6_LPSMAP90_sars$simul_summary[,1]),
                      mean(sim6_LPSMAP90_sars$simul_summary[,3]),
                      mean(sim6_LPSMAP90_sars$simul_summary[,5]))
Scenario6_sars[1,4] <- mean(sim6_LPSMAP95_sars$simul_summary[,5])
Scenario6_sars[1,5] <- mean(sim6_LPSMAP90_sars$ciwidthepilps)
Scenario6_sars[1,6] <- mean(sim6_LPSMAP95_sars$ciwidthepilps)

# LPSMALA
Scenario6_sars[2,1:3] <- c(mean(sim6_LPSMALA90_sars$simul_summary[,1]),
                      mean(sim6_LPSMALA90_sars$simul_summary[,3]),
                      mean(sim6_LPSMALA90_sars$simul_summary[,5]))
Scenario6_sars[2,4] <- mean(sim6_LPSMALA95_sars$simul_summary[,5])
Scenario6_sars[2,5] <- mean(sim6_LPSMALA90_sars$ciwidthepilps)
Scenario6_sars[2,6] <- mean(sim6_LPSMALA95_sars$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario6_sars[3,1:3] <- c(mean(sim6_LPSMAP90_sars$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim6_LPSMAP90_sars$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim6_LPSMAP90_sars$simul_summary[,6],
                           na.rm = TRUE))
Scenario6_sars[3,4] <- mean(sim6_LPSMAP95_sars$simul_summary[,6],
                            na.rm = TRUE)
Scenario6_sars[3,5] <- mean(sim6_LPSMAP90_sars$ciwidthepiestim,
                            na.rm = TRUE)
Scenario6_sars[3,6] <- mean(sim6_LPSMAP95_sars$ciwidthepiestim,
                            na.rm = TRUE)

# EpiEstim (3 days sliding window)
Scenario6_sars[4,1:3] <- c(mean(sim6_3d90_sars$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim6_3d90_sars$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim6_3d90_sars$simul_summary[,6],
                           na.rm = TRUE))
Scenario6_sars[4,4] <- mean(sim6_3d95_sars$simul_summary[,6],
                            na.rm = TRUE)
Scenario6_sars[4,5] <- mean(sim6_3d90_sars$ciwidthepiestim,
                            na.rm = TRUE)
Scenario6_sars[4,6] <- mean(sim6_3d95_sars$ciwidthepiestim,
                            na.rm = TRUE)

# Extract tables
S6sars <- round(Scenario6_sars, 3)
write.table(S6sars, file="S6sarsNegBin.txt", row.names = TRUE, col.names = TRUE)
write.xlsx(S6sars, file="S6sarsNegBin.xls")
S6sars

# Save environment
save.image(file="Scenario6-CW.RData")





































