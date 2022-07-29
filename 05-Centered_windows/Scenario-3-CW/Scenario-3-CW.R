#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 3)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")
library("xlsx")

S <- 100 # Total number of epidemics considered
seedval <- 1314

#------------------------------------------------------------------------------#
#               Scenario 3 "Wiggly R(t)" SI-FLU (NegBin DGP)                   #
#------------------------------------------------------------------------------#

si_flu <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001)

# Above si obtained with discr_si function of EpiEstim package
# si_alter <- round(EpiEstim::discr_si(k = seq(1,11), mu = 2.6, sigma=1.5),3)
# si_alter
# sum(si_alter)
simcheck <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                   dist="negbin", overdisp = 1000, verbose = TRUE,
                   plotsim = TRUE)

# LPSMAP, LPSMALA and EpiEstim with 7 days window
pdf(file = "Figures/LPSMAP_S3_7d.pdf", width = 15, height = 4.6)
sim3_LPSMAP90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 3,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90,
                           seed = seedval, themetype = "gray", epidays = 40,
                           midwindow = TRUE)
dev.off()

sim3_LPSMALA90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 3,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 1000, seed = seedval,
                                themetype = "gray", epidays = 40,
                                midwindow = TRUE)

sim3_LPSMAP95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 3,
                               K = 40, method = "LPSMAP", dist = "negbin",
                               overdisp = 1000, ci_level = 0.95,
                               seed = seedval, themetype = "gray", epidays = 40,
                               midwindow = TRUE)

sim3_LPSMALA95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 3,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.95, dist = "negbin",
                                overdisp = 1000, seed = seedval,
                                themetype = "gray", epidays = 40,
                                midwindow = TRUE)

# EpiEstim with 3 days window
sim3_3d90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 3,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90,
                           seed = seedval, themetype = "gray", epidays = 40,
                           slidewindow = 2, midwindow = TRUE)

sim3_3d95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 3,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.95,
                           seed = seedval, themetype = "gray", epidays = 40,
                           slidewindow = 2, midwindow = TRUE)


#---- Extracting plots
pdf(file = "Figures/S3_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim3_LPSMAP90_flu$inciplot+ggplot2::ggtitle("Incidence (Scenario 3)")
dev.off()
pdf(file = "Figures/S3_LPSMAP_flu.pdf", width = 6, height = 4.5) 
sim3_LPSMAP90_flu$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories")
dev.off()
pdf(file = "Figures/S3_LPSMALA_flu.pdf", width = 6, height = 4.5) 
sim3_LPSMALA90_flu$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories")
dev.off()
pdf(file = "Figures/S3_EpiEstim7d_flu.pdf", width = 6, height = 4.5) 
sim3_LPSMAP90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories")
dev.off()
pdf(file = "Figures/S3_EpiEstim3d_flu.pdf", width = 6, height = 4.5) 
sim3_3d90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories")
dev.off()


png(file = "Figures/Scenario3_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim3_LPSMAP90_flu$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 3)"),
                        sim3_LPSMAP90_flu$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories"),
                        sim3_LPSMALA90_flu$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories"),
                        sim3_LPSMAP90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories"),
                        sim3_3d90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for flu SI
Scenario3_flu <- matrix(0, nrow = 4, ncol = 6)
rownames(Scenario3_flu) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d")
colnames(Scenario3_flu) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario3_flu[1,1:3] <- c(mean(sim3_LPSMAP90_flu$simul_summary[,1]),
                      mean(sim3_LPSMAP90_flu$simul_summary[,3]),
                      mean(sim3_LPSMAP90_flu$simul_summary[,5]))
Scenario3_flu[1,4] <- mean(sim3_LPSMAP95_flu$simul_summary[,5])
Scenario3_flu[1,5] <- mean(sim3_LPSMAP90_flu$ciwidthepilps)
Scenario3_flu[1,6] <- mean(sim3_LPSMAP95_flu$ciwidthepilps)

# LPSMALA
Scenario3_flu[2,1:3] <- c(mean(sim3_LPSMALA90_flu$simul_summary[,1]),
                      mean(sim3_LPSMALA90_flu$simul_summary[,3]),
                      mean(sim3_LPSMALA90_flu$simul_summary[,5]))
Scenario3_flu[2,4] <- mean(sim3_LPSMALA95_flu$simul_summary[,5])
Scenario3_flu[2,5] <- mean(sim3_LPSMALA90_flu$ciwidthepilps)
Scenario3_flu[2,6] <- mean(sim3_LPSMALA95_flu$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario3_flu[3,1:3] <- c(mean(sim3_LPSMAP90_flu$simul_summary[,2]),
                      mean(sim3_LPSMAP90_flu$simul_summary[,4]),
                      mean(sim3_LPSMAP90_flu$simul_summary[,6]))
Scenario3_flu[3,4] <- mean(sim3_LPSMAP95_flu$simul_summary[,6])
Scenario3_flu[3,5] <- mean(sim3_LPSMAP90_flu$ciwidthepiestim)
Scenario3_flu[3,6] <- mean(sim3_LPSMAP95_flu$ciwidthepiestim)

# EpiEstim (3 days sliding window)
Scenario3_flu[4,1:3] <- c(mean(sim3_3d90_flu$simul_summary[,2]),
                      mean(sim3_3d90_flu$simul_summary[,4]),
                      mean(sim3_3d90_flu$simul_summary[,6]))
Scenario3_flu[4,4] <- mean(sim3_3d95_flu$simul_summary[,6])
Scenario3_flu[4,5] <- mean(sim3_3d90_flu$ciwidthepiestim)
Scenario3_flu[4,6] <- mean(sim3_3d95_flu$ciwidthepiestim)


# Extract tables
S3flu  <- round(Scenario3_flu, 3)
write.table(S3flu, file="S3fluNegBin.txt", row.names = TRUE, col.names = TRUE)
write.xlsx(S3flu, file="S3fluNegBin.xls")
S3flu

# Save environment
save.image(file="Scenario3-CW.RData")







































