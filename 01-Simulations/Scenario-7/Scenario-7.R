#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 7)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")

S <- 100 # Total number of epidemics considered
seedval <- 115 

#------------------------------------------------------------------------------#
#               Scenario 7 "Wiggly R(t)" SI-SARS (NegBin DGP)                  #
#------------------------------------------------------------------------------#

si_sars <- c(0.001,0.012,0.043,0.078,0.104,0.117,0.116,0.108,0.094,0.078,0.063,
             0.049,0.038,0.028,0.021,0.015,0.011,0.008,0.005,0.004,0.003,0.002,
             0.001,0.001)

# Above si obtained with discr_si function of EpiEstim package
# si_alter <- round(EpiEstim::discr_si(k = seq(1,24), mu = 8.4, sigma=3.8),3)
# si_alter
# sum(si_alter)
simcheck <- episim(serial_interval = si_sars, endepi = 40, Rpattern = 3,
                   dist="negbin", overdisp = 5, verbose = TRUE, plotsim = TRUE)


# LPSMAP, LPSMALA and EpiEstim with 7 days window
sim7_LPSMAP90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 5, ci_level = 0.90,
                           seed = seedval, themetype = "gray", epidays = 40)

sim7_LPSMALA90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 5, seed = seedval,
                                themetype = "gray", epidays = 40)

sim7_LPSMAP95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                                K = 40, method = "LPSMAP", dist = "negbin",
                                overdisp = 5, ci_level = 0.95, 
                                seed = seedval, themetype = "gray",
                                epidays = 40)

sim7_LPSMALA95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.95, dist = "negbin",
                                overdisp = 5, seed = seedval,
                                themetype = "gray", epidays = 40)

# EpiEstim with 3 days window
sim7_3d90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.90,
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 2)

sim7_3d95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.95, 
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 2)
  
# EpiEstim with 1 day window
sim7_1d90_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.90, 
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 0)

sim7_1d95_sars <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 5, ci_level = 0.95, 
                            seed = seedval, themetype = "gray", epidays = 40,
                            slidewindow = 0)

#---- Extracting plots
pdf(file = "Figures/S7_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim7_LPSMAP90_sars$inciplot+ggplot2::ggtitle("Incidence (Scenario 6)")
dev.off()
pdf(file = "Figures/S7_LPSMAP_sars.pdf", width = 6, height = 4.5) 
sim7_LPSMAP90_sars$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories")
dev.off()
pdf(file = "Figures/S7_LPSMALA_sars.pdf", width = 6, height = 4.5) 
sim7_LPSMALA90_sars$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories")
dev.off()
pdf(file = "Figures/S7_EpiEstim7d_sars.pdf", width = 6, height = 4.5) 
sim7_LPSMAP90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories")
dev.off()
pdf(file = "Figures/S7_EpiEstim3d_sars.pdf", width = 6, height = 4.5) 
sim7_3d90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories")
dev.off()
pdf(file = "Figures/S7_EpiEstim1d_sars.pdf", width = 6, height = 4.5) 
sim7_1d90_sars$Repiesplot+ggplot2::ggtitle("EpiEstim 1d windows trajectories")
dev.off()

png(file = "Figures/Scenario7_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim7_LPSMAP90_sars$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 6)"),
                        sim7_LPSMAP90_sars$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories"),
                        sim7_LPSMALA90_sars$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories"),
                        sim7_LPSMAP90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories"),
                        sim7_3d90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories"),
                        sim7_1d90_sars$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 1d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for SARS SI
Scenario7_sars <- matrix(0, nrow = 5, ncol = 6)
rownames(Scenario7_sars) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d",
                             "EpiEstim 1d")
colnames(Scenario7_sars) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario7_sars[1,1:3] <- c(mean(sim7_LPSMAP90_sars$simul_summary[,1]),
                      mean(sim7_LPSMAP90_sars$simul_summary[,3]),
                      mean(sim7_LPSMAP90_sars$simul_summary[,5]))
Scenario7_sars[1,4] <- mean(sim7_LPSMAP95_sars$simul_summary[,5])
Scenario7_sars[1,5] <- mean(sim7_LPSMAP90_sars$ciwidthepilps)
Scenario7_sars[1,6] <- mean(sim7_LPSMAP95_sars$ciwidthepilps)

# LPSMALA
Scenario7_sars[2,1:3] <- c(mean(sim7_LPSMALA90_sars$simul_summary[,1]),
                      mean(sim7_LPSMALA90_sars$simul_summary[,3]),
                      mean(sim7_LPSMALA90_sars$simul_summary[,5]))
Scenario7_sars[2,4] <- mean(sim7_LPSMALA95_sars$simul_summary[,5])
Scenario7_sars[2,5] <- mean(sim7_LPSMALA90_sars$ciwidthepilps)
Scenario7_sars[2,6] <- mean(sim7_LPSMALA95_sars$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario7_sars[3,1:3] <- c(mean(sim7_LPSMAP90_sars$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim7_LPSMAP90_sars$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim7_LPSMAP90_sars$simul_summary[,6],
                           na.rm = TRUE))
Scenario7_sars[3,4] <- mean(sim7_LPSMAP95_sars$simul_summary[,6],
                            na.rm = TRUE)
Scenario7_sars[3,5] <- mean(sim7_LPSMAP90_sars$ciwidthepiestim,
                            na.rm = TRUE)
Scenario7_sars[3,6] <- mean(sim7_LPSMAP95_sars$ciwidthepiestim,
                            na.rm = TRUE)

# EpiEstim (3 days sliding window)
Scenario7_sars[4,1:3] <- c(mean(sim7_3d90_sars$simul_summary[,2],
                                na.rm = TRUE),
                      mean(sim7_3d90_sars$simul_summary[,4],
                           na.rm = TRUE),
                      mean(sim7_3d90_sars$simul_summary[,6],
                           na.rm = TRUE))
Scenario7_sars[4,4] <- mean(sim7_3d95_sars$simul_summary[,6],
                            na.rm = TRUE)
Scenario7_sars[4,5] <- mean(sim7_3d90_sars$ciwidthepiestim,
                            na.rm = TRUE)
Scenario7_sars[4,6] <- mean(sim7_3d95_sars$ciwidthepiestim,
                            na.rm = TRUE)

# EpiEstim (1 day sliding window)
Scenario7_sars[5,1:3] <- c(mean(sim7_1d90_sars$simul_summary[,2],
                                na.rm = TRUE),
                          mean(sim7_1d90_sars$simul_summary[,4],
                               na.rm = TRUE),
                          mean(sim7_1d90_sars$simul_summary[,6],
                               na.rm = TRUE))
Scenario7_sars[5,4] <- mean(sim7_1d95_sars$simul_summary[,6],
                            na.rm = TRUE)
Scenario7_sars[5,5] <- mean(sim7_1d90_sars$ciwidthepiestim,
                            na.rm = TRUE)
Scenario7_sars[5,6] <- mean(sim7_1d95_sars$ciwidthepiestim,
                            na.rm = TRUE)

# Extract tables
S7sars <- round(Scenario7_sars, 3)
write.table(S7sars, file="S7sarsNegBin.txt", row.names = TRUE, col.names = TRUE)
S7sars

# Save environment
save.image(file="Scenario7.RData")


# Plot mean variance (of yt) relationship for the Scenario
mu_gen <- matrix(0,nrow = S, ncol = 40)
var_gen <- matrix(0,nrow = S, ncol = 40)
set.seed(seedval)
for(s in 1:S){
  simepidemic <- episim(serial_interval = si_sars, endepi = 40, Rpattern = 3,
                        dist="negbin", overdisp = 5, 
                        verbose = FALSE,plotsim = FALSE)
  mu_gen[s, ] <- simepidemic$mu_y
  var_gen[s,] <- simepidemic$mu_y*(1+simepidemic$mu_y/1000) 
}

# Plot the means and variances
png(file = "Figures/Scenario7_Mean_Variance.png", width = 1500, height = 900)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 40), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 7",
     cex.lab= 1.8, cex.axis = 1.8,
     cex.main = 1.8)
lines(var_gen[1, ], type = "p", pch = 19, col = "red")
for (s in 2:S) {
  # Plot the means of y
  lines(mu_gen[s,], type = "p", pch = 25, col = "blue")
  # Plot the variances of y
  lines(var_gen[s,], type = "p", pch = 19, col = "red")
}
legend("topleft", c("Mean of yt","Variance of yt"),
       col=c("blue","red"), pch = c(25,19), bty = "n",cex=1.8)
dev.off()


# Plot the means and variances
pdf(file = "Figures/Scenario7_Mean_Variance.pdf", width = 13.5, height = 7)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 40), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 7",
     cex.lab= 1.8, cex.axis = 1.8,
     cex.main = 1.8)
lines(var_gen[1, ], type = "p", pch = 19, col = "red")
for (s in 2:S) {
  # Plot the means of y
  lines(mu_gen[s,], type = "p", pch = 25, col = "blue")
  # Plot the variances of y
  lines(var_gen[s,], type = "p", pch = 19, col = "red")
}
legend("topleft", c("Mean of yt","Variance of yt"),
       col=c("blue","red"), pch = c(25,19), bty = "n",cex=1.8)
dev.off()





































