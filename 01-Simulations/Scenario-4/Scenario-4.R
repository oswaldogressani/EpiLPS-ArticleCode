#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 4)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")

S <- 100 # Total number of epidemics considered
seedval <- 1966

#------------------------------------------------------------------------------#
#               Scenario 4 "Decaying R(t)" SI-FLU (NegBin DGP)                 #
#------------------------------------------------------------------------------#

si_flu <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001)

# Above si obtained with discr_si function of EpiEstim package
# si_alter <- round(EpiEstim::discr_si(k = seq(1,11), mu = 2.6, sigma=1.5),3)
# si_alter
# sum(si_alter)
simcheck <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 4,
                   dist="negbin", overdisp = 1000, verbose = TRUE,
                   plotsim = TRUE)


# LPSMAP, LPSMALA and EpiEstim with 7 days window
sim4_LPSMAP90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90, seed = seedval, 
                           themetype = "gray", epidays = 40)

sim4_LPSMALA90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 1000, seed = seedval,
                                themetype = "gray", epidays = 40)

sim4_LPSMAP95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                               K = 40, method = "LPSMAP", dist = "negbin",
                               overdisp = 1000, ci_level = 0.95, seed = seedval, 
                               themetype = "gray", epidays = 40)

sim4_LPSMALA95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.95, dist = "negbin",
                                overdisp = 1000, seed = seedval,
                                themetype = "gray", epidays = 40)

# EpiEstim with 3 days window
sim4_3d90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90, seed = seedval, 
                           themetype = "gray", epidays = 40, slidewindow = 2)

sim4_3d95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.95, seed = seedval, 
                           themetype = "gray", epidays = 40, slidewindow = 2)

# EpiEstim with 1 day window
sim4_1d90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90, seed = seedval, 
                           themetype = "gray", epidays = 40, slidewindow = 0)

sim4_1d95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 4,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.95, seed = seedval, 
                           themetype = "gray", epidays = 40, slidewindow = 0)

#---- Extracting plots
pdf(file = "Figures/S4_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim4_LPSMAP90_flu$inciplot+ggplot2::ggtitle("Incidence (Scenario 4)")
dev.off()
pdf(file = "Figures/S4_LPSMAP_flu.pdf", width = 6, height = 4.5) 
sim4_LPSMAP90_flu$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories")
dev.off()
pdf(file = "Figures/S4_LPSMALA_flu.pdf", width = 6, height = 4.5) 
sim4_LPSMALA90_flu$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories")
dev.off()
pdf(file = "Figures/S4_EpiEstim7d_flu.pdf", width = 6, height = 4.5) 
sim4_LPSMAP90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories")
dev.off()
pdf(file = "Figures/S4_EpiEstim3d_flu.pdf", width = 6, height = 4.5) 
sim4_3d90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories")
dev.off()
pdf(file = "Figures/S4_EpiEstim1d_flu.pdf", width = 6, height = 4.5) 
sim4_1d90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 1d windows trajectories")
dev.off()

png(file = "Figures/Scenario4_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim4_LPSMAP90_flu$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 4)"),
                        sim4_LPSMAP90_flu$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories"),
                        sim4_LPSMALA90_flu$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories"),
                        sim4_LPSMAP90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories"),
                        sim4_3d90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories"),
                        sim4_1d90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 1d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for flu SI

Scenario4_flu <- matrix(0, nrow = 5, ncol = 6)
rownames(Scenario4_flu) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d",
                             "EpiEstim 1d")
colnames(Scenario4_flu) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario4_flu[1,1:3] <- c(mean(sim4_LPSMAP90_flu$simul_summary[,1]),
                      mean(sim4_LPSMAP90_flu$simul_summary[,3]),
                      mean(sim4_LPSMAP90_flu$simul_summary[,5]))
Scenario4_flu[1,4] <- mean(sim4_LPSMAP95_flu$simul_summary[,5])
Scenario4_flu[1,5] <- mean(sim4_LPSMAP90_flu$ciwidthepilps)
Scenario4_flu[1,6] <- mean(sim4_LPSMAP95_flu$ciwidthepilps)

# LPSMALA
Scenario4_flu[2,1:3] <- c(mean(sim4_LPSMALA90_flu$simul_summary[,1]),
                      mean(sim4_LPSMALA90_flu$simul_summary[,3]),
                      mean(sim4_LPSMALA90_flu$simul_summary[,5]))
Scenario4_flu[2,4] <- mean(sim4_LPSMALA95_flu$simul_summary[,5])
Scenario4_flu[2,5] <- mean(sim4_LPSMALA90_flu$ciwidthepilps)
Scenario4_flu[2,6] <- mean(sim4_LPSMALA95_flu$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario4_flu[3,1:3] <- c(mean(sim4_LPSMAP90_flu$simul_summary[,2]),
                      mean(sim4_LPSMAP90_flu$simul_summary[,4]),
                      mean(sim4_LPSMAP90_flu$simul_summary[,6]))
Scenario4_flu[3,4] <- mean(sim4_LPSMAP95_flu$simul_summary[,6])
Scenario4_flu[3,5] <- mean(sim4_LPSMAP90_flu$ciwidthepiestim)
Scenario4_flu[3,6] <- mean(sim4_LPSMAP95_flu$ciwidthepiestim)

# EpiEstim (3 days sliding window)
Scenario4_flu[4,1:3] <- c(mean(sim4_3d90_flu$simul_summary[,2]),
                      mean(sim4_3d90_flu$simul_summary[,4]),
                      mean(sim4_3d90_flu$simul_summary[,6]))
Scenario4_flu[4,4] <- mean(sim4_3d95_flu$simul_summary[,6])
Scenario4_flu[4,5] <- mean(sim4_3d90_flu$ciwidthepiestim)
Scenario4_flu[4,6] <- mean(sim4_3d95_flu$ciwidthepiestim)

# EpiEstim (1 day sliding window)
Scenario4_flu[5,1:3] <- c(mean(sim4_1d90_flu$simul_summary[,2]),
                          mean(sim4_1d90_flu$simul_summary[,4]),
                          mean(sim4_1d90_flu$simul_summary[,6]))
Scenario4_flu[5,4] <- mean(sim4_1d95_flu$simul_summary[,6])
Scenario4_flu[5,5] <- mean(sim4_1d90_flu$ciwidthepiestim)
Scenario4_flu[5,6] <- mean(sim4_1d95_flu$ciwidthepiestim)

# Extract tables
S4flu  <- round(Scenario4_flu, 3)
write.table(S4flu, file="S4fluNegBin.txt", row.names = TRUE, col.names = TRUE)
S4flu

# Save environment
save.image(file="Scenario4.RData")


# Plot mean variance (of yt) relationship for the Scenario
mu_gen <- matrix(0,nrow = S, ncol = 40)
var_gen <- matrix(0,nrow = S, ncol = 40)
set.seed(seedval)
for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 4,
                        dist="negbin", overdisp = 1000, verbose = FALSE, 
                        plotsim = FALSE)
  mu_gen[s, ] <- simepidemic$mu_y
  var_gen[s,] <- simepidemic$mu_y*(1+simepidemic$mu_y/1000) 
}

# Plot the means and variances
png(file = "Figures/Scenario4_Mean_Variance.png", width = 1500, height = 900)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 40), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 4",
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
       col=c("blue","red"), pch = c(25,19), bty = "n", cex=1.8)
dev.off()


# Plot the means and variances
pdf(file = "Figures/Scenario4_Mean_Variance.pdf", width = 13.5, height = 7)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 40), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 4",
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
       col=c("blue","red"), pch = c(25,19), bty = "n", cex=1.8)
dev.off()



































