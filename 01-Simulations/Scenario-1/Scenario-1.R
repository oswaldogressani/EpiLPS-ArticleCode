#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 1)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")

S <- 100 # Total number of epidemics considered
seedval <- 1325

#------------------------------------------------------------------------------#
#               Scenario 1 "Constant R(t)" SI-FLU (NegBin DGP)                 #
#------------------------------------------------------------------------------#

si_flu <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001)

# Above si obtained with discr_si function of EpiEstim package
# si_alter <- round(EpiEstim::discr_si(k = seq(1,11), mu = 2.6, sigma=1.5),3)
# si_alter
# sum(si_alter)
# simcheck <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 1,
#                    Rconst = 1.3, dist="negbin",
#                    overdisp = 1000, verbose = TRUE, plotsim = TRUE)


# LPSMAP, LPSMALA and EpiEstim with 7 days window
sim1_LPSMAP90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40)

sim1_LPSMALA90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 1000, Rconst = 1.3, seed = seedval,
                                themetype = "gray", epidays = 40)

sim1_LPSMAP95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                               K = 40, method = "LPSMAP", dist = "negbin",
                               overdisp = 1000, ci_level = 0.95, Rconst = 1.3,
                               seed = seedval, themetype = "gray", epidays = 40)

sim1_LPSMALA95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.95, dist = "negbin",
                                overdisp = 1000, Rconst = 1.3, seed = seedval,
                                themetype = "gray", epidays = 40)

# EpiEstim with 3 days window
sim1_3d90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40,
                           slidewindow = 2)

sim1_3d95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.95, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40,
                           slidewindow = 2)

# EpiEstim with 1 day window
sim1_1d90_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.90, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40,
                           slidewindow = 0)

sim1_1d95_flu <- perfcheck(S = S, serial_interval = si_flu, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 1000, ci_level = 0.95, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40,
                           slidewindow = 0)

#---- Extracting plots
pdf(file = "Figures/S1_LPSMAP_Incidence.pdf", width = 6, height = 4.5) 
sim1_LPSMAP90_flu$inciplot+ggplot2::ggtitle("Incidence (Scenario 1)")
dev.off()
pdf(file = "Figures/S1_LPSMAP_flu.pdf", width = 6, height = 4.5) 
sim1_LPSMAP90_flu$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories")
dev.off()
pdf(file = "Figures/S1_LPSMALA_flu.pdf", width = 6, height = 4.5) 
sim1_LPSMALA90_flu$Rlpsplot+ggplot2::ggtitle("LPSMALA trajectories")
dev.off()
pdf(file = "Figures/S1_EpiEstim7d_flu.pdf", width = 6, height = 4.5) 
sim1_LPSMAP90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories")
dev.off()
pdf(file = "Figures/S1_EpiEstim3d_flu.pdf", width = 6, height = 4.5) 
sim1_3d90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 3d windows trajectories")
dev.off()
pdf(file = "Figures/S1_EpiEstim1d_flu.pdf", width = 6, height = 4.5) 
sim1_1d90_flu$Repiesplot+ggplot2::ggtitle("EpiEstim 1d windows trajectories")
dev.off()

png(file = "Figures/Scenario1_Summary_plots.png", width = 1000, height = 1100)
gridExtra::grid.arrange(sim1_LPSMAP90_flu$inciplot+
                          ggplot2::ggtitle("Incidence (Scenario 1)"),
                        sim1_LPSMAP90_flu$Rlpsplot+
                          ggplot2::ggtitle("LPSMAP trajectories"),
                        sim1_LPSMALA90_flu$Rlpsplot+
                          ggplot2::ggtitle("LPSMALA trajectories"),
                        sim1_LPSMAP90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 7d windows trajectories"),
                        sim1_3d90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 3d windows trajectories"),
                        sim1_1d90_flu$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 1d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()
                        
#------ Populating table for flu SI

Scenario1_flu <- matrix(0, nrow = 5, ncol = 6)
rownames(Scenario1_flu) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d",
                             "EpiEstim 1d")
colnames(Scenario1_flu) <- c("Bias", "MSE", "CP90%", "CP95%", "90% CI width",
                             "95% CI width")

# LPSMAP
Scenario1_flu[1,1:3] <- c(mean(sim1_LPSMAP90_flu$simul_summary[,1]),
                      mean(sim1_LPSMAP90_flu$simul_summary[,3]),
                      mean(sim1_LPSMAP90_flu$simul_summary[,5]))
Scenario1_flu[1,4] <- mean(sim1_LPSMAP95_flu$simul_summary[,5])
Scenario1_flu[1,5] <- mean(sim1_LPSMAP90_flu$ciwidthepilps)
Scenario1_flu[1,6] <- mean(sim1_LPSMAP95_flu$ciwidthepilps)

# LPSMALA
Scenario1_flu[2,1:3] <- c(mean(sim1_LPSMALA90_flu$simul_summary[,1]),
                      mean(sim1_LPSMALA90_flu$simul_summary[,3]),
                      mean(sim1_LPSMALA90_flu$simul_summary[,5]))
Scenario1_flu[2,4] <- mean(sim1_LPSMALA95_flu$simul_summary[,5])
Scenario1_flu[2,5] <- mean(sim1_LPSMALA90_flu$ciwidthepilps)
Scenario1_flu[2,6] <- mean(sim1_LPSMALA95_flu$ciwidthepilps)

# EpiEstim (7 days sliding window)
Scenario1_flu[3,1:3] <- c(mean(sim1_LPSMAP90_flu$simul_summary[,2]),
                      mean(sim1_LPSMAP90_flu$simul_summary[,4]),
                      mean(sim1_LPSMAP90_flu$simul_summary[,6]))
Scenario1_flu[3,4] <- mean(sim1_LPSMAP95_flu$simul_summary[,6])
Scenario1_flu[3,5] <- mean(sim1_LPSMAP90_flu$ciwidthepiestim)
Scenario1_flu[3,6] <- mean(sim1_LPSMAP95_flu$ciwidthepiestim)

# EpiEstim (3 days sliding window)
Scenario1_flu[4,1:3] <- c(mean(sim1_3d90_flu$simul_summary[,2]),
                      mean(sim1_3d90_flu$simul_summary[,4]),
                      mean(sim1_3d90_flu$simul_summary[,6]))
Scenario1_flu[4,4] <- mean(sim1_3d95_flu$simul_summary[,6])
Scenario1_flu[4,5] <- mean(sim1_3d90_flu$ciwidthepiestim)
Scenario1_flu[4,6] <- mean(sim1_3d95_flu$ciwidthepiestim)

# EpiEstim (1 day sliding window)
Scenario1_flu[5,1:3] <- c(mean(sim1_1d90_flu$simul_summary[,2]),
                          mean(sim1_1d90_flu$simul_summary[,4]),
                          mean(sim1_1d90_flu$simul_summary[,6]))
Scenario1_flu[5,4] <- mean(sim1_1d95_flu$simul_summary[,6])
Scenario1_flu[5,5] <- mean(sim1_1d90_flu$ciwidthepiestim)
Scenario1_flu[5,6] <- mean(sim1_1d95_flu$ciwidthepiestim)

# Extract tables
S1flu  <- round(Scenario1_flu, 3)
write.table(S1flu, file="S1fluNegBin.txt", row.names = TRUE, col.names = TRUE)
S1flu

# Save environment
save.image(file = "Scenario1.RData")

# Plot mean variance (of yt) relationship for the Scenario
mu_gen <- matrix(0,nrow = S, ncol = 40)
var_gen <- matrix(0,nrow = S, ncol = 40)
set.seed(seedval)
for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 1,
                        Rconst = 1.3, dist = "negbin",
                        overdisp = 1000, verbose = FALSE, plotsim = FALSE)
  mu_gen[s, ] <- simepidemic$mu_y
  var_gen[s,] <- simepidemic$mu_y*(1+simepidemic$mu_y/1000) 
}

# Print the plots 
png(file = "Figures/Scenario1_Mean_Variance.png", width = 1500, height = 900)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 40), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 1",
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
       col=c("blue","red"), pch = c(25,19), bty = "n",
       cex=1.8)
dev.off()

pdf(file = "Figures/Scenario1_Mean_Variance.pdf", width = 13.5, height = 7)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 40), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 1",
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
       col=c("blue","red"), pch = c(25,19), bty = "n",
       cex=1.8)
dev.off()



































