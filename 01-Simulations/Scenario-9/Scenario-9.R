#------------------------------------------------------------------------------#
#                   Simulations comparing EpiLPS and EpiEstim                  #
#                     Negative Binomial setting (Scenario 9)                   #
#                           Oswaldo Gressani, 2022                             #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
library("gridExtra")

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
                           seed = seedval, themetype = "gray", epidays = 60)
dev.off()

sim9_LPSMALA90_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                                K = 40, method = "LPSMALA", chain_length = 3000,
                                burn = 1000, ci_level = 0.90, dist = "negbin",
                                overdisp = 50, seed = seedval,
                                themetype = "gray", epidays = 60)

sim9_LPSMAP95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                                K = 40, method = "LPSMAP", dist = "negbin",
                                overdisp = 50, ci_level = 0.95,
                                seed = seedval, themetype = "gray", epidays = 60)

sim9_LPSMALA95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                                 K = 40, method = "LPSMALA", chain_length = 3000,
                                 burn = 1000, ci_level = 0.95, dist = "negbin",
                                 overdisp = 50, seed = seedval,
                                 themetype = "gray", epidays = 60)

# EpiEstim with 3 days window
sim9_3d90_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 50, ci_level = 0.90,
                            seed = seedval, themetype = "gray", epidays = 60,
                            slidewindow = 2)

sim9_3d95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                            K = 40, method = "LPSMAP", dist = "negbin",
                            overdisp = 50, ci_level = 0.95,
                            seed = seedval, themetype = "gray", epidays = 60,
                            slidewindow = 2)
  
# EpiEstim with 1 day window
sim9_1d90_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                             K = 40, method = "LPSMAP", dist = "negbin",
                             overdisp = 50, ci_level = 0.90,
                             seed = seedval, themetype = "gray", epidays = 60,
                             slidewindow = 0)

sim9_1d95_mers <- perfcheck(S = S, serial_interval = si_mers, scenario = 5,
                             K = 40, method = "LPSMAP", dist = "negbin",
                             overdisp = 50, ci_level = 0.95,
                             seed = seedval, themetype = "gray", epidays = 60,
                             slidewindow = 0)

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
pdf(file = "Figures/S9_EpiEstim1d_mers.pdf", width = 6, height = 4.5) 
sim9_1d90_mers$Repiesplot+ggplot2::ggtitle("EpiEstim 1d windows trajectories")
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
                        sim9_1d90_mers$Repiesplot+
                          ggplot2::ggtitle("EpiEstim 1d windows trajectories"),
                        nrow = 3, ncol = 2)
dev.off()

# LPSMAP, LPSMALA and EpiEstim with 7 days window
pdf(file = "Figures/LPSMAP_S9_7d.pdf", width = 15, height = 4.6)
gridExtra::grid.arrange(
  sim9_LPSMAP90_mers$inciplot+ggplot2::ggtitle("Incidence (Scenario 9)"), 
  sim9_LPSMAP90_mers$Rlpsplot+ggplot2::ggtitle("LPSMAP trajectories") +
    ggplot2::xlim(c(8,62)),
  sim9_LPSMAP90_mers$Repiesplot+ggplot2::ggtitle("EpiEstim 7d windows trajectories") +
    ggplot2::xlim(c(8,62)),
  nrow = 1, ncol = 3)
dev.off()
                        
#------ Populating table for MERS SI
Scenario9_mers <- matrix(0, nrow = 5, ncol = 6)
rownames(Scenario9_mers) <- c("LPSMAP", "LPSMALA", "EpiEstim 7d", "EpiEstim 3d",
                             "EpiEstim 1d")
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

# EpiEstim (1 day sliding window)
Scenario9_mers[5,1:3] <- c(mean(sim9_1d90_mers$simul_summary[,2],
                                na.rm = TRUE),
                          mean(sim9_1d90_mers$simul_summary[,4],
                               na.rm = TRUE),
                          mean(sim9_1d90_mers$simul_summary[,6],
                               na.rm = TRUE))
Scenario9_mers[5,4] <- mean(sim9_1d95_mers$simul_summary[,6],
                            na.rm = TRUE)
Scenario9_mers[5,5] <- mean(sim9_1d90_mers$ciwidthepiestim,
                            na.rm = TRUE)
Scenario9_mers[5,6] <- mean(sim9_1d95_mers$ciwidthepiestim,
                            na.rm = TRUE)

# Extract tables
S9mers <- round(Scenario9_mers, 3)
write.table(S9mers, file="S9mersNegBin.txt", row.names = TRUE, col.names = TRUE)
S9mers

# Save environment
save.image(file="Scenario9.RData")


# Plot mean variance (of yt) relationship for the Scenario
mu_gen <- matrix(0,nrow = S, ncol = 60)
var_gen <- matrix(0,nrow = S, ncol = 60)
set.seed(seedval)
for(s in 1:S){
  simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                        dist="negbin", overdisp = 50, 
                        verbose = FALSE,plotsim = FALSE)
  mu_gen[s, ] <- simepidemic$mu_y
  var_gen[s,] <- simepidemic$mu_y*(1+simepidemic$mu_y/1000) 
}

# Plot the means and variances
png(file = "Figures/Scenario9_Mean_Variance.png", width = 1500, height = 900)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 60), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 9",
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
pdf(file = "Figures/Scenario9_Mean_Variance.pdf", width = 13.5, height = 7)
plot(mu_gen[1, ], type = "p", pch = 25, 
     ylim = c(0, max(cbind(mu_gen, var_gen))),
     xlim = c(0, 60), ylab = "", xlab = "Time (in days)",
     col = "blue", main = "Scenario 9",
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






































