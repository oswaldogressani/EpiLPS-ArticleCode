#------------------------------------------------------------------------------#
#                    Simulations to assess the performance of the              #
#                      estimate for the overdisperion parameter                #
#                             Oswaldo Gressani, 2022                           #
#------------------------------------------------------------------------------#

library("EpiLPS")

S <- 100 # Total number of epidemics considered
seedval <- 2025

si_sars <- c(0.001,0.012,0.043,0.078,0.104,0.117,0.116,0.108,0.094,0.078,0.063,
             0.049,0.038,0.028,0.021,0.015,0.011,0.008,0.005,0.004,0.003,0.002,
             0.001,0.001)

#------------------------------------------------------------------------------#
#               Scenario 5 "Constant R(t)" SI-SARS (NegBin DGP)                #
#------------------------------------------------------------------------------#

sim5_LPSMAP <- perfcheck(S = S, serial_interval = si_sars, scenario = 1,
                           K = 40, method = "LPSMAP", dist = "negbin",
                           overdisp = 5, ci_level = 0.90, Rconst = 1.3,
                           seed = seedval, themetype = "gray", epidays = 40,
                           midwindow = TRUE)
sim5_logrho_LPSMAP <- log(sim5_LPSMAP$dispvec)

#------------------------------------------------------------------------------#
#               Scenario 6 "Step R(t)" SI-SARS (NegBin DGP)                    #
#------------------------------------------------------------------------------#

sim6_LPSMAP <- perfcheck(S = S, serial_interval = si_sars, scenario = 2,
                          K = 40, method = "LPSMAP", dist = "negbin",
                          overdisp = 5, ci_level = 0.90,
                          seed = seedval, themetype = "gray", epidays = 40,
                          midwindow = TRUE)
sim6_logrho_LPSMAP <- log(sim6_LPSMAP$dispvec)


#------------------------------------------------------------------------------#
#               Scenario 7 "Wiggly R(t)" SI-SARS (NegBin DGP)                  #
#------------------------------------------------------------------------------#

sim7_LPSMAP <- perfcheck(S = S, serial_interval = si_sars, scenario = 3,
                                K = 40, method = "LPSMAP", dist = "negbin",
                                overdisp = 5, ci_level = 0.90,
                                seed = seedval, themetype = "gray", epidays = 40,
                                midwindow = TRUE)
sim7_logrho_LPSMAP <- log(sim7_LPSMAP$dispvec)


#------------------------------------------------------------------------------#
#               Scenario 8 "Decaying R(t)" SI-SARS (NegBin DGP)                #
#------------------------------------------------------------------------------#

sim8_LPSMAP <- perfcheck(S = S, serial_interval = si_sars, scenario = 4,
                              K = 40, method = "LPSMAP", dist = "negbin",
                              overdisp = 5, ci_level = 0.90,
                              seed = seedval, themetype = "gray", epidays = 40,
                              midwindow = TRUE)
sim8_logrho_LPSMAP <- log(sim8_LPSMAP$dispvec)


#------------------------------------------------------------------------------#
#                                  Boxplot                                     #
#------------------------------------------------------------------------------#

pdf(file = "Boxplot_overdispersion.pdf", width = 12, height = 7.5)
boxplot(sim5_logrho_LPSMAP, sim6_logrho_LPSMAP, sim7_logrho_LPSMAP,
        sim8_logrho_LPSMAP, 
        names = c("Scenario 5", "Scenario 6","Scenario 7", "Scenario 8"),
        col=c("cornflowerblue","coral2","cyan2","yellow"))
abline(h = log(5), col = "red", lwd = 2, lty = 2)
legend("topright", col="red", "Target", lwd = 2, lty = 2, bty = "n")
dev.off()



# Save environment
save.image(file="Simulations_overdispersion.RData")



































