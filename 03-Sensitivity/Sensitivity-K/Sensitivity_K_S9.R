#------------------------------------------------------------------------------#
#                        Sensitivity analysis for EpiLPS                       #
#                        Check effect of free parameters                       #
#                             Oswaldo Gressani, 2022                           #
#------------------------------------------------------------------------------#

library("EpiLPS")

S <- 50 # Total number of epidemics considered

#------------------------------------------------------------------------------#
#               Scenario 9 "Wiggly then flat R(t)" SI-MERS (NegBin DGP)        #
#------------------------------------------------------------------------------#

si_mers <- round(EpiEstim::discr_si(k = seq(1,20), mu = 6.8, sigma=4.1),3)
si_mers <- si_mers/sum(si_mers)


simcheck <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                dist="negbin", overdisp = 50, verbose = TRUE, plotsim = TRUE)


#--------- Sensitivity with respect to the number of B-splines K 
RK<- matrix(0, nrow = 5, ncol = 53)

# K1: K=20
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 53)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                        dist="negbin", overdisp = 50)
  incidence <- simepidemic$y


  epilps_fit <- epilps(incidence = incidence, K = 20, method = "LPSMAP",
                       serial_interval = si_mers, penorder = 2,
                       ci_level = 0.90, hyperprior = c(10,10),
                       verbose = FALSE)

Rhat_mat[s, ] <- epilps_fit$epifit[8:60, 2]
}

RK[1, ] <- apply(Rhat_mat, 2, "median")

# K2: K=30
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 53)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                        dist="negbin", overdisp = 50)
  incidence <- simepidemic$y
  
  
  epilps_fit <- epilps(incidence = incidence, K = 30, method = "LPSMAP",
                       serial_interval = si_mers, penorder = 2,
                       ci_level = 0.90, hyperprior = c(10,10),
                       verbose = FALSE)
  
  Rhat_mat[s, ] <- epilps_fit$epifit[8:60, 2]
}

RK[2, ] <- apply(Rhat_mat, 2, "median")

# K3: K=40
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 53)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                        dist="negbin", overdisp = 50)
  incidence <- simepidemic$y
  
  
  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                       serial_interval = si_mers, penorder = 2,
                       ci_level = 0.90, hyperprior = c(10,10),
                       verbose = FALSE)
  
  Rhat_mat[s, ] <- epilps_fit$epifit[8:60, 2]
}

RK[3, ] <- apply(Rhat_mat, 2, "median")

# K4: K=50
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 53)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                        dist="negbin", overdisp = 50)
  incidence <- simepidemic$y
  
  
  epilps_fit <- epilps(incidence = incidence, K = 50, method = "LPSMAP",
                       serial_interval = si_mers, penorder = 2,
                       ci_level = 0.90, hyperprior = c(10,10),
                       verbose = FALSE)
  
  Rhat_mat[s, ] <- epilps_fit$epifit[8:60, 2]
}

RK[4, ] <- apply(Rhat_mat, 2, "median")


# K5: K=60
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 53)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                        dist="negbin", overdisp = 50)
  incidence <- simepidemic$y
  
  
  epilps_fit <- epilps(incidence = incidence, K = 60, method = "LPSMAP",
                       serial_interval = si_mers, penorder = 2,
                       ci_level = 0.90, hyperprior = c(10,10),
                       verbose = FALSE)
  
  Rhat_mat[s, ] <- epilps_fit$epifit[8:60, 2]
}

RK[5, ] <- apply(Rhat_mat, 2, "median")


#------ Plot sensitivity of R(t) estimated with respect to hyperparameter delta

#Png 
tdom <- seq(1, 60, by = 0.01)
tdiscr <- seq(8, 60)
Rtarget <- sapply(tdom, simepidemic$Rtrue)
png(file = "Sensitivity_K_Scenario9.png", width = 1300, height = 900)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.3)
plot(tdom, Rtarget, type = "l", col = "black", ylab="R(t)", xlab="Time (days)",
     ylim = c(0,4), lty = 2, main = "Scenario 9", cex.lab= 1.4, cex.axis = 1.5,
     cex.main = 1.8)
grid(nx = 10, ny = 10)
abline(h=1, col="purple", lty = 3)
lines(tdiscr, RK[1,], col="blue", lwd = 2)            # K=20
lines(tdiscr, RK[2,], col="orange", lwd = 2)          # K=30
lines(tdiscr, RK[3,], col="firebrick", lwd = 2)       # K=40
lines(tdiscr, RK[4,], col="darkolivegreen4", lwd = 2) # K=50
lines(tdiscr, RK[5,], col="gold4", lwd = 2)           # K=60
legend("topright",
       c("Target", "K=20", "K=30", "K=40", "K=50", "K=60"),
       col=c("black","blue","orange","firebrick","darkolivegreen4","gold4"),
       lty=c(2,1,1,1,1,1),
       cex = 1.5, bty = "n", y.intersp = 1)
dev.off()

#Pdf 
tdom <- seq(1, 60, by = 0.01)
tdiscr <- seq(8, 60)
Rtarget <- sapply(tdom, simepidemic$Rtrue)
pdf(file = "Sensitivity_K_Scenario9.pdf", width = 13.5, height = 7)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.3)
plot(tdom, Rtarget, type = "l", col = "black", ylab="R(t)", xlab="Time (days)",
     ylim = c(0,4), lty = 2, main = "Scenario 9", cex.lab= 1.5, cex.axis = 1.5,
     cex.main = 1.8)
grid(nx = 10, ny = 10)
abline(h=1, col="purple", lty = 3)
lines(tdiscr, RK[1,], col="blue", lwd = 2)            # K=20
lines(tdiscr, RK[2,], col="orange", lwd = 2)          # K=30
lines(tdiscr, RK[3,], col="firebrick", lwd = 2)       # K=40
lines(tdiscr, RK[4,], col="darkolivegreen4", lwd = 2) # K=50
lines(tdiscr, RK[5,], col="gold4", lwd = 2)           # K=60
legend("topright",
       c("Target", "K=20", "K=30", "K=40", "K=50", "K=60"),
       col=c("black","blue","orange","firebrick","darkolivegreen4","gold4"),
       lty=c(2,1,1,1,1,1),
       cex = 1.5, bty = "n", y.intersp = 0.9)
dev.off()































