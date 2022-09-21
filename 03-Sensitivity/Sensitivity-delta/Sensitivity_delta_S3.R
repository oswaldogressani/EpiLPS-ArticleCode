#------------------------------------------------------------------------------#
#                        Sensitivity analysis for EpiLPS                       #
#                        Check effect of free parameters                       #
#                             Oswaldo Gressani, 2022                           #
#------------------------------------------------------------------------------#

library("EpiLPS")

S <- 50 # Total number of epidemics considered

#------------------------------------------------------------------------------#
#               Scenario 3 "Wiggly R(t)" SI-FLU (NegBin DGP)                   #
#------------------------------------------------------------------------------#

si_flu <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001)

simcheck <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                   dist="negbin", overdisp = 1000, verbose = TRUE,
                   plotsim = TRUE)

#------ Sensitivity with respect to Gamma parameters of hyperparameter delta
Rpriors  <- matrix(0, nrow = 6, ncol = 33)

# P1: a_delta=5; b_delta=5
set.seed(123)

Rhat_mat <- matrix(0, nrow = S, ncol = 33)
RCIlowa5 <- matrix(0, nrow = S, ncol = 33)
RCIupa5 <- matrix(0, nrow = S, ncol = 33)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                  dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y

  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                     serial_interval = si_flu, penorder = 2,
                     ci_level = 0.90, hyperprior = c(5,5),
                     verbose = FALSE)

  Rhat_mat[s, ] <- epilps_fit$epifit[8:40, 2]
  RCIlowa5[s, ] <- epilps_fit$epifit[8:40, 3]
  RCIupa5[s, ] <- epilps_fit$epifit[8:40, 4]
}


Rpriors[1, ] <- apply(Rhat_mat, 2, "median")

# P2: a_delta=10; b_delta=10
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 33)
RCIlowa10 <- matrix(0, nrow = S, ncol = 33)
RCIupa10 <- matrix(0, nrow = S, ncol = 33)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                      dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y
  
  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                    serial_interval = si_flu, penorder = 2,
                    ci_level = 0.90, hyperprior = c(10,10),
                    verbose = FALSE)
  Rhat_mat[s, ] <- epilps_fit$epifit[8:40, 2]
  RCIlowa10[s, ] <- epilps_fit$epifit[8:40, 3]
  RCIupa10[s, ] <- epilps_fit$epifit[8:40, 4]

}

Rpriors[2, ] <- apply(Rhat_mat, 2, "median")


# P3: a_delta=20; b_delta=20
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 33)
RCIlowa20 <- matrix(0, nrow = S, ncol = 33)
RCIupa20 <- matrix(0, nrow = S, ncol = 33)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                        dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y
  
  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                       serial_interval = si_flu, penorder = 2,
                       ci_level = 0.90, hyperprior = c(20,20),
                       verbose = FALSE)
  Rhat_mat[s, ] <- epilps_fit$epifit[8:40, 2]
  RCIlowa20[s, ] <- epilps_fit$epifit[8:40, 3]
  RCIupa20[s, ] <- epilps_fit$epifit[8:40, 4]
}

Rpriors[3, ] <- apply(Rhat_mat, 2, "median")


# P4: a_delta=30; b_delta=30
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 33)
RCIlowa30 <- matrix(0, nrow = S, ncol = 33)
RCIupa30 <- matrix(0, nrow = S, ncol = 33)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                        dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y
  
  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                       serial_interval = si_flu, penorder = 2,
                       ci_level = 0.90, hyperprior = c(30,30),
                       verbose = FALSE)
  Rhat_mat[s, ] <- epilps_fit$epifit[8:40, 2]
  RCIlowa30[s, ] <- epilps_fit$epifit[8:40, 3]
  RCIupa30[s, ] <- epilps_fit$epifit[8:40, 4]
}

Rpriors[4, ] <- apply(Rhat_mat, 2, "median")


# P5: a_delta=50; b_delta=50
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 33)
RCIlowa50 <- matrix(0, nrow = S, ncol = 33)
RCIupa50 <- matrix(0, nrow = S, ncol = 33)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                        dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y
  
  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                       serial_interval = si_flu, penorder = 2,
                       ci_level = 0.90, hyperprior = c(50,50),
                       verbose = FALSE)
  
  Rhat_mat[s, ] <- epilps_fit$epifit[8:40, 2]
  RCIlowa50[s, ] <- epilps_fit$epifit[8:40, 3]
  RCIupa50[s, ] <- epilps_fit$epifit[8:40, 4]
}
Rpriors[5, ] <- apply(Rhat_mat, 2, "median")


# P6: a_delta=60; b_delta=60
set.seed(123)
Rhat_mat <- matrix(0, nrow = S, ncol = 33)
RCIlowa60 <- matrix(0, nrow = S, ncol = 33)
RCIupa60 <- matrix(0, nrow = S, ncol = 33)


for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                        dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y
  
  epilps_fit <- epilps(incidence = incidence, K = 40, method = "LPSMAP",
                       serial_interval = si_flu, penorder = 2,
                       ci_level = 0.90, hyperprior = c(60,60),
                       verbose = FALSE)
  RCIlowa60[s, ] <- epilps_fit$epifit[8:40, 3]
  RCIupa60[s, ] <- epilps_fit$epifit[8:40, 4]
  Rhat_mat[s, ] <- epilps_fit$epifit[8:40, 2]
}

Rpriors[6, ] <- apply(Rhat_mat, 2, "median")

#------ Plot sensitivity of estimated R(t) with respect to hyperparameter delta

# Png plot
tdom <- seq(1, 40, by = 0.01)
tdiscr <- seq(8, 40)
Rtarget <- sapply(tdom, simepidemic$Rtrue)
png(file = "Sensitivity_delta_Scenario3.png", width = 1300, height = 900)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.3)
plot(tdom, Rtarget, type = "l", col = "black", ylab="R(t)", xlab="Time (days)",
     ylim = c(0,4), lty = 2, main = "Scenario 3", cex.lab= 1.4, cex.axis = 1.5,
     cex.main = 1.8)
grid(nx = 10, ny = 10)
abline(h = 1, col = "purple", lty = 3)
lines(tdiscr, Rpriors[1,], col="blue", lwd = 2)     # a = 5; b = 5
lines(tdiscr, Rpriors[2,], col="orange", lwd = 2)   # a = 10; b = 10
lines(tdiscr, Rpriors[3,], col="firebrick", lwd = 2)# a = 20; b = 20
lines(tdiscr, Rpriors[4,], col="darkolivegreen4", lwd = 2)  # a = 30; b = 30
lines(tdiscr, Rpriors[5,], col="mediumorchid", lwd = 2)     # a = 50; b = 50
lines(tdiscr, Rpriors[6,], col="gray", lwd = 2)             # a = 60; b = 60
legend("topright", c("Target", "a_delta=5, b_delta=5","a_delta=10, b_delta=10",
                     "a_delta=20, b_delta=20", "a_delta=30, b_delta=30",
                     "a_delta=50, b_delta=50","a_delta=60, b_delta=60"),
    col=c("black","blue","orange","firebrick","darkolivegreen4","mediumorchid",
          "gray"), lty=c(2,1,1,1,1,1,1), cex = 1.5, bty = "n", y.intersp = 1)
dev.off()


# Pdf plot
tdom <- seq(1, 40, by = 0.01)
tdiscr <- seq(8, 40)
Rtarget <- sapply(tdom, simepidemic$Rtrue)
pdf(file = "Sensitivity_delta_Scenario3.pdf", width = 13.5, height = 7)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.3)
plot(tdom, Rtarget, type = "l", col = "black", ylab="R(t)", xlab="Time (days)",
     ylim = c(0,4), lty = 2, main = "Scenario 3", cex.lab= 1.5, cex.axis = 1.5,
     cex.main = 1.8)
grid(nx = 10, ny = 10)
abline(h = 1, col = "purple", lty = 3)
lines(tdiscr, Rpriors[1,], col="blue", lwd = 2)        # a = 5; b = 5
lines(tdiscr, Rpriors[2,], col="orange", lwd = 2)      # a = 10; b = 10
lines(tdiscr, Rpriors[3,], col="firebrick", lwd = 2)   # a = 20; b = 20
lines(tdiscr, Rpriors[4,], col="darkolivegreen4", lwd = 2)  # a = 30; b = 30
lines(tdiscr, Rpriors[5,], col="mediumorchid", lwd = 2)     # a = 50; b = 50
lines(tdiscr, Rpriors[6,], col="gray", lwd = 2)             # a = 60; b = 60
legend("topright", c("Target", "a_delta=5, b_delta=5","a_delta=10, b_delta=10",
                     "a_delta=20, b_delta=20", "a_delta=30, b_delta=30",
                     "a_delta=50, b_delta=50","a_delta=60, b_delta=60"),
       col=c("black","blue","orange","firebrick","darkolivegreen4","mediumorchid",
             "gray"), lty=c(2,1,1,1,1,1,1), cex = 1.3, bty = "n", 
       y.intersp = 0.9)
dev.off()


deltapriors <- matrix(0,nrow = 6, ncol=4)
rownames(deltapriors) <- c("P1", "P2", "P3", "P4", "P5", "P6")
colnames(deltapriors) <- c("a_delta","b_delta","Mean","Variance")
deltapriors[1,] <- c(5,5,5/5,5/(5^2))
deltapriors[2,] <- c(10,10,10/10,10/(10^2))
deltapriors[3,] <- c(20,20,20/20,20/(20^2))
deltapriors[4,] <- c(30,30,30/30,30/(30^2))
deltapriors[5,] <- c(50,50,50/50,50/(50^2))
deltapriors[6,] <- c(60,60,60/60,60/(60^2))
round(deltapriors, 3)

#--------------------------- Sensitivity of CI 


# Plot the CI lower bounds
pdf(file = "Sensitivity_CI_Scenario3.pdf", width = 15, height = 21)

par(mfrow = c(2,1))

# P1: a_delta=5; b_delta=5
plot(seq(1,33),RCIlowa5[1,], type="p", ylim = c(0,5), col="darkolivegreen",
     pch = 0, main = "Sensitivity of 90% CI lower bound (Scenario 3)", 
     xlab = "Time (days)", ylab = "CI lower bound")
for(s in 2:S){
  lines(seq(1,33),RCIlowa5[s,], type = "p", col = "darkolivegreen", pch = 0)
}
# P2: a_delta=10; b_delta=10
for(s in 1:S){
  lines(seq(1,33),RCIlowa10[s,], type = "p", col = "coral", pch = 1)
}
# P3: a_delta=20; b_delta=20
for(s in 1:S){
  lines(seq(1,33),RCIlowa20[s,], type = "p", col = "cornflowerblue", pch = 2)
}
# P4: a_delta=30; b_delta=30
for(s in 1:S){
  lines(seq(1,33),RCIlowa20[s,], type = "p", col = "gray84", pch = 5)
}
# P5: a_delta=50; b_delta=50
for(s in 1:S){
  lines(seq(1,33),RCIlowa50[s,], type = "p", col = "lightgoldenrod2", pch = 3)
}
# P6: a_delta=60; b_delta=60
for(s in 1:S){
  lines(seq(1,33),RCIlowa60[s,], type = "p", col = "darksalmon", pch = 6)
}
legend("top", horiz = TRUE, c("5","10","20","30","50","60"),
       col=c("darkolivegreen","coral","cornflowerblue","gray84",
             "lightgoldenrod2","darksalmon"), pch = c(0,1,2,5,3,6), 
       bty = "n", text.width = 2)
legend("topleft", "(a_delta=b_delta)", bty = "n")



# Plot the CI upper bounds

# P1: a_delta=5; b_delta=5
plot(seq(1,33),RCIupa5[1,], type="p", ylim = c(0,5), col="darkolivegreen", 
     pch = 0, main = "Sensitivity of 90% CI upper bound (Scenario 3)", 
     xlab = "Time (days)", ylab = "CI upper bound")
for(s in 2:S){
  lines(seq(1,33),RCIupa5[s,], type = "p", col = "darkolivegreen", pch = 0)
}
# P2: a_delta=10; b_delta=10
for(s in 1:S){
  lines(seq(1,33),RCIupa10[s,], type = "p", col = "coral", pch = 1)
}
# P3: a_delta=20; b_delta=20
for(s in 1:S){
  lines(seq(1,33),RCIupa20[s,], type = "p", col = "cornflowerblue", pch = 2)
}
# P4: a_delta=30; b_delta=30
for(s in 1:S){
  lines(seq(1,33),RCIupa20[s,], type = "p", col = "gray84", pch = 5)
}
# P5: a_delta=50; b_delta=50
for(s in 1:S){
  lines(seq(1,33),RCIupa50[s,], type = "p", col = "lightgoldenrod2", pch = 3)
}
# P6: a_delta=60; b_delta=60
for(s in 1:S){
  lines(seq(1,33),RCIupa60[s,], type = "p", col = "darksalmon", pch = 6)
}
legend("top", horiz = TRUE, c("5","10","20","30","50","60"),
       col=c("darkolivegreen","coral","cornflowerblue","gray84",
             "lightgoldenrod2","darksalmon"), pch = c(0,1,2,5,3,6), 
       bty = "n", text.width = 2)
legend("topleft", "(a_delta=b_delta)", bty = "n")


dev.off()
































