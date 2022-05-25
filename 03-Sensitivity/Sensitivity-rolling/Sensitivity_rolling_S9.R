#------------------------------------------------------------------------------#
#                        Sensitivity analysis for EpiLPS                       #
#                        Check effect of free parameters                       #
#                             Oswaldo Gressani, 2022                           #
#------------------------------------------------------------------------------#

library("EpiLPS")

#------------------------------------------------------------------------------#
#               Scenario 9 "Wiggly then flat R(t)" SI-MERS (NegBin DGP)        #
#------------------------------------------------------------------------------#

si_mers <- round(EpiEstim::discr_si(k = seq(1,20), mu = 6.8, sigma=4.1),3)
si_mers <- si_mers/sum(si_mers)


#--------- Sensitivity of EpiLPS (LPSMAP) estimate as time passes... 
set.seed(2022)
simepidemic <- episim(serial_interval = si_mers, endepi = 60, Rpattern = 5,
                  dist="negbin", overdisp = 50, verbose = TRUE)
incidence <- simepidemic$y
slidewindow <- 6 # Weekly sliding windows for EpiEstim

Rhatlps <- list()
Rhatepiestim <- list()

for(j in 1:51){
  
# Estimation with EpiLPS
epilps_fit <- epilps(incidence = incidence[1:(9+j)], K = 20, method = "LPSMAP",
                     serial_interval = si_mers, penorder = 2,
                     ci_level = 0.90, hyperprior = c(10,10),
                     verbose = FALSE)
Rhatlps[[j]] <- epilps_fit$epifit$R_estim[8:(9+j)]

# Estimation with EpiEstim
t_start <- seq(2, (9+j) - slidewindow)
t_end <- t_start + slidewindow
epiestim_fit <- suppressMessages(suppressWarnings(
  EpiEstim::estimate_R(incidence[1:(9+j)], method = "non_parametric_si",
                  config = EpiEstim::make_config(list(si_distr = c(0, si_mers),
                                                           t_start = t_start,
                                                           t_end = t_end)))))
Rhatepiestim[[j]] <- epiestim_fit$R$`Mean(R)`

}

# Plot rolling estimates (png)
tdom <- seq(1, 60, by = 0.01)
tdiscr <- seq(8, 60)
Rtarget <- sapply(tdom, simepidemic$Rtrue)
png(file = "Sensitivity_rolling.png", width = 1500, height = 900)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.3)
plot(tdom, Rtarget, type = "l", col = "black", ylab="R(t)", xlab="Time (days)",
     ylim = c(0,4), lty = 2, main = "Scenario 9", cex.lab= 1.4, cex.axis = 1.5,
     cex.main = 1.8)
grid(nx = 10, ny = 10)
for(j in 1:51){
  lines(seq(8,8+length(Rhatlps[[j]])-1), Rhatlps[[j]], type="p",
        col="cornflowerblue", pch = 2)
  lines(seq(8,8+length(Rhatepiestim[[j]])-1), Rhatepiestim[[j]], type="p",
        col="brown3", pch = 5)
}
lines(seq(8,8+length(Rhatlps[[11]])-1), Rhatlps[[11]], type="p",
      col="darkorange2", pch = 17) # t = 20
lines(seq(8,8+length(Rhatlps[[31]])-1), Rhatlps[[31]], type="p",
      col="gold1", pch = 17) # t = 40
lines(seq(8,8+length(Rhatlps[[51]])-1), Rhatlps[[51]], type="p",
      col="darkblue", pch = 17) # t = 60
legend("topright", c("Target", "Rolling estimate with LPSMAP",
        "Rolling estimate with EpiEstim", "Estimate with T=20 using LPSMAP",
        "Estimate with T=40 using LPSMAP","Estimate with T=60 using LPSMAP"),
       lty=c(2,NA,NA,NA,NA,NA),
       pch=c(NA,2,5,17,17,17),
       col=c("black","cornflowerblue","brown3","darkorange2",
             "gold1","darkblue"), bty="n",cex= 1.5)
dev.off()


# Plot rolling estimates (pdf)
tdom <- seq(1, 60, by = 0.01)
tdiscr <- seq(8, 60)
Rtarget <- sapply(tdom, simepidemic$Rtrue)
pdf(file = "Sensitivity_rolling.pdf", width = 13.5, height = 7)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.3)
plot(tdom, Rtarget, type = "l", col = "black", ylab="R(t)", xlab="Time (days)",
     ylim = c(0,4), lty = 2, main = "Scenario 9", cex.lab= 1.5, cex.axis = 1.5,
     cex.main = 1.8)
grid(nx = 10, ny = 10)
for(j in 1:51){
  lines(seq(8,8+length(Rhatlps[[j]])-1), Rhatlps[[j]], type="p",
        col="cornflowerblue", pch = 2)
  lines(seq(8,8+length(Rhatepiestim[[j]])-1), Rhatepiestim[[j]], type="p",
        col="brown3", pch = 5)
}
lines(seq(8,8+length(Rhatlps[[11]])-1), Rhatlps[[11]], type="p",
      col="darkorange2", pch = 17) # t = 20
lines(seq(8,8+length(Rhatlps[[31]])-1), Rhatlps[[31]], type="p",
      col="gold1", pch = 17) # t = 40
lines(seq(8,8+length(Rhatlps[[51]])-1), Rhatlps[[51]], type="p",
      col="darkblue", pch = 17) # t = 60
legend("topright", c("Target", "Rolling estimate with LPSMAP",
       "Rolling estimate with EpiEstim", "Estimate with T=20 using LPSMAP",
        "Estimate with T=40 using LPSMAP","Estimate with T=60 using LPSMAP"),
       lty=c(2,NA,NA,NA,NA,NA), pch=c(NA,2,5,17,17,17),
       col=c("black","cornflowerblue","brown3","darkorange2",
             "gold1","darkblue"), bty="n",cex= 1.2)
dev.off()


# Mean Absolute error of the difference between estimated R(t) and true R(t)

# t=20
MAE_EpiLPS_t20 <- mean(abs(Rhatlps[[11]]-simepidemic$Rtrue(seq(8,20))))
MAE_EpiEstim_t20 <- mean(abs(Rhatepiestim[[11]]-simepidemic$Rtrue(seq(8,20))))

# t=25
MAE_EpiLPS_t25 <- mean(abs(Rhatlps[[16]]-simepidemic$Rtrue(seq(8,25))))
MAE_EpiEstim_t25 <- mean(abs(Rhatepiestim[[16]]-simepidemic$Rtrue(seq(8,25))))

# t=30
MAE_EpiLPS_t30 <- mean(abs(Rhatlps[[21]]-simepidemic$Rtrue(seq(8,30))))
MAE_EpiEstim_t30 <- mean(abs(Rhatepiestim[[21]]-simepidemic$Rtrue(seq(8,30))))

# t=35
MAE_EpiLPS_t35 <- mean(abs(Rhatlps[[26]]-simepidemic$Rtrue(seq(8,35))))
MAE_EpiEstim_t35 <- mean(abs(Rhatepiestim[[26]]-simepidemic$Rtrue(seq(8,35))))

# t=40
MAE_EpiLPS_t40 <- mean(abs(Rhatlps[[31]]-simepidemic$Rtrue(seq(8,40))))
MAE_EpiEstim_t40 <- mean(abs(Rhatepiestim[[31]]-simepidemic$Rtrue(seq(8,40))))

# t=45
MAE_EpiLPS_t45 <- mean(abs(Rhatlps[[36]]-simepidemic$Rtrue(seq(8,45))))
MAE_EpiEstim_t45 <- mean(abs(Rhatepiestim[[36]]-simepidemic$Rtrue(seq(8,45))))

# t=50
MAE_EpiLPS_t50 <- mean(abs(Rhatlps[[41]]-simepidemic$Rtrue(seq(8,50))))
MAE_EpiEstim_t50 <- mean(abs(Rhatepiestim[[41]]-simepidemic$Rtrue(seq(8,50))))

# t=55
MAE_EpiLPS_t55 <- mean(abs(Rhatlps[[46]]-simepidemic$Rtrue(seq(8,55))))
MAE_EpiEstim_t55 <- mean(abs(Rhatepiestim[[46]]-simepidemic$Rtrue(seq(8,55))))

# t=60
MAE_EpiLPS_t60 <- mean(abs(Rhatlps[[51]]-simepidemic$Rtrue(seq(8,60))))
MAE_EpiEstim_t60 <- mean(abs(Rhatepiestim[[51]]-simepidemic$Rtrue(seq(8,60))))

MAE_summary <- matrix(0, nrow = 9, ncol = 2)
colnames(MAE_summary) <- c("MAE_EpiLPS","MAE_EpiEstim")
rownames(MAE_summary) <- c("t=20","t=25","t=30","t=35","t=40",
                           "t=45", "t=50", "t=55","t=60")
MAE_summary[,1] <- c(MAE_EpiLPS_t20, MAE_EpiLPS_t25, MAE_EpiLPS_t30,
                     MAE_EpiLPS_t35,MAE_EpiLPS_t40,MAE_EpiLPS_t45,
                     MAE_EpiLPS_t50,MAE_EpiLPS_t55,MAE_EpiLPS_t60)
MAE_summary[,2] <- c(MAE_EpiEstim_t20, MAE_EpiEstim_t25, MAE_EpiEstim_t30,
                     MAE_EpiEstim_t35,MAE_EpiEstim_t40,MAE_EpiEstim_t45,
                     MAE_EpiEstim_t50,MAE_EpiEstim_t55,MAE_EpiEstim_t60)
round(MAE_summary,3)














