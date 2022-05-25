#------------------------------------------------------------------------------#
#                       Timing EpiLPS (LPSMAP/LPSMALA)                         #
#                          Oswaldo Gressani, 2022                              #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("progress")
library("crayon")
library("xlsx")

S <- 10 # Total number of epidemics considered

#------------------------------------------------------------------------------#
#               Scenario 3 "Wiggly R(t)" SI-FLU (NegBin DGP)                   #
#------------------------------------------------------------------------------#

si_flu <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001)

simcheck <- episim(serial_interval = si_flu, endepi = 40, Rpattern = 3,
                   dist="negbin", overdisp = 1000, verbose = TRUE,
                   plotsim = TRUE)

###############################################################################
# Timing LPSMAP 
###############################################################################
LPSMAP_TIME <- matrix(0, nrow = 5, ncol = 5)
colnames(LPSMAP_TIME) <- c("K=20","K=30","K=40","K=50","K=60")
rownames(LPSMAP_TIME) <- c("T=20","T=30","T=40","T=50","T=60")

set.seed(2489)
epitime <- c(20,30,40,50,60)

progbar <- progress::progress_bar$new(
  format = crayon::white$green("Timing LPSMAP [:elapsed :spin] [:bar] :percent"),
  total = 5,
  clear = FALSE
)

for(timeidx in 1:5){

  timingsT <- matrix(0, nrow = S, ncol = 5)

for(s in 1:S){
  simepidemic <- episim(serial_interval = si_flu, endepi = epitime[timeidx],
                        Rpattern = 3, dist="negbin", overdisp = 1000)
  incidence <- simepidemic$y

  sim3_K20 <- epilps_fit <- epilps(incidence = incidence, K = 20, 
                                   method = "LPSMAP", serial_interval = si_flu,
                                   verbose = FALSE, tictoc = TRUE)
  sim3_K30 <- epilps_fit <- epilps(incidence = incidence, K = 30, 
                                   method = "LPSMAP", serial_interval = si_flu,
                                   verbose = FALSE, tictoc = TRUE)
  sim3_K40 <- epilps_fit <- epilps(incidence = incidence, K = 40, 
                                   method = "LPSMAP", serial_interval = si_flu,
                                   verbose = FALSE, tictoc = TRUE)
  sim3_K50 <- epilps_fit <- epilps(incidence = incidence, K = 50, 
                                   method = "LPSMAP", serial_interval = si_flu,
                                   verbose = FALSE, tictoc = TRUE)
  sim3_K60 <- epilps_fit <- epilps(incidence = incidence, K = 60, 
                                   method = "LPSMAP", serial_interval = si_flu,
                                   verbose = FALSE, tictoc = TRUE)
  
  timingsT[s, ] <- c(sim3_K20$elapsed, sim3_K30$elapsed, sim3_K40$elapsed,
                      sim3_K50$elapsed, sim3_K60$elapsed)
  
}
  LPSMAP_TIME[timeidx, ] <- colMeans(timingsT)
  progbar$tick()
}

###############################################################################
# Timing LPSMALA 
###############################################################################
LPSMALA_TIME <- matrix(0, nrow = 5, ncol = 5)
colnames(LPSMALA_TIME) <- c("K=20","K=30","K=40","K=50","K=60")
rownames(LPSMALA_TIME) <- c("T=20","T=30","T=40","T=50","T=60")

set.seed(2489)
epitime <- c(20,30,40,50,60)

progbar <- progress::progress_bar$new(
  format = crayon::white$green("Timing LPSMALA [:elapsed :spin] [:bar] :percent"),
  total = 5,
  clear = FALSE
)

for(timeidx in 1:5){
  
  timingsT <- matrix(0, nrow = S, ncol = 5)
  
  for(s in 1:S){
    simepidemic <- episim(serial_interval = si_flu, endepi = epitime[timeidx],
                          Rpattern = 3, dist="negbin", overdisp = 1000)
    incidence <- simepidemic$y
    
    sim3_K20 <- epilps_fit <- epilps(incidence = incidence, K = 20, 
                                     method = "LPSMALA", 
                                     chain_length = 3000, burn = 1000,
                                     serial_interval = si_flu,
                                     verbose = FALSE, tictoc = TRUE)
    sim3_K30 <- epilps_fit <- epilps(incidence = incidence, K = 30, 
                                     method = "LPSMALA", 
                                     chain_length = 3000, burn = 1000,
                                     serial_interval = si_flu,
                                     verbose = FALSE, tictoc = TRUE)
    sim3_K40 <- epilps_fit <- epilps(incidence = incidence, K = 40, 
                                     method = "LPSMALA", 
                                     chain_length = 3000, burn = 1000,
                                     serial_interval = si_flu,
                                     verbose = FALSE, tictoc = TRUE)
    sim3_K50 <- epilps_fit <- epilps(incidence = incidence, K = 50, 
                                     method = "LPSMALA", 
                                     chain_length = 3000, burn = 1000,
                                     serial_interval = si_flu,
                                     verbose = FALSE, tictoc = TRUE)
    sim3_K60 <- epilps_fit <- epilps(incidence = incidence, K = 60, 
                                     method = "LPSMALA", 
                                     chain_length = 3000, burn = 1000,
                                     serial_interval = si_flu,
                                     verbose = FALSE, tictoc = TRUE)
    
    timingsT[s, ] <- c(sim3_K20$elapsed, sim3_K30$elapsed, sim3_K40$elapsed,
                       sim3_K50$elapsed, sim3_K60$elapsed)
    
  }
  LPSMALA_TIME[timeidx, ] <- colMeans(timingsT)
  progbar$tick()
}


# write.xlsx(as.data.frame(LPSMAP_TIME), "LPSMAP_time.xlsx",
#            sheetName = "LPSMAP", row.names = TRUE)
# write.xlsx(as.data.frame(LPSMALA_TIME), "LPSMALA_time.xlsx",
#            sheetName = "LPSMALA", row.names = TRUE)



















