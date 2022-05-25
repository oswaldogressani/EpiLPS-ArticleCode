#------------------------------------------------------------------------------#
#                         EpiLPS real data applications                        #
#                            Oswaldo Gressani, 2022                            #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("EpiEstim")

#------------------------------------------------------------------------------#
#                   1. Application to SARS epidemic in Hong Kong               #
#------------------------------------------------------------------------------#
data("SARS2003")
incidence <- SARS2003$incidence
si <- SARS2003$si_distr[2:25]

set.seed(15478)
epifit_sars <- epilps(incidence = incidence, K = 40, serial_interval = si,
                 method = "LPSMALA", chain_length = 25000, burn = 12000,
                 etainit = c(0,0))

round(epifit_sars$epifit[22,],2) # First peak
round(epifit_sars$epifit[39,],2) # Second peak

#------------------------------------------------------------------------------#
#        2. Application to the pandemic influenza in Pennsylvania, USA         #
#------------------------------------------------------------------------------#
data("Flu2009")
incidence <- Flu2009$incidence[,2]
si <- Flu2009$si_distr[-1]

set.seed(31415)
epifit_flu <- epilps(incidence = incidence, K = 40, serial_interval = si,
                     method = "LPSMALA", chain_length = 25000, burn = 15000,
                     etainit=c(0,0))

head(epifit_flu$epifit,28)
round(epifit_flu$epifit[13,2:4],2)

#----- Plot results
pdf(file = "SARS-Influenza.pdf", width = 13.5, height = 7)
suppressWarnings(gridExtra::grid.arrange(plot(epifit_sars, plotout = "epicurve", 
                             epicol = "chocolate4",
                             epititle = "SARS, Hong Kong, 2003",
                             incibars = TRUE, themetype = "light"),
                      plot(epifit_flu, plotout = "epicurve", 
                      epicol = "chocolate4",
                      epititle = "Pandemic influenza, Pennsylvania (USA), 2009",
                             incibars = TRUE, themetype = "light"),
                        plot(epifit_sars, plotout = "rt",
                             rtcol = "deepskyblue4", themetype = "light",
                             overlayEpiestim = FALSE),
                        plot(epifit_flu, plotout = "rt",
                             rtcol = "deepskyblue4", themetype = "light",
                             overlayEpiestim = FALSE),
                        nrow = 2, ncol = 2))
dev.off()

#------- Overlay EpiEstim
suppressWarnings(gridExtra::grid.arrange(plot(epifit_sars, plotout = "epicurve", 
                                        epicol = "chocolate4",
                                        epititle = "SARS, Hong Kong, 2003",
                                        incibars = TRUE, themetype = "light"),
                                        plot(epifit_flu, plotout = "epicurve", 
                                        epicol = "chocolate4",
                      epititle = "Pandemic influenza, Pennsylvania (USA), 2009",
                                    incibars = TRUE, themetype = "light"),
                                    plot(epifit_sars, plotout = "rt",
                                    rtcol = "deepskyblue4", themetype = "light",
                                    overlayEpiestim = TRUE),
                                    plot(epifit_flu, plotout = "rt",
                                    rtcol = "deepskyblue4", themetype = "light",
                                    overlayEpiestim = TRUE),
                                    nrow = 2, ncol = 2))














