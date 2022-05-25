#------------------------------------------------------------------------------#
#                       EpiLPS Covid-19 data application                       #
#                            Oswaldo Gressani, 2022                            #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("EpiEstim")

#library("COVID19")
load("Covid19Data.Rda")

# Uniform serial interval over 5 days
si <- discr_si(k = seq(1, 5), mu = 3, sigma = 2.48)
si <- si/sum(si)

# Belgium
# Belgium <- covid19(country = "BEL", level = 1)
dateBEL <- Covid19DAT$dateBEL
inciBEL <- Covid19DAT$inciBEL

# Denmark
# Denmark<- covid19(country = "DNK", level = 1)
dateDNK <- Covid19DAT$dateDNK
inciDNK <- Covid19DAT$inciDNK

# Portugal
# Portugal<- covid19(country = "PRT", level = 1)
datePRT <- Covid19DAT$datePRT
inciPRT <- Covid19DAT$inciPRT

# France
# France <- covid19(country = "FRA", level = 1)
dateFRA <- Covid19DAT$dateFRA
inciFRA <- Covid19DAT$inciFRA

tic <- proc.time()
# Fit with EpiLPS
epiBEL <- epilps(incidence = inciBEL, serial_interval = si)
epiDNK <- epilps(incidence = inciDNK, serial_interval = si)
epiPRT <- epilps(incidence = inciPRT, serial_interval = si)
epiFRA <- epilps(incidence = inciFRA, serial_interval = si)
toc <- proc.time() - tic
toc

pdf(file = "SARS-CoV-2.pdf", width = 15, height = 12)
gridExtra::grid.arrange(
plot(epiBEL, dates = dateBEL, datelab = "3m",
     rtcol = "steelblue", Rtitle = "Estimated R Belgium",
     themetype = "classic", overlayEpiestim = TRUE, tcut = 2),
plot(epiDNK, dates = dateDNK, datelab = "3m",
      rtcol = "chartreuse4", Rtitle = "Estimated R Denmark",
      themetype = "classic", overlayEpiestim = TRUE, tcut = 1),
plot(epiPRT, dates = datePRT, datelab = "3m",
     rtcol = "brown2", Rtitle = "Estimated R Portugal",
     themetype = "classic", overlayEpiestim = TRUE, tcut = 3),
plot(epiFRA, dates = dateFRA, datelab = "3m",
     rtcol = "darkorchid1", Rtitle = "Estimated R France",
     themetype = "classic", overlayEpiestim = TRUE, tcut = 3),
nrow = 4, ncol = 1)
dev.off()





