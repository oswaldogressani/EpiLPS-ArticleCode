
#------------------- Obtaining the dataset Covid19Data.Rda from the COVID19 pkg

#---- Date retrieved: 2022-07-14

library("COVID19")

#-------------------------   Belgium hospitalization data
Belgium <- covid19(country = "BEL", level = 1, 
                   start = "2020-04-05",
                   end = "2021-10-31")

inciBEL <- Belgium$hosp
dateBEL <- Belgium$date

#-------------------------   Denmark hospitalization data
Denmark <- covid19(country = "DNK", level = 1, 
                   start = "2020-04-05",
                   end = "2021-10-31")

inciDNK <- Denmark$hosp
dateDNK <- Denmark$date

#-------------------------   Portugal hospitalization data
Portugal <- covid19(country = "PRT", level = 1, 
                   start = "2020-04-05",
                   end = "2021-10-31")

inciPRT <- Portugal$hosp
datePRT <- Portugal$date

#-------------------------   Portugal hospitalization data
France <- covid19(country = "FRA", level = 1, 
                    start = "2020-04-05",
                    end = "2021-10-31")

inciFRA <- France$hosp
dateFRA <- France$date




Covid19DAT <- data.frame(dateBEL, inciBEL, dateDNK, inciDNK, datePRT, inciPRT, 
                          dateFRA, inciFRA)

save(Covid19DAT, file="Covid19Data.Rda")