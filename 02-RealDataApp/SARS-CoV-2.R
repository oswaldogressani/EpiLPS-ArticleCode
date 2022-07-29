#------------------------------------------------------------------------------#
#                       EpiLPS Covid-19 data application                       #
#                            Oswaldo Gressani, 2022                            #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("EpiEstim")
library("reshape2")

# epicurve function
epicurve <- function(epicounts, smoothfit = NULL, dates = NULL,
                     datelab = c("7d", "1m", "3m", "6m"),
                     seriesname = c("Series-1","Series-2","Series-3"),
                     themetype = c("gray","classic","light","dark"),
                     col1 = "deepskyblue4",
                     title = "Epidemic curve",
                     poslegend = "right",
                     dirlegend = "vertical",
                     xlabel = "Time (days)",
                     ylabel = "Incidence"){
  
  if(!is.data.frame(epicounts)) stop("Count data should be a data frame")
  numSeries <- ncol(epicounts)
  if(numSeries > 3) stop("Not more than 3 time series allowed")
  totDays <- nrow(epicounts)
  if (!is.null(dates)) {
    epicounts$Time <- dates
  } else{
    epicounts$Time <- seq_len(totDays)
  }
  if (!is.null(smoothfit)) {
    colnames(epicounts) <- c(paste0(seriesname, " (smooth)")[1:numSeries], "Time")
  } else{
    colnames(epicounts) <- c(seriesname[1:numSeries], "Time")
  }
  epiData <- reshape2::melt(epicounts, id = "Time")
  colnames(epiData) <- c("Time","Series","Incidence")
  if(!is.null(smoothfit)){# augment epiData
    epiData$muhat <- (unlist(sapply(smoothfit, "[[", 4)[5,]))
    epiData$muCIlow <- (unlist(sapply(smoothfit, "[[", 4)[6,]))
    epiData$muCIup <- (unlist(sapply(smoothfit, "[[", 4)[7,]))
  }
  
  # Adapt x-axis labels if dates are provided
  if(!is.null(dates)) {
    datelab <- match.arg(datelab)
    if (datelab == "7d") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '7 days')"))
    } else if (datelab == "1m") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '1 month')"))
    } else if (datelab == "3m") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '3 months')"))
    } else if (datelab == "6m") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '6 months')"))
    }
  } else{
    xlabtype <- NULL
  }
  
  # Theme choice
  themetype <- match.arg(themetype)
  if (themetype == "classic") {
    themeval <- eval(parse(text = "ggplot2::theme_classic()"))
  } else if (themetype == "gray") {
    themeval <- eval(parse(text = "ggplot2::theme_gray()"))
  } else if (themetype == "light") {
    themeval <- eval(parse(text = "ggplot2::theme_light()"))
  } else if (themetype == "dark") {
    themeval <- eval(parse(text = "ggplot2::theme_dark()"))
  }
  
  colpalette <- c(col1, "brown1", "darkolivegreen4")
  epicurvePlot <- 
    ggplot2::ggplot(data = epiData, ggplot2::aes(x = Time))+
    ggplot2::geom_bar(stat = "identity", width = 1, 
                      ggplot2::aes(y = Incidence, fill = Series),
                      position = ggplot2::position_dodge(1)) +
    ggplot2::scale_fill_manual(values=colpalette[1:numSeries])+
    ggplot2::labs(fill="Daily incidence")+
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(title) +
    themeval +
    xlabtype +
    ggplot2::theme(legend.position = poslegend,
                   legend.direction = dirlegend,
                   legend.title = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 13),
                   plot.title = ggplot2::element_text(size = 17),
                   axis.title.x = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_text(size = 14),
                   axis.text.x = ggplot2::element_text(size = 14),
                   axis.text.y = ggplot2::element_text(size = 14)) 
  
  if(!is.null(smoothfit)){
    epicurvePlot <- epicurvePlot + 
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = muCIlow, ymax = muCIup, fill = Series), alpha = 0.3) +
      ggplot2::geom_line(inherit.aes = FALSE, ggplot2::aes(x = Time, y = muhat,
                                                           colour = Series), show.legend = FALSE, size = 1.1) +
      ggplot2::scale_color_manual(values=colpalette[1:numSeries])
  }
  
  
  
  return(epicurvePlot)
}

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

# Plot the epidemic curve for each country

pdf(file = "SARS-CoV-2-Epicurves.pdf", width = 22, height = 9)
gridExtra::grid.arrange(
epicurve(as.data.frame(inciBEL), dates = dateBEL, datelab = "3m", 
           seriesname = "Belgium", title = "Epidemic curve Belgium"),  
epicurve(as.data.frame(inciDNK), dates = dateDNK, datelab = "3m", 
         seriesname = "Denmark", title = "Epidemic curve Denmark", 
         col1 = "orange"),
epicurve(as.data.frame(inciPRT), dates = datePRT, datelab = "3m", 
         seriesname = "Portugal", title = "Epidemic curve Portugal",
         col1 = "darkolivegreen"), 
epicurve(as.data.frame(inciFRA), dates = dateFRA, datelab = "3m", 
         seriesname = "France", title = "Epidemic curve France",
         col1 = "darkorchid"),
nrow = 2, ncol = 2)
dev.off()



























