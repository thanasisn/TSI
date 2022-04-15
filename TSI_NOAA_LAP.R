#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisthanasis@gmail.com> */

#' ---
#' title:  "TSI NOAA data preparation."
#' author: "Natsis Athanasios"
#' institute: "AUTH"
#' affiliation: "Laboratory of Atmospheric Physics"
#' abstract: "Read original total solar irradiance data from NOAA
#'            last measurement as the Solar Constant."
#' output:
#'   html_document:
#'     toc: true
#'     fig_width:  9
#'     fig_height: 5
#'   pdf_document:
#' date: "`r format(Sys.time(), '%F')`"
#' ---

#+ include=FALSE

rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
options("width" = 130)
tic <- Sys.time()
Script.Name <- tryCatch({ funr::sys.script() },
                        error = function(e) { cat(paste("\nUnresolved script name: ", e),"\n")
                            return("Undefined R script name!!") })
if(!interactive()) {
    pdf(  file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink( file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".lock", Script.Name))), timeout = 0)
}




library(knitr)
opts_chunk$set( comment    = NA )
opts_chunk$set( fig.width  = 8,
                fig.height = 5  )

library(pander)
library(caTools)
library(RAerosols)
library(data.table)


####    Variables    ####
source("~/TSI/DEFINITIONS.R")


ASTROPY_data <- data.table(readRDS(ASTROPYdb))
NOAA_data    <- data.table(readRDS(OUTPUT_NOAA))


ASTROPY_data <- ASTROPY_data[ Date > TSI_START ]
NOAA_data    <- NOAA_data[    time > TSI_START ]

NOAA_data[, file     := NULL ]
NOAA_data[, time_low := NULL ]
NOAA_data[, time_upp := NULL ]



#+ include=FALSE

stop()

#' ##### Test if the sun is on (zero measurements may occur) #####
####  Subset to project dates  ####
data  <- data[data$time > TSI_START,]

####  Remove invalid records
lower_tsi_limit <- min(mean(data$TSI), median(data$TSI)) * 0.70
data            <- data[ data$TSI > lower_tsi_limit, ]


#+ include=FALSE



####  Extend measurements until today  ####
## use last measurement date and time step from origin
# timestep        = tsi_df$nominal_date_jdn[100] - tsi_df$nominal_date_jdn[99]
# last_measu_date = max(tsi_df$nominal_date_jdn) + timestep

#' #### Interpolate TSI measurements to our data ####
#' Make functions from TSI measurements to match out data.
#' Interpolate between measurements and extend tails as closest constant.
tsi_fun <- approxfun(x      = data$time,
                     y      = data$TSI,
                     method = "linear",
                     rule   = 2,
                     f      = 0,
                     ties   = mean )

#+ include=FALSE
cat("Interpolate between measurements and extend tail with closest constant\n")

#+ include=TRUE
unc_fuc <- approxfun(x      = data$time,
                     y      = data$TSI_UNC,
                     method = "linear",
                     rule   = 2,
                     f      = 0,
                     ties   = mean )

#' Interpolate the data, we have assumed that dates from Astropy are complete (we made them).
tsi_all      <- tsi_fun( ASTROPY_data$Date )  # SORCE dates to project dates
unc_all      <- unc_fuc( ASTROPY_data$Date )

#' Compute TSI on earth using TSI at 1 au
tsi_astropy  <- tsi_all / ( ASTROPY_data$Dist ^ 2 )

#+ include=FALSE
#### Constructed TSI data for output --------------------------------
tsi_comb <- data.frame(
                nominal_dates         = ASTROPY_data$Date, # Dates from SORCE extended to today
                sun_dist              = ASTROPY_data$Dist, # Astropy sun distance not optimal
                # tsi_true_earth_compex = tsi_astropy,       # Original data and extension with Astropy distance
                TSIextEARTH_comb      = tsi_astropy,       # Original data and extension with Astropy distance
                measur_error          = unc_all,           # Original data and extension of last value
                tsi_1au               = tsi_all            # TSI at 1 au
            )

#+ include=FALSE
## Some test graphs
# pdf(file = tsigraphs)
#     par(mar = c(2,4,2,1))
#     pdf.options(height=5)
#     plot(ASTROPY_data$Date, tsi_all,           "l", ylab = "TSI (Interpolated) watt/m^2")
#     plot(ASTROPY_data$Date, tsi_astropy,       "l", ylab = "TSI True Earth (Astropy) watt/m^2")
#     plot(ASTROPY_data$Date, unc_all,           "l", ylab = "TSI Uncertainty (Interpolated) watt/m^2")
#     plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l", ylab = "Sun-Thessaloniki Distance AU")
# dev.off()


#+ include=TRUE, echo=FALSE
par(mar = c(2,4,2,1))
#' ### Solar irradiance at 1 au ###
#+ include=TRUE, echo=FALSE
ylim = range(c( tsi_all + 0.1, tsi_all - 0.1 ), na.rm = T)
plot( ASTROPY_data$Date, tsi_all, "l", xlab = "", ylab = "TSI (Interpolated) watt/m^2", ylim = ylim)
lines(ASTROPY_data$Date, runmean(tsi_all, 15000), col = 5, lwd = 3)
qq <- quantile(tsi_all, na.rm = T)

abline( h = qq[3], col = 'orange',lwd = 3 )
abline( h = qq[2], col = 'green', lwd = 2 )
abline( h = qq[4], col = 'green', lwd = 2 )
abline( h = qq[1], col = 'red',   lwd = 1 )
abline( h = qq[5], col = 'red',   lwd = 1 )

text(ASTROPY_data$Date[1], y = qq[3],
     as.character(round(qq[3],1)), adj = c(0,1.5), col = 'orange', lwd = 4 )

text(ASTROPY_data$Date[1], y = qq[2],
     as.character(round(qq[2],1)), adj = c(0,1.5), col = 'green', lwd = 4 )
text(ASTROPY_data$Date[1], y = qq[4],
     as.character(round(qq[4],1)), adj = c(0,1.5), col = 'green', lwd = 4 )

text(ASTROPY_data$Date[1], y = qq[1],
     as.character(round(qq[1],1)), adj = c(0,1.5), col = 'red', lwd = 4 )
text(ASTROPY_data$Date[1], y = qq[5],
     as.character(round(qq[5],1)), adj = c(0,1.5), col = 'red', lwd = 4 )

pander(summary(tsi_all),digits = 8)

#' ### Solar irradiance at TOA
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, tsi_astropy,       "l", xlab = "", ylab = "TSI True Earth (Astropy) watt/m^2")
pander(summary(tsi_astropy),digits = 8)

#' ### Solar irradiance uncertainty
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, unc_all,           "l", xlab = "", ylab = "TSI Uncertainty (Interpolated) watt/m^2")
pander(summary(unc_all))

#' ### Earth - Sun distance
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l", xlab = "", ylab = "Sun-Thessaloniki Distance AU")
pander(summary(ASTROPY_data$Dist))

#+ include=FALSE
## Export last day of SORCE measurements to file
# LAST_TSI_DATE = max(tsi_df$nominal_date_yyyymmdd)
# save( LAST_TSI_DATE, file = GLOBAL_file)


#' #### Output data for use ####
#+ include=TRUE, echo=FALSE
myRtools::write_RDS(object = tsi_comb,
                    file   = tsi_Rds_noaa  )

#+ include=FALSE


panderOptions("table.style", 'rmarkdown')
panderOptions("table.split.table", 100 )
panderOptions("table.continues", '')
panderOptions("graph.fontsize", 9)
panderOptions("table.alignment.default", "right")


#' ### Statistics on input data ###
#+ include=TRUE, echo=FALSE
pander(summary(data,   digits = 5))

#' ### Statistics on output data ###
#+ include=TRUE, echo=FALSE
pander(summary(tsi_comb, digits = 5))



#+ include=FALSE, echo=FALSE
tac = Sys.time();
cat(paste("\n  --  ",  Script.Name, " DONE  --  \n\n"))
cat(sprintf("%s %-10s %-10s %-20s  %f mins\n\n",Sys.time(),Sys.info()["nodename"],Sys.info()["login"],Script.Name,difftime(tac,tic,units="mins")))
write(sprintf("%s %-10s %-10s %-50s %f",Sys.time(),Sys.info()["nodename"],Sys.info()["login"],Script.Name,difftime(tac,tic,units="mins")),"~/Aerosols/run.log",append=T)
