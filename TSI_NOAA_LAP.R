#!/usr/bin/env Rscript
# /* Copyright (C) 2022-2023 Athanasios Natsis <natsisphysicist@gmail.com> */
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

#+ include=T, echo=F

## __ Set environment  ---------------------------------------------------------
rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
options("width" = 130)
tic <- Sys.time()
Script.Name <- "~/TSI/TSI_NOAA_LAP.R"

if (!interactive()) {
    pdf( file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink(file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".lock", Script.Name))), timeout = 0)
}

library(knitr)
library(pander)
library(caTools)
# library(RAerosols)
library(data.table)

opts_chunk$set(dev        = "png")
opts_chunk$set(comment    = NA   )
opts_chunk$set(fig.width  = 8,
               fig.height = 5    )

source("~/TSI/DEFINITIONS.R")

## __ Load data  ---------------------------------------------------------------
ASTROPY_data <- data.table(readRDS(ASTROPYdb))
NOAA_data    <- data.table(readRDS(OUTPUT_NOAA))

## NOAA data does not extend to now
ASTROPY_data <- ASTROPY_data[Date > TSI_START & Date <= max(NOAA_data$time)]
NOAA_data    <- NOAA_data[   time > TSI_START]

## remove columns
NOAA_data[, file     := NULL ]
NOAA_data[, time_low := NULL ]
NOAA_data[, time_upp := NULL ]

plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l",
     xlab = "", ylab = "Distance [au]")

#'
#' ## Interpolate TSI measurements to our dates
#'
#' Make functions from TSI measurements to match out data.
#'
#' Interpolate between measurements only.
#'
#' These data are the main TSI data we use.
#' It is completed by TSIS for the current year.
#'
#+ include = TRUE, echo = FALSE

## TSI
tsi_fun <- approxfun(x      = NOAA_data$time,
                     y      = NOAA_data$TSI,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )
## TSI uncertainty
unc_fuc <- approxfun(x      = NOAA_data$time,
                     y      = NOAA_data$TSI_UNC,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )

#'
#' Interpolate the data, assuming that dates from Astropy are complete (we made them).
#'
#+ include = TRUE, echo = FALSE
tsi_all      <- tsi_fun(ASTROPY_data$Date)
unc_all      <- unc_fuc(ASTROPY_data$Date)

#'
#' Compute TSI on earth using TSI at 1 au
#'
#+ include = TRUE, echo = FALSE
tsi_astropy  <- tsi_all / (ASTROPY_data$Dist^2)

##  Constructed TSI data for output  -------------------------------------------
tsi_comb <- data.frame(
                nominal_dates         = ASTROPY_data$Date, # Dates from SORCE extended to today
                sun_dist              = ASTROPY_data$Dist, # Astropy sun distance not optimal
                # tsi_true_earth_compex = tsi_astropy,       # Original data and extension with Astropy distance
                TSIextEARTH_comb      = tsi_astropy,       # Original data and extension with Astropy distance
                measur_error          = unc_all,           # Original data and extension of last value
                tsi_1au               = tsi_all            # TSI at 1 au
            )


par(mar = c(2, 4, 2, 1))

#'
#' ### Solar irradiance at 1 au
#'
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

pander(summary(tsi_all), digits = 8)

#'
#' ### Solar irradiance at TOA
#'
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, tsi_astropy,       "l", xlab = "", ylab = "TSI True Earth (Astropy) watt/m^2")
pander(summary(tsi_astropy),digits = 8)

#'
#' ### Solar irradiance uncertainty
#'
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, unc_all,           "l", xlab = "", ylab = "TSI Uncertainty (Interpolated) watt/m^2")
pander(summary(unc_all))

#'
#' ### Earth - Sun distance
#'
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l", xlab = "", ylab = "Sun-Thessaloniki Distance AU")
pander(summary(ASTROPY_data$Dist))


#'
#' ## Output data for use
#'
#+ include=TRUE, echo=FALSE
myRtools::write_RDS(object = tsi_comb,
                    file   = OUTPUT_NOAA_LAP)

panderOptions("table.style", 'rmarkdown')
panderOptions("table.split.table", 100 )
panderOptions("table.continues", '')
panderOptions("graph.fontsize", 9)
panderOptions("table.alignment.default", "right")

#'
#' ### Statistics on input data
#'
#+ include=TRUE, echo=FALSE
pander(summary(NOAA_data, digits = 5))

#'
#' ### Statistics on output data
#'
#+ include=TRUE, echo=FALSE
pander(summary(tsi_comb, digits = 5))

#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
