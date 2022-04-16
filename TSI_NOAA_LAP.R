







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


####    Check if need to run    ####


if ( !file.exists(OUTPUT_NOAA_LAP) |
     file.mtime(ASTROPYdb)   > file.mtime(OUTPUT_NOAA_LAP) |
     file.mtime(OUTPUT_NOAA) > file.mtime(OUTPUT_NOAA_LAP) ) {
    cat("New data to parse\n\n")
} else {
    cat("NO new data to parse\n\n")
}


ASTROPY_data <- data.table(readRDS(ASTROPYdb))
NOAA_data    <- data.table(readRDS(OUTPUT_NOAA))


ASTROPY_data <- ASTROPY_data[ Date > TSI_START ]
NOAA_data    <- NOAA_data[    time > TSI_START ]

NOAA_data[, file     := NULL ]
NOAA_data[, time_low := NULL ]
NOAA_data[, time_upp := NULL ]


plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l")


#+ include=FALSE



#' #### Interpolate TSI measurements to our dates ####
#' Make functions from TSI measurements to match out data.
#' Interpolate between measurements only.
tsi_fun <- approxfun(x      = NOAA_data$time,
                     y      = NOAA_data$TSI,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )

#+ include=FALSE
cat("Interpolate between measurements\n")

#+ include=TRUE
unc_fuc <- approxfun(x      = NOAA_data$time,
                     y      = NOAA_data$TSI_UNC,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )



#' Interpolate the data, we have assumed that dates from Astropy are complete (we made them).
tsi_all      <- tsi_fun( ASTROPY_data$Date )
unc_all      <- unc_fuc( ASTROPY_data$Date )

#' Compute TSI on earth using TSI at 1 au
tsi_astropy  <- tsi_all / ( ASTROPY_data$Dist ^ 2 )

#+ include=FALSE
#### Constructed TSI data for output --------------------------------
tsi_comb <- data.frame(
    nominal_dates = ASTROPY_data$Date, # Dates from SORCE extended to today
    sun_dist      = ASTROPY_data$Dist, # Astropy sun distance not optimal
    TSIextEARTH   = tsi_astropy,       # Original data with Astropy distance
    measur_error  = unc_all,           # Original data
    tsi_1au       = tsi_all            # TSI at 1 au
)


#' #### Output data for use ####
#+ include=TRUE, echo=FALSE
myRtools::write_RDS(object = tsi_comb,
                    file   = OUTPUT_NOAA_LAP  )

#+ include=FALSE


panderOptions("table.style", 'rmarkdown')
panderOptions("table.split.table", 100 )
panderOptions("table.continues", '')
panderOptions("graph.fontsize", 9)
panderOptions("table.alignment.default", "right")



#' ### Statistics on output data ###
#+ include=TRUE, echo=FALSE
pander(summary(tsi_comb, digits = 5))




#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
