#!/usr/bin/env Rscript
# /* Copyright (C) 2022-2023 Athanasios Natsis <natsisphysicist@gmail.com> */

#### Prepare TSI data from TSIS

#'
#' Process TSI from TSIS.
#'
#' This is used in operational, to fill the most recent days,
#' will be replaced by new NOAA data.
#'
#+ include=T, echo=F

## __ Set environment  ---------------------------------------------------------
rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
options("width" = 130)
tic <- Sys.time()
Script.Name <- "~/TSI/TSI_TSIS_LAP.R"

if(!interactive()) {
    pdf( file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink(file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".lock", Script.Name))), timeout = 0)
}

library(knitr)
library(pander)
library(data.table)

opts_chunk$set(comment    = NA )
opts_chunk$set(fig.width  = 8,
               fig.height = 5  )

source("~/TSI/DEFINITIONS.R")


## __ Check if need to run  ----------------------------------------------------
if (!file.exists(OUTPUT_TSIS_LAP) |
    file.mtime(ASTROPYdb)   > file.mtime(OUTPUT_TSIS_LAP) |
    file.mtime(OUTPUT_TSIS) > file.mtime(OUTPUT_TSIS_LAP) ) {
    cat("\nNew data to parse\n\n")
} else {
    stop("\nNO new data to parse\n\n")
}

## Load data
ASTROPY_data <- data.table(readRDS(ASTROPYdb))
TSIS_data    <- data.table(readRDS(OUTPUT_TSIS))

## we want the most recent data from TSIS
ASTROPY_data <- ASTROPY_data[ Date > TSI_START & Date > min(TSIS_data$Date)]
TSIS_data    <- TSIS_data[    Date > TSI_START ]

plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l",
     xlab = "", ylab = "Distance [au]")

#'
#' ## Interpolate TSI measurements to our dates
#'
#' Make functions from TSI measurements to match out data.
#' Interpolate between measurements only.
#'
#+ include=T, echo=F
tsi_fun <- approxfun(x      = TSIS_data$Date,
                     y      = TSIS_data$tsi_1au,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )

unc_fuc <- approxfun(x      = TSIS_data$Date,
                     y      = TSIS_data$measurement_uncertainty_1au,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )

#'
#' Interpolate the data, assuming that dates from Astropy are complete.
#'
#+ include=T, echo=F
tsi_all      <- tsi_fun(ASTROPY_data$Date)
unc_all      <- unc_fuc(ASTROPY_data$Date)

#'
#' Compute TSI on earth using TSI at 1 au
#'
#+ include=T, echo=F
tsi_astropy  <- tsi_all / (ASTROPY_data$Dist ^ 2)


##  Constructed TSI data for output  -------------------------------------------
tsi_comb <- data.frame(
    Date          = ASTROPY_data$Date, # Dates from SORCE extended to today
    sun_dist      = ASTROPY_data$Dist, # Astropy sun distance not optimal
    TSIextEARTH   = tsi_astropy,       # Original data with Astropy distance
    measur_error  = unc_all,           # Original data
    tsi_1au       = tsi_all            # TSI at 1 au
)
tsi_comb <- tsi_comb[ !is.na(tsi_comb$tsi_1au), ]


##  Output TSIS data for use  --------------------------------------------------
myRtools::write_RDS(object = tsi_comb,
                    file   = OUTPUT_TSIS_LAP  )


panderOptions("table.style", 'rmarkdown')
panderOptions("table.split.table", 100 )
panderOptions("table.continues", '')
panderOptions("graph.fontsize", 9)
panderOptions("table.alignment.default", "right")

#'
#' ### Statistics on output data
#'
#+ include=T, echo=F
pander(summary(tsi_comb, digits = 5))


#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
