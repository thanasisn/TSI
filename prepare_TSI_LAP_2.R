#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisthanasis@gmail.com> */

#' ---
#' title:  "TSI data preparation."
#' author: "Natsis Athanasios"
#' institute: "AUTH"
#' affiliation: "Laboratory of Atmospheric Physics"
#' abstract: "Read original total solar irradiance data from SORCE
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
####  Clean environment  -------------------------------------------------------
rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
tic <- Sys.time()
Script.Name <- tryCatch({ funr::sys.script() },
                        error = function(e) { cat(paste("\nUnresolved script name: ", e),"\n")
                            return("Undefined R script name!!") })
if(!interactive()) {
    pdf( file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink(file = paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/", basename(sub("\\.R$",".lock", Script.Name))), timeout = 0)
}




library(knitr)
opts_chunk$set( comment    = NA    )
opts_chunk$set( fig.width  = 8,
                fig.height = 5     )
opts_chunk$set( dev        = "png" )

library(pander)    # just to plot tables
library(caTools)
library(data.table)


####    Variables    ####
source("~/TSI/DEFINITIONS.R")


NOAA <- readRDS(OUTPUT_NOAA)
PMOD <- readRDS(OUTPUT_PMOD)
SORC <- readRDS(OUTPUT_SORCE)
TSIS <- readRDS(OUTPUT_TSIS)



ylim <- range(NOAA$TSI,  PMOD$tsi_1au, SORC$tsi_1au, TSIS$tsi_1au)
xlim <- range(NOAA$time, PMOD$Date, SORC$Date, TSIS$Date)
xlim[1] <- TSI_START


plot(  TSIS$Date, TSIS$tsi_1au, ylim = ylim, xlim = xlim, pch = ".")
points(NOAA$time, NOAA$TSI,               pch = ".", col = 2 )
points(SORC$Date, SORC$tsi_1au,           pch = ".", col = 3 )
points(PMOD$Date, PMOD$tsi_1au,           pch = ".", col = 4 )
# points(PMOD$Date, PMOD$tsi_1au_old_VIRGO, pch = ".", col = 5 )


## USE TSIS and NOAA



####    Load LAP Sun data    ###################################################
NOAA         <- data.table(readRDS(OUTPUT_NOAA_LAP))
TSIS         <- data.table(readRDS(OUTPUT_TSIS_LAP))
ASTROPY_data <- data.table(readRDS(ASTROPYdb))


names(NOAA)[names(NOAA) == "time"] <- "Date"

NOAA         <- NOAA[ Date >= TSI_START, ]
TSIS         <- TSIS[ Date >= TSI_START, ]
ASTROPY_data <- ASTROPY_data[ Date >= TSI_START, ]


names(NOAA)[names(NOAA) == "TSIextEARTH" ] <- "TSIextEARTH_NOAA"
names(TSIS)[names(TSIS) == "TSIextEARTH" ] <- "TSIextEARTH_TSIS"
names(NOAA)[names(NOAA) == "measur_error"] <- "measur_error_NOAA"
names(TSIS)[names(TSIS) == "measur_error"] <- "measur_error_TSIS"
names(NOAA)[names(NOAA) == "tsi_1au"     ] <- "tsi_1au_NOAA"
names(TSIS)[names(TSIS) == "tsi_1au"     ] <- "tsi_1au_TSIS"


#' combine data and get mean diff
tsi_merge <- merge(NOAA, TSIS, all = T )
rm(NOAA, TSIS, PMOD, SORC)



meandiff  <- mean(  tsi_merge$TSIextEARTH_NOAA - tsi_merge$TSIextEARTH_TSIS, na.rm = T)
meaddiff  <- median(tsi_merge$TSIextEARTH_NOAA - tsi_merge$TSIextEARTH_TSIS, na.rm = T)

ylim <- range(tsi_merge$tsi_1au_NOAA, tsi_merge$tsi_1au_TSIS, na.rm = T)
plot(  tsi_merge$Date, tsi_merge$tsi_1au_NOAA, pch = ".", ylim = ylim, col = 2)
points(tsi_merge$Date, tsi_merge$tsi_1au_TSIS, pch = ".", col = 3)
title(main = "Original data")


plot(  tsi_merge$Date, tsi_merge$tsi_1au_NOAA,            pch = ".", ylim = ylim, col = 2)
points(tsi_merge$Date, tsi_merge$tsi_1au_TSIS + meandiff, pch = ".", col = 3)
title(main = "Moved data 1au")


ylim <- range(tsi_merge$TSIextEARTH_NOAA, tsi_merge$TSIextEARTH_TSIS + meandiff, na.rm = T)
plot(  tsi_merge$Date, tsi_merge$TSIextEARTH_NOAA,            pch = ".", ylim = ylim, col = 2)
points(tsi_merge$Date, tsi_merge$TSIextEARTH_TSIS + meandiff, pch = ".", col = 3)
title(main = "Moved data")



####   Create composite output   ###############################################
tsi_merge$tsi_1au_TSIS     <- tsi_merge$tsi_1au_TSIS     + meandiff
tsi_merge$TSIextEARTH_TSIS <- tsi_merge$TSIextEARTH_TSIS + meandiff



#' Fill with both data
vec        <- !is.na(tsi_merge$TSIextEARTH_NOAA)
tsi_merge[ vec, TSIextEARTH_comb  := TSIextEARTH_NOAA  ]
tsi_merge[ vec, measur_error_comb := measur_error_NOAA ]
tsi_merge[ vec, tsi_1au_comb      := tsi_1au_NOAA      ]
tsi_merge[ vec, Source := "NOAA"                       ]
vec        <- is.na(tsi_merge$TSIextEARTH_comb)
tsi_merge[ vec, Source := "TSIS_adjusted"              ]
tsi_merge[ vec, TSIextEARTH_comb  := TSIextEARTH_TSIS  ]
tsi_merge[ vec, tsi_1au_comb      := tsi_1au_TSIS      ]
tsi_merge[ vec, measur_error_comb := measur_error_TSIS ]


#+ include=FALSE
tsi_merge[, TSIextEARTH_NOAA   := NULL ]
tsi_merge[, measur_error_NOAA  := NULL ]
tsi_merge[, tsi_1au_NOAA       := NULL ]
tsi_merge[, TSIextEARTH_TSIS   := NULL ]
tsi_merge[, measur_error_TSIS  := NULL ]
tsi_merge[, tsi_1au_TSIS       := NULL ]


pander(summary(tsi_merge))

## TODO fill last data

plot(tsi_merge$Date, tsi_merge$sun_dist,          pch = ".")
plot(tsi_merge$Date, tsi_merge$TSIextEARTH_comb,  pch = ".")
plot(tsi_merge$Date, tsi_merge$measur_error_comb, pch = ".")
plot(tsi_merge$Date, tsi_merge$tsi_1au_comb,      pch = ".")

#' #### Output data for use ####
#+ include=TRUE, echo=FALSE
myRtools::write_RDS(object = tsi_merge,
                    file   = COMP_TSI  )
#+ include=FALSE


#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
