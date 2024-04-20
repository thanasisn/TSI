#!/usr/bin/env Rscript
# /* Copyright (C) 2022-2023 Athanasios Natsis <natsisphysicist@gmail.com> */
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

#+ include=T, echo=F

## __ Set environment  ---------------------------------------------------------
rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
tic <- Sys.time()
Script.Name <- "~/TSI/prepare_TSI_LAP.R"

if (!interactive()) {
    pdf( file = paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink(file = paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".lock", Script.Name))))
}

library(pander    , quietly = TRUE, warn.conflicts = FALSE)
library(caTools   , quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(knitr     , quietly = TRUE, warn.conflicts = FALSE)
library(arrow, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)


opts_chunk$set(comment    = NA    )
opts_chunk$set(fig.width  = 8,
               fig.height = 5     )
opts_chunk$set(dev        = "png" )

source("~/TSI/DEFINITIONS.R")


## __ Load all data  -----------------------------------------------------------

#'
#' ## Available TSI data
#'
#+ include=T, echo=F
NOAA <- readRDS(OUTPUT_NOAA)
PMOD <- readRDS(OUTPUT_PMOD)
SORC <- readRDS(OUTPUT_SORCE)
TSIS <- readRDS(OUTPUT_TSIS)

ylim <- range(NOAA$TSI,  PMOD$tsi_1au, SORC$tsi_1au, TSIS$tsi_1au)
xlim <- range(NOAA$time, PMOD$Date, SORC$Date, TSIS$Date)
xlim[1] <- TSI_START

plot(  TSIS$Date, TSIS$tsi_1au, ylim = ylim, xlim = xlim, pch = ".",
       xlab = "", ylab = "TSI at 1au")
points(NOAA$time, NOAA$TSI,               pch = ".", col = 2 )
points(SORC$Date, SORC$tsi_1au,           pch = ".", col = 3 )
points(PMOD$Date, PMOD$tsi_1au,           pch = ".", col = 4 )
# points(PMOD$Date, PMOD$tsi_1au_old_VIRGO, pch = ".", col = 5 )

legend("bottom", ncol = 4,
       legend = c("TSIS", "NOAA", "SORC", "PMOD"),
       col    = c(1:4),
       pch    = 19, bty = "n", cex = .7)
rm(NOAA, PMOD, SORC, TSIS)



#'
#' ## Use TSI interpolated data from NOAA and TSIS
#'
#+ include=T, echo=F

## __  Load LAP Sun data  ------------------------------------------------------
NOAA    <- data.table(readRDS(OUTPUT_NOAA_LAP))
TSIS    <- data.table(readRDS(OUTPUT_TSIS_LAP))

names(NOAA)[names(NOAA) == "nominal_dates"] <- "Date"
names(NOAA)[names(NOAA) == "time"         ] <- "Date"
names(NOAA)[names(NOAA) == "TSIextEARTH"  ] <- "TSIextEARTH_NOAA"
names(TSIS)[names(TSIS) == "TSIextEARTH"  ] <- "TSIextEARTH_TSIS"
names(NOAA)[names(NOAA) == "measur_error" ] <- "measur_error_NOAA"
names(TSIS)[names(TSIS) == "measur_error" ] <- "measur_error_TSIS"
names(NOAA)[names(NOAA) == "tsi_1au"      ] <- "tsi_1au_NOAA"
names(TSIS)[names(TSIS) == "tsi_1au"      ] <- "tsi_1au_TSIS"


cat("\nNOAA TSI range:", format(range(NOAA$Date)), "\n")
cat("\nTSIS TSI range:", format(range(TSIS$Date)), "\n")


#'
#' Combine data and get mean difference for common dates
#'
#+ include=T, echo=F
tsi_merge <- merge(NOAA, TSIS, all = TRUE )
rm(NOAA, TSIS)

meandiff  <- mean(  tsi_merge$TSIextEARTH_NOAA - tsi_merge$TSIextEARTH_TSIS, na.rm = T)
medidiff  <- median(tsi_merge$TSIextEARTH_NOAA - tsi_merge$TSIextEARTH_TSIS, na.rm = T)

stopifnot(is.na(medidiff) == FALSE)
stopifnot(is.na(meandiff) == FALSE)

ylim <- range(tsi_merge$tsi_1au_NOAA, tsi_merge$tsi_1au_TSIS, na.rm = T)
plot(  tsi_merge$Date, tsi_merge$tsi_1au_NOAA, pch = ".", ylim = ylim, col = 2)
points(tsi_merge$Date, tsi_merge$tsi_1au_TSIS, pch = ".", col = 3)
title(main = "Original data")
cat("\n\n")

plot(  tsi_merge$Date, tsi_merge$tsi_1au_NOAA,            pch = ".", ylim = ylim, col = 2)
points(tsi_merge$Date, tsi_merge$tsi_1au_TSIS + medidiff, pch = ".", col = 3)
title(main = "Moved data with median differance")
cat("\n\n")

plot(  tsi_merge$Date, tsi_merge$tsi_1au_NOAA,            pch = ".", ylim = ylim, col = 2)
points(tsi_merge$Date, tsi_merge$tsi_1au_TSIS + meandiff, pch = ".", col = 3)
title(main = "Moved data with mean differance")
cat("\n\n")

# ylim <- range(tsi_merge$TSIextEARTH_NOAA, tsi_merge$TSIextEARTH_TSIS + meandiff, na.rm = T)
# plot(  tsi_merge$Date, tsi_merge$TSIextEARTH_NOAA,            pch = ".", ylim = ylim, col = 2)
# points(tsi_merge$Date, tsi_merge$TSIextEARTH_TSIS + meandiff, pch = ".", col = 3)
# title(main = "Moved data with mean difference")



##  Create composite output  ---------------------------------------------------
## move TSIS to NOAA values
tsi_merge$tsi_1au_TSIS     <- tsi_merge$tsi_1au_TSIS     + meandiff
tsi_merge$TSIextEARTH_TSIS <- tsi_merge$TSIextEARTH_TSIS + meandiff

## Complete missing from NOAA with TSIS
vec        <- !is.na(tsi_merge$TSIextEARTH_NOAA)
tsi_merge[ vec, TSIextEARTH_comb  := TSIextEARTH_NOAA ]
tsi_merge[ vec, measur_error_comb := measur_error_NOAA]
tsi_merge[ vec, tsi_1au_comb      := tsi_1au_NOAA     ]
tsi_merge[ vec, Source := "NOAA"                      ]
vec        <- is.na(tsi_merge$TSIextEARTH_comb)
tsi_merge[ vec, Source := "TSIS_adjusted"             ]
tsi_merge[ vec, TSIextEARTH_comb  := TSIextEARTH_TSIS ]
tsi_merge[ vec, tsi_1au_comb      := tsi_1au_TSIS     ]
tsi_merge[ vec, measur_error_comb := measur_error_TSIS]

## clean columns
tsi_merge[, TSIextEARTH_NOAA   := NULL]
tsi_merge[, measur_error_NOAA  := NULL]
tsi_merge[, tsi_1au_NOAA       := NULL]
tsi_merge[, TSIextEARTH_TSIS   := NULL]
tsi_merge[, measur_error_TSIS  := NULL]
tsi_merge[, tsi_1au_TSIS       := NULL]


cat("\n\n")
pander(summary(tsi_merge))
cat("\n\n")

xlim <- range(tsi_merge$Date,         na.rm = TRUE)
ylim <- range(tsi_merge$tsi_1au_comb, na.rm = TRUE)
plot(tsi_merge[Source == "NOAA", Date], tsi_merge[Source == "NOAA", tsi_1au_comb],
     col = 2,  pch = ".",
     xlim = xlim, ylim = ylim,
     xlab = "", ylab = "TSI composite")
points(tsi_merge[Source == "TSIS_adjusted", Date], tsi_merge[Source == "TSIS_adjusted", tsi_1au_comb],
       col = 3,  pch = ".")
title(main = "Composite TSI from NOAA and TSIS")
cat("\n\n")


xlim <- range(tsi_merge$Date,             na.rm = TRUE)
ylim <- range(tsi_merge$TSIextEARTH_comb, na.rm = TRUE)
plot(tsi_merge[Source == "NOAA", Date], tsi_merge[Source == "NOAA", TSIextEARTH_comb],
     col = 2,  pch = ".",
     xlim = xlim, ylim = ylim,
     xlab = "", ylab = "TOA TSI composite")
points(tsi_merge[Source == "TSIS_adjusted", Date], tsi_merge[Source == "TSIS_adjusted", TSIextEARTH_comb],
       col = 3,  pch = ".")
title(main = "Composite TSI at TOA from NOAA and TSIS")
cat("\n\n")


plot(tsi_merge$Date, tsi_merge$measur_error_comb, pch = ".")
cat("\n\n")


##  Output composite data for use  ---------------------------------------------
myRtools::write_RDS(object = tsi_merge,
                    file   = COMP_TSI)


## Parquet data set output  ----------------------------------------------------

# # tsi_merge <- readRDS(COMP_TSI)
#
# ## copied from BB not tried to test efficiency
# DB_compress_codec <- "brotli"
# DB_compress_level <- 5
#
# tsi_merge[, year := year(Date)]
#
# write_dataset(dataset           = tsi_merge,
#               path              = COMP_TSI_dataset,
#               format            = "parquet",
#               compression       = DB_compress_codec,
#               compression_level = DB_compress_level,
#               partitioning      = c("year"),
#               hive_style        = FALSE)
# cat("Written dataset: ", COMP_TSI_dataset, "\n")



#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
