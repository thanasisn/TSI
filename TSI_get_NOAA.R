#!/usr/bin/env Rscript
# /* Copyright (C) 2022-2023 Athanasios Natsis <natsisphysicist@gmail.com> */

#### Prepare TSI data from NOAA

#'
#' Get and parse TSI from NOAA.
#'
#' Data is stored asis.
#'
#' This is used in operational as the main TSI data.
#'

## __ Set environment  ---------------------------------------------------------
rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
options("width" = 130)
tic <- Sys.time()
Script.Name <- "~/TSI/TSI_get_NOAA.R"

if (!interactive()) {
    pdf( file = paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink(file = paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".lock", Script.Name))), timeout = 0)
}

library(data.table)
library(RNetCDF)

source("~/TSI/DEFINITIONS.R")


##  Get data  ------------------------------------------------------------------
system(paste("wget -N -r -np -nd -nH -A .nc -P", DEST_NOAA, FROM_NOAA))


## __ Check if we have to parse  -----------------------------------------------
ncfiles <- list.files(path       = DEST_NOAA,
                      pattern    = "_daily_s.*.nc",
                      recursive  = T,
                      full.names = T )

if (!file.exists(OUTPUT_NOAA) | file.mtime(OUTPUT_NOAA) < max(file.mtime(ncfiles))) {
    cat("\nNew data to parse!\n")
} else {
    stop("NO new data to parse!\n")
}


##  Parse data  ----------------------------------------------------------------
gather <- data.table()
for (af in ncfiles) {
    anc  <- open.nc(af, write  = FALSE)
    data <- read.nc(anc, unpack = TRUE )
    # print.nc(anc)

    data <- data.table( TSI      = as.vector(data$TSI),
                        time     = as.vector(data$time),
                        TSI_UNC  = as.vector(data$TSI_UNC),
                        time_low = as.vector(data$time_bnds[1,]),
                        time_upp = as.vector(data$time_bnds[2,]) )
    ## check sanity
    dateorigin <- "1610-01-01 00:00:00"
    stopifnot(grepl(dateorigin, att.get.nc(anc, "time",      "units")))
    stopifnot(grepl(dateorigin, att.get.nc(anc, "time",      "units")))
    stopifnot(grepl(dateorigin, att.get.nc(anc, "time_bnds", "units")))
    ## format data
    data$time     <- as.Date(    data$time,     origin = dateorigin)
    data$time     <- as.POSIXct( data$time ) + 12 * 3600
    data$time_low <- as.Date(    data$time_low, origin = dateorigin)
    data$time_upp <- as.Date(    data$time_upp, origin = dateorigin)
    data$file     <- af

    gather <- rbind(gather, data)
}
setorder(gather, time)

## __ Inspect data  ------------------------------------------------------------
hist(gather$TSI)

plot(gather$time, gather$TSI, pch = ".",
     xlab = "", ylab = "NOAA TSI")

cat("\nNOAA TSI range:", format(range(gather$time)), "\n")


##  Save TSI data  -------------------------------------------------------------
myRtools::write_RDS(object = gather,
                    file   = OUTPUT_NOAA)

tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
