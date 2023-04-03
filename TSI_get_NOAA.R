#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisphysicist@gmail.com> */

#### Prepare TSI data form NOAA


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


library(data.table)
library(RNetCDF)

####    Variables    ####
source("~/TSI/DEFINITIONS.R")


####    Get data    ####
system(paste("wget -N -r -np -nd -nH -A .nc -P", DEST_NOAA, FROM_NOAA))


####    Check if we have to parse    ####
ncfiles <- list.files(path       = DEST_NOAA,
                      pattern    = "*.nc",
                      recursive  = T,
                      full.names = T )

if (!file.exists(OUTPUT_NOAA) | file.mtime(OUTPUT_NOAA) < max(file.mtime(ncfiles))) {
    cat("New data to parse!\n")
} else {
    cat("NO new data to parse!\n")
    stop("NO new data to parse!\n")
}


####    Parse data    ####
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

    dateorigin <- "1610-01-01 00:00:00"
    stopifnot( grepl( dateorigin, att.get.nc(anc, "time",      "units") ) )
    stopifnot( grepl( dateorigin, att.get.nc(anc, "time",      "units") ) )
    stopifnot( grepl( dateorigin, att.get.nc(anc, "time_bnds", "units") ) )

    data$time     <- as.Date(    data$time,     origin = "1610-01-01 00:00:00" )
    data$time     <- as.POSIXct( data$time ) + 12 * 3600
    data$time_low <- as.Date(    data$time_low, origin = "1610-01-01 00:00:00" )
    data$time_upp <- as.Date(    data$time_upp, origin = "1610-01-01 00:00:00" )
    data$file     <- af

    gather <- rbind(gather, data)
}
setorder(gather, time)

hist( gather$TSI )

plot(gather$time, gather$TSI, pch = ".")

myRtools::write_RDS(object = gather,
                    file   = OUTPUT_NOAA  )

tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
