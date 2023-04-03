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


####    Variables    ####
source("~/TSI/DEFINITIONS.R")


####    Get data    ####
system( paste0("curl \"", FROM_TSIS, "\" > ", DEST_TSIS) )



####    Parse data    ####
tsis_data <- fread(DEST_TSIS)

## fix names
names(tsis_data)[grep("time",names(tsis_data))]    <- "Date"
names(tsis_data)[grep( " \\(W/m\\^2\\)", names(tsis_data))] <-
    sub(  " \\(W/m\\^2\\)", "", grep( " \\(W/m\\^2\\)", names(tsis_data), value = T))


#' ignore zeros
tsis_data <- tsis_data[ tsi_1au >= 1 ]

setorder(tsis_data, Date)

tsis_data[ , provisional_flag := NULL ]

wecare    <- grep( "true_earth" ,names(tsis_data), value = T, invert = T)
tsis_data <- tsis_data[, ..wecare ]


hist( tsis_data$tsi_1au )

plot( tsis_data$Date, tsis_data$tsi_1au, pch = ".")

myRtools::write_RDS(object = tsis_data,
                    file   = OUTPUT_TSIS  )

tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
