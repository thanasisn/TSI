#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisthanasis@gmail.com> */

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
system( paste0("curl \"", FROM_SORCE, "\" > ", DEST_SORCE) )



####    Parse data    ####
SORCE_data <- fread(DEST_SORCE)

## fix names
names(SORCE_data)[grep("time",names(SORCE_data))]    <- "Date"
names(SORCE_data)[grep( " \\(W/m\\^2\\)", names(SORCE_data))] <-
    sub(  " \\(W/m\\^2\\)", "", grep( " \\(W/m\\^2\\)", names(SORCE_data), value = T))


#' ignore zeros
SORCE_data <- SORCE_data[ tsi_1au >= 1 ]

setorder(SORCE_data, Date)


wecare    <- grep( "true_earth" ,names(SORCE_data), value = T, invert = T)
SORCE_data <- SORCE_data[, ..wecare ]


hist( SORCE_data$tsi_1au )

plot( SORCE_data$Date, SORCE_data$tsi_1au, pch = ".")

myRtools::write_RDS(object = SORCE_data,
                    file   = OUTPUT_SORCE  )

tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
