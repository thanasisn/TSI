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
system(paste("wget -N -r -np -nd -nH -O", DEST_PMOD, FROM_PMOD))


####    Parse data    ####
PMOD_data <- read.table(file         = DEST_PMOD,
                        comment.char = ";",
                        skip         = 1 )

plot(PMOD_data$V2)


stop()
paste()
strptime(PMOD_data$V1,"%y%m%d")

stop()
PMOD_data$Date <- as.POSIXct( strptime(PMOD_data$V1,"%y%m%d") ) + 12*3600
PMOD_data <- data.table(PMOD_data)
PMOD_data[ , V1 := NULL ]
PMOD_data[ , V2 := NULL ]
PMOD_data <- PMOD_data[ V3 > 100 & V4 > 100  ]


names(PMOD_data)[names(PMOD_data) == "V3"] <- "tsi_1au"
names(PMOD_data)[names(PMOD_data) == "V4"] <- "tsi_1au_old_VIRGO"

setorder(PMOD_data, Date)




hist( PMOD_data$tsi_1au )

plot( PMOD_data$Date, PMOD_data$tsi_1au, pch = ".")

myRtools::write_RDS(object = PMOD_data,
                    file   = OUTPUT_PMOD  )

tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
