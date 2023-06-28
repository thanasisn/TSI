#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisphysicist@gmail.com> */

#### Prepare TSI data

# https://www.pmodwrc.ch/en/research-development/solar-physics/tsi-composite/


rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
options("width" = 130)
tic <- Sys.time()
Script.Name <- "~/TSI/retreive_data/TSI_get_PMOD.R"

if (!interactive()) {
    pdf( file = paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".pdf", Script.Name))))
    sink(file = paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".out", Script.Name))), split = TRUE)
    filelock::lock(paste0("~/TSI/REPORTS/RUNTIME/", basename(sub("\\.R$",".lock", Script.Name))), timeout = 0)
}


library(data.table)


##  Variables  -----------------------------------------------------------------
source("~/TSI/DEFINITIONS.R")


##  Get data  ------------------------------------------------------------------
system(paste("wget -N -r -np -nd -nH -O", DEST_PMOD, FROM_PMOD))


##  Parse data  ----------------------------------------------------------------
PMOD_data <- read.table(file         = DEST_PMOD,
                        comment.char = ";", colClasses = c("character"),
                        skip         = 1 )

## Convert variables
PMOD_data$Date <- as.POSIXct(strptime(PMOD_data$V1, "%y%m%d")) + 12*3600
PMOD_data$V1   <- NULL
PMOD_data$V2   <- NULL
PMOD_data$V3   <- as.numeric(PMOD_data$V3)
PMOD_data$V4   <- as.numeric(PMOD_data$V4)
PMOD_data      <- data.table(PMOD_data)
## Filter some bad values
PMOD_data      <- PMOD_data[ V3 > 100 & V4 > 100  ]

## Set names
names(PMOD_data)[names(PMOD_data) == "V3"] <- "tsi_1au"
names(PMOD_data)[names(PMOD_data) == "V4"] <- "tsi_1au_old_VIRGO"
setorder(PMOD_data, Date)

## Save data
myRtools::write_RDS(object = PMOD_data,
                    file   = OUTPUT_PMOD)


##  Do some data investigation  ------------------------------------------------

hist( PMOD_data$tsi_1au )

ylim <- range(PMOD_data$tsi_1au, PMOD_data$tsi_1au_old_VIRGO)

plot(  PMOD_data$Date, PMOD_data$tsi_1au,           pch = ".", ylim = ylim)
points(PMOD_data$Date, PMOD_data$tsi_1au_old_VIRGO, pch = ".", col = "yellow")


tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
