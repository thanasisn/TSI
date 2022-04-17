#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisthanasis@gmail.com> */


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
library(data.table)


####    Variables    ####
source("~/TSI/DEFINITIONS.R")


####    Check if need to run    ####


if ( !file.exists(OUTPUT_PMOD_LAP) |
     file.mtime(ASTROPYdb)   > file.mtime(OUTPUT_PMOD_LAP) |
     file.mtime(OUTPUT_PMOD) > file.mtime(OUTPUT_PMOD_LAP) ) {
    cat("New data to parse\n\n")
} else {
    cat("NO new data to parse\n\n")
    stop("NO new data to parse\n\n")
}



## Load data
ASTROPY_data <- data.table(readRDS(ASTROPYdb))
PMOD_data    <- data.table(readRDS(OUTPUT_PMOD))

ASTROPY_data <- ASTROPY_data[ Date > TSI_START ]
PMOD_data    <- PMOD_data[    Date > TSI_START ]


plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l")


#' #### Interpolate TSI measurements to our dates ####
#' Make functions from TSI measurements to match out data.
#' Interpolate between measurements only.
tsi_fun <- approxfun(x      = PMOD_data$Date,
                     y      = PMOD_data$tsi_1au,
                     method = "linear",
                     rule   = 1,
                     ties   = mean )

tsi_fun2 <- approxfun(x      = PMOD_data$Date,
                      y      = PMOD_data$tsi_1au_old_VIRGO,
                      method = "linear",
                      rule   = 1,
                      ties   = mean )




#' Interpolate the data, we have assumed that dates from Astropy are complete.
tsi_all           <- tsi_fun(  ASTROPY_data$Date )
tsi_all_old_VIRGO <- tsi_fun2( ASTROPY_data$Date )

#' Compute TSI on earth using TSI at 1 au
tsi_astropy            <- tsi_all           / ( ASTROPY_data$Dist ^ 2 )
tsi_astropy_old_VIRGO  <- tsi_all_old_VIRGO / ( ASTROPY_data$Dist ^ 2 )


#### Constructed TSI data for output --------------------------------
tsi_comb <- data.frame(
    Date                  = ASTROPY_data$Date,     # Dates from SORCE extended to today
    sun_dist              = ASTROPY_data$Dist,     # Astropy sun distance not optimal
    TSIextEARTH           = tsi_astropy,           # Original data with Astropy distance
    TSIextEARTH_old_VIRGO = tsi_astropy_old_VIRGO, # Original data with Astropy distance
    tsi_1au               = tsi_all                # TSI at 1 au
)
tsi_comb <- tsi_comb[ !is.na(tsi_comb$tsi_1au), ]


#' #### Output data for use ####
myRtools::write_RDS(object = tsi_comb,
                    file   = OUTPUT_PMOD_LAP  )


panderOptions("table.style", 'rmarkdown')
panderOptions("table.split.table", 100 )
panderOptions("table.continues", '')
panderOptions("graph.fontsize", 9)
panderOptions("table.alignment.default", "right")



#' ### Statistics on output data ###
pander(summary(tsi_comb, digits = 5))


#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
