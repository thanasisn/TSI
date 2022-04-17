#!/usr/bin/env Rscript
# /* Copyright (C) 2022 Athanasios Natsis <natsisthanasis@gmail.com> */

#' ---
#' title:  "TSI data preparation."
#' author: "Natsis Athanasios"
#' institute: "AUTH"
#' affiliation: "Laboratory of Atmospheric Physics"
#' abstract: "Read original total solar irradiance data from SORCE
#'            use sorce_tsi_L3_c06h file to create a R data.frame
#'            corresponding to the originals and sub-setting to CHP1 data.
#'            TSI measurements are extended to present day by using
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
# rm(list = (ls()[ls() != ""]))
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


#' #### Load Astropy prebuild sun distance data ####
ASTROPY_data <- readRDS(paste0(astropy_db,".Rds"))
#+ include=FALSE

####  Function to convert Julian date number to POSIXct  -----------------------
## Accepts a vector of jdn dates and return POSIXct vector
## it is not a general Julian date converter due to different
## Julian date specs
jd2posixct <- function(jdn) {
    times_hrs = (( jdn - 0.5 ) %% 1 ) * 24
    times_min = ( times_hrs    %% 1 ) * 60
    times_sec = ( times_min    %% 1 ) * 60
    str_dates = format.jd(jdn, "%F %T")
    str_times = sprintf("%02d:%02d:%02d", as.integer(times_hrs), as.integer(times_min), as.integer(times_sec))
    dates     = as.POSIXct(paste(str_dates, str_times),tz = "UTC")
    return(dates)
}


#+ include=FALSE
##    SORCE data preparation ####

#'
#' ## SORCE data preparation
#'
#' #### Get raw TSI data from file
#'
#' sorce has been discontinued !!
#'
sun_data  <- read.table(sorce_tsi, comment.char = ";")
#+ include=FALSE


## variables names from original input file
names(sun_data) <- c("nominal_date_yyyymmdd",
                     "nominal_date_jdn",
                     "avg_measurement_date_jdn",
                     "std_dev_measurement_date",
                     "tsi_1au",
                     "instrument_accuracy_1au",
                     "instrument_precision_1au",
                     "solar_standard_deviation_1au",
                     "measurement_uncertainty_1au",
                     "tsi_true_earth",
                     "instrument_accuracy_true_earth",
                     "instrument_precision_true_earth",
                     "solar_standard_deviation_true_earth",
                     "measurement_uncertainty_true_earth",
                     "provisional_flag")




####  Build a R data.frame  ----------------------------------------------------
## sorce
tsi_df <- data.frame( nominal_date_yyyymmdd               = sun_data$nominal_date_yyyymmdd ,
                      # nominal_date_jdn                    = jd2posixct(sun_data$nominal_date_jdn),
                      # avg_measurement_date_jdn            = jd2posixct(sun_data$avg_measurement_date_jdn),
                      std_dev_measurement_date            = sun_data$std_dev_measurement_date,
                      tsi_1au                             = sun_data$tsi_1au,
                      # instrument_accuracy_1au             = sun_data$instrument_accuracy_1au,
                      # instrument_precision_1au            = sun_data$instrument_precision_1au,
                      # solar_standard_deviation_1au        = sun_data$solar_standard_deviation_1au,
                      measurement_uncertainty_1au         = sun_data$measurement_uncertainty_1au,
                      tsi_true_earth                      = sun_data$tsi_true_earth,
                      instrument_accuracy_true_earth      = sun_data$instrument_accuracy_true_earth,
                      instrument_precision_true_earth     = sun_data$instrument_precision_true_earth,
                      solar_standard_deviation_true_earth = sun_data$solar_standard_deviation_true_earth,
                      measurement_uncertainty_true_earth  = sun_data$measurement_uncertainty_true_earth #,
                      # provisional_flag                    = sun_data$provisional_flag
)

strdt      <- paste(tsi_df$nominal_date_yyyymmdd)
dayseconds <- ( tsi_df$nominal_date_yyyymmdd %% 1 ) * 24 * 3600
hours      <-     dayseconds %/% 3600
minutes    <-   ( dayseconds - hours * 3600 ) %/% 60
seconds    <- ((( dayseconds - hours * 3600 ) %/% 60 ) * 60) %/% 60
times      <- sprintf("%02d:%02d:%02d", hours, minutes, seconds)

tsi_df$nominal_date_yyyymmdd <- strptime( paste(substr(strdt,1,8), times), "%Y%m%d %H:%M:%S" )


#' ##### Test if the sun is on (zero measurements may occur) #####
tsi_df  <- tsi_df[tsi_df$tsi_true_earth > 500,]


#+ include=FALSE
####  Subset to project dates  ####
# tsi_df  <- tsi_df[tsi_df$nominal_date_jdn > TSI_START,]
tsi_df          <- tsi_df[tsi_df$nominal_date_yyyymmdd > TSI_START,]

####  Remove invalid records
lower_tsi_limit <- min(mean(tsi_df$tsi_1au), median(tsi_df$tsi_1au)) * 0.70
tsi_df          <- tsi_df[ tsi_df$tsi_1au > lower_tsi_limit, ]



####  Extend measurements until today  ####
## use last measurement date and time step from origin
# timestep        = tsi_df$nominal_date_jdn[100] - tsi_df$nominal_date_jdn[99]
# last_measu_date = max(tsi_df$nominal_date_jdn) + timestep

#' #### Interpolate TSI measurements to our data ####
#' Make functions from TSI measurements to match out data.
#' Interpolate between measurements and extend tails as closest constant.
tsi_fun <- approxfun(x      = tsi_df$nominal_date_yyyymmdd,
                     y      = tsi_df$tsi_1au,
                     method = "linear",
                     rule   = 2:1,
                     f      = 0,
                     ties   = mean )

#+ include=FALSE
cat("Interpolate between measurements\n")

#+ include=TRUE
unc_fuc <- approxfun(x      = tsi_df$nominal_date_yyyymmdd,
                     y      = tsi_df$measurement_uncertainty_1au,
                     method = "linear",
                     rule   = 2:1,
                     f      = 0,
                     ties   = mean )

#' Interpolate the data, we have assumed that dates from Astropy are complete (we made them).
tsi_all      <- tsi_fun( ASTROPY_data$Date )  # SORCE dates to project dates
unc_all      <- unc_fuc( ASTROPY_data$Date )

#' Compute TSI on earth using TSI at 1 au
tsi_astropy  <- tsi_all / ( ASTROPY_data$Dist ^ 2 )

#+ include=FALSE
#### Constructed TSI data for output --------------------------------
tsi_comb <- data.frame(
                nominal_dates         = ASTROPY_data$Date, # Dates from SORCE extended to today
                sun_dist              = ASTROPY_data$Dist, # Astropy sun distance not optimal
                # tsi_true_earth_compex = tsi_astropy,       # Original data and extension with Astropy distance
                TSIextEARTH_comb      = tsi_astropy,       # Original data and extension with Astropy distance
                measur_error          = unc_all,           # Original data and extension of last value
                tsi_1au               = tsi_all            # TSI at 1 au
            )




#+ include=TRUE, echo=FALSE
par(mar = c(2,4,2,1))
#' ### Solar irradiance at 1 au ###
#+ include=TRUE, echo=FALSE
ylim = range(c( tsi_all + 0.1, tsi_all - 0.1 ), na.rm = T)
plot( ASTROPY_data$Date, tsi_all, "l", xlab = "", ylab = "SORCE TSI (Interpolated) watt/m^2", ylim = ylim)
lines(ASTROPY_data$Date, runmean(tsi_all, 15000), col = 5, lwd = 3)
qq <- quantile(tsi_all, na.rm = T)

abline( h = qq[3], col = 'orange',lwd = 3 )
abline( h = qq[2], col = 'green', lwd = 2 )
abline( h = qq[4], col = 'green', lwd = 2 )
abline( h = qq[1], col = 'red',   lwd = 1 )
abline( h = qq[5], col = 'red',   lwd = 1 )

text(ASTROPY_data$Date[1], y = qq[3],
     as.character(round(qq[3],1)), adj = c(0,1.5), col = 'orange', lwd = 4 )

text(ASTROPY_data$Date[1], y = qq[2],
     as.character(round(qq[2],1)), adj = c(0,1.5), col = 'green', lwd = 4 )
text(ASTROPY_data$Date[1], y = qq[4],
     as.character(round(qq[4],1)), adj = c(0,1.5), col = 'green', lwd = 4 )

text(ASTROPY_data$Date[1], y = qq[1],
     as.character(round(qq[1],1)), adj = c(0,1.5), col = 'red', lwd = 4 )
text(ASTROPY_data$Date[1], y = qq[5],
     as.character(round(qq[5],1)), adj = c(0,1.5), col = 'red', lwd = 4 )

pander(summary(tsi_all),digits = 8)

#' ### Solar irradiance at TOA
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, tsi_astropy,       "l", xlab = "", ylab = "SORCE TSI True Earth (Astropy) watt/m^2")
pander(summary(tsi_astropy),digits = 8)

#' ### Solar irradiance uncertainty
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, unc_all,           "l", xlab = "", ylab = "SORCE TSI Uncertainty (Interpolated) watt/m^2")
pander(summary(unc_all))

#' ### Earth - Sun distance
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, ASTROPY_data$Dist, "l", xlab = "", ylab = "Sun-Thessaloniki Distance AU")
pander(summary(ASTROPY_data$Dist))



#+ include=FALSE
## Export last day of SORCE measurements to file
LAST_TSI_DATE = max(tsi_df$nominal_date_yyyymmdd)
save( LAST_TSI_DATE, file = GLOBAL_file)

#' #### Output old data for use ####
#+ include=TRUE, echo=FALSE
myRtools::write_RDS(object = tsi_comb,
                    file   = tsi_Rds  )
#+ include=FALSE



panderOptions("table.style", 'rmarkdown')
panderOptions("table.split.table", 100 )
panderOptions("table.continues", '')
panderOptions("graph.fontsize", 9)
panderOptions("table.alignment.default", "right")


#' ### Statistics on input data ###
#+ include=TRUE, echo=FALSE
pander(summary(tsi_df,   digits = 5))

#' ### Statistics on output data ###
#+ include=TRUE, echo=FALSE
pander(summary(tsi_comb, digits = 5))

#' ### Raw Input data specifications ###

#' |     Field name                      | type |format | Col. #, description                                                     |
#' |-------------------------------------|------|-------|-------------------------------------------------------------------------|
#' | nominal_date_yyyymmdd               | R8   | f12.3 | Column  1: Nominal Data Time, YYYYMMDD                                  |
#' | nominal_date_jdn                    | R8   | f12.3 | Column  2: Nominal Data Time, Julian Day Number                         |
#' | avg_measurement_date_jdn            | R8   | f15.6 | Column  3: Average Data Time, Julian Day Number                         |
#' | std_dev_measurement_date            | R4   |  f7.4 | Column  4: Stdev of Average Data Time, days, 1 sigma                    |
#' | tsi_1au                             | R8   | f10.4 | Column  5: Total Solar Irradiance (TSI) at 1-AU, $W/m^2$                |
#' | instrument_accuracy_1au             | R4   | e10.3 | Column  6: Instrument Accuracy in 1-AU TSI, $W/m^2$, 1 sigma            |
#' | instrument_precision_1au            | R4   | e10.3 | Column  7: Instrument Precision in TSI at 1-AU, $W/m^2$, 1 sigma        |
#' | solar_standard_deviation_1au        | R4   | e10.3 | Column  8: Solar Standard Deviation in 1-AU TSI, $W/m^2$, 1 sigma       |
#' | measurement_uncertainty_1au         | R4   | e10.3 | Column  9: Total Uncertainty in TSI at 1-AU, $W/m^2$, 1 sigma           |
#' | tsi_true_earth                      | R8   | f10.4 | Column 10: Total Solar Irradiance at Earth distance, $W/m^2$            |
#' | instrument_accuracy_true_earth      | R4   | e10.3 | Column 11: Instrument Accuracy at Earth distance, $W/m^2$, 1 sigma      |
#' | instrument_precision_true_earth     | R4   | e10.3 | Column 12: Instrument Precision at Earth distance, $W/m^2$, 1 sigma     |
#' | solar_standard_deviation_true_earth | R4   | e10.3 | Column 13: Solar Standard Deviation in TSI at Earth, $W/m^2$, 1 sigma   |
#' | measurement_uncertainty_true_earth  | R4   | e10.3 | Column 14: Total Uncertainty in TSI at Earth distance, $W/m^2$, 1 sigma |
#' | provisional_flag                    | I2   |    i2 | Column 15: Provisional Flag, 1=provisional data, 0=final data           |



#+ include=FALSE
#### lisird/latis data preparation ####

#'
#' ## lisird/latis data preparation
#'
#' #### Get raw TSI data from file
#'
#'
## from lasp.colorado.edu/lisird/latis
tsis_data <- data.table::fread(tsis_tsi)

#' ignore zeros
tsis_data <- tsis_data[ `tsi_1au (W/m^2)` >= 1 ]
#+ include=FALSE




names(tsis_data)[grep("time",names(tsis_data))]    <- "Date"
names(tsis_data)[grep("tsi_1au",names(tsis_data))] <- "TSI_lasp"


tsis_data <- tsis_data[ Date >= TSI_START]



#' #### Interpolate TSI measurements to our data ####
#' Make functions from TSI measurements to match out data.
#' Interpolate between measurements and extend tails as closest constant.
tsi_fun2 <- approxfun(x      = tsis_data$Date,
                      y      = tsis_data$TSI_lasp,
                      method = "linear",
                      rule   = 1:2,
                      f      = 0,
                      ties   = mean )
#+ include=FALSE
cat("Interpolate between measurements and extend tail with closest constant\n")

#+ include=TRUE
unc_fuc2 <- approxfun(x      = tsis_data$Date,
                      y      = tsis_data$`measurement_uncertainty_1au (W/m^2)`,
                      method = "linear",
                      rule   = 1:2,
                      f      = 0,
                      ties   = mean )

#' Interpolate the data, we have assumed that dates from Astropy are complete (we made them).
tsis_all      <- tsi_fun2( ASTROPY_data$Date )  # SORCE dates to project dates
uncs_all      <- unc_fuc2( ASTROPY_data$Date )

#' Compute TSI on earth using TSI at 1 au
tsis_astropy  <- tsis_all / ( ASTROPY_data$Dist ^ 2 )

#+ include=FALSE
#### Constructed TSI data for output --------------------------------
tsis_comb <- data.frame(
    Date                  = ASTROPY_data$Date, # Dates from SORCE extended to today
    sun_dist              = ASTROPY_data$Dist, # Astropy sun distance not optimal
    TSIextEARTH_comb      = tsis_astropy,       # Original data and extension with Astropy distance
    measur_error          = uncs_all,           # Original data and extension of last value
    tsi_1au               = tsis_all            # TSI at 1 au
)




#+ include=TRUE, echo=FALSE
par(mar = c(2,4,2,1))
#' ### Solar irradiance at 1 au ###
#+ include=TRUE, echo=FALSE
ylim = range(c( tsis_all + 0.1, tsis_all - 0.1 ), na.rm = T)
plot( ASTROPY_data$Date, tsis_all, "l", xlab = "", ylab = "LASP TSI (Interpolated) watt/m^2", ylim = ylim)
lines(ASTROPY_data$Date, runmean(tsis_all, 15000), col = 5, lwd = 3)
qq <- quantile(tsis_all, na.rm = T)

abline( h = qq[3], col = 'orange',lwd = 3 )
abline( h = qq[2], col = 'green', lwd = 2 )
abline( h = qq[4], col = 'green', lwd = 2 )
abline( h = qq[1], col = 'red',   lwd = 1 )
abline( h = qq[5], col = 'red',   lwd = 1 )

text(ASTROPY_data$Date[1], y = qq[3],
     as.character(round(qq[3],1)), adj = c(0,1.5), col = 'orange', lwd = 4 )

text(ASTROPY_data$Date[1], y = qq[2],
     as.character(round(qq[2],1)), adj = c(0,1.5), col = 'green', lwd = 4 )
text(ASTROPY_data$Date[1], y = qq[4],
     as.character(round(qq[4],1)), adj = c(0,1.5), col = 'green', lwd = 4 )

text(ASTROPY_data$Date[1], y = qq[1],
     as.character(round(qq[1],1)), adj = c(0,1.5), col = 'red', lwd = 4 )
text(ASTROPY_data$Date[1], y = qq[5],
     as.character(round(qq[5],1)), adj = c(0,1.5), col = 'red', lwd = 4 )

pander(summary(tsis_all),digits = 8)


#' ### Solar irradiance at TOA
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, tsis_astropy,   "l", xlab = "", ylab = "LASP TSI True Earth (Astropy) watt/m^2")
pander(summary(tsi_astropy),digits = 8)

#' ### Solar irradiance uncertainty
#+ include=TRUE, echo=FALSE
plot(ASTROPY_data$Date, uncs_all,       "l", xlab = "", ylab = "LASP TSI Uncertainty (Interpolated) watt/m^2")
pander(summary(unc_all))





#' ## Combine time series
#+ include=FALSE

names(tsi_comb)[grep("dates",names(tsi_comb))]              <- "Date"
names(tsi_comb)[grep("TSIextEARTH_comb",names(tsi_comb))]   <- "TSIextEARTH_comb_SORCE"
names(tsi_comb)[grep("measur_error",names(tsi_comb))]       <- "measur_error_SORCE"
names(tsi_comb)[grep("tsi_1au",names(tsi_comb))]            <- "tsi_1au_SORCE"

names(tsis_comb)[grep("TSIextEARTH_comb",names(tsis_comb))] <- "TSIextEARTH_comb_LASP"
names(tsis_comb)[grep("measur_error",names(tsis_comb))]     <- "measur_error_LASP"
names(tsis_comb)[grep("tsi_1au",names(tsis_comb))]          <- "tsi_1au_LASP"



#' combine data and get mean diff
tsi_merge <- data.table(merge(tsi_comb, tsis_comb, all = T))

meandiff  <- mean(  tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP, na.rm = T)
meaddiff  <- median(tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP, na.rm = T)

ylim <- range(tsi_merge$tsi_1au_SORCE, tsi_merge$tsi_1au_LASP, na.rm = T)
plot(  tsi_merge$Date, tsi_merge$tsi_1au_SORCE, "l", ylim = ylim)
lines( tsi_merge$Date, tsi_merge$tsi_1au_LASP  )


plot( tsi_merge$TSIextEARTH_comb_SORCE / tsi_merge$TSIextEARTH_comb_LASP )
plot( tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP )

plot( - meandiff + tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP )
plot( - meaddiff + tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP )

plot( (tsi_merge$TSIextEARTH_comb_SORCE - meandiff) / tsi_merge$TSIextEARTH_comb_LASP )
plot( (tsi_merge$TSIextEARTH_comb_SORCE - meaddiff) / tsi_merge$TSIextEARTH_comb_LASP )

plot( (tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP)/tsi_merge$TSIextEARTH_comb_LASP  )

summary( tsi_merge$TSIextEARTH_comb_SORCE / tsi_merge$TSIextEARTH_comb_LASP )
summary( tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP )
summary((tsi_merge$TSIextEARTH_comb_SORCE - tsi_merge$TSIextEARTH_comb_LASP)/tsi_merge$TSIextEARTH_comb_LASP )


#' apply correction on SORCE data
tsi_merge$TSIextEARTH_comb_SORCE <- tsi_merge$TSIextEARTH_comb_SORCE - meandiff
tsi_merge$tsi_1au_SORCE          <- tsi_merge$tsi_1au_SORCE          - meandiff


#' fill with both data
vec        <- !is.na(tsi_merge$TSIextEARTH_comb_LASP)
tsi_merge[ vec, TSIextEARTH_comb  := TSIextEARTH_comb_LASP ]
tsi_merge[ vec, measur_error_comb := measur_error_LASP     ]
tsi_merge[ vec, tsi_1au_comb      := tsi_1au_LASP          ]
tsi_merge[ vec, Source := "LASP"  ]
vec        <- is.na(tsi_merge$TSIextEARTH_comb)
tsi_merge[ vec, Source := "SORCE" ]
tsi_merge[ vec, TSIextEARTH_comb  := TSIextEARTH_comb_SORCE ]
tsi_merge[ vec, tsi_1au_comb      := tsi_1au_SORCE ]
tsi_merge[ vec, measur_error_comb := measur_error_SORCE ]


#+ include=FALSE
tsi_merge[, TSIextEARTH_comb_LASP  := NULL ]
tsi_merge[, measur_error_LASP      := NULL ]
tsi_merge[, tsi_1au_LASP           := NULL ]
tsi_merge[, TSIextEARTH_comb_SORCE := NULL ]
tsi_merge[, measur_error_SORCE     := NULL ]
tsi_merge[, tsi_1au_SORCE          := NULL ]



#+ include=TRUE, echo=FALSE
par(mar = c(2,4,2,1))
#' ### Solar irradiance at 1 au ###
#+ include=TRUE, echo=FALSE
ylim = range(c( tsi_merge$tsi_1au_comb + 0.1, tsi_merge$tsi_1au_comb - 0.1 ), na.rm = T)
plot( tsi_merge$Date, tsi_merge$tsi_1au_comb, "l", xlab = "", ylab = "LASP TSI (Interpolated) watt/m^2", ylim = ylim)
lines(tsi_merge$Date, runmean(tsi_merge$tsi_1au_comb, 15000), col = 5, lwd = 3)
qq <- quantile(tsi_merge$tsi_1au_comb, na.rm = T)

abline( h = qq[3], col = 'orange',lwd = 3 )
abline( h = qq[2], col = 'green', lwd = 2 )
abline( h = qq[4], col = 'green', lwd = 2 )
abline( h = qq[1], col = 'red',   lwd = 1 )
abline( h = qq[5], col = 'red',   lwd = 1 )

text(tsi_merge$Date[1], y = qq[3],
     as.character(round(qq[3],1)), adj = c(0,1.5), col = 'orange', lwd = 4 )

text(tsi_merge$Date[1], y = qq[2],
     as.character(round(qq[2],1)), adj = c(0,1.5), col = 'green', lwd = 4 )
text(tsi_merge$Date[1], y = qq[4],
     as.character(round(qq[4],1)), adj = c(0,1.5), col = 'green', lwd = 4 )

text(tsi_merge$Date[1], y = qq[1],
     as.character(round(qq[1],1)), adj = c(0,1.5), col = 'red', lwd = 4 )
text(tsi_merge$Date[1], y = qq[5],
     as.character(round(qq[5],1)), adj = c(0,1.5), col = 'red', lwd = 4 )

pander(summary(tsi_merge$TSIextEARTH_comb),digits = 8)


#' ### Solar irradiance at TOA
#+ include=TRUE, echo=FALSE
plot(tsi_merge$Date,  tsi_merge$sun_dist,  "l", xlab = "", ylab = "LASP TSI True Earth (Astropy) watt/m^2")
pander(summary(tsi_astropy),digits = 8)

#' ### Solar irradiance uncertainty
#+ include=TRUE, echo=FALSE
plot(tsi_merge$Date, tsi_merge$measur_error_comb,      "l", xlab = "", ylab = "LASP TSI Uncertainty (Interpolated) watt/m^2")
pander(summary(unc_all))




##### TODO plot line by source

#' #### Output data for use ####
#+ include=TRUE, echo=FALSE
myRtools::write_RDS(object = tsi_merge,
                    file   = tsi_comp_out  )
#+ include=FALSE




#' **END**
#+ include=T, echo=F
tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
