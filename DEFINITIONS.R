#'
#' Here are some global or common variables for this project.
#'

TSI_START   <- as.POSIXct("1993-01-01 00:00", tz = "UTC")
ASTROPYdb   <- "/home/athan/DATA_RAW/SUN/Astropy_LAP.Rds"

####    Modeled TSI from NOAA    ###############################################

FROM_NOAA       <- "https://www.ncei.noaa.gov/data/total-solar-irradiance/access/daily/"
DEST_NOAA       <- "~/DATA/SUN/TSI_model_NOAA/"
OUTPUT_NOAA     <- "~/DATA/SUN/TSI_model_NOAA.Rds"
OUTPUT_NOAA_LAP <- "~/DATA/SUN/TSI_model_NOAA_LAP.Rds"







