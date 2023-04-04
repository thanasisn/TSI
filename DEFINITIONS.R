
#'
#' Global or common variables for this project.
#'

TSI_START        <- as.POSIXct("1993-01-01 00:00", tz = "UTC")
ASTROPYdb        <- "/home/athan/DATA_RAW/SUN/Astropy_LAP.Rds"

## TSI data to use by other ---------------------------------------------------
COMP_TSI         <- "~/DATA/SUN/TSI_COMPOSITE.Rds"


##  NOAA TSI  -----------------------------------------------------------------
FROM_NOAA        <- "https://www.ncei.noaa.gov/data/total-solar-irradiance/access/daily/"
DEST_NOAA        <- "~/DATA/SUN/TSI_model_NOAA/"
OUTPUT_NOAA      <- "~/DATA/SUN/TSI_model_NOAA.Rds"
OUTPUT_NOAA_LAP  <- "~/DATA/SUN/TSI_model_NOAA_LAP.Rds"


##  TSIS TSI from LISIRD  -----------------------------------------------------
FROM_TSIS        <- "https://lasp.colorado.edu/lisird/latis/dap/tsis_tsi_6hr.csv?&format_time(yyyy-MM-dd'T'HH:mm:ss.SSSZ)"
DEST_TSIS        <- "~/DATA/SUN/TSI_tsis_6hr.csv"
OUTPUT_TSIS      <- "~/DATA/SUN/TSI_tsis_6hr.Rds"
OUTPUT_TSIS_LAP  <- "~/DATA/SUN/TSI_tsis_6hr_LAP.Rds"


##  PMOD TSI from pmodwrc.ch  ------------------------------------------------
FROM_PMOD        <- "ftp://ftp.pmodwrc.ch/pub/data/irradiance/composite/DataPlots/ext_composite_42_65_1805.dat"
DEST_PMOD        <- "~/DATA/SUN/ext_composite_42_65_1805.dat"
OUTPUT_PMOD      <- "~/DATA/SUN/TSI_PMOD_ext_composite_42_65_1805.Rds"
OUTPUT_PMOD_LAP  <- "~/DATA/SUN/TSI_PMOD_ext_composite_42_65_1805_LAP.Rds"


##  SORCE TSI from LISIRD  ---------------------------------------------------
FROM_SORCE       <- "https://lasp.colorado.edu/lisird/latis/dap/sorce_tsi_6hr_l3.csv?&format_time(yyyy-MM-dd'T'HH:mm:ss.SSS)"
DEST_SORCE       <- "~/DATA/SUN/TSI_sorce_tsi_6hr_l3.csv"
OUTPUT_SORCE     <- "~/DATA/SUN/TSI_sorce_tsi_6hr_l3.Rds"
OUTPUT_SORCE_LAP <- "~/DATA/SUN/TSI_sorce_tsi_6hr_l3_LAP.Rds"


