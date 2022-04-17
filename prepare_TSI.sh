#!/bin/bash
## created on 2022-04-15

#### This is for daily usage.
## For building TSI data


info() { echo ; echo "$(date +%F_%T) :: $* " >&1; }
mkdir -p "$(dirname "$0")/LOGs/"
LOG_FILE="$(dirname "$0")/LOGs/$(basename "$0")_$(date +%F_%T).log"
ERR_FILE="$(dirname "$0")/LOGs/$(basename "$0")_$(date +%F_%T).err"
exec  > >(tee -i "${LOG_FILE}")
exec 2> >(tee -i "${ERR_FILE}" >&2)
info "START :: $0 :: $* ::"



# info "Get TSI model from NOAA"
# Rscript "$(dirname "$0")/TSI_get_NOAA.R"

# info "Prepare NOAA TSI for LAP"
# Rscript "$(dirname "$0")/TSI_NOAA_LAP.R"

info "Get TSI from TSIS"
Rscript "$(dirname "$0")/TSI_get_TSIS.R"

info "Prepare TSIS TSI for LAP"
Rscript "$(dirname "$0")/TSI_TSIS_LAP.R"

## Frozen data
# info "Get TSI from PMOD"
# Rscript "$(dirname "$0")/TSI_get_PMOD.R"

info "Prepare PMOD TSI for LAP"
Rscript "$(dirname "$0")/TSI_PMOD_LAP.R"

## Frozen data
# info "Get TSI from SORCE"
# Rscript "$(dirname "$0")/TSI_get_SORCE.R"

info "Prepare TSIS SORCE for LAP"
Rscript "$(dirname "$0")/TSI_SORCE_LAP.R"








## end coding
printf "%s %-10s %-10s %-50s %f\n" "$(date +"%F %H:%M:%S")" "$HOSTNAME" "$USER" "$(basename $0)" "$SECONDS"
exit 0
