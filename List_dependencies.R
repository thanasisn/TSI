#!/opt/R/4.2.3/bin/Rscript
# /* Copyright (C) 2022-2023 Athanasios Natsis <natsisphysicist@gmail.com> */


## __ Set environment  ---------------------------------------------------------
rm(list = (ls()[ls() != ""]))
Sys.setenv(TZ = "UTC")
tic <- Sys.time()
Script.Name <- "~/TSI/List_dependencies.R"


## __ Describe environment -----------------------------------------------------
env_info <- "~/TSI/Dependencies.md"

cat("\n## TSI Current Running Environment\n", file = env_info)
cat("\n", R.version.string, "\n", file = env_info, append = TRUE)

pkgs <- renv::dependencies("~/TSI")
for (ap in unique(pkgs$Package)) {
    cat(sprintf("%16s:  %8s\n", ap, packageVersion(ap)), file = env_info, append = TRUE)
}


tac <- Sys.time()
cat(sprintf("%s %s@%s %s %f mins\n\n",Sys.time(),Sys.info()["login"],Sys.info()["nodename"],Script.Name,difftime(tac,tic,units="mins")))
