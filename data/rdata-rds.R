## code to prepare `rdata.rds` dataset goes here

rdata <- readRDS(paste(getwd(),"/rdata.rds",sep=""))
usethis::use_data(rdata, overwrite = TRUE)

vdata <- readRDS(paste(getwd(),"/vdata.rds",sep=""))
usethis::use_data(vdata, overwrite = TRUE)

# usethis::load_all()
usethis::use_data_raw()
