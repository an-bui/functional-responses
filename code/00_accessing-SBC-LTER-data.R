##########################################################################-
# Getting SBC LTER data from EDI
# last modified: 2024-04-09
##########################################################################-

##########################################################################-
# 1. libraries ------------------------------------------------------------
##########################################################################-

# file path organization
library(here)

# pulling from the EDI
library(EDIutils)

# general use
library(tidyverse)

##########################################################################-
# 2. SBC LTER data on EDI -------------------------------------------------
##########################################################################-

# SBC LTE package ID
packageID <- "knb-lter-sbc.119.10"

# downloads zipped file
read_data_package_archive(packageID, path = here("data", "SBC-LTE"))

# creates an intermediate object
LTE_All_Species_Biomass_at_transect_20230530 <- read_csv(
  here("data", 
       "SBC-LTE", 
       "raw-data", 
       "knb-lter-sbc.119.10", 
       "LTE_All_Species_Biomass_at_transect_20230530.csv")
)

# creates an RDS file
saveRDS(object = LTE_All_Species_Biomass_at_transect_20230530,
        file = here("data",
                    "SBC-LTE",
                    "rds",
                    "LTE_All_Species_Biomass_at_transect_20230530.rds"))



