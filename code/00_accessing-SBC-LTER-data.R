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

# unzips file
unzip(zipfile = here("data", "SBC-LTE", "knb-lter-sbc.119.10.zip"),
      exdir = here("data", "SBC-LTE", "knb-lter-sbc.119.10"))




