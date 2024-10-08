##########################################################################-
# Getting SBC LTER data from EDI
# last modified: 2024-04-09
# This only needs to be run to get data from EDI; otherwise, start with
# script `01_source.R`.
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
packageID <- "knb-lter-sbc.119.11"

# downloads zipped file
read_data_package_archive(packageID, path = here("data", "SBC-LTE"))

# unzips file
unzip(zipfile = here("data", "SBC-LTE", "knb-lter-sbc.119.11.zip"),
      exdir = here("data", "SBC-LTE", "knb-lter-sbc.119.11"))

# SBC LTE NPP dataset
nppID <- "knb-lter-sbc.58.18"

# downloads zipped file
read_data_package_archive(nppID, path = here("data", "SBC-LTER"))

# unzips file
unzip(zipfile = here("data", "SBC-LTE", "knb-lter-sbc.58.18.zip"),
      exdir = here("data", "SBC-LTE", "knb-lter-sbc.58.18"))
