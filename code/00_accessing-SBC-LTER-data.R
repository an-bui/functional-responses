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
packageID <- "knb-lter-sbc.119.12"

# downloads zipped file
read_data_package_archive(packageID, path = here("data", "SBC-LTE"))

# unzips file
unzip(zipfile = here("data", "SBC-LTE", "knb-lter-sbc.119.12.zip"),
      exdir = here("data", "SBC-LTE", "knb-lter-sbc.119.12"))

# SBC LTE NPP dataset
nppID <- "knb-lter-sbc.58.19"

# downloads zipped file
read_data_package_archive(nppID, path = here("data", "SBC-LTE"))

# unzips file
unzip(zipfile = here("data", "SBC-LTE", "knb-lter-sbc.58.19.zip"),
      exdir = here("data", "SBC-LTE", "knb-lter-sbc.58.19"))

# annual surveys
annualID <- "knb-lter-sbc.50.17"

# downloads zipped file
read_data_package_archive(annualID, path = here("data", "SBC-LTER"))

# unzips file
unzip(zipfile = here("data", "SBC-LTER", "knb-lter-sbc.50.17.zip"),
      exdir = here("data", "SBC-LTER", "knb-lter-sbc.50.17"))

# substrate from annual surveys
substrateID <- "knb-lter-sbc.138.5"

# downloads zipped file
read_data_package_archive(substrateID, path = here("data", "SBC-LTER"))

# unzips file
unzip(zipfile = here("data", "SBC-LTER", "knb-lter-sbc.138.5.zip"),
      exdir = here("data", "SBC-LTER", "knb-lter-sbc.138.5"))

# photosynthesis parameters
photosynthesisID <- "knb-lter-sbc.127.3"

# downloads zipped file
read_data_package_archive(photosynthesisID, path = here("data", "SBC-LTER"))

# unzips file
unzip(zipfile = here("data", "SBC-LTER", "knb-lter-sbc.127.3.zip"),
      exdir = here("data", "SBC-LTER", "knb-lter-sbc.127.3"))
