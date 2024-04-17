##########################################################################-
# Source script
# last modified: 2024-04-09

# This script contains all packages, data, and other objects for
# downstream use.
##########################################################################-

##########################################################################-
# 1. libraries ------------------------------------------------------------
##########################################################################-

# general use
library(here)
library(tidyverse)
library(janitor)
library(readxl)

# analysis
library(vegan)
library(ade4)
library(FD)
library(cluster)
library(NbClust)

# visualization
library(factoextra)

##########################################################################-
# 2. start and end dates --------------------------------------------------
##########################################################################-

# ⟞ a. Arroyo Quemado (AQUE) ----------------------------------------------

aque_start_dates <- c("AQUE_CONTROL_2008-01-30", 
                      "AQUE_ANNUAL_2008-01-30", 
                      "AQUE_CONTINUAL_2010-04-26")

aque_start_date <- as_date("2008-01-30")

aque_after_date <- as_date("2017-03-02")

aque_after_date_annual <- as_date("2018-05-10")

aque_after_date_continual <- as_date("2017-08-16")

# ⟞ b. Naples (NAPL) ------------------------------------------------------

napl_start_dates <- c("NAPL_CONTROL_2008-01-10", 
                      "NAPL_ANNUAL_2008-01-10", 
                      "NAPL_CONTINUAL_2010-04-27")

napl_start_date <- as_date("2008-01-10")

napl_after_date <- as_date("2016-02-19") # wrong in methods? 

napl_after_date_annual <- as_date("2017-05-16")

napl_after_date_continual <- as_date("2016-08-14")


# ⟞ c. Isla Vista (IVEE) --------------------------------------------------

ivee_start_dates <- c("IVEE_CONTROL_2011-10-26", 
                      "IVEE_ANNUAL_2011-10-26")

ivee_start_date <- as_date("2011-10-26")

ivee_after_date <- as_date("2016-02-18")

ivee_after_date_annual <- as_date("2017-05-15")

# ⟞ d. Mohawk (MOHK) ------------------------------------------------------

mohk_start_dates <- c("MOHK_ANNUAL_2008-01-17", 
                      "MOHK_CONTROL_2008-01-17", 
                      "MOHK_CONTINUAL_2010-05-05")

mohk_start_date <- as_date("2008-01-17")

mohk_after_date <- as_date("2017-02-13")

mohk_after_date_annual <- as_date("2018-05-15")

mohk_after_date_continual <- as_date("2017-08-11")

# ⟞ e. Carpintera (CARP) --------------------------------------------------


carp_start_dates <- c("CARP_CONTROL_2008-02-12", 
                      "CARP_ANNUAL_2008-02-12", 
                      "CARP_CONTINUAL_2010-04-23")

carp_start_date <- as_date("2008-02-12")

carp_after_date <- as_date("2017-02-15")

carp_after_date_annual <- as_date("2018-05-22")

carp_after_date_continual <- as_date("2017-08-10")

##########################################################################-
# 3. wrangling functions --------------------------------------------------
##########################################################################-

# create a column for "after" experimental removal
after_dates_column <- function(df) {
  df %>% 
    mutate(after_dates = case_when(
      site == "aque" & treatment == "annual" ~ aque_after_date_annual,
      site == "aque" & treatment == "continual" ~ aque_after_date_continual,
      site == "aque" & treatment == "control" ~ aque_after_date_continual,
      site == "napl" & treatment == "annual" ~ napl_after_date_annual,
      site == "napl" & treatment == "continual" ~ napl_after_date_continual,
      site == "napl" & treatment == "control" ~ napl_after_date_continual,
      site == "ivee" & treatment == "annual" ~ ivee_after_date_annual,
      site == "mohk" & treatment == "annual" ~ mohk_after_date_annual,
      site == "mohk" & treatment == "continual" ~ mohk_after_date_continual,
      site == "mohk" & treatment == "control" ~ mohk_after_date_continual,
      site == "carp" & treatment == "annual" ~ carp_after_date_annual,
      site == "carp" & treatment == "continual" ~ carp_after_date_continual,
      site == "carp" & treatment == "control" ~ carp_after_date_continual
    ))
}

# make a new column for during and after and set factor levels
exp_dates_column <- function(df) {
  df %>% 
    mutate(exp_dates = case_when(
      # after for annual removal:
      site == "aque" & treatment == "annual" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "annual" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "annual" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "annual" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "annual" & date > carp_after_date_annual ~ "after",
      # after for continual removal:
      site == "aque" & treatment == "continual" & date > aque_after_date_continual ~ "after",
      site == "napl" & treatment == "continual" & date > napl_after_date_continual ~ "after",
      site == "mohk" & treatment == "continual" & date > mohk_after_date_continual ~ "after",
      site == "carp" & treatment == "continual" & date > carp_after_date_continual ~ "after",
      # after for control:
      site == "aque" & treatment == "control" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "control" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "control" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "control" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "control" & date > carp_after_date_annual ~ "after",
      # everything else is "during" the experiment
      TRUE ~ "during"
    ),
    exp_dates = fct_relevel(exp_dates, c("during", "after")))  
}

# annual removal: make a new column for during and after and set factor levels
exp_dates_column_annual <- function(df) {
  df %>% 
    mutate(exp_dates = case_when(
      # after for annual removal:
      site == "aque" & treatment == "annual" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "annual" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "annual" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "annual" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "annual" & date > carp_after_date_annual ~ "after",
      # after for control:
      site == "aque" & treatment == "control" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "control" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "control" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "control" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "control" & date > carp_after_date_annual ~ "after",
      # everything else is "during" the experiment
      TRUE ~ "during"
    ),
    exp_dates = fct_relevel(exp_dates, c("during", "after")))  
}

# continual removal: make a new column for during and after and set factor levels
exp_dates_column_continual <- function(df) {
  df %>% 
    mutate(exp_dates = case_when(
      # after for continual removal:
      site == "aque" & treatment == "continual" & date > aque_after_date_continual ~ "after",
      site == "napl" & treatment == "continual" & date > napl_after_date_continual ~ "after",
      site == "mohk" & treatment == "continual" & date > mohk_after_date_continual ~ "after",
      site == "carp" & treatment == "continual" & date > carp_after_date_continual ~ "after",
      # after for control:
      site == "aque" & treatment == "control" & date > aque_after_date_continual ~ "after",
      site == "napl" & treatment == "control" & date > napl_after_date_continual ~ "after",
      # site == "ivee" & treatment == "control" & date > ivee_after_date_continual ~ "after",
      site == "mohk" & treatment == "control" & date > mohk_after_date_continual ~ "after",
      site == "carp" & treatment == "control" & date > carp_after_date_continual ~ "after",
      # everything else is "during" the experiment
      TRUE ~ "during"
    ),
    exp_dates = fct_relevel(exp_dates, c("during", "after")))  
}

# create a new column for season and set factor levels
season_column <- function(df) {
  df %>% 
    mutate(season = case_when(
      month %in% c(12, 1, 2, 3) ~ "winter",
      month %in% c(4, 5) ~ "spring",
      month %in% c (6, 7, 8) ~ "summer",
      month %in% c(9, 10, 11) ~ "fall"
    ),
    season = fct_relevel(season, "spring", "summer", "fall", "winter")) 
}

# create a new column for new groupings
new_group_column <- function(df) {
  df %>% 
    mutate(new_group = case_when(
      group == "algae" ~ "algae",
      group == "fish" ~ "fish",
      group == "invert" & taxon_family != "Pholadidae" & mobility == "sessile" ~ "epi_inverts",
      taxon_family == "Pholadidae" ~ "endo_inverts",
      sp_code %in% c("OPSP", "SPL", "LIGL", "SFL", "MECR") ~ "herb_inverts",
      sp_code %in% c("COCA", "AML", "PGL", "PAIN", "PHL", "LOGR", "DIL", "KEKE", "PBL") ~ "carn_inverts"
    ))
}

# function to calculate standard error
se <- function(x,...){
  sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
}

# time since start and time since end
# expects data frame with the following columns:
# site, year, month, date, control, annual, delta_annual
time_since_columns_annual <- function(df) {
  df %>% 
    # create a column for quarter
    mutate(quarter = case_when(
      month <= 3 ~ "Q1",
      month <= 6 ~ "Q2",
      month <= 9 ~ "Q3",
      TRUE ~ "Q4"
    )) %>% 
    # calculate time since start of experiment
    mutate(time_yrs = case_when(
      # AQUE, NAPL, MOHK, CARP: control and annual started in 2008
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q1" ~ year + 0.125 - 2008,
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q2" ~ year + 0.375 - 2008,
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q3" ~ year + 0.625 - 2008,
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q4" ~ year + 0.875 - 2008, 
      # IVEE control and annual started in 2011
      site == "ivee" & quarter == "Q1" ~ year + 0.125 - 2011,
      site == "ivee" & quarter == "Q2" ~ year + 0.375 - 2011,
      site == "ivee" & quarter == "Q3" ~ year + 0.625 - 2011,
      site == "ivee" & quarter == "Q4" ~ year + 0.875 - 2011
    )) %>% 
    group_by(site) %>% 
    mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
    ungroup() %>% 
    # calculate time since end of experiment
    group_by(site, exp_dates) %>% 
    # if "after", then simple: the time in years - the minimum time in years
    # if "during", then more complex: take the max time in years and add 0.25, then subtract the time in years
    mutate(time_since_end = case_when(
      exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
      exp_dates == "after" ~ time_yrs - min(time_yrs)
    )) %>% 
    ungroup()
}

time_since_columns_continual <- function(df) {
  df %>% 
    # create a column for quarter
    mutate(quarter = case_when(
      month <= 3 ~ "Q1",
      month >= 4 & month <= 6 ~ "Q2",
      month >= 7 & month <= 9 ~ "Q3",
      month >= 10 & month <= 12 ~ "Q4"
    )) %>% 
    # mutate(quarter = case_when(
    #   month <= 3 ~ "Q1",
    #   month <= 6 ~ "Q2",
    #   month <= 9 ~ "Q3",
    #   TRUE ~ "Q4"
    # )) %>% 
    # calculate time since start of experiment
    mutate(time_yrs = case_when(
      # AQUE, NAPL, MOHK, CARP: continual started in 2010
      quarter == "Q1" ~ year + 0.125 - 2010,
      quarter == "Q2" ~ year + 0.375 - 2010,
      quarter == "Q3" ~ year + 0.625 - 2010,
      quarter == "Q4" ~ year + 0.875 - 2010
    )) %>% 
    group_by(site) %>% 
    mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
    ungroup() %>% 
    # calculate time since end of experiment
    group_by(site, exp_dates) %>% 
    mutate(time_since_end = case_when(
      exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
      exp_dates == "after" ~ time_yrs - min(time_yrs)
    )) %>% 
    mutate(test_min_time_yrs = min(time_yrs)) %>% 
    ungroup()
}

# kelp year
# expects data frame that already has quarter from time_since_columns functions
kelp_year_column <- function(df) {
  df %>% 
    # create a new column for "kelp year"
    mutate(quarter = fct_relevel(quarter, "Q2", "Q3", "Q4", "Q1")) %>% 
    mutate(kelp_year = case_when(
      quarter %in% c("Q2", "Q3", "Q4") ~ year,
      quarter == "Q1" ~ year - 1
    )) %>% 
    mutate(kelp_year = paste("kelp_", kelp_year, "-", kelp_year + 1, sep = "")) 
}

# comparison column for annual removal
# first 2 or 3 years of start, last 2 or 3 years during experimental removal, most recent 2 or 3 years
comparison_column_annual <- function(df) {
  df %>% 
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022") ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    # create a column for the points to compare for "3 year interval"
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010", "kelp_2010-2011") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2019-2020", "kelp_2020-2021", "kelp_2021-2022") ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    # add in full names of sites
    left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
    rename(site_full = value) %>% 
    mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Isla Vista", "Mohawk", "Carpinteria"))
}

comparison_column_annual_new <- function(df) {
  df %>% 
    mutate(comp_1yrs = case_when(
      site %in% c("aque", "carp") & between(time_since_end, -7.25, -6.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -5.25) ~ "start",
      site == "mohk" & between(time_since_end, -7, -6) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -1.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 4.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 5.75, 6.75) ~ "after"
      
    )) %>% 
    mutate(comp_1yrs = fct_relevel(comp_1yrs, "start", "during", "after")) %>% 
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "carp") & between(time_since_end, -7.25, -5.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -4.25) ~ "start",
      site == "mohk" & between(time_since_end, -7, -5) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -2.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 3.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 4.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "carp") & between(time_since_end, -7.25, -4.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -3.25) ~ "start",
      site == "mohk" & between(time_since_end, -7, -4) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -3.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 2.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 3.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
    unite("sample_ID_short", site, date, remove = FALSE)
}

comparison_column_continual <- function(df) {
  df %>% 
    # create a column for the points to compare for "1 year interval"
    mutate(comp_1yr = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2015-2016") ~ "during",
      kelp_year %in% c("kelp_2022-2023") ~ "after"
    )) %>% 
    mutate(comp_1yr = fct_relevel(comp_1yr, "start", "during", "after")) %>% 
    # create a column for the points to compare for "2 year interval"
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2021-2022", "kelp_2022-2023") ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    # create a column for the points to compare for "3 year interval"
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012", "kelp_2012-2013") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022", "kelp_2022-2023") ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    # create a new sample ID that is site, date, quarter
    unite("sample_ID", site, date, quarter, remove = FALSE)
}

comparison_column_continual_new <- function(df) {
  df %>% 
    mutate(comp_1yrs = case_when(
      site %in% c("aque", "carp") & between(time_since_end, -7.25, -6.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -5.25) ~ "start",
      site == "mohk" & between(time_since_end, -7, -6) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -1.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 4.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 5.75, 6.75) ~ "after"
      
    )) %>% 
    mutate(comp_1yrs = fct_relevel(comp_1yrs, "start", "during", "after")) %>% 
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "carp") & between(time_since_end, -7.25, -5.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -4.25) ~ "start",
      site == "mohk" & between(time_since_end, -7, -5) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -2.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 3.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 4.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "carp") & between(time_since_end, -7.25, -4.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -3.25) ~ "start",
      site == "mohk" & between(time_since_end, -7, -4) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -3.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 2.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 3.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
    unite("sample_ID_short", site, date, remove = FALSE)
}

# from AlexisW: https://forum.posit.co/t/rstudio-error-with-fviz-nbclust-for-nbclust-results/176083
my_fviz_nbclust <- function(x, print.summary = TRUE, barfill = "steelblue", barcolor = "steelblue"){
  best_nc <- x$Best.nc
  best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
  best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
  
  ss <- summary(best_nc$Number_clusters)
  cat("Among all indices: \n===================\n")
  for (i in 1:length(ss)) {
    cat("*", ss[i], "proposed ", names(ss)[i], "as the best number of clusters\n")
  }
  cat("\nConclusion\n=========================\n")
  cat("* According to the majority rule, the best number of clusters is ", 
      names(which.max(ss)), ".\n\n")
  
  df <- data.frame(Number_clusters = names(ss), freq = ss, 
                   stringsAsFactors = TRUE)
  p <- ggpubr::ggbarplot(df, x = "Number_clusters", y = "freq", 
                         fill = "steelblue", color = "steelblue") +
    ggplot2::labs(x = "Number of clusters k", 
                  y = "Frequency among all indices",
                  title = paste0("Optimal number of clusters - k = ", 
                                 names(which.max(ss))))
  p
}

##########################################################################-
# 3. data -----------------------------------------------------------------
##########################################################################-

guilds <- read_csv(here::here("data", "SBC-LTE", "LTE_guild_data.csv")) %>% 
  mutate(sp.code = replace_na(sp.code, "Nandersoniana")) %>% 
  rename("new_group" = biomass.guild) %>% 
  # for whatever reason Yellowtail Rockfish are not in the guild csv
  add_row(sp.code = "SFLA", new_group = "fish", diversity.guild = "fish")

biomass <- read_csv(
  here("data",
       "SBC-LTE",
       "knb-lter-sbc.119.10",
       "LTE_All_Species_Biomass_at_transect_20230530.csv")
) %>% 
  clean_names() %>%
  # ANOB is incorrectly coded as having "SESSILE" mobility
  mutate(mobility = replace(mobility, sp_code == "ANOB", "MOBILE")) %>%
  # UEC is incorrectly coded as Amphiroa beauvoisii
  mutate(scientific_name = replace(scientific_name, sp_code == "UEC", "Unidentified erect coralline spp.")) %>%
  # replace NA sp_code with Nandersoniana
  mutate(sp_code = case_when(
    scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
    TRUE ~ sp_code
  )) %>%
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         density = replace(density, density < 0, NA)) %>%
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>%
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "treatment", "site"), str_to_lower) %>%
  # make a new column for during and after and set factor levels
  exp_dates_column() %>%
  # create a new column for season and set factor levels
  season_column() %>%
  # new group column %>%
  left_join(., guilds, by = c("sp_code" = "sp.code")) %>%
  # take out all the first dates
  filter(!(sample_ID %in% c(aque_start_dates, napl_start_dates, ivee_start_dates, mohk_start_dates, carp_start_dates))) %>%
  # dangling controls (from annual plot surveys) makes things harder
  filter(!(sample_ID %in% c("NAPL_CONTROL_2010-04-27", "CARP_CONTROL_2010-04-23",
                            "AQUE_CONTROL_2010-04-26", "MOHK_CONTROL_2010-05-05"))) %>%
  # calculating average biomass (across sampling dates for 2010-2012, when sampling was done 8x per year)
  time_since_columns_continual() %>%
  group_by(site, year, treatment, quarter, sp_code) %>%
  mutate(dry_gm2 = mean(dry_gm2),
         wm_gm2 = mean(wm_gm2),
         density = mean(density)) %>%
  # take out the "duplicates": only one sampling date per quarter in the dataframe, with values averaged across the two sampling dates
  slice(1L) %>%
  ungroup() %>%
  # take out extraneous columns from time_since_columns_continual()
  select(!quarter:test_min_time_yrs)


traits <- read_csv(here("data", 
                         "functional-traits",
                         "joe-traits-lter_2024-04-16.csv"))


##########################################################################-
# 4. objects --------------------------------------------------------------
##########################################################################-

# function to calculate standard error
se <- function(x,...){
  sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
}

sites_full <- setNames(c("Arroyo Quemado",
                         "Naples",
                         "Isla Vista",
                         "Mohawk",
                         "Carpinteria"),
                       c("aque",
                         "napl",
                         "ivee",
                         "mohk",
                         "carp"))

site_quality <- tribble(
  ~site, ~quality,
  "mohk", "high",
  "ivee", "high",
  "aque", "medium",
  "napl", "medium",
  "carp", "low"
) %>% 
  mutate(quality = fct_relevel(quality, c("low", "medium", "high")))

algae_spp <- biomass %>% 
  filter(new_group == "algae") %>% 
  select(scientific_name, sp_code) %>% 
  unique()

excluded_spp <- c(
    "Halymenia spp.; Schizymenia pacifica",
    "crustose coralline algae spp.",
    "Ectocarpaceae spp.",
    "Ulva spp.; Sponogomorpha spp.",
    "Rhodophyta",
    "Neoptilota spp.; Ptilota spp.; Rhodoptilum spp.",
    "Unidentifiable Branching Red Alga",
    "Unidentifiable juvenile kelp",
    "small Ceramiaceae spp.",
    "Unidentifiable small brown blade",
    "Unidentified erect coralline spp."
  )

algae_taxa <- biomass %>% 
  filter(new_group == "algae") %>% 
  filter(!(scientific_name %in% excluded_spp)) %>% 
  select(scientific_name, taxon_phylum:taxon_genus) %>% 
  unique()

##########################################################################-
# 5. plot themes ----------------------------------------------------------
##########################################################################-

chloro_col <- "#8AA789"

ochro_col <- "#985E5C"

rhodo_col <- "#4CA2B0"

theme_set(theme_bw() +
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12)))






