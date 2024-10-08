##########################################################################-
# Source script
# last modified: 2024-09-27

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
library(repmis)

# tables
library(flextable)
library(gtsummary)

# analysis
library(vegan)
library(ade4)
library(FD)
library(mFD)
library(cluster)
library(NbClust)
library(glmmTMB)
library(DHARMa)
library(performance)
library(RVAideMemoire)
library(modelsummary)
library(MuMIn)
library(car)

# visualization
library(factoextra)
library(ggdendro)
library(ggeffects)
library(funspace)
library(ggdist)
library(ggConvexHull)

##########################################################################-
# 2. start and end dates --------------------------------------------------
##########################################################################-

# ⟞ a. Arroyo Quemado (AQUE) ---------------------------------------------

# date of first removal: 2010-01-29
# first survey after first removal: 2010-04-26
aque_start_date_continual <- as_date("2010-04-26")
# date of last removal: 2017-03-02
# first survey after last removal: 2017-05-18
# start of recovery period: 2017-08-16
aque_after_date_continual <- as_date("2017-08-16")


# ⟞ b. Naples (NAPL) -----------------------------------------------------

# date of first removal: 2010-01-28
# first survey after first removal: 2010-04-27
napl_start_date_continual <- as_date("2010-04-27")

# date of last removal: 2016-02-09
# first survey after last removal: 2016-05-17
# start of recovery period: 2016-08-16
napl_after_date_continual <- as_date("2016-08-16")

# ⟞ c. Mohawk (MOHK) ------------------------------------------------------

# date of first removal: 2010-05-05
# first survey after first removal: 2010-06-14
mohk_start_date_continual <- as_date("2010-07-23")

# date of last removal: 2017-02-13
# first survey after last removal: 2017-05-17
# start of recovery period: 2017-08-11
mohk_after_date_continual <- as_date("2017-08-11")

# ⟞ d. Carpinteria (CARP)  ------------------------------------------------

# date of first removal: 2010-02-04 
# first survey after first removal: 2010-04-23
carp_start_date_continual <- as_date("2010-04-23")
# first survey after first removal: 2010-03-11

# date of last removal: 2017-02-15
# first survey after last removal: 2017-05-19
# start of recovery period: 2017-08-10
carp_after_date_continual <- as_date("2017-08-10")


##########################################################################-
# 3. wrangling functions --------------------------------------------------
##########################################################################-

# create a column for "after" experimental removal
# after_dates_column <- function(df) {
#   df %>% 
#     mutate(after_dates = case_when(
#       site == "aque" & treatment == "annual" ~ aque_after_date_annual,
#       site == "aque" & treatment == "continual" ~ aque_after_date_continual,
#       site == "aque" & treatment == "control" ~ aque_after_date_continual,
#       site == "napl" & treatment == "annual" ~ napl_after_date_annual,
#       site == "napl" & treatment == "continual" ~ napl_after_date_continual,
#       site == "napl" & treatment == "control" ~ napl_after_date_continual,
#       site == "ivee" & treatment == "annual" ~ ivee_after_date_annual,
#       site == "mohk" & treatment == "annual" ~ mohk_after_date_annual,
#       site == "mohk" & treatment == "continual" ~ mohk_after_date_continual,
#       site == "mohk" & treatment == "control" ~ mohk_after_date_continual,
#       site == "carp" & treatment == "annual" ~ carp_after_date_annual,
#       site == "carp" & treatment == "continual" ~ carp_after_date_continual,
#       site == "carp" & treatment == "control" ~ carp_after_date_continual
#     ))
# }

# make a new column for during and after and set factor levels
exp_dates_column <- function(df) {
  df %>%
    mutate(exp_dates = case_when(
      # after for annual removal:
      # site == "aque" & treatment == "annual" & date > aque_after_date_annual ~ "after",
      # site == "napl" & treatment == "annual" & date > napl_after_date_annual ~ "after",
      # site == "ivee" & treatment == "annual" & date > ivee_after_date_annual ~ "after",
      # site == "mohk" & treatment == "annual" & date > mohk_after_date_annual ~ "after",
      # site == "carp" & treatment == "annual" & date > carp_after_date_annual ~ "after",
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
# exp_dates_column_annual <- function(df) {
#   df %>% 
#     mutate(exp_dates = case_when(
#       # after for annual removal:
#       site == "aque" & treatment == "annual" & date > aque_after_date_annual ~ "after",
#       site == "napl" & treatment == "annual" & date > napl_after_date_annual ~ "after",
#       site == "ivee" & treatment == "annual" & date > ivee_after_date_annual ~ "after",
#       site == "mohk" & treatment == "annual" & date > mohk_after_date_annual ~ "after",
#       site == "carp" & treatment == "annual" & date > carp_after_date_annual ~ "after",
#       # after for control:
#       site == "aque" & treatment == "control" & date > aque_after_date_annual ~ "after",
#       site == "napl" & treatment == "control" & date > napl_after_date_annual ~ "after",
#       site == "ivee" & treatment == "control" & date > ivee_after_date_annual ~ "after",
#       site == "mohk" & treatment == "control" & date > mohk_after_date_annual ~ "after",
#       site == "carp" & treatment == "control" & date > carp_after_date_annual ~ "after",
#       # everything else is "during" the experiment
#       TRUE ~ "during"
#     ),
#     exp_dates = fct_relevel(exp_dates, c("during", "after")))  
# }

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
# time_since_columns_annual <- function(df) {
#   df %>% 
#     # create a column for quarter
#     mutate(quarter = case_when(
#       month <= 3 ~ "Q1",
#       month <= 6 ~ "Q2",
#       month <= 9 ~ "Q3",
#       TRUE ~ "Q4"
#     )) %>% 
#     # calculate time since start of experiment
#     mutate(time_yrs = case_when(
#       # AQUE, NAPL, MOHK, CARP: control and annual started in 2008
#       site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q1" ~ year + 0.125 - 2008,
#       site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q2" ~ year + 0.375 - 2008,
#       site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q3" ~ year + 0.625 - 2008,
#       site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q4" ~ year + 0.875 - 2008, 
#       # IVEE control and annual started in 2011
#       site == "ivee" & quarter == "Q1" ~ year + 0.125 - 2011,
#       site == "ivee" & quarter == "Q2" ~ year + 0.375 - 2011,
#       site == "ivee" & quarter == "Q3" ~ year + 0.625 - 2011,
#       site == "ivee" & quarter == "Q4" ~ year + 0.875 - 2011
#     )) %>% 
#     group_by(site) %>% 
#     mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
#     ungroup() %>% 
#     # calculate time since end of experiment
#     group_by(site, exp_dates) %>% 
#     # if "after", then simple: the time in years - the minimum time in years
#     # if "during", then more complex: take the max time in years and add 0.25, then subtract the time in years
#     mutate(time_since_end = case_when(
#       exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
#       exp_dates == "after" ~ time_yrs - min(time_yrs)
#     )) %>% 
#     ungroup()
# }

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
# comparison_column_annual <- function(df) {
#   df %>% 
#   mutate(comp_2yrs = case_when(
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010") ~ "start",
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during", 
#     site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
#     site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022") ~ "after"
#   )) %>% 
#     mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
#     # create a column for the points to compare for "3 year interval"
#     mutate(comp_3yrs = case_when(
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010", "kelp_2010-2011") ~ "start",
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during", 
#       site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
#       site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
#       kelp_year %in% c("kelp_2019-2020", "kelp_2020-2021", "kelp_2021-2022") ~ "after"
#     )) %>% 
#     mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
#     # add in full names of sites
#     left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
#     rename(site_full = value) %>% 
#     mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Isla Vista", "Mohawk", "Carpinteria"))
# }

# comparison_column_continual <- function(df) {
#   df %>%
#     # create a column for the points to compare for "1 year interval"
#     mutate(comp_1yr = case_when(
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011") ~ "start",
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016") ~ "during",
#       site == "ivee" & kelp_year %in% c("kelp_2011-2012") ~ "start",
#       site == "ivee" & kelp_year %in% c("kelp_2015-2016") ~ "during",
#       kelp_year %in% c("kelp_2022-2023") ~ "after"
#     )) %>%
#     mutate(comp_1yr = fct_relevel(comp_1yr, "start", "during", "after")) %>%
#     # create a column for the points to compare for "2 year interval"
#     mutate(comp_2yrs = case_when(
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012") ~ "start",
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
#       site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
#       site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
#       kelp_year %in% c("kelp_2021-2022", "kelp_2022-2023") ~ "after"
#     )) %>%
#     mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>%
#     # create a column for the points to compare for "3 year interval"
#     mutate(comp_3yrs = case_when(
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012", "kelp_2012-2013") ~ "start",
#       site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
#       site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
#       site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
#       kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022", "kelp_2022-2023") ~ "after"
#     )) %>%
#     mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>%
#     # create a new sample ID that is site, date, quarter
#     unite("sample_ID", site, date, quarter, remove = FALSE)
# }

comparison_column_continual_new <- function(df) {
  df %>% 
    mutate(comp_1yrs = case_when(
      site %in% c("aque", "mohk", "carp") & between(time_since_end, -7.25, -6.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -5.25) ~ "start",
      # site == "mohk" & between(time_since_end, -7, -6) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -1.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 4.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 5.75, 6.75) ~ "after"
      
    )) %>% 
    mutate(comp_1yrs = fct_relevel(comp_1yrs, "start", "during", "after")) %>% 
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "mohk", "carp") & between(time_since_end, -7.25, -5.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -4.25) ~ "start",
      # site == "mohk" & between(time_since_end, -7, -5) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -2.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 3.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 4.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "mohk", "carp") & between(time_since_end, -7.25, -4.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -3.25) ~ "start",
      # site == "mohk" & between(time_since_end, -7, -4) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -3.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 2.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 3.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
    unite("sample_ID_short", site, date, remove = FALSE)
}

# anova_summary_fxn <- function(adonis2.obj) {
#   # turn object name into string
#   name <- deparse(substitute(adonis2.obj))
#   
#   adonis2.obj %>% 
#     # turn adonis2 result into data frame
#     as.data.frame() %>% 
#     # make rownames "variables"
#     rownames_to_column("variables") %>% 
#     # rename Pr(>F) column into something intelligible
#     rename(p = `Pr(>F)`) %>% 
#     # round values to 2 decimal points
#     mutate(across(SumOfSqs:p, ~ round(.x, digits = 3))) %>% 
#     # replace comp_.yrs with time period
#     mutate(variables = str_replace(variables, "comp_.yrs", "time period")) %>% 
#     mutate(variables = str_replace(variables, "comp_.yr", "time period")) %>% 
#     # make object name column
#     mutate(model = name) 
# }



# anova_summary_fxn <- function(adonis2.obj) {
#   # turn object name into string
#   name <- deparse(substitute(adonis2.obj))
#   
#   adonis2.obj %>% 
#     # turn adonis2 result into data frame
#     as.data.frame() %>% 
#     # make rownames "variables"
#     rownames_to_column("variables") %>% 
#     # rename Pr(>F) column into something intelligible
#     rename(p = `Pr(>F)`) %>% 
#     # round values to 2 decimal points
#     mutate(across(SumOfSqs:p, ~ round(.x, digits = 3))) %>% 
#     # replace comp_.yrs with time period
#     mutate(variables = str_replace(variables, "comp_.yrs", "time period")) %>% 
#     mutate(variables = str_replace(variables, "comp_.yr", "time period")) %>% 
#     # make object name column
#     mutate(model = name) 
# }

# difflsmeans_summary_fxn <- function(anova.obj) {
#   anova.obj %>% 
#     difflsmeans(test.effs = "Group", ddf = "Kenward-Roger") %>% 
#     as.data.frame() %>% 
#     clean_names() %>% 
#     rownames_to_column("rowname") %>% 
#     select(levels, estimate, std_error, df, t_value, pr_t) %>% 
#     mutate(across(c(estimate, std_error, t_value, pr_t), ~round(., digits = 3))) %>% 
#     mutate(levels = fct_relevel(levels, c("start - during", "during - after", "start - after"))) %>% 
#     arrange(levels)
# }

# function to extract model summaries
model_summary_fxn <- function(model) {
  model %>% 
    # use tidy to get model summary and calculate 95% CI
    tidy(conf.int = TRUE) %>% 
    # only include fixed conditional effects
    filter(effect == "fixed" & component == "cond") %>%
    select(term, estimate, p.value, conf.low, conf.high) %>% 
    # create a new column that indicates whether an effect is significant
    mutate(signif = case_when(
      p.value <= 0.05 ~ "yes",
      TRUE ~ "no"
    )) %>% 
    # create a p-value column that converts very small values to < 0.001
    # and rounds all other values to relevant digits
    mutate(p.value = case_when(
      between(p.value, 0, 0.001) ~ "<0.001",
      between(p.value, 0.001, 0.01) ~ as.character(round(p.value, digits = 3)),
      between(p.value, 0.01, 1) ~ as.character(round(p.value, digits = 2))
    )) %>%
    # round other numeric values to two digits
    mutate(across(where(is.numeric), ~ round(., digits = 2))) %>%
    # create a confidence interval column
    unite(ci_interval, conf.low, conf.high, sep = ", ") %>%
    # rename the terms to be neater
    mutate(term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "time_since_end" ~ "Time since end",
      term == "treatmentremoval" ~ "Treatment (removal)",
      term == "time_since_end:treatmentremoval" ~ "Time since end × treatment (removal)"
    ))
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
       "knb-lter-sbc.119.11",
       "LTE_All_Species_Biomass_at_transect_20240501.csv")
) %>%
  clean_names() %>%
  # replace NA sp_code with Nandersoniana
  mutate(sp_code = case_when(
    scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
    TRUE ~ sp_code
  )) %>%
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         density = replace(density, density < 0, NA)) %>%
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "treatment", "site"), str_to_lower) %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # filter to only include continual removal plots and control plots
  filter(treatment %in% c("continual", "control")) %>% 
  left_join(., guilds, by = c("sp_code" = "sp.code")) %>% 
  mutate(exp_dates = case_when(
    site == "aque" & date >= aque_start_date_continual & date < aque_after_date_continual ~ "during",
    site == "aque" & date >= aque_after_date_continual ~ "after",
    site == "napl" & date >= napl_start_date_continual & date < napl_after_date_continual ~ "during",
    site == "napl" & date >= napl_after_date_continual ~ "after",
    site == "mohk" & date >= mohk_start_date_continual & date < mohk_after_date_continual ~ "during",
    site == "mohk" & date >= mohk_after_date_continual ~ "after",
    site == "carp" & date >= carp_start_date_continual & date < carp_after_date_continual ~ "during",
    site == "carp" & date >= carp_after_date_continual ~ "after"
  ),
  exp_dates = fct_relevel(exp_dates, "during", "after")) %>% 
  # take out all surveys that were before the removal experiment started
  drop_na(exp_dates) %>% 
  time_since_columns_continual() %>%
  group_by(site, year, treatment, quarter, sp_code) %>%
  mutate(dry_gm2 = mean(dry_gm2),
         wm_gm2 = mean(wm_gm2),
         density = mean(density)) %>%
  # take out the "duplicates": only one sampling date per quarter in the dataframe, with values averaged across the two sampling dates
  slice(1L) %>%
  ungroup() %>%
  # take out extraneous columns from time_since_columns_continual()
  select(!test_min_time_yrs) %>% 
  comparison_column_continual_new() %>% 
  kelp_year_column() %>% 
  season_column() %>% 
  unite("season_ID", year, season, site, treatment, remove = FALSE)

npp <- read_csv(
  here("data",
       "SBC-LTE",
       "knb-lter-sbc.58.18",
       "Understory_NPP_All_Year_season_20240501.csv")
) %>%
  clean_names() %>%
  # replace NA sp_code with Nandersoniana
  mutate(sp_code = case_when(
    scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
    TRUE ~ sp_code
  )) %>%
  # change to lower case
  mutate_at(c("treatment", "site"), str_to_lower) %>% 
  mutate(season = case_match(
    season,
    "1-WINTER" ~ "winter",
    "2-SPRING" ~ "spring",
    "3-SUMMER" ~ "summer",
    "4-AUTUMN" ~ "fall"
  ),
  season = fct_relevel(season, "spring", "summer", "fall", "winter")) %>% 
  filter(!(sp_code == "MAPY")) %>% 
  group_by(year, season, site, treatment) %>% 
  summarize(total_npp = sum(npp_season_g_c_m2_day)) %>% 
  unite("season_ID", year, season, site, treatment)
  # # filter to only include continual removal plots and control plots
  # filter(treatment %in% c("continual", "control")) %>% 
  # mutate(exp_dates = case_when(
  #   site == "aque" & date >= aque_start_date_continual & date < aque_after_date_continual ~ "during",
  #   site == "aque" & date >= aque_after_date_continual ~ "after",
  #   site == "napl" & date >= napl_start_date_continual & date < napl_after_date_continual ~ "during",
  #   site == "napl" & date >= napl_after_date_continual ~ "after",
  #   site == "mohk" & date >= mohk_start_date_continual & date < mohk_after_date_continual ~ "during",
  #   site == "mohk" & date >= mohk_after_date_continual ~ "after",
  #   site == "carp" & date >= carp_start_date_continual & date < carp_after_date_continual ~ "during",
  #   site == "carp" & date >= carp_after_date_continual ~ "after"
  # ),
  # exp_dates = fct_relevel(exp_dates, "during", "after")) %>% 
  # # take out all surveys that were before the removal experiment started
  # drop_na(exp_dates) %>% 
  # time_since_columns_continual() %>%
  # group_by(site, year, treatment, quarter, sp_code) %>%
  # mutate(dry_gm2 = mean(dry_gm2),
  #        wm_gm2 = mean(wm_gm2),
  #        density = mean(density)) %>%
  # # take out the "duplicates": only one sampling date per quarter in the dataframe, with values averaged across the two sampling dates
  # slice(1L) %>%
  # ungroup() %>%
  # # take out extraneous columns from time_since_columns_continual()
  # select(!test_min_time_yrs) %>% 
  # comparison_column_continual_new() %>% 
  # kelp_year_column()

# categorical traits 
traits <- read_csv(here("data", 
                        "functional-traits",
                        "joe-traits-lter_2024-09-27.csv"))

# Steneck and Dethier and Littler and Littler traits
coarse_traits <- repmis::source_data("https://www.dropbox.com/scl/fi/s1gb2f3f13ry1oqtzoyyj/00-coarse_traits.csv?rlkey=tgd6j5q3y7bfsdxz7sgs77klb&st=moyfpk93&dl=1") %>% 
  mutate(scientific_name = str_replace(scientific_name, "_", " ")) %>% 
  select(scientific_name, sd_growth_form, ll_func_form)


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
# 58 total "species" from LTE dataset

excluded_spp <- tribble(
  ~ "scientific_name",
  # not present in LTE surveys
  "Amphiroa beauvoisii",
  "Eisenia arborea",
  "Mazzaella spp.",
  "Agardhiella subulata",
  "Scytosiphon lomentaria",
  # species not specific enough
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
    "Unidentified Erect Coralline spp."
  ) %>% 
  left_join(., algae_spp, by = "scientific_name")
# 16 excluded species

algae_taxa <- biomass %>% 
  filter(new_group == "algae") %>% 
  filter(!(scientific_name %in% excluded_spp)) %>% 
  select(sp_code, scientific_name, taxon_phylum:taxon_genus) %>% 
  unique()

# 42 species ultimately in dataset (72% of species)

##########################################################################-
# 5. plot themes ----------------------------------------------------------
##########################################################################-

chloro_col <- "#8AA789"

ochro_col <- "#985E5C"

rhodo_col <- "#4CA2B0"

continual_col <- "#CC7540"

control_col <- "#6D5A18"

low_col <- "#C70000"

medium_col <- "#54662C"

high_col <- "#114C54"

# cluster colors
cluster1 <- "#DE7424" 
cluster2 <- "#EDAD30" 
cluster3 <- "#DDB531" 
cluster4 <- "#AD8D26" 
cluster5 <- "#6A743D" 
cluster6 <- "#525D5C" 
cluster7 <- "#654783"

# Littler & Litter functional form colors
coa_bra_col <- "#C70000"
cru_ff_col <- "#DDB531"
fil_ff_col <- "#DDB531"
joi_cal_col <- "#54662C"
sheet_col <- "#009BB0"
thi_lea_col <- "#114C54"

# Steneck and Dethier growth form colors
art_cal_col <- "#1D457F"
cor_fol_col <- "#61599D"
cor_mac_col <- "#C36377"
lea_mac_col <- "#EB7F54"
fol_col <- "#F2AF4A"

art_cal_col <- "#84A6A2"
cor_fol_col <- "#BE5A47"
cor_mac_col <- "#604A76"
lea_mac_col <- "#C2607F"
fol_col <- "#5D8FBC"

cru_col <- "#DDB531"
fil_gf_col <- "#DDB531"



cluster_cols <- calecopal::cal_palette(name = "superbloom2", n = 7, type = "continuous")
ff_cols <- c("jointed_calcareous" = joi_cal_col,
             "coarsely_branched" = coa_bra_col,
             "thick_leathery" = thi_lea_col,
             "filamentous" = fil_ff_col,
             "sheet" = sheet_col,
             "crustose" = cru_ff_col)
gf_cols <- c("articulated_calcareous" = art_cal_col, 
             "corticated_foliose" = cor_fol_col, 
             "corticated_macrophytes" = cor_mac_col, 
             "leathery_macrophyte" = lea_mac_col, 
             "crustose" = cru_col,
             "filamentous_algae" = fil_gf_col, 
             "foliose" = fol_col)

theme_set(theme_bw() +
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12),
                  panel.grid = element_blank()))

model_preds_aesthetics <- list(
  scale_color_manual(values = c(control = control_col, 
                                continual = continual_col),
                     labels = c(control = "Reference", 
                                continual = "Removal")),
  scale_linetype_manual(values = c(control = "22", 
                                   continual = "solid"),
                        labels = c(control = "Reference", 
                                   continual = "Removal")) 
)

model_preds_theme <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          panel.grid = element_blank())
}




