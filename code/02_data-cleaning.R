##########################################################################-
# Data cleaning
# last modified: 2024-06-20

# This is a script to get data files in order for downstream analysis.
##########################################################################-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 1. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to do this once per session
source(here::here("code", "01_source.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 2. traits -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

photosynthesis_clean <- photosynthesis %>% 
  # take out the juveniles of large brown species that use MAPY as proxy
  filter(pe_spp != "MPJ") %>% 
  # take out CYOS reproductive parts
  filter(pe_spp != "CYOS_R") %>% 
  # create new column to match species
  # if the values for Pmax, alpha, etc. are not for that specific species, remove
  mutate(pe_matched_spp = case_when(
    sp_code == pe_spp ~ "match",
    TRUE ~ "no match"
  )) %>% 
  mutate(chn_matched_spp = case_when(
    sp_code == chn_spp ~ "match",
    TRUE ~ "no match"
  )) %>% 
  filter(pe_matched_spp == "match" | chn_matched_spp == "match") %>% 
  mutate(pmax_new = case_when(
    pe_matched_spp == "match" ~ pmax,
    pe_matched_spp == "no match" ~ NA
  ),
  pe_alpha_new = case_when(
    pe_matched_spp == "match" ~ pe_alpha,
    pe_matched_spp == "no match" ~ NA
  ),
  cn_new = case_when(
    chn_matched_spp == "match" ~ cn,
    chn_matched_spp == "no match" ~ NA
  )) %>% 
  select(scientific_name, pmax_new, pe_alpha_new, cn_new) %>% 
  unique()

traits_clean <- traits %>% 
  # subset of traits
  # select(scientific_name,
  #        thickness, position_to_benthos, stipe, branching, 
  #        blades, blade_category, attachment, calcification) %>% 
  # all traits
  select(scientific_name, 
         size_cm, thickness, position_to_benthos, articulated, stipe,
         midrib, branching, branch_shape, blades, blade_category, longevity,
         coenocyte, attachment, tissue_complexity, growth, calcification) %>%
  filter(!(scientific_name %in% pull(excluded_spp, scientific_name))) %>% 
  left_join(., coarse_traits, by = "scientific_name") %>% 
  left_join(., photosynthesis_clean, by = "scientific_name") 

# turns the cleaned trait data frame into a matrix for dissimilarity stuff
trait_matrix <- traits_clean %>% 
  filter(scientific_name != "Macrocystis pyrifera") %>% 
  column_to_rownames("scientific_name") %>% 
  mutate(across(.cols = thickness:ll_func_form, as_factor)) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 3. community matrices -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. LTE ----------------------------------------------------------------

# broad community matrix
comm_df <- biomass %>% 
  # join with site quality data frame and data frame of full names of site
  left_join(., site_quality, 
            by = "site") %>% 
  left_join(., enframe(sites_full), 
            by = c("site" = "name")) %>% 
  rename(site_full = value) %>% 
  mutate(site_full = fct_relevel(
    site_full, 
    "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>%
  # create new sample ID with treatment
  unite("sample_ID", site, treatment, date, remove = FALSE) # %>% 
  # mutate(across(where(is.character), as.factor))

# metadata for all plots
comm_meta <- comm_df %>% 
  select(sample_ID, season_ID, site, date, year, month, 
         treatment, exp_dates, quarter, time_yrs, 
         time_since_start, time_since_end, kelp_year, 
         comp_1yrs, comp_2yrs, comp_3yrs, quality, site_full, season) %>% 
  unique() %>% 
  mutate(time_since_zero = case_when(
    time_since_end < 0 ~ time_since_start,
    TRUE ~ time_since_end
  )) %>% 
  mutate(treatment = case_match(
    treatment,
    "continual" ~ "Removal",
    "control" ~ "Reference"
  ))


# widening function
widen <- function(df) {
  df %>% 
    # select columns of interest
    select(sample_ID, scientific_name, dry_gm2) %>% 
    # get into wide format for community analysis 
    pivot_wider(names_from = scientific_name, values_from = dry_gm2) %>% 
    # make sample_ID column row names
    column_to_rownames("sample_ID") %>% 
    replace(is.na(.), 0)
}

# getting a list of sites to exclude because they have 0 species

sites_to_exclude <- comm_df %>% 
  # algae only
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  select(sample_ID, scientific_name, dry_gm2) %>% 
  filter(!(scientific_name %in% pull(excluded_spp, scientific_name))) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = scientific_name, values_from = dry_gm2) %>% 
  column_to_rownames("sample_ID") %>% 
  mutate(sum = rowSums(., na.rm = TRUE)) %>% 
  filter(!(sum > 0)) %>% 
  rownames()

comm_mat_algae <- comm_df %>% 
  # algae only
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  # filter(exp_dates == "during") %>% 
  select(sample_ID, scientific_name, dry_gm2) %>% 
  filter(!(scientific_name %in% pull(excluded_spp, scientific_name))) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = scientific_name, values_from = dry_gm2) %>% 
  filter(!(sample_ID %in% sites_to_exclude)) %>% 
  # make the sample_ID column row names
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0) %>% 
  # reordering species to be in the same order as in the trait data frame
  select(rownames(trait_matrix)) %>% 
  as.matrix()

comm_meta_algae <- comm_meta %>% 
  filter(sample_ID %in% rownames(comm_mat_algae))

# comm_mat_continual_algae <- comm_df %>% 
#   filter(new_group == "algae") %>% 
#   # only include sampling from removal plots
#   filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
#   widen()
# 
# comm_mat_control_algae <- comm_df %>% 
#   filter(new_group == "algae") %>% 
#   # only include sampling from control plots
#   filter(sample_ID %in% comm_meta_control$sample_ID) %>% 
#   widen()
  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------------- 4. total biomass ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# total biomass
algae_biomass <- biomass %>% 
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  filter(!(scientific_name %in% c(
    "Halymenia spp.; Schizymenia pacifica",
    "crustose coralline algae spp.",
    "Ectocarpaceae spp.",
    "Ulva spp.; Sponogomorpha spp.",
    "Rhodophyta",
    "Neoptilota spp.; Ptilota spp.; Rhodoptilum spp."
  )))

# delta biomass
# delta_algae_biomass <- algae_biomass %>% 
#   dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
#   group_by(site, year, month, treatment, date) %>% 
#   summarize(total_dry = sum(dry_gm2, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = treatment, values_from = total_dry) %>% 
#   mutate(delta_annual = annual - control,
#          delta_continual = continual - control) 

# long format
# algae_continual_long <- delta_algae_biomass %>% 
#   dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
#   # take out years where continual removal hadn't happened yet
#   drop_na(delta_continual) %>% 
#   select(!delta_continual) %>% 
#   mutate(exp_dates = case_when(
#     # after for annual removal:
#     site == "aque" & date >= aque_after_date_continual ~ "after",
#     site == "napl" & date >= napl_after_date_continual ~ "after",
#     site == "mohk" & date >= mohk_after_date_continual ~ "after",
#     site == "carp" & date >= carp_after_date_continual ~ "after",
#     # everything else is "during" the experiment
#     TRUE ~ "during"
#   ),
#   exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
#   time_since_columns_continual() %>% 
#   kelp_year_column() %>% 
#   comparison_column_continual_new() %>% 
#   # make it longer
#   pivot_longer(cols = c(control, continual)) %>% 
#   # rename columns
#   rename(treatment = name, algae_biomass = value) %>% 
#   # change treatment names
#   mutate(treatment = case_match(treatment, "control" ~ "reference", "continual" ~ "removal")) %>% 
#   # create a new sample ID that is site, year, quarter, treatment
#   unite("sample_ID", site, date, quarter, treatment, remove = FALSE)

# algae_annual_long <- delta_algae_biomass %>% 
#   dplyr::select(site, year, month, date, control, annual, delta_annual) %>% 
#   mutate(exp_dates = case_when(
#     # after for annual removal:
#     site == "aque" & date >= aque_after_date_annual ~ "after",
#     site == "napl" & date >= napl_after_date_annual ~ "after",
#     site == "ivee" & date >= napl_after_date_annual ~ "after",
#     site == "mohk" & date >= mohk_after_date_annual ~ "after",
#     site == "carp" & date >= carp_after_date_annual ~ "after",
#     # everything else is "during" the experiment
#     TRUE ~ "during"
#   ),
#   exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
#   time_since_columns_annual() %>% 
#   kelp_year_column() %>% 
#   comparison_column_annual_new() %>% 
#   # make it longer
#   pivot_longer(cols = c(control, annual)) %>% 
#   # rename columns
#   rename(treatment = name, algae_biomass = value) %>% 
#   # change treatment names
#   mutate(treatment = case_match(treatment, "control" ~ "reference", "continual" ~ "removal")) %>% 
#   # create a new sample ID that is site, year, quarter, treatment
#   unite("sample_ID", site, date, quarter, treatment, remove = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------- 5. excluded sampling events ----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

zero_table <- enframe(sites_to_exclude) %>% 
  separate_wider_delim(cols = value, 
                       delim = "_",
                       names = c("site", "treatment", "date")) %>% 
  select(-name) %>% 
  flextable() %>% 
  autofit()

# save_as_docx(zero_table,
#              path = here("tables", "misc", "zero-sites.docx"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 6. trait matrices ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# categories for each trait (needed to use `mFD` functions)
algae_traits_cat <- colnames(trait_matrix) %>% 
  enframe() %>% 
  select(value) %>% 
  rename(trait_name = value) %>% 
  mutate(trait_type = case_match(
    trait_name,
    "size_cm" ~ "Q", # quantitative
    "thickness" ~ "N", 
    "position_to_benthos" ~ "N", # nominal
    "articulated" ~ "N",
    "stipe" ~ "N",
    "midrib" ~ "N",
    "branching" ~ "N",
    "branch_shape" ~ "N",
    "blades" ~ "N",
    "blade_category" ~ "N",
    "coenocyte" ~ "N",
    "attachment" ~ "N",
    "tissue_complexity" ~ "N",
    "growth" ~ "N",
    "calcification" ~ "N",
    "longevity" ~ "N"
  ))

# reduced trait matrix
trait_matrix_reduced <- trait_matrix %>% 
  select(size_cm, position_to_benthos,
         stipe, branching, branch_shape, blade_category,
         calcification, longevity, attachment, pmax_new, pe_alpha_new, 
         cn_new
         )

# reduced trait categories
algae_traits_cat_reduced <- algae_traits_cat %>% 
  filter(trait_name %in% colnames(trait_matrix_reduced))

