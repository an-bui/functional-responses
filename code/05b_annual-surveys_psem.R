# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 1. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# general use
library(here)
library(tidyverse)
library(janitor)
library(readxl)
library(repmis)

# analysis
library(FD)
library(mFD)
library(glmmTMB)
library(lmerTest)
library(DHARMa)
library(MuMIn)
library(car)


guilds <- read_csv(here::here("data", "SBC-LTE", "LTE_guild_data.csv")) %>% 
  mutate(sp.code = replace_na(sp.code, "Nandersoniana")) %>% 
  rename("new_group" = biomass.guild) %>% 
  # for whatever reason Yellowtail Rockfish are not in the guild csv
  add_row(sp.code = "SFLA", new_group = "fish", diversity.guild = "fish")

photosynthesis <- read_csv(here(
  "data", 
  "SBC-LTER", 
  "knb-lter-sbc.127.3",
  "Algal_Biomass_Relationships_NPP_calculation_20210113.csv"
)) %>% 
  clean_names()

benthics <- read_csv(
  here("data", 
       "SBC-LTER",
       "knb-lter-sbc.50.17",
       "Annual_All_Species_Biomass_at_transect_20240823.csv")
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
  mutate_at(c("group", "mobility", "growth_morph", "site"), str_to_lower) %>% 
  left_join(., guilds, by = c("sp_code" = "sp.code")) %>% 
  # create a sample ID
  unite("sample_ID", site, year, sep = "_", remove = FALSE)

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

# categorical traits 
traits <- read_csv(here("data", 
                        "functional-traits",
                        "joe-traits-lter_2024-11-14.csv"))

# Steneck and Dethier and Littler and Littler traits
coarse_traits <- read_rds(here("data",
                               "functional-traits",
                               "coarse-traits.RDS"))

algae_spp <- benthics %>% 
  filter(new_group == "algae") %>% 
  select(scientific_name, sp_code) %>% 
  unique()

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
  "Unidentified Erect Coralline spp.",
  "Cryptopleura spp."
) %>% 
  left_join(., benthics, by = "scientific_name")

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
         calcification, longevity, attachment, # pmax_new, pe_alpha_new, 
         cn_new
  )

# reduced trait categories
algae_traits_cat_reduced <- algae_traits_cat %>% 
  filter(trait_name %in% colnames(trait_matrix_reduced))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 2. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# trait matrix for annual survey
# Anisocladella never actually shows up in the annual surveys, surprisingly
trait_matrix_benthics <- trait_matrix_reduced %>% 
  filter(rownames(.) != "Anisocladella pacifica")

# benthics is created in the source script
benthics_comm_df <- benthics %>%
  # only include species that are in the trait matrix
  filter(scientific_name %in% rownames(trait_matrix_benthics)) %>%
  # create a sample ID
  unite("sample_ID", site, year, remove = FALSE) %>% 
  # sum biomass for all transects at a site in a given year
  group_by(sample_ID, scientific_name) %>% 
  summarize(dry_gm2 = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup()

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

benthics_comm_df_wide <- benthics_comm_df %>% 
  widen() %>%   
  # reordering species to be in the same order as in the trait data frame
  select(rownames(trait_matrix_benthics)) %>% 
  # filter out surveys where there were no species
  # only 2 surveys excluded: Arroyo Hondo in 2017 and in 2019
  mutate(zero_spp = rowSums(across(where(is.numeric)))) %>% 
  filter(zero_spp > 0) %>% 
  select(!zero_spp)

# calculating total understory biomass
benthics_biomass <- benthics %>%
  filter(new_group == "algae") %>% 
  filter(sp_code != "MAPY") %>% 
  group_by(sample_ID) %>% 
  summarize(total_biomass = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() 

# kelp biomass
benthics_kelp <- benthics %>% 
  filter(sp_code == "MAPY") %>% 
  group_by(sample_ID) %>% 
  summarize(total_kelp_biomass = sum(dry_gm2, na.rm = TRUE),
            total_kelp_density = sum(density, na.rm = TRUE)) %>% 
  ungroup()

# metadata
benthics_comm_meta <- benthics %>%
  select(sample_ID,
         site, year) %>% 
  unique() %>% 
  filter(sample_ID %in% rownames(benthics_comm_df_wide))


benthics_fd <- dbFD(x = trait_matrix_benthics,
                    a = benthics_comm_df_wide,
                    corr = "none",
                    print.pco = TRUE)

benthics_div <- vegan::diversity(x = benthics_comm_df_wide,
                                 index = "simpson") %>% 
  enframe() %>% 
  rename(simpson = value)

benthics_fd_metrics <- benthics_fd$nbsp %>% 
  enframe() %>% 
  rename(sample_ID = name,
         spp_rich = value) %>% 
  left_join(., enframe(benthics_fd$FRic), by = c("sample_ID" = "name")) %>% 
  rename(fric = value) %>% 
  left_join(., enframe(benthics_fd$RaoQ), by = c("sample_ID" = "name")) %>% 
  rename(raoq = value) %>% 
  left_join(., enframe(benthics_fd$FDis), by = c("sample_ID" = "name")) %>% 
  rename(fdis = value) %>%  
  left_join(., enframe(benthics_fd$FEve), by = c("sample_ID" = "name")) %>% 
  rename(feve = value) %>% 
  left_join(., benthics_div, by = c("sample_ID" = "name")) %>% 
  mutate(redund = simpson - raoq) %>% 
  left_join(., benthics_comm_meta, by = "sample_ID") %>% 
  left_join(., benthics_biomass, by = "sample_ID") %>% 
  # rough estimate of NPP from Shannon's paper
  mutate(npp_estimate = 2.77*total_biomass + 0.14) %>% 
  left_join(., benthics_kelp, by = "sample_ID") %>% 
  # left_join(., benthics_substrate, by = "sample_ID") %>% 
  # filter out the 4 surveys that have too few species to calculate FRic
  # Hondo 2020, Hondo 2018, SCTW 2005, SCTW 2023
  drop_na(fric) %>% 
  # calculating variation in kelp
  group_by(site) %>% 
  mutate(mean_kelp = mean(total_kelp_biomass),
         diff_from_mean = total_kelp_biomass - mean_kelp) %>% 
  ungroup()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------------- 3. SEM --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(piecewiseSEM)

algae_psem <- psem(
  
  # species richness as function of kelp
  glmmTMB(spp_rich ~ total_kelp_biomass + 
            (1|site) + (1|year),
          family = nbinom2(link = "log"),
          na.action = na.omit,
          data = benthics_fd_metrics),
  
  # functional richness as function of kelp and species richness and hard substrate
  glmmTMB(fric ~ total_kelp_biomass + spp_rich +
            (1|site) + (1|year),
          family = beta_family(link = "logit"),
          na.action = na.omit,
          data = benthics_fd_metrics),
  
  # NPP as function of functional richness and species richness and hard substrate
  glmmTMB(npp_estimate ~ spp_rich + fric +
            (1|site) + (1|year),
          family = Gamma(link = "log"),
          na.action = na.omit,
          data = benthics_fd_metrics)
  
)

algae_psem_gaussian <- psem(
  
  # species richness as function of kelp
  lmer(spp_rich ~ total_kelp_biomass + 
         (1|site) + (1|year),
       na.action = na.omit,
       data = benthics_fd_metrics),
  
  # functional richness as function of kelp and species richness
  lmer(fric ~ total_kelp_biomass + spp_rich + 
         (1|site) + (1|year),
       na.action = na.omit,
       data = benthics_fd_metrics),
  
  # NPP as function of functional richness and species richness
  lmer(npp_estimate ~ spp_rich + fric +
         (1|site) + (1|year),
       na.action = na.omit,
       data = benthics_fd_metrics)
  
)

summary(algae_psem, conserve = TRUE)
summary(algae_psem_gaussian, conserve = TRUE)
gaussian_sem_plot <- plot(algae_psem_gaussian)
gaussian_sem_plot
