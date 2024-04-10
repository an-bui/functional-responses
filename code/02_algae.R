
traits_clean <- traits %>% 
  select(scientific_name, 
         size_cm, thickness, position_to_benthos, articulated, stipe, 
         midrib, branching, branch_shape, blades, blade_category,
         coenocyte, attachment, tissue_complexity, growth, calcification) %>% 
  filter(!(scientific_name %in% c(
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
  )))

trait_matrix <- traits_clean %>% 
  filter(scientific_name != "Macrocystis pyrifera") %>% 
  column_to_rownames("scientific_name") 

trait_gower <- gowdis(trait_matrix)
trait_pcoa <- wcmdscale(d = trait_gower)
trait_pcoa_scores <- scores(trait_pcoa, choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  left_join(., algae_taxa, by = "scientific_name")

trait_pcoa_plot <- ggplot(trait_pcoa_scores,
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = taxon_phylum, shape = taxon_phylum),
             size = 3,
             alpha = 0.9) +
  scale_color_manual(values = c("Chlorophyta" = chloro_col,
                                "Ochrophyta" = ochro_col,
                                "Rhodophyta" = rhodo_col)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.4, 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.4, 0.5)) +
  guides(color = guide_legend(position = "inside")) +
  theme(legend.title = element_blank(),
        legend.position.inside = c(0.85, 0.1),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid = element_blank())

trait_pcoa_plot

# ggsave(filename = here("figures", 
#                        "trait-ordination", 
#                        paste0("gower-traits_", today(), ".jpg")),
#        plot = trait_pcoa_plot,
#        dpi = 300,
#        width = 6, 
#        height = 6)

spp_trait_data <- scores(trait_pcoa, choices = c(1, 2, 3, 4)) %>% 
  as_tibble(rownames = "scientific_name") %>% 
  # mutate(Dim1_new = Dim1 + 1,
  #        Dim2_new = Dim2 + 2) %>% 
  # select(scientific_name, Dim1_new, Dim2_new) %>% 
  column_to_rownames("scientific_name") %>% 
  as.matrix()

comm_mat_algae_matrix <- comm_mat_algae %>% 
  # putting the columns in the right order
  select(rownames(spp_trait_data)) %>% 
  as.matrix()

site_by_trait <- comm_mat_algae_matrix %*% spp_trait_data

nmds <- metaMDS(comm = site_by_trait, distance = "gower", zerodist = "ignore")
nmds_scores <- scores(nmds, choices = c(1, 2), display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_ID") %>% 
  left_join(., comm_meta_algae, by = "sample_ID") %>% 
  filter(!(sample_ID %in% c("napl_control_2023-05-18", "mohk_continual_2012-11-15")))

with(comm_meta_algae, adonis2(site_by_trait ~ treatment*time_since_end, data = comm_meta_algae, strata = "site"))
adonis2(site_by_trait ~ treatment*time_since_end, data = comm_meta_algae)

ggplot(nmds_scores,
       aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = treatment))
# weird point is napl_control_2023-05-18

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
delta_algae_biomass <- algae_biomass %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  group_by(site, year, month, treatment, date) %>% 
  summarize(total_dry = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = treatment, values_from = total_dry) %>% 
  mutate(delta_annual = annual - control,
         delta_continual = continual - control) 

# long format
algae_continual_long <- delta_algae_biomass %>% 
  dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  select(!delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for annual removal:
    site == "aque" & date >= aque_after_date_continual ~ "after",
    site == "napl" & date >= napl_after_date_continual ~ "after",
    site == "mohk" & date >= mohk_after_date_continual ~ "after",
    site == "carp" & date >= carp_after_date_continual ~ "after",
    # everything else is "during" the experiment
    TRUE ~ "during"
  ),
  exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual_new() %>% 
  # make it longer
  pivot_longer(cols = c(control, continual)) %>% 
  # rename columns
  rename(treatment = name, algae_biomass = value) %>% 
  # change treatment names
  mutate(treatment = case_match(treatment, "control" ~ "reference", "continual" ~ "removal")) %>% 
  # create a new sample ID that is site, year, quarter, treatment
  unite("sample_ID", site, date, quarter, treatment, remove = FALSE)

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



comm_df <- biomass %>% 
  filter(treatment %in% c("control", "continual")) %>% 
  # select columns of interest 
  dplyr::select(site, year, month, treatment, date, new_group, sp_code, dry_gm2) %>% 
  unite("sample_ID_short", site, date, remove = FALSE) %>% 
  # filtered from kelp delta data frame created in upstream script
  # filter(sample_ID_short %in% (delta_continual$sample_ID_short)) %>% 
  exp_dates_column_continual() %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual_new() %>% 
  full_join(., site_quality, by = "site") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename(site_full = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>%
  # create new sample ID with treatment
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # only include 3 year sampling sites
  drop_na(comp_3yrs)

# metadata for all plots
comm_meta <- comm_df %>% 
  select(sample_ID, site, date, year, month, treatment, exp_dates, quarter, time_yrs, time_since_start, time_since_end, kelp_year, comp_1yrs, comp_2yrs, comp_3yrs, quality, site_full) %>% 
  unique()

# metadata for continual removal plots
comm_meta_continual <- comm_meta %>% 
  filter(treatment == "continual")

# metadata for control plots
comm_meta_control <- comm_meta %>% 
  filter(treatment == "control")



widen <- function(df) {
  df %>% 
    # select columns of interest
    select(sample_ID, sp_code, dry_gm2) %>% 
    # get into wide format for community analysis 
    pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
    # make sample_ID column row names
    column_to_rownames("sample_ID") %>% 
    replace(is.na(.), 0)
}





comm_mat_algae <- comm_df %>% 
  # algae only
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  left_join(., algae_spp, by = "sp_code") %>% 
  select(sample_ID, scientific_name, dry_gm2) %>% 
  filter(!(scientific_name %in% c(
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
  ))) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = scientific_name, values_from = dry_gm2) %>% 
  filter(!(sample_ID %in% sites_to_exclude)) %>% 
  # make the sample_ID column row names
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)

sites_to_exclude <- comm_mat_algae %>% 
  mutate(sum = rowSums(., na.rm=TRUE)) %>% 
  filter(sum == 0) %>% 
  select(sum) %>% 
  rownames_to_column("sample_ID") %>% 
  pull(sample_ID)

comm_meta_algae <- comm_meta %>% 
  filter(sample_ID %in% rownames(comm_mat_algae))

comm_mat_continual_algae <- comm_df %>% 
  filter(new_group == "algae") %>% 
  # only include sampling from removal plots
  filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
  widen()

comm_mat_control_algae <- comm_df %>% 
  filter(new_group == "algae") %>% 
  # only include sampling from control plots
  filter(sample_ID %in% comm_meta_control$sample_ID) %>% 
  widen()




