# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 1. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


source(here::here("code", "02_data-cleaning.R"))


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

benthics_comm_df_wide <- benthic_comm_df %>% 
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
  summarize(total_kelp = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup()

# metadata
benthics_comm_meta <- benthics %>%
  select(sample_ID,
         site, year) %>% 
  unique() %>% 
  filter(sample_ID %in% rownames(benthic_comm_df_wide))

# substrate
benthics_substrate <- substrate 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 3. functional diversity metrics --------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. calculating metrics ------------------------------------------------

benthics_fd <- dbFD(x = trait_matrix_benthics,
                    a = benthic_comm_df_wide,
                    corr = "none",
                    print.pco = TRUE)

# messages:
# Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept. 
# FEVe: Could not be calculated for communities with <3 functionally singular species. 
# FDis: Equals 0 in communities with only one functionally singular species. 
# FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
# FRic: Dimensionality reduction was required. The last 14 PCoA axes (out of 16 in total) were removed. 
# FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) = 0.2428839 
# FDiv: Could not be calculated for communities with <3 functionally singular species. 

benthics_div <- vegan::diversity(x = benthics_comm_df_wide,
                                 index = "simpson") %>% 
  enframe() %>% 
  rename(simpson = value)


# ⟞ b. wrangling ----------------------------------------------------------

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
  # filter out the 4 surveys that have too few species to calculate FRic
  # Hondo 2020, Hondo 2018, SCTW 2005, SCTW 2023
  drop_na(fric) %>% 
  # calculating variation in kelp
  group_by(site) %>% 
  mutate(mean_kelp = mean(total_kelp),
         diff_from_mean = total_kelp - mean_kelp) %>% 
  ungroup()

# ⟞ c. exploratory visualization ------------------------------------------

ggplot(data = benthics_fd_metrics,
       aes(x = fric,
           y = npp_estimate, 
           size = total_kelp,
           color = total_kelp)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich,
           y = npp_estimate, 
           size = total_kelp,
           color = total_kelp)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich,
           y = fric)) +
  geom_point() 

ggplot(data = benthics_fd_metrics,
       aes(x = npp_estimate)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = fric)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 4. modeling ------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. NPP ~ diversity ----------------------------------------------------

spp_rich_model <- glmmTMB(npp_estimate ~ spp_rich + (1|site) + (1|year),
                          family = gaussian(link = "log"),
                          data = benthics_fd_metrics)

plot(simulateResiduals(spp_rich_model))

summary(spp_rich_model)

ggpredict(spp_rich_model,
          terms = "spp_rich") %>% 
  plot(show_data = TRUE) +
  labs(x = "Species richness",
       y = "NPP (estimate)",
       title = "Species richness predicts understory NPP") +
  theme(panel.grid = element_blank())

fric_model <- glmmTMB(npp_estimate ~ fric + (1|site) + (1|year),
                          family = gaussian(link = "log"),
                          data = benthics_fd_metrics)

plot(simulateResiduals(fric_model))

summary(fric_model)

ggpredict(fric_model,
          terms = c("fric[0.01:0.42 by = 0.01]")) %>% 
  plot(show_data = TRUE) +
  labs(x = "Functional richness",
       y = "NPP (estimate)",
       title = "Functional richness predicts understory NPP") +
  theme(panel.grid = element_blank())

AICc(spp_rich_model, fric_model)

r.squaredGLMM(spp_rich_model)
r.squaredGLMM(fric_model)


# ⟞ b. diversity ~ kelp ---------------------------------------------------

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp,
           y = spp_rich)) +
  geom_point() +
  labs(x = "Total kelp biomass",
       y = "Species richness")

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp,
           y = fric)) +
  geom_point()

asy_model <- nls(spp_rich ~ SSasymp(total_kelp, Asym, R0, lrc),
                 data = benthics_fd_metrics)
summary(asy_model)
ggpredict(asy_model,
          terms = c("total_kelp[0:10000, by = 100]")) %>% 
  plot(show_data = TRUE)

spp_rich_kelp_model <- glmmTMB(spp_rich ~ total_kelp,
                               family = poisson(link = "log"),
                               data = benthics_fd_metrics)

plot(simulateResiduals(spp_rich_kelp_model))

summary(spp_rich_kelp_model)

ggpredict(spp_rich_kelp_model,
          terms = "total_kelp[0:10000, by = 100]") %>% 
  plot(show_data = TRUE)

asy_model <- nls(fric ~ SSasymp(total_kelp, Asym, R0, lrc),
                 data = benthics_fd_metrics)
hist(resid(asy_model))
summary(asy_model)
ggpredict(asy_model,
          terms = c("total_kelp[0:10000, by = 100]")) %>% 
  plot(show_data = TRUE) +
  labs(x = "Total kelp biomass",
       y = "Functional richness",
       title = "Total kelp biomass has an asymptotic relationship with functional richness") +
  theme(panel.grid = element_blank())


fric_kelp_model <- glmmTMB(fric ~ total_kelp,
                           family = beta_family(link = "logit"),
                           ziformula = ~1,
                           data = benthics_fd_metrics)
hist(resid(fric_kelp_model))

plot(simulateResiduals(fric_kelp_model))

summary(fric_kelp_model)

ggpredict(fric_kelp_model,
          terms = "total_kelp") %>% 
  plot(show_data = TRUE)

AICc(asy_model, fric_kelp_model)

