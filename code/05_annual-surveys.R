# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 1. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


source(here::here("code", "02_data-cleaning.R"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 2. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# trait matrix for annual survey
# Anisocladella never actually shows up in the annual surveys, surprisingly
trait_matrix_benthics <- trait_matrix_reduced |> 
  filter(row.names(trait_matrix_benthics) != "Anisocladella pacifica")

# benthics is created in the source script
benthics_comm_df <- benthics |>
  # filter out SCI sites
  filter(!(site %in% c("sctw", "scdi"))) |> 
  # only include species that are in the trait matrix
  filter(scientific_name %in% rownames(trait_matrix_benthics)) |>
  # create a sample ID
  unite("sample_ID", site, year, transect, remove = FALSE)

benthics_comm_df_wide <- benthics_comm_df |> 
  widen() |>   
  # reordering species to be in the same order as in the trait data frame
  select(rownames(trait_matrix_benthics)) |> 
  # filter out surveys where there were no species
  # only 2 surveys excluded: Arroyo Hondo in 2017 and in 2019
  mutate(zero_spp = rowSums(across(where(is.numeric)))) |> 
  filter(zero_spp > 0) |> 
  select(!zero_spp)

# calculating total understory biomass
benthics_biomass <- benthics |>
  filter(new_group == "algae") |> 
  filter(sp_code != "MAPY") |> 
  group_by(sample_ID) |> 
  summarize(total_biomass = sum(dry_gm2, na.rm = TRUE)) |> 
  ungroup() 

# kelp biomass
benthics_kelp <- benthics |> 
  filter(sp_code == "MAPY") |> 
  select(sample_ID, dry_gm2) |> 
  rename("total_kelp_biomass" = "dry_gm2")

# metadata
benthics_comm_meta <- benthics |>
  select(sample_ID,
         site, year, transect) |> 
  unique() |> 
  filter(sample_ID %in% rownames(benthics_comm_df_wide))

# substrate
benthics_substrate <- substrate |> 
  # create a sample ID
  unite("sample_ID", site, year, transect, remove = FALSE) |> 
  # create a new column with "hard" and "soft" substrate
  # these are based on Miller et al. 2018 groupings
  mutate(substrate_new = case_match(
    common_name,
    c("bedrock", "boulder large", "boulder medium", "boulder small", "cobble") ~ "hard",
    c("sand", "shell debris", "shallow sand") ~ "soft"
  )) |> 
  filter(substrate_new == "hard") |> 
  group_by(sample_ID) |> 
  summarize(mean_hard_pc = mean(percent_cover, na.rm = TRUE)) |> 
  ungroup() 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 3. functional diversity metrics --------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. calculating metrics ------------------------------------------------

benthics_fd <- dbFD(x = trait_matrix_benthics,
                    a = benthics_comm_df_wide,
                    corr = "none",
                    print.pco = TRUE)
# some help with distances: https://stat.ethz.ch/pipermail/r-sig-ecology/2016-January/005264.html

# messages:
# Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept. 
# FEVe: Could not be calculated for communities with <3 functionally singular species. 
# FDis: Equals 0 in communities with only one functionally singular species. 
# FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
# FRic: Dimensionality reduction was required. The last 14 PCoA axes (out of 16 in total) were removed. 
# FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) = 0.2428839 
# FDiv: Could not be calculated for communities with <3 functionally singular species. 

benthics_div <- vegan::diversity(x = benthics_comm_df_wide,
                                 index = "simpson") |> 
  enframe() |> 
  rename(simpson = value)

benthics_shan <- vegan::diversity(x = benthics_comm_df_wide,
                                 index = "shannon") |> 
  enframe() |> 
  rename(shannon = value)


# ⟞ b. wrangling ----------------------------------------------------------

benthics_fd_metrics <- benthics_fd$nbsp |> 
  enframe() |> 
  rename(sample_ID = name,
         spp_rich = value) |> 
  left_join(enframe(benthics_fd$FRic), by = c("sample_ID" = "name")) |> 
  rename(fric = value) |> 
  left_join(enframe(benthics_fd$RaoQ), by = c("sample_ID" = "name")) |> 
  rename(raoq = value) |> 
  left_join(enframe(benthics_fd$FDis), by = c("sample_ID" = "name")) |> 
  rename(fdis = value) |>  
  left_join(enframe(benthics_fd$FEve), by = c("sample_ID" = "name")) |> 
  rename(feve = value) |> 
  left_join(benthics_div, by = c("sample_ID" = "name")) |> 
  mutate(redund = simpson - raoq) |> 
  left_join(benthics_shan, by = c("sample_ID" = "name")) |> 
  left_join(benthics_comm_meta, by = "sample_ID") |> 
  left_join(benthics_biomass, by = "sample_ID") |> 
  # rough estimate of NPP from Shannon's paper
  mutate(npp_estimate = 2.77*total_biomass + 0.14) |> 
  left_join(benthics_kelp, by = "sample_ID") |> 
  left_join(benthics_substrate, by = "sample_ID") |> 
  # filter out the 4 surveys that have too few species to calculate FRic
  # Hondo 2020, Hondo 2018, SCTW 2005, SCTW 2023
  drop_na(fric) |> 
  # calculating variation in kelp
  group_by(site) |> 
  mutate(mean_kelp = mean(total_kelp_biomass),
         diff_from_mean = total_kelp_biomass - mean_kelp) |> 
  ungroup() |> 
  drop_na(spp_rich, fric, total_kelp_biomass, npp_estimate) |> 
  mutate(site = as_factor(site),
         transect = as_factor(transect))
  

# ⟞ c. exploratory visualization ------------------------------------------

# ⟞ ⟞ scatter plots -------------------------------------------------------

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp_biomass,
           y = fric)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp_biomass,
           y = feve)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = mean_hard_pc,
           y = spp_rich)) +
  geom_point()


ggplot(data = benthics_fd_metrics,
       aes(x = mean_hard_pc,
           y = npp_estimate)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = fric,
           y = npp_estimate, 
           size = total_kelp_biomass,
           color = total_kelp_biomass)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp_biomass,
           y = spp_rich)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich,
           y = npp_estimate, 
           size = total_kelp_biomass,
           color = total_kelp_biomass)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = shannon,
           y = npp_estimate, 
           size = total_kelp_biomass,
           color = total_kelp_biomass)) +
  geom_point()

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich,
           y = fric)) +
  geom_point() 

ggplot(data = benthics_fd_metrics,
       aes(x = feve,
           y = npp_estimate)) +
  geom_point() 

# ⟞ ⟞ histograms ----------------------------------------------------------

ggplot(data = benthics_fd_metrics,
       aes(x = total_biomass)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = npp_estimate)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = total_kelp_biomass)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = fric)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = fric)) +
  geom_density(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(sample = fric)) +
  stat_qq_line() +
  stat_qq()

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black",
                 bins = 12) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich)) +
  geom_density(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = shannon)) +
  geom_density(fill = "cornflowerblue",
               color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = feve)) +
  geom_density(fill = "cornflowerblue",
               color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 4. modeling ------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. NPP ~ diversity ----------------------------------------------------


# ⟞ ⟞ species richness ----------------------------------------------------

spp_rich_model <- glmmTMB(
  scale(total_biomass, center = FALSE) ~ scale(spp_rich, center = FALSE) + (1|site/transect) + (1|year),
  family = Gamma(link = "log"),
  data = benthics_fd_metrics
)

plot(simulateResiduals(spp_rich_model))

r.squaredGLMM(spp_rich_model)

summary(spp_rich_model)

predict_response(spp_rich_model,
          terms = "spp_rich") |> 
  plot(show_data = TRUE) +
  labs(x = "Species richness",
       y = "Understory biomass",
       title = "Species richness predicts understory biomass") +
  theme(panel.grid = element_blank())


# ⟞ ⟞ functional richness -------------------------------------------------

fric_model <- glmmTMB(
  total_biomass ~ fric + (1|site/transect) + (1|year),
  family = Gamma(link = "log"),
  data = benthics_fd_metrics
)

plot(simulateResiduals(fric_model))

summary(fric_model)

predict_response(fric_model,
          terms = c("fric[0.01:0.4 by = 0.001")) |> 
  plot(show_data = TRUE) +
  labs(x = "Functional richness",
       y = "Understory biomass",
       title = "Functional richness predicts understory biomass") +
  theme(panel.grid = element_blank())


# ⟞ ⟞ kelp biomass --------------------------------------------------------

kelp_model <- glmmTMB(
  total_biomass ~ total_kelp_biomass + (1|site/transect) + (1|year),
  family = Gamma(link = "log"),
  data = benthics_fd_metrics
)

plot(simulateResiduals(kelp_model))

summary(kelp_model)


# ⟞ ⟞ all predictors ------------------------------------------------------

spp_rich_and_fric <- glmmTMB(
  total_biomass ~ spp_rich + fric + total_kelp_biomass +  
    (1|site/transect) + (1|year),
  family = Gamma(link = "log"),
  data = benthics_fd_metrics
)

plot(simulateResiduals(spp_rich_and_fric))

summary(spp_rich_and_fric)

ggpredict(spp_rich_and_fric,
          terms = c("spp_rich")) |> 
  plot(show_data = TRUE)

ggpredict(spp_rich_and_fric,
          terms = c("fric")) |> 
  plot(show_data = TRUE)

ggpredict(spp_rich_and_fric,
          terms = c("total_kelp_biomass")) |> 
  plot(show_data = TRUE)

# ⟞ ⟞ model selection -----------------------------------------------------

AICc(spp_rich_model, 
     fric_model, 
     spp_rich_and_fric, 
     kelp_model) |> 
  arrange(AICc)

# ⟞ b. diversity ~ kelp ---------------------------------------------------

# ggplot(data = benthics_fd_metrics,
#        aes(x = total_kelp,
#            y = spp_rich)) +
#   geom_point() +
#   labs(x = "Total kelp biomass",
#        y = "Species richness")

# asy_model <- nls(spp_rich ~ SSasymp(total_kelp, Asym, R0, lrc),
#                  data = benthics_fd_metrics)
# summary(asy_model)
# ggpredict(asy_model,
#           terms = c("total_kelp[0:10000, by = 100]")) |> 
#   plot(show_data = TRUE)
#   

# ⟞ ⟞ species richness ----------------------------------------------------

# no convergence problems with poisson
spp_rich_kelp_model <- glmmTMB(
  spp_rich ~ total_kelp_biomass + (1|site/transect) + (1|year),
  family = poisson(link = "log"),
  # family = nbinom2(link = "log"),
  data = benthics_fd_metrics)

plot(simulateResiduals(spp_rich_kelp_model))

summary(spp_rich_kelp_model)

ggpredict(spp_rich_kelp_model,
          terms = "total_kelp_biomass[0:3000, by = 100]") |> 
  plot(show_data = TRUE)


# ⟞ ⟞ functional richness -------------------------------------------------

fric_kelp_model <- glmmTMB(
  fric ~ total_kelp_biomass + (1|site/transect) + (1|year),
  # family = beta_family(link = "logit"),
  data = benthics_fd_metrics
  )

# lmerTest::lmer(fric ~ total_kelp_biomass + spp_rich +
#                  (1|site/transect) + (1|year),
#                # family = beta_family(link = "logit"),
#                data = benthics_fd_metrics) |> summary()

plot(simulateResiduals(fric_kelp_model))

summary(fric_kelp_model)

ggpredict(fric_kelp_model,
          terms = "total_kelp_biomass[0:3000, by = 100]") |> 
  plot(show_data = TRUE)

fric_spp_rich_model <- glmmTMB(
  fric ~ spp_rich + (1|site/transect) + (1|year),
  # family = beta_family(link = "logit"),
  data = benthics_fd_metrics
  )

plot(simulateResiduals(fric_spp_rich_model))

summary(fric_spp_rich_model)

ggpredict(fric_spp_rich_model,
          terms = "spp_rich") |> 
  plot(show_data = TRUE)

fric_kelp_spp_rich_model <- glmmTMB(
  fric ~ total_kelp_biomass + spp_rich + (1|site/transect) + (1|year),
  # family = tweedie(link = "log"),
  # family = beta_family(link = "logit"),
  data = benthics_fd_metrics
)

plot(simulateResiduals(fric_kelp_spp_rich_model))

summary(fric_kelp_spp_rich_model)

ggpredict(fric_kelp_spp_rich_model,
          terms = "total_kelp_biomass") |> 
  plot(show_data = TRUE)

ggpredict(fric_kelp_spp_rich_model,
          terms = "spp_rich") |> 
  plot(show_data = TRUE)

# asy_model <- nls(fric ~ SSasymp(total_kelp, Asym, R0, lrc),
#                  data = benthics_fd_metrics)
# hist(resid(asy_model))
# summary(asy_model)
# ggpredict(asy_model,
#           terms = c("total_kelp[0:10000, by = 100]")) |> 
#   plot(show_data = TRUE) +
#   labs(x = "Total kelp biomass",
#        y = "Functional richness",
#        title = "Total kelp biomass has an asymptotic relationship with functional richness") +
#   theme(panel.grid = element_blank())




# benthics_fd_metrics$kelp_01 <- ifelse(
#   benthics_fd_metrics$total_kelp > 0, "present", "absent"
# )
# 
# benthics_fd_metrics$total_kelp_corrected <- ifelse(
#   benthics_fd_metrics$total_kelp == 0, NA, benthics_fd_metrics$total_kelp
# )
# 
# 
# fric_kelp_model <- glmmTMB(log(fric) ~ log(total_kelp_corrected),
#                            # family = gaussian(link = "log"),
#                            # family = beta_family(link = "logit"),
#                            # ziformula = ~1,
#                            na.action = na.exclude,
#                            data = benthics_fd_metrics)
# hist(resid(fric_kelp_model))
# 
# plot(simulateResiduals(fric_kelp_model))
# 
# summary(fric_kelp_model)
# 
# ggplot(data = benthics_fd_metrics,
#        aes(x = log(total_kelp_corrected),
#            y = log(fric))) +
#   geom_point()
# 
# ggpredict(fric_kelp_model,
#           terms = "total_kelp_corrected") |> 
#   plot(show_data = TRUE)
# 
# AICc(asy_model, fric_kelp_model)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------------- 5. piecewise SEM ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

algae_psem <- psem(
  
  # species richness as function of kelp
  glmmTMB(
    spp_rich ~ total_kelp_biomass + (1|site/transect) + (1|year),
    family = poisson(link = "log"),
    na.action = na.omit,
    data = benthics_fd_metrics
  ),
  
  # functional richness as function of kelp and species richness
  glmmTMB(
    fric ~ total_kelp_biomass + spp_rich + (1|site/transect) + (1|year),
    # family = beta_family(link = "logit"),
    na.action = na.omit,
    data = benthics_fd_metrics
  ),

  # NPP as function of functional richness and species richness and hard substrate
  glmmTMB(total_biomass ~ spp_rich + fric + total_kelp_biomass +
            (1|site/transect) + (1|year),
          family = Gamma(link = "log"),
          na.action = na.omit,
          data = benthics_fd_metrics)
  
)

algae_psem_gaussian <- psem(
  
  # species richness as function of kelp
  lmer(spp_rich ~ total_kelp_biomass + 
         (1|site/transect) + (1|year),
          na.action = na.omit,
          data = benthics_fd_metrics),
  
  # functional richness as function of kelp and species richness
  lmer(fric ~ total_kelp_biomass + spp_rich + 
         (1|site/transect) + (1|year),
          na.action = na.omit,
          data = benthics_fd_metrics),
  
  # NPP as function of functional richness and species richness
  lmer(total_biomass ~ spp_rich + fric + total_kelp_biomass +
         (1|site/transect) + (1|year),
          na.action = na.omit,
          data = benthics_fd_metrics)
  
)

summary(algae_psem, conserve = TRUE)
coefs(algae_psem)
summary(algae_psem_gaussian, conserve = TRUE)
coefs(algae_psem_gaussian, standardize.type = "latent.linear")

gaussian_sem_plot <- plot(algae_psem_gaussian)
gaussian_sem_plot
performance::check_singularity(algae_psem_gaussian)

boot <- bootEff(algae_psem_gaussian, R = 100, seed = 666, parallel = "no", ran.eff = "site")
effs <- semEff(boot)
effs$Summary$npp.estimate$Effect
effs$Summary$npp.estimate |> 
  select(c(1, 2, 4, 10, 11))
summary(effs, response = "npp.estimate")

library(lavaan)
library(piecewiseSEM)

set.seed(6)

data <- data.frame(y = runif(100), x = runif(100))

xy_model <- lm(y ~ x, data = data)

spp_rich_kelp_model <- glmmTMB(
  spp_rich ~ total_kelp_biomass + (1|site/transect) + (1|year),
  family = poisson(link = "log"),
  # family = nbinom2(link = "log"),
  data = benthics_fd_metrics)

# perform manual standardization
beta <- summary(spp_rich_kelp_model)$coefficients$cond[2, 1]
(beta_std <- beta * (sd(benthics_fd_metrics$total_kelp_biomass)/sd(benthics_fd_metrics$spp_rich)))



# model <- lmer(scale(fric) ~ scale(total_kelp) + scale(spp_rich) +
#                 (1|site) + (1|year),
#               # family = Gamma(link = "log"),
#               na.action = na.omit,
#               data = benthics_fd_metrics)
# 
# summary(model)
# 
# model <- glmmTMB(fric ~ total_kelp + spp_rich +
#                      (1|site) + (1|year),
#                    family = beta_family(link = "logit"),
#                    na.action = na.omit,
#                    data = benthics_fd_metrics)
# 
# plot(simulateResiduals(model))
# 
# summary(model)
