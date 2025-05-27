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
  rownames_to_column() |> 
  filter(rowname != "Anisocladella pacifica") |> 
  column_to_rownames("rowname")

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
         transect = as_factor(transect)) |> 
  mutate(fric_stand = standardize(fric),
         spp_rich_stand = standardize(spp_rich),
         kelp_stand = standardize(total_kelp_biomass),
         macroalgae_stand = standardize(total_biomass)) |> 
  mutate(total_biomass_log = log(total_biomass))
  

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

ggplot(data = benthics_fd_metrics,
       aes(x = kelp_stand,
           y = macroalgae_stand)) +
  geom_point() 

ggplot(data = benthics_fd_metrics,
       aes(x = kelp_stand,
           y = spp_rich_stand)) +
  geom_point() 

ggplot(data = benthics_fd_metrics,
       aes(x = kelp_stand,
           y = fric_stand)) +
  geom_point() 

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich_stand,
           y = fric_stand)) +
  geom_point() 

# ⟞ ⟞ histograms ----------------------------------------------------------

ggplot(data = benthics_fd_metrics,
       aes(x = total_biomass)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = log(total_biomass))) +
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
       aes(sample = total_biomass)) +
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

ggplot(data = benthics_fd_metrics,
       aes(x = macroalgae_stand)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black",
                 bins = 15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = kelp_stand)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black",
                 bins = 15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = fric_stand)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black",
                 bins = 15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggplot(data = benthics_fd_metrics,
       aes(x = spp_rich_stand)) +
  geom_histogram(fill = "cornflowerblue",
                 color = "black",
                 bins = 15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------- 4. individual models -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. complementarity ----------------------------------------------------

# These models represent paths from kelp --> diversity (functional and species)
# --> understory biomass. The individual models that will go into the SEM are:
# 1. functional richness as a function of giant kelp biomass and species richness
# 2. species richness as a function of giant kelp biomass
# 3. understory biomass (natural log transformed) as a function of functional 
#    richness, species richness, and giant kelp biomass
# We added the pathway between giant kelp biomass and understory biomass because
# the missing pathway was identified in the SEM as being significant.

# ⟞ ⟞ model 1. functional richness ----------------------------------------

# model
comp_fric_model <- lmer(
  fric ~ total_kelp_biomass + spp_rich + (1|site/transect) + (1|year),
  data = benthics_fd_metrics
)

# DHARMa residuals
plot(simulateResiduals(comp_fric_model))

# model summary
# no effect of giant kelp biomass, positive effect of species richness
summary(comp_fric_model)

# model predictions given giant kelp biomass
predict_response(comp_fric_model,
                 terms = c("total_kelp_biomass")) |> 
  plot(show_data = TRUE) +
  labs(x = "Giant kelp biomass",
       y = "Understory functional richness") +
  theme(panel.grid = element_blank())

# model predictions given species richness
predict_response(comp_fric_model,
                 terms = c("spp_rich")) |> 
  plot(show_data = TRUE) +
  labs(x = "Understory species richness",
       y = "Understory functional richness") +
  theme(panel.grid = element_blank())

# ⟞ ⟞ model 2. species richness -------------------------------------------

# no convergence problems with poisson
# poisson better than regular gaussian
comp_spp_rich_model <- glmmTMB(
  spp_rich ~ total_kelp_biomass + (1|site/transect) + (1|year),
  family = poisson(link = "log"),
  data = benthics_fd_metrics
)

# DHARMa residuals
plot(simulateResiduals(comp_spp_rich_model))

# model summary
# negative effect of giant kelp biomass
summary(comp_spp_rich_model)

# model predictions given giant kelp biomass
ggpredict(comp_spp_rich_model,
          terms = "total_kelp_biomass[0:3000, by = 100]") |> 
  plot(show_data = TRUE) +
  labs(x = "Giant kelp biomass",
       y = "Understory species richness")

# ⟞ ⟞ model 3. understory biomass -----------------------------------------

# model
comp_understory_biomass_model <- lmer(
  log(total_biomass) ~ total_kelp_biomass + spp_rich + fric +
    (1|site/transect) + (1|year),
  data = benthics_fd_metrics
)

# DHARMa residuals
plot(simulateResiduals(comp_understory_biomass_model))

# model summary
# negative effect of giant kelp biomass (barely)
# positive effect of species richness
# positive effect of functional richness
summary(comp_understory_biomass_model)

# model predictions given species richness
ggpredict(comp_understory_biomass_model,
          terms = c("spp_rich")) |> 
  plot(show_data = TRUE) +
  labs(x = "Understory species richness",
       y = "Understory biomass")

# model predictions given functional richness
ggpredict(comp_understory_biomass_model,
          terms = c("fric")) |> 
  plot(show_data = TRUE) +
  labs(x = "Understory functional richness",
       y = "Understory biomass")

# model predictions given giant kelp biomass
ggpredict(comp_understory_biomass_model,
          terms = c("total_kelp_biomass")) |> 
  plot(show_data = TRUE) +
  labs(x = "Giant kelp biomass",
       y = "Understory biomass")

# ⟞ b. selection effect ---------------------------------------------------

# These models represent paths from kelp --> diversity (functional and species)
# and understory biomass. There are direct pathways between understory biomass
# and diversity whereby understory biomass influences species and functional
# diversity, which is the opposite of the models in section 4a. The individual
# models that will go into the SEM are:
# 1. functional richness as a function of giant kelp biomass, understory 
#    biomass (log transformed), and species richness
# 2. species richness as a function of giant kelp biomass and understory biomass
#    (log transformed)
# 3. understory biomass (log transformed) as a function of giant kelp biomass
# For any models using understory biomass as predictors, I used the 
# pre-transformed variable instead of transforming within the model call. In 
# the SEM, I transformed within the model.

# ⟞ ⟞ model 1. functional richness ----------------------------------------

# better with gaussian than beta family
sele_fric_model <- lmer(
  fric ~ total_kelp_biomass + total_biomass_log + spp_rich +
    (1|site/transect) + (1|year),
  data = benthics_fd_metrics
)

# DHARMa residuals
plot(simulateResiduals(sele_fric_model))

# model summary
# no effect of giant kelp biomass
# positive effects of understory biomass (natural log transformed) and species richness
summary(sele_fric_model)

# model predictions given giant kelp biomass
ggpredict(sele_fric_model,
          terms = "total_kelp_biomass") |> 
  plot(show_data = TRUE) +
  labs(x = "Giant kelp biomass",
       y = "Understory functional richness",
       title = "No effect of giant kelp biomass on \n understory functional richness")

# model predictions given log transformed understory biomass
ggpredict(sele_fric_model,
          terms = "total_biomass_log") |> 
  plot(show_data = TRUE) +
  labs(x = "Understory biomass (natural log transformed)",
       y = "Understory functional richness")

# model predictions given species richness
ggpredict(sele_fric_model,
          terms = "spp_rich") |> 
  plot(show_data = TRUE) +
  labs(x = "Understory species richness",
       y = "Understory functional richness")


# ⟞ ⟞ model 2. species richness -------------------------------------------

# model
sele_spp_rich_model <- glmmTMB(
  spp_rich ~ total_kelp_biomass + total_biomass_log +
    (1|site/transect) + (1|year),
  family = poisson(link = "log"),
  data = benthics_fd_metrics
)

# DHARMa residuals
plot(simulateResiduals(sele_spp_rich_model))

# model summary
# no effect of giant kelp biomass
# positive effect of log understory biomass
summary(sele_spp_rich_model)

# model predictions given giant kelp biomass
ggpredict(sele_spp_rich_model,
          terms = "total_kelp_biomass") |> 
  plot(show_data = TRUE) +
  labs(x = "Giant kelp biomass",
       y = "Understory species richness",
       title = "No effect of giant kelp biomass on \n 
       understory species richness")

# model predictions given understory biomass
ggpredict(sele_spp_rich_model,
          terms = "total_biomass_log") |> 
  plot(show_data = TRUE) +
  labs(x = "Understory biomass (natural log transformed)",
       y = "Understory species richness")

# ⟞ ⟞ model 3. understory biomass -----------------------------------------

# model
sele_understory_biomass_model <- lmer(
  log(total_biomass) ~ total_kelp_biomass +
    (1|site/transect) + (1|year),
  data = benthics_fd_metrics
)

# DHARMa residuals
plot(simulateResiduals(sele_understory_biomass_model))

# model summary
# negative effect of giant kelp biomass
summary(sele_understory_biomass_model)

# model predictions given giant kelp biomass
ggpredict(sele_understory_biomass_model,
          terms = "total_kelp_biomass") |> 
  plot(show_data = TRUE) +
  labs(x = "Giant kelp biomass",
       y = "Understory biomass") 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------------- 5. piecewise SEM ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. complementarity ----------------------------------------------------

comp_model <- psem(
  
  # 1. functional richness as a function of giant kelp biomass and species richness
  lmer(
    fric ~ total_kelp_biomass + spp_rich + (1|site/transect) + (1|year),
    data = benthics_fd_metrics
  ),
  
  # 2. species richness as a function of giant kelp biomass
  glmmTMB(
    spp_rich ~ total_kelp_biomass + (1|site/transect) + (1|year),
    family = poisson(link = "log"),
    data = benthics_fd_metrics
  ),
  
  # 3. understory biomass (natural log transformed) as a function of functional 
  #    richness, species richness, and giant kelp biomass
  lmer(
    log(total_biomass) ~ total_kelp_biomass + spp_rich + fric +
      (1|site/transect) + (1|year),
    data = benthics_fd_metrics
  )
)

summary(comp_model)

coefs(comp_model)

plot(sele_model)

# ⟞ b. selection effect ---------------------------------------------------

sele_model <- psem(
  # 1. functional richness as a function of giant kelp biomass, understory 
  #    biomass (log transformed), and species richness
  lmer(
    fric ~ total_kelp_biomass + log(total_biomass) + spp_rich +
      (1|site/transect) + (1|year),
    data = benthics_fd_metrics
  ),
  
  # 2. species richness as a function of giant kelp biomass and understory biomass
  #    (log transformed)
  glmmTMB(
    spp_rich ~ total_kelp_biomass + log(total_biomass) +
      (1|site/transect) + (1|year),
    family = poisson(link = "log"),
    data = benthics_fd_metrics
  ),
  
  # 3. understory biomass (log transformed) as a function of giant kelp biomass
  lmer(
    log(total_biomass) ~ total_kelp_biomass +
      (1|site/transect) + (1|year),
    data = benthics_fd_metrics
  )
)

summary(sele_model)

coefs(sele_model)

plot(sele_model)

# ⟞ c. saving coefficients ------------------------------------------------

coefs_all <- bind_rows(
  coefs(comp_model) |> clean_names() |> mutate(model = "complementarity"),
  coefs(sele_model) |> clean_names() |> mutate(model = "selection effect")
) |> 
  rename("sig_stars" = "x") |> 
  relocate(model, .before = response) |> 
 # mutate(across(where(is.numeric),~round(., digits = 2))) |> 
  as_flextable() |> 
  merge_v(j = ~ model) |> 
  colformat_num(j = 4:9, 
                digits = 3) |> 
  autofit() |> 
  fit_to_width(10)
coefs_all

save_as_docx(coefs_all,
             path = here("tables", "SEM", paste0("SEM-coefs_", today(), ".docx")),
             pr_section = prop_section(
               page_size = page_size(
                 orient = "landscape"
             )))
