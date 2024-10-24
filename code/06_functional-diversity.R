# this is junk!!!

# Functional diversity revealed by removal experiments Diaz 2003

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(here::here("code", "02_data-cleaning.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 1. excluded species and surveys --------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to figure out what proportion of the data set is
# excluded for purposes of identification (observations not identified to 
# fine taxonomic resolution to be assigned trait values) or because there were
# too few species observed in the survey.

# ⟞ a. excluded species ---------------------------------------------------

# This relies on the `excluded_spp` object created in `01_source.R`.

# Double checked that calculations of biomass in excluded and included groups
# made sense; they did.

excluded_biomass <- comm_df %>% 
  filter(new_group == "algae") %>% 
  mutate(included = case_when(
    sp_code %in% pull(excluded_spp, sp_code) ~ "no",
    TRUE ~ "yes"
  )) %>% 
  select(sample_ID, sp_code, dry_gm2, included) %>% 
  group_by(sample_ID, included) %>% 
  summarize(biomass = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(sample_ID) %>% 
  mutate(total_biomass = sum(biomass),
         prop = biomass/total_biomass)

# histogram of excluded biomass proportion
ggplot(data = excluded_biomass %>% filter(included == "no"),
       aes(x = prop)) +
  geom_histogram(bins = 14,
                 fill = "cornflowerblue",
                 color = "black") +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125)) +
  labs(x = "Proportion of excluded biomass \n(excluded biomass/total biomass)",
       title = "For most surveys, the biomass of excluded species is less than 
       25% of total survey biomass")

# histogram of included biomass proportion (mirror image of previous histogram)
ggplot(data = excluded_biomass %>% filter(included == "yes"),
       aes(x = prop)) +
  geom_histogram(bins = 14,
                 fill = "darkgreen",
                 color = "black") +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125)) +
  labs(x = "Proportion of included biomass \n(included biomass/total biomass)",
       title = "For most surveys, the biomass of included surveys is more than
       75% of total survey biomass")

# what is the proportion of surveys where > 25% of species are excluded?
# (and thus, traits do not capture the majority of community biomass that
# is performing some ecosystem function)
excluded_count <- excluded_biomass %>% 
  filter(included == "no") %>% 
  mutate(more_than_25 = case_when(
    prop > 0.25 ~ "more",
    TRUE ~ "less"
  )) 

excluded_count %>% 
  group_by(more_than_25) %>% 
  count()
# for 104 out of 412 surveys, more than 25% of species are excluded. This is 
# ~ 25% of surveys.

# what is the distribution of sites where more than 25% of species are excluded?
# what sites does this occur at? during what period of the experimental removal?

more_than_25_sites <- excluded_count %>% 
  # filter only for sites where > 25% of species are excluded
  filter(more_than_25 == "more") %>% 
  select(sample_ID) %>% 
  left_join(., comm_meta, by = "sample_ID") %>% 
  group_by(site_full, exp_dates, treatment) %>% 
  count() %>% 
  ungroup() %>% 
  complete(site_full, exp_dates, treatment,
           fill = list(n = 0)) %>% 
  mutate(exp_dates = case_match(
    exp_dates,
    "during" ~ "During removal",
    "after" ~ "After removal"
  ),
  treatment = case_match(
    treatment,
    "continual" ~ "Removal",
    "control" ~ "Reference"
  )) %>% 
  flextable() %>% 
  merge_v(j = c("site_full", "exp_dates")) %>% 
  valign(j = c("site_full", "exp_dates"),
         part = "body",
         valign = "top") %>% 
  set_header_labels(values = c(
    site_full = "Site",
    exp_dates = "Removal",
    treatment = "Plot",
    n = "Count"
  )) %>% 
  autofit()

more_than_25_sites

# ⟞ b. excluded surveys ---------------------------------------------------

# What is the distribution of species number across surveys?

survey_spp_count <- comm_df %>% 
  # algae only
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  # filter(exp_dates == "during") %>% 
  select(sample_ID, scientific_name, dry_gm2) %>% 
  # exclude "excluded species"
  filter(!(scientific_name %in% pull(excluded_spp, scientific_name))) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = scientific_name, values_from = dry_gm2) %>% 
  # replace NAs with 0
  replace(is.na(.), 0) %>% 
  # changing to occurrences
  mutate_if(is.numeric, ~1 * (. > 0)) %>% 
  mutate(spp_sum = rowSums(across(where(is.numeric)))) 

survey_count <- survey_spp_count %>% 
  group_by(spp_sum) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(total_surveys = sum(n),
         prop_surveys = n/total_surveys)

ggplot(data = survey_count,
       aes(x = spp_sum,
           y = n)) +
  geom_col(fill = "orange",
           color = "black") +
  scale_x_continuous(breaks = seq(from = 0, to = 20, by = 1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
  labs(x = "Number of species per survey",
       y = "Count",
       title = "Most surveys have between 8-9 species",
       subtitle = "~11% of surveys have 0-1 species")

few_spp_surveys <- survey_spp_count %>% 
  filter(spp_sum < 3) %>% 
  select(sample_ID, spp_sum) %>% 
  left_join(., comm_meta, by = "sample_ID")
# 47 surveys with 0-1 species
# 28 surveys with 2 species
# 23 surveys with 3 species

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 2. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# community matrix transformed to occurrences (1 = species is present)
comm_mat_bin <- comm_mat_algae %>% 
  as.data.frame() %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) %>% 
  as.matrix()

# reduced surveys
comm_mat_algae_reduced <- comm_mat_algae %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_ID") %>% 
  filter(!(sample_ID %in% pull(few_spp_surveys, sample_ID))) %>% 
  column_to_rownames("sample_ID") %>% 
  as.matrix()

comm_mat_algae_reduced_bin <- comm_mat_algae %>% 
  as.data.frame() %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) %>% 
  rownames_to_column("sample_ID") %>% 
  filter(!(sample_ID %in% pull(few_spp_surveys, sample_ID))) %>% 
  column_to_rownames("sample_ID") %>% 
  as.matrix()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------- 3. diversity metrics -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. Gower + PCoA -------------------------------------------------------

# This section includes code to generate the Gower dissimilarity matrix from 
# the trait matrix. The Gower matrix then goes into a Principal Coordinates
# Analysis (PCoA) to generate synthetic traits from the Gower matrix.
# Ultimately, the species values along PCoA axes will be the "traits" in 
# downstream analyses.

# generate Gower matrix
trait_gower <- gowdis(trait_matrix_reduced)

# doing PCoA to get dimensions
trait_pcoa <- ape::pcoa(D = trait_gower)

# extracting proportion of inertia explained by each axis
proportion_inertia <- trait_pcoa$values %>% 
  rownames_to_column("axis") %>% 
  mutate(axis = as_factor(axis))

# plotting for visualization
prop_inertia_plot <- ggplot(data = proportion_inertia,
                            aes(x = axis,
                                y = Relative_eig)) +
  geom_point()
# prop_inertia_plot

# extracting scores
spp_pcoa_scores <- trait_pcoa$vectors %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  left_join(., algae_taxa, by = "scientific_name") %>% 
  left_join(., traits_clean, by = "scientific_name")

axes_12_attachment <- ggplot(data = spp_pcoa_scores,
                  aes(x = Axis.1,
                      y = Axis.2,
                      color = attachment,
                      shape = attachment)) +
  geom_point(size = 3)
# axes_12_attachment

axes_12_calcification <- ggplot(data = spp_pcoa_scores,
                                aes(x = Axis.1,
                                    y = Axis.2,
                                    color = calcification,
                                    shape = calcification)) +
  geom_point(size = 3)
# axes_12_calcification

axes_12_position <- ggplot(data = spp_pcoa_scores,
                           aes(x = Axis.1,
                               y = Axis.2,
                               color = position_to_benthos,
                               shape = position_to_benthos)) +
  geom_point(size = 3)
# axes_12_position

axes_12_phyla <- ggplot(data = spp_pcoa_scores,
                        aes(x = Axis.1,
                            y = Axis.2,
                            color = taxon_phylum,
                            shape = taxon_phylum)) +
  geom_point(size = 3)
# axes_12_phyla

# ⟞ a. `FD` ---------------------------------------------------------------

# ⟞ ⟞ i. calculating metrics ----------------------------------------------

# This section includes code to calculate functional diversity metrics using 
# `FD::dbFD()`. It takes the original trait matrix and goes through the Gower
# and PCoA steps internally, choosing 2 axes.

algae_fd <- dbFD(x = trait_matrix,
                 a = comm_mat_algae,
                 corr = "none",
                 print.pco = TRUE)

algae_fd_bin <- dbFD(x = trait_matrix,
                 a = comm_mat_bin,
                 corr = "none",
                 print.pco = TRUE)

algae_fd_reduced <- dbFD(x = trait_matrix_reduced,
                     a = comm_mat_algae,
                     corr = "none",
                     print.pco = TRUE)

# messages:
# Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept. 
# FEVe: Could not be calculated for communities with <3 functionally singular species. 
# FDis: Equals 0 in communities with only one functionally singular species. 
# FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
# FRic: Dimensionality reduction was required. The last 16 PCoA axes (out of 18 in total) were removed. 
# FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) = 0.2849965 
# FDiv: Could not be calculated for communities with <3 functionally singular species. 

algae_div <- vegan::diversity(x = comm_mat_algae,
                              index = "simpson") %>% 
  enframe() %>% 
  rename(simpson = value)

ll_group_biomass <- biomass %>% 
  filter(new_group == "algae") %>% 
  filter(!(sp_code %in% pull(excluded_spp, sp_code))) %>%
  left_join(., coarse_traits, by = "scientific_name") %>% 
  group_by(year, season, site, treatment, ll_func_form) %>% 
  summarize(group_bio = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite("season_ID", year, season, site, treatment, sep = "_") %>% 
  pivot_wider(names_from = ll_func_form,
              values_from = group_bio)

sd_group_biomass <- biomass %>% 
  filter(new_group == "algae") %>% 
  filter(!(sp_code %in% pull(excluded_spp, sp_code))) %>%
  left_join(., coarse_traits, by = "scientific_name") %>% 
  group_by(year, season, site, treatment, sd_growth_form) %>% 
  summarize(group_bio = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite("season_ID", year, season, site, treatment, sep = "_") %>% 
  pivot_wider(names_from = sd_growth_form,
              values_from = group_bio)

fd_metrics <- algae_fd$nbsp %>% 
  enframe() %>% 
  rename(sample_ID = name,
         spp_rich = value) %>% 
  left_join(., enframe(algae_fd$FRic), by = c("sample_ID" = "name")) %>% 
  rename(fric = value) %>% 
  left_join(., enframe(algae_fd$RaoQ), by = c("sample_ID" = "name")) %>% 
  rename(raoq = value) %>% 
  left_join(., enframe(algae_fd$FDis), by = c("sample_ID" = "name")) %>% 
  rename(fdis = value) %>%  
  left_join(., enframe(algae_fd$FEve), by = c("sample_ID" = "name")) %>% 
  rename(feve = value) %>% 
  left_join(., algae_div, by = c("sample_ID" = "name")) %>% 
  mutate(redund = simpson - raoq) %>% 
  left_join(., comm_meta, by = "sample_ID") %>% 
  left_join(., npp, by = "season_ID") %>% 
  left_join(., ll_group_biomass, by = "season_ID") %>% 
  left_join(., sd_group_biomass, by = "season_ID")

fd_metrics_reduced <- algae_fd_reduced$nbsp %>% 
  enframe() %>% 
  rename(sample_ID = name,
         spp_rich = value) %>% 
  left_join(., enframe(algae_fd_reduced$FRic), by = c("sample_ID" = "name")) %>% 
  rename(fric = value) %>% 
  left_join(., enframe(algae_fd_reduced$RaoQ), by = c("sample_ID" = "name")) %>% 
  rename(raoq = value) %>% 
  left_join(., enframe(algae_fd_reduced$FDis), by = c("sample_ID" = "name")) %>% 
  rename(fdis = value) %>%  
  left_join(., enframe(algae_fd_reduced$FEve), by = c("sample_ID" = "name")) %>% 
  rename(feve = value) %>% 
  left_join(., algae_div, by = c("sample_ID" = "name")) %>% 
  mutate(redund = simpson - raoq) %>% 
  left_join(., comm_meta, by = "sample_ID") %>% 
  left_join(., npp, by = "season_ID") %>% 
  left_join(., ll_group_biomass, by = "season_ID") %>% 
  left_join(., sd_group_biomass, by = "season_ID")

# ⟞ ⟞ ii. models ----------------------------------------------------------


# ⟞ ⟞ ⟞ richness ----------------------------------------------------------

spp_rich_during <- glmmTMB(
  spp_rich ~ treatment*quality + (1|site) + (1|year),
  family = "poisson",
  data = fd_metrics_reduced %>% filter(exp_dates == "during")
)

# spp_rich_during <- glmmTMB(
#   spp_rich ~ time_since_end*treatment*quality + (1|year),
#   family = "poisson",
#   data = fd_metrics_reduced_reduced %>% filter(exp_dates == "during")
# )

plot(simulateResiduals(spp_rich_during))

spp_rich_after <- glmmTMB(
  spp_rich ~ treatment*quality + (1|site) + (1|year),
  family = poisson(link = "log"),
  data = fd_metrics_reduced %>% filter(exp_dates == "after")
)

# spp_rich_after <- glmmTMB(
#   spp_rich ~ time_since_end*treatment*quality + (1|year),
#   family = "poisson",
#   data = fd_metrics_reduced_reduced %>% filter(exp_dates == "after")
# )

plot(simulateResiduals(spp_rich_after))

fric_during <- glmmTMB(
  fric ~ treatment*quality + (1|site) + (1|year),
  family = beta_family(link = "logit"),
  data = fd_metrics_reduced %>% filter(exp_dates == "during")
)

# fric_during <- glmmTMB(
#   fric ~ time_since_end*treatment*quality + (1|year),
#   family = beta_family(link = "logit"),
#   data = fd_metrics_reduced %>% filter(exp_dates == "during")
# )

plot(simulateResiduals(fric_during))

fric_after <- glmmTMB(
  fric ~ treatment*quality + (1|site) + (1|year),
  # family = beta_family(link = "logit"),
  data = fd_metrics_reduced %>% filter(exp_dates == "after")
)

# fric_after <- glmmTMB(
#   fric ~ time_since_end*treatment*quality + (1|site) + (1|year),
#   family = beta_family(link = "logit"),
#   data = fd_metrics_reduced %>% filter(exp_dates == "after")
# )

plot(simulateResiduals(fric_after))
plotResiduals(fric_after, form = (fd_metrics_reduced %>% filter(exp_dates == "after"))$treatment)

fd_metrics_reduced %>% 
  filter(exp_dates == "during") %>% 
  ggplot(aes(x = redund)) +
  geom_histogram(bins = 12)

feve_during <- glmmTMB(
  feve ~ treatment*quality + (1|site) + (1|year),
  family = beta_family(link = "logit"),
  ziformula = ~ 1,
  data = fd_metrics_reduced %>% filter(exp_dates == "during")
)

# redund_during <- glmmTMB(
#   redund ~ time_since_end*treatment*quality + (1|year),
#   family = beta_family(link = "logit"),
#   ziformula = ~ 1,
#   data = fd_metrics_reduced %>% filter(exp_dates == "during")
# )

plot(simulateResiduals(feve_during))

feve_after <- glmmTMB(
  feve ~ treatment*quality + (1|site) + (1|year),
  family = beta_family(link = "logit"),
  ziformula = ~1,
  data = fd_metrics_reduced %>% filter(exp_dates == "after")
)

# redund_after <- glmmTMB(
#   redund ~ time_since_end*treatment*quality + (1|year),
#   family = beta_family(link = "logit"),
#   ziformula = ~1,
#   data = fd_metrics_reduced %>% filter(exp_dates == "after")
# )

plot(simulateResiduals(feve_after)) 

summary(spp_rich_during)
Anova(spp_rich_during, type = "II")
pairs(emmeans(object = spp_rich_during, 
              specs = c("treatment", "quality")),
      adjust = "bh")
ggpredict(spp_rich_during,
          terms = c("quality", "treatment")) %>% 
  plot(show_data = TRUE, jitter = TRUE)
# significant interaction of treatment and quality

summary(spp_rich_after)
Anova(spp_rich_after, type = "II")
pairs(emmeans(spp_rich_after, c("treatment", "quality")),
      adjust = "bh")
ggpredict(spp_rich_after,
          terms = c("quality", "treatment")) %>% 
  plot(show_data = TRUE, jitter = TRUE)
# effect of treatment but not of quality

summary(fric_during)
Anova(fric_during, type = "II")
pairs(emmeans(fric_during, c("treatment", "quality")),
      adjust = "bh")
ggpredict(fric_during,
          terms = c("quality", "treatment")) %>% 
  plot(show_data = TRUE, jitter = TRUE)
# significant effect of treatment and quality
# fric different in high quality habitats between treatment and removal, not in low quality habitats

summary(fric_after)
Anova(fric_after, type = "II")
pairs(emmeans(fric_after, c("treatment", "quality")),
      adjust = "bh")
ggpredict(fric_after,
          terms = c("quality", "treatment")) %>% 
  plot(show_data = TRUE, jitter = TRUE)
# significant effect of treatment only

summary(feve_during)
Anova(feve_during, type = "II")
pairs(emmeans(feve_during, c("treatment", "quality")),
      adjust = "bh")
ggpredict(feve_during,
          terms = c("quality", "treatment")) %>% 
  plot(show_data = TRUE, jitter = TRUE)
# significant effect of treatment and quality
# feve different in high quality habitats between treatment and removal, not in low quality habitats

summary(feve_after)
Anova(feve_after, type = "II")
pairs(emmeans(feve_after, c("treatment", "quality")),
      adjust = "bh")
ggpredict(feve_after,
          terms = c("quality", "treatment")) %>% 
  plot(show_data = TRUE, jitter = TRUE)
# interactions are not significant, going to type II
# significant interaction between treatment and quality, time and quality

spp_rich_pred_during <- ggpredict(spp_rich_during,
                              terms = c("treatment", "quality")) %>% 
  rename(treatment = x,
         quality = group)

ggemmeans(spp_rich_during, 
          terms = c("treatment", "quality")) %>% plot()

spp_rich_pred_after <- ggpredict(spp_rich_after,
                             terms = c("treatment", "quality")) %>% 
  rename(treatment = x,
         quality = group)

ggemmeans(spp_rich_after, 
          terms = c("treatment", "quality")) %>% plot()

fric_pred_during <- ggpredict(fric_during,
                                terms = c("treatment", "quality")) %>% 
  rename(treatment = x,
         quality = group)

ggemmeans(fric_during,
          terms = c("treatment", "quality")) %>% 
  plot() +
  coord_cartesian(ylim = c(0, 0.35))

fric_pred_after <- ggpredict(fric_after,
                               terms = c("treatment", "quality")) %>% 
  rename(treatment = x,
         quality = group)

ggemmeans(fric_after,
          terms = c("treatment", "quality")) %>% 
  plot()

feve_pred_during <- ggpredict(feve_during,
                                terms = c("treatment", "quality")) %>% 
  rename(treatment = x,
         quality = group)

ggemmeans(feve_during,
          terms = c("treatment", "quality")) %>% plot()

feve_pred_after <- ggpredict(feve_after,
                               terms = c("treatment", "quality")) %>% 
  rename(treatment = x,
         quality = group)

ggemmeans(feve_after,
          terms = c("treatment", "quality")) %>% plot()

spp_rich_during_plot <- ggplot() +
  geom_point(data = fd_metrics_reduced %>% filter(exp_dates == "during"),
             aes(x = treatment,
                 y = spp_rich,
                 color = treatment),
             alpha = 0.2,
             shape = 21,
             position = position_jitter(width = 0.1, 
                                        seed = 666)) +
  geom_pointrange(data = spp_rich_pred_during,
                  aes(x = treatment,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high,
                      color = treatment)) +
  model_preds_theme() +
  model_preds_aesthetics +
  labs(y = "Species richness",
       title = "During kelp removal") +
  facet_grid(cols = vars(quality)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 20))

spp_rich_after_plot <- ggplot() +
  geom_point(data = fd_metrics_reduced %>% filter(exp_dates == "after"),
             aes(x = treatment,
                 y = spp_rich,
                 color = treatment),
             alpha = 0.2,
             shape = 21,
             position = position_jitter(width = 0.1, 
                                        seed = 666)) +
  geom_pointrange(data = spp_rich_pred_after,
                  aes(x = treatment,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high,
                      color = treatment)) +
  model_preds_theme() +
  model_preds_aesthetics +
  labs(y = "Species richness",
       title = "Kelp recovery") +
  facet_grid(cols = vars(quality)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 20))

spp_rich_during_plot / spp_rich_after_plot

fric_during_plot <- ggplot() +
  geom_point(data = fd_metrics_reduced %>% filter(exp_dates == "during"),
             aes(x = treatment,
                 y = fric,
                 color = treatment),
             alpha = 0.2,
             shape = 21,
             position = position_jitter(width = 0.1, 
                                        seed = 666)) +
  geom_pointrange(data = fric_pred_during,
                  aes(x = treatment,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high,
                      color = treatment)) +
  model_preds_theme() +
  model_preds_aesthetics +
  labs(y = "Functional richness",
       title = "During kelp removal") +
  facet_grid(cols = vars(quality)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.42))

fric_after_plot <- ggplot() +
  geom_point(data = fd_metrics_reduced %>% filter(exp_dates == "after"),
             aes(x = treatment,
                 y = fric,
                 color = treatment),
             alpha = 0.2,
             shape = 21,
             position = position_jitter(width = 0.1, 
                                        seed = 666)) +
  geom_pointrange(data = fric_pred_after,
                  aes(x = treatment,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high,
                      color = treatment)) +
  model_preds_theme() +
  model_preds_aesthetics +
  labs(y = "Functional richness",
       title = "Kelp recovery") +
  facet_grid(cols = vars(quality)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.42))

fric_during_plot / fric_after_plot

feve_during_plot <- ggplot() +
  geom_point(data = fd_metrics_reduced %>% filter(exp_dates == "during"),
             aes(x = treatment,
                 y = feve,
                 color = treatment),
             alpha = 0.2,
             shape = 21,
             position = position_jitter(width = 0.1, 
                                        seed = 666)) +
  geom_pointrange(data = feve_pred_during,
                  aes(x = treatment,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high,
                      color = treatment)) +
  model_preds_theme() +
  model_preds_aesthetics +
  labs(y = "Functional evenness",
       title = "During kelp removal") +
  facet_grid(cols = vars(quality)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1))

feve_after_plot <- ggplot() +
  geom_point(data = fd_metrics_reduced %>% filter(exp_dates == "after"),
             aes(x = treatment,
                 y = feve,
                 color = treatment),
             alpha = 0.2,
             shape = 21,
             position = position_jitter(width = 0.1, 
                                        seed = 666)) +
  geom_pointrange(data = feve_pred_after,
                  aes(x = treatment,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high,
                      color = treatment)) +
  model_preds_theme() +
  model_preds_aesthetics +
  labs(y = "Functional evenness",
       title = "Kelp recovery") +
  facet_grid(cols = vars(quality)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1))

feve_during_plot / feve_after_plot



# spp_rich_time <- ggplot() +
#   coord_cartesian(ylim = c(-0.01, 18)) +
#   geom_point(data = fd_metrics_reduced,
#              aes(x = time_since_end,
#                  y = spp_rich,
#                  color = treatment),
#              alpha = 0.2, 
#              shape = 21) +
#   geom_ribbon(data = spp_rich_pred_during,
#               aes(x = time_since_end,
#                   ymin = conf.low,
#                   ymax = conf.high,
#                   group = treatment),
#               alpha = 0.05) +
#   geom_ribbon(data = spp_rich_pred_after,
#               aes(x = time_since_end,
#                   ymin = conf.low,
#                   ymax = conf.high,
#                   group = treatment),
#               alpha = 0.05) +
#   geom_line(data = spp_rich_pred_during,
#             aes(x = time_since_end,
#                 y = predicted,
#                 group = treatment,
#                 color = treatment,
#                 linetype = treatment),
#             linewidth = 1) +
#   geom_line(data = spp_rich_pred_after,
#             aes(x = time_since_end,
#                 y = predicted,
#                 group = treatment,
#                 color = treatment,
#                 linetype = treatment),
#             linewidth = 1) +
#   model_preds_aesthetics +
#   model_preds_theme() +
#   labs(title = "Species richness") +
#   facet_wrap(~ quality)
# spp_rich_time
# 
# redund_time <- ggplot() +
#   geom_point(data = fd_metrics_reduced,
#              aes(x = time_since_end,
#                  y = redund,
#                  color = treatment),
#              alpha = 0.2, 
#              shape = 21) +
#   geom_ribbon(data = redund_pred_during,
#               aes(x = time_since_end,
#                   ymin = conf.low,
#                   ymax = conf.high,
#                   group = treatment),
#               alpha = 0.1) +
#   geom_ribbon(data = redund_pred_after,
#               aes(x = time_since_end,
#                   ymin = conf.low,
#                   ymax = conf.high,
#                   group = treatment),
#               alpha = 0.05) +
#   geom_line(data = redund_pred_during,
#             aes(x = time_since_end,
#                 y = predicted,
#                 group = treatment,
#                 color = treatment,
#                 linetype = treatment),
#             linewidth = 1) +
#   geom_line(data = redund_pred_after,
#             aes(x = time_since_end,
#                 y = predicted,
#                 group = treatment,
#                 color = treatment,
#                 linetype = treatment),
#             linewidth = 1) +
#   model_preds_aesthetics + 
#   model_preds_theme() +
#   labs(title = "Functional redundancy") +
#   facet_wrap(~ quality)
# redund_time
# 
# fric_time <- ggplot() +
#   coord_cartesian(ylim = c(-0.01, 0.32)) +
#   geom_point(data = fd_metrics_reduced,
#              aes(x = time_since_end,
#                  y = fric,
#                  color = treatment),
#              alpha = 0.2, 
#              shape = 21) +
#   geom_ribbon(data = fric_pred_during,
#               aes(x = time_since_end,
#                   ymin = conf.low,
#                   ymax = conf.high,
#                   group = treatment),
#               alpha = 0.05) +
#   geom_ribbon(data = fric_pred_after,
#               aes(x = time_since_end,
#                   ymin = conf.low,
#                   ymax = conf.high,
#                   group = treatment),
#               alpha = 0.05) +
#   geom_line(data = fric_pred_during,
#             aes(x = time_since_end,
#                 y = predicted,
#                 group = treatment,
#                 color = treatment,
#                 linetype = treatment),
#             linewidth = 1) +
#   geom_line(data = fric_pred_after,
#             aes(x = time_since_end,
#                 y = predicted,
#                 group = treatment,
#                 color = treatment,
#                 linetype = treatment),
#             linewidth = 1) +
#   model_preds_aesthetics + 
#   model_preds_theme() +
#   labs(title = "Functional richness") +
#   facet_wrap(~ quality)
# fric_time
# 
# plots_together <- rich_time / fric_time / redund_time

ggsave(filename = here::here(
  "figures",
  "model-predictions",
  paste0("div-models_", today(), ".jpg")),
  plots_together,
  height = 14,
  width = 16,
  units = "cm",
  dpi = 300)

# ⟞ ⟞ ⟞  redundancy -------------------------------------------------------

rich_mod <- glmmTMB(fric ~ spp_rich,
                    data = fd_metrics,
                    family = beta_family(link = "logit"))

plot(simulateResiduals(rich_mod))

summary(rich_mod)

ggpredict(rich_mod,
          terms = "spp_rich") %>% 
  plot(show_data = TRUE) 

ggplot(data = fd_metrics,
       aes(x = spp_rich,
           y = fric)) +
  geom_point() +
  geom_smooth()

asy_model <- nls(fric ~ SSasymp(spp_rich, Asym, R0, lrc),
    data = fd_metrics)
summary(asy_model)
ggpredict(asy_model,
          terms = c("spp_rich")) %>% 
  plot(show_data = TRUE)

AICc(asy_model, rich_mod)

asy_preds <- propagate::predictNLS(asy_model,
                                   newdata = data.frame(spp_rich = seq(3, 19, by = 1)))

asy_preds_df <- enframe(asy_preds)

# ⟞ ⟞ ⟞ NPP ---------------------------------------------------------------

ggplot(data = fd_metrics %>% filter(exp_dates == "during"),
       aes(x = (total_npp))) +
  geom_histogram(bins = 12,
                 fill = "cornflowerblue",
                 color = "black") # log normal?

ggplot(data = fd_metrics %>% filter(exp_dates == "after"),
       aes(x = (total_npp))) +
  geom_histogram(bins = 12,
                 fill = "cornflowerblue",
                 color = "black") # log normal?

npp_df <- fd_metrics_reduced %>% 
  drop_na(fric) %>% 
  filter(sample_ID != "napl_control_2023-05-18")

npp_fdis <- glmmTMB(total_npp ~ fdis + (1|site) + (1|year) + (1|season),
                    family = lognormal(link = "log"),
                    na.action = na.fail,
                    data = npp_df)
plot(simulateResiduals(npp_fdis))
summary(npp_fdis)
r.squaredGLMM(npp_fdis)
ggpredict(npp_fdis,
          terms = c("fdis")) %>% 
  plot(show_data = TRUE) +
  coord_cartesian(ylim = c(0, 50))

npp_spp_rich <- glmmTMB(total_npp ~ spp_rich + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   na.action = na.fail,
                   data = npp_df)
plot(simulateResiduals(npp_spp_rich))
summary(npp_spp_rich)
r.squaredGLMM(npp_spp_rich)
ggpredict(npp_spp_rich,
          terms = c("spp_rich")) %>% 
  plot(show_data = TRUE) +
  coord_cartesian(ylim = c(0, 50))

ggplot(data = fd_metrics,
       aes(x = spp_rich,
           y = total_npp)) +
  geom_point() +
  geom_smooth()

npp_mod_fric <- glmmTMB(total_npp ~ fric + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   data = npp_df)
plot(simulateResiduals(npp_mod_fric))
summary(npp_mod_fric)
Anova(npp_mod_fric, type = "III")
ggpredict(npp_mod_fric,
          terms = c("fric")) %>% 
  plot(show_data = TRUE) +
  coord_cartesian(ylim = c(0, 50))

npp_mod_redund <- glmmTMB(total_npp ~ redund + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   data = npp_df)
plot(simulateResiduals(npp_mod_redund))
summary(npp_mod_redund)
Anova(npp_mod_redund, type = "III")
ggpredict(npp_mod_redund,
          terms = c("redund")) %>% 
  plot(show_data = TRUE) +
  coord_cartesian(ylim = c(0, 50))

npp_mod_thick_leathery <- glmmTMB(total_npp ~ thick_leathery + (1|site) + (1|year) + (1|season),
                                  family = lognormal(link = "log"),
                                  na.action = na.fail,
                                  data = npp_df)
plot(simulateResiduals(npp_mod_thick_leathery))
summary(npp_mod_thick_leathery)

npp_mod_leathery_macrophyte <- glmmTMB(total_npp ~ leathery_macrophyte + (1|site) + (1|year) + (1|season),
                                       family = lognormal(link = "log"),
                                       na.action = na.fail,
                                       data = npp_df)

plot(simulateResiduals(npp_mod_leathery_macrophyte))
summary(npp_mod_leathery_macrophyte)

ggplot(data = npp_df %>% 
         filter(sample_ID != "napl_control_2023-05-18"), 
       aes(x = thick_leathery,
           y = total_npp)) +
  geom_point()

ggpredict(npp_mod_thick_leathery,
          terms = "thick_leathery") %>% 
  plot(show_data = TRUE)

AICc(npp_mod_fric, npp_spp_rich, npp_mod_thick_leathery, npp_mod_leathery_macrophyte)

ggpredict(npp_spp_rich,
          terms = c("spp_rich")) %>% 
  plot(show_data = TRUE)

ggplot(data = fd_metrics,
       aes(x = fric,
           y = total_npp)) +
  geom_point(alpha = 0.5, 
             shape = 21,
             aes(size = thick_leathery))

ggpredict(npp_mod_fric,
          terms = c("fric")) %>% 
  plot(show_data = TRUE)

npp_mod <- glmmTMB(total_npp ~ spp_rich + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   data = fd_metrics %>% filter(exp_dates == "after"))
plot(simulateResiduals(npp_mod))
summary(npp_mod)
Anova(npp_mod, type = "III")
ggpredict(npp_mod,
          terms = c("spp_rich")) %>% 
  plot(show_data = TRUE)

npp_mod <- glmmTMB(total_npp ~ fric*quality*treatment + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   data = fd_metrics %>% filter(exp_dates == "after"))
plot(simulateResiduals(npp_mod))
summary(npp_mod)
Anova(npp_mod, type = "III")

ggpredict(npp_mod,
          terms = c("fric", "treatment", "quality")) %>% 
  plot(show_data = TRUE)

npp_mod <- glmmTMB(total_npp ~ redund*quality*treatment + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   data = fd_metrics %>% filter(exp_dates == "after"))
plot(simulateResiduals(npp_mod))
summary(npp_mod)
Anova(npp_mod, type = "III")

ggpredict(npp_mod,
          terms = c("redund", "treatment", "quality")) %>% 
  plot(show_data = TRUE)

npp_mod <- glmmTMB(total_npp ~ fric*quality*treatment + (1|site) + (1|year) + (1|season),
                   family = lognormal(link = "log"),
                   data = fd_metrics %>% filter(exp_dates == "after"))
plot(simulateResiduals(npp_mod))
plotResiduals(npp_mod, form = fd_metrics$redund)

summary(npp_mod)
Anova(npp_mod, type = "II")

ggpredict(npp_mod,
          terms = c("fric", "treatment", "quality")) %>% 
  plot(show_data = TRUE)


# ⟞ ⟞ ⟞ change in richness ------------------------------------------------

cont_redund <- fd_metrics %>% 
  filter(treatment == "continual") %>% 
  unite("yss", year, quarter, site, sep = "_", remove = FALSE)

delta_richness <- fd_metrics %>% 
  select(year, quarter, site, date, treatment, fric) %>% 
  pivot_wider(names_from = treatment,
              values_from = fric) %>% 
  replace_na(list(continual = 0,
                  control = 0)) %>% 
  mutate(delta_fric = continual - control) %>% 
  unite("yss", year, quarter, site, sep = "_") %>% 
  left_join(., cont_redund, by = "yss")

ggplot(data = delta_richness,
       aes(x = redund,
           y = delta_fric)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_smooth(method = "lm")

delta_mod <- glmmTMB(delta_fric ~ redund + (1|site) + (1|year),
                     data = delta_richness)
plot(simulateResiduals(delta_mod))
summary(delta_mod)
ggpredict(delta_mod,
          terms = c("redund")) %>% 
  plot(show_data = TRUE) +
  coord_cartesian(xlim = c(0, 0.8))


# ⟞ ⟞ iii. other plots ----------------------------------------------------

spric_plot <- ggplot(fd_metrics %>% filter(treatment == "continual"),
                     aes(x = quality,
                         y = spp_rich)) +
  geom_violin() +
  geom_point(alpha = 0.05, 
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Species richness")

fric_plot <- ggplot(fd_metrics %>% filter(treatment == "continual"),
                    aes(x = quality,
                        y = fric)) +
  geom_violin() +
  geom_point(alpha = 0.05, 
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates) +
  labs(title = "Functional richness (convex hull volume)")

redund_plot <- ggplot(fd_metrics %>% filter(treatment == "continual"),
                      aes(x = quality,
                          y = redund)) +
  geom_violin() +
  geom_point(alpha = 0.05, 
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates) +
  labs(title = "Functional redundancy")

cowplot::plot_grid(spric_plot, redund_plot, fric_plot, ncol = 3)


# ⟞ b. `mFD` --------------------------------------------------------------

# ⟞ ⟞ i. functional entities ----------------------------------------------

algae_fe <- sp.to.fe(
  sp_tr = trait_matrix,
  tr_cat = algae_traits_cat,
  fe_nm_type = "fe_rank",
  check_input = TRUE
)

# functional entity names
algae_fe$fe_nm

# how each species is distributed into a functional entity
algae_fe$sp_fe

# traits in each functional entity
algae_fe$fe_tr

# number of species per functional entity
algae_fe$fe_nb_sp

alphas <- alpha.fd.fe(
  asb_sp_occ = comm_mat_bin,
  sp_to_fe = algae_fe,
  check_input = TRUE,
  details_returned = TRUE
)

df <- alphas$asb_fdfe %>% 
  as_tibble(rownames = "sample_ID") %>% 
  left_join(., comm_meta, by = "sample_ID")

fred <- ggplot(data = df %>% filter(treatment == "continual"),
       aes(x = quality,
           y = fred)) +
  geom_point(alpha = 0.1,
             position = position_jitter(width = 0.05, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Functional redundancy (Mouillot et al. 2014)")

frich <- ggplot(data = df %>% filter(treatment == "continual"),
       aes(x = quality,
           y = nb_sp)) +
  geom_point(alpha = 0.1,
             position = position_jitter(width = 0.05, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Species richness")

cowplot::plot_grid(frich, fred, ncol = 2)

comm_mat_fe <- comm_mat_algae %>% 
  as_tibble(rownames = "sample_ID") %>% 
  pivot_longer(cols = 2:43,
               names_to = "scientific_name",
               values_to = "dry_gm2") %>% 
  left_join(., enframe(algae_fe$sp_fe), by = c("scientific_name" = "name")) %>% 
  rename(fe = value) %>% 
  group_by(sample_ID, fe) %>% 
  summarize(dry_gm2 = sum(dry_gm2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = fe,
              values_from = dry_gm2) %>% 
  filter(!(sample_ID %in% pull(few_spp_surveys, sample_ID))) %>% 
  filter(!(sample_ID == "carp_control_2016-02-22")) %>% 
  # filter(!(sample_ID %in% c(
  #   "aque_continual_2010-06-15",
  #   "aque_continual_2010-07-22",
  #   "aque_continual_2011-01-11",
  #   "aque_continual_2011-04-13"
  #   "aque_continual_2014-05-15",
  #   "aque_control_2013-02-13",
  #   "aque_control_2013-05-21",
  #   "aque_control_2020-08-14",
  #   "aque_control_2020-11-13",
  #   "aque_control_2021-02-25",
  #   "aque_control_2021-08-13",
  #   "aque_control_2022-05-23",
  #   "carp_continual_2010-07-19",
  #   "carp_continual_2010-10-14",
  #   "carp_continual_2011-07-18",
  #   "carp_continual_2012-05-10",
  #   "carp_continual_2012-08-15",
  #   "carp_continual_2013-05-20",
  #   "carp_continual_2014-05-20",
  #   "carp_continual_2014-11-17",
  #   "carp_continual_2015-08-18",
  #   "carp_continual_2015-11-23",
  #   "carp_continual_2016-02-22",
  #   "carp_continual_2018-02-20",
  #   "carp_continual_2019-02-25",
  #   "carp_continual_2019-05-13",
  #   "carp_continual_2019-08-12",
  #   "carp_continual_2019-11-30",
  #   "carp_control_2010-07-19",
  #   "carp_control_2010-10-14"
  # ))) %>%
  column_to_rownames("sample_ID") %>% 
  as.matrix()

# ⟞ ⟞ ii. diversity metrics -----------------------------------------------

# computing distances
algae_trait_distance <- mFD::funct.dist(
  sp_tr         = trait_matrix,
  tr_cat        = algae_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = FALSE)
# output is the same as trait_gower

# this basically does a PCoA
algae_quality <- mFD::quality.fspaces(
  sp_dist             = algae_trait_distance,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

# this gives the same axis values as trait_pcoa$vectors
sp_faxes_coord_algae <- algae_quality$"details_fspaces"$"sp_pc_coord"

alpha_fd_indices_algae <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_algae[ , c("PC1", "PC2")],
  asb_sp_w         = comm_mat_algae_reduced,
  ind_vect         = c('fide', 'fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 'fori', 'fspe'),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)



fun_div_df <- alpha_fd_indices_algae$functional_diversity_indices %>% 
  rownames_to_column("sample_ID") %>% 
  left_join(., comm_meta_algae, by = "sample_ID") %>% 
  filter(treatment == "continual")

sprich_mfd <- ggplot(data = fun_div_df,
       aes(x = quality,
           y = sp_richn)) +
  geom_point(position = position_jitter(width = 0.2, seed = 666),
             alpha = 0.1) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  facet_wrap(~exp_dates, ncol = 2) +
  labs(title = "Species richness")

fric_mfd <- ggplot(data = fun_div_df,
       aes(x = quality,
           y = fric)) +
  geom_point(position = position_jitter(width = 0.2, seed = 666),
             alpha = 0.1) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  facet_wrap(~exp_dates, ncol = 2) +
  labs(title = "Functional richness (convex hull volume)")

cowplot::plot_grid(sprich_mfd, fric_mfd, redundancy, nrow = 3)


# ⟞ c. `fundiversity` -----------------------------------------------------

library(fundiversity)

# same as results from FD
fric_fundiversity <- fd_fric(traits = trait_pcoa$vectors[, 1:2], 
        sp_com = comm_mat_algae) %>% 
  rename(sample_ID = site) %>% 
  left_join(., comm_meta, by = c("sample_ID"))

# different from FD!!!!
raoq_fundiversity <- fd_raoq(traits = trait_pcoa$vectors[, 1:2], 
                             sp_com = comm_mat_algae) %>% 
  rename(sample_ID = site) %>% 
  left_join(., comm_meta, by = c("sample_ID"))

ggplot(data = fric_fundiversity %>% filter(treatment == "continual"),
       aes(x = quality,
           y = FRic)) +
  geom_point(alpha = 0.1) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  facet_wrap(~exp_dates, ncol = 2)


# ⟞ ⟞ d. SYNCSA -----------------------------------------------------------

library(SYNCSA)

# different raoQ, different redundancy
syncsa_metrics <- rao.diversity(comm = comm_mat_algae,
                                traits = trait_matrix) 

syncsa_metrics_df <- bind_cols(
  enframe(syncsa_metrics$Simpson),
  enframe(syncsa_metrics$FunRao),
  enframe(syncsa_metrics$FunRedundancy)
) %>% 
  rename(sample_ID = name...1,
         simpson = value...2,
         raoq = value...4,
         redund = value...6) %>% 
  select(sample_ID, simpson, raoq, redund) %>% 
  left_join(., comm_meta, by = "sample_ID")

syncsa_metrics_df %>% 
  filter(treatment == "continual") %>% 
  ggplot(aes(x = quality,
             y = redund)) +
  geom_point(alpha = 0.1,
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  facet_wrap(~exp_dates)



# ⟞ d. functional space plots ---------------------------------------------

# ⟞ ⟞ i. global species pool ----------------------------------------------

# get the PCoA axes for axis
m <- trait_pcoa$vectors

# get the set of triangles that defines the hull
tr <- tri.mesh(m[,1], m[,2])

# get the points on the outside of the hull
ch <- convex.hull(tr) %>% 
  data.frame() %>% 
  rownames_to_column("scientific_name")
# points are species in PCoA axes (as in Teixido et al.)

# ⟞ ⟞ ii. survey trait space ----------------------------------------------

species_presence <- comm_mat_algae %>% 
  # turn the matrix into a data frame
  as_tibble(rownames = "sample_ID") %>% 
  # if the species biomass is greater than 0, replace that number with "yes"
  # if not, then "no"
  mutate(across(where(is.numeric), 
                ~ case_when(. > 0 ~ "yes", TRUE ~ "no"))) %>% 
  # make the data frame longer
  pivot_longer(cols = 2:43,
               names_to = "scientific_name",
               values_to = "presence") %>% 
  # join with PCoA scroes from above
  left_join(., spp_pcoa_scores, by = "scientific_name") %>% 
  # nest the data frame by sample_ID to get convex hulls for each survey
  nest(.by = sample_ID, data = everything()) %>% 
  mutate(filtered_comm = map2(
    data, sample_ID,
    ~ filter(.x, sample_ID == .y & presence == "yes")
  )) %>% 
  mutate(coords = map(
    filtered_comm,
    ~ select(.x, 
             scientific_name, Axis.1, Axis.2) %>% 
      column_to_rownames("scientific_name") %>% 
      as.matrix()
  )) %>% 
  mutate(hull_outside_species = map(
    filtered_comm,
    ~ nrow(.x)
  )) %>% 
  filter(hull_outside_species > 2) %>% 
  mutate(trimesh = map(
    coords,
    ~ tri.mesh(.x[, 1], .x[, 2])
  )) %>% 
  mutate(convex_hull = map(
    trimesh,
    ~ convex.hull(.x) %>% 
      data.frame() %>% 
      rownames_to_column("scientific_name")
  )) %>% 
  left_join(., comm_meta, by = "sample_ID")

# surveys that have too few species to plot a convex hull
too_few <- comm_mat_algae %>% 
  as_tibble(rownames = "sample_ID") %>% 
  mutate(across(where(is.numeric), 
                ~ case_when(. > 0 ~ "yes", TRUE ~ "no"))) %>% 
  pivot_longer(cols = 2:43,
               names_to = "scientific_name",
               values_to = "presence") %>% 
  left_join(., spp_pcoa_scores, by = "scientific_name") %>% 
  nest(.by = sample_ID, data = everything()) %>% 
  mutate(filtered_comm = map2(
    data, sample_ID,
    ~ filter(.x, sample_ID == .y & presence == "yes")
  )) %>% 
  mutate(coords = map(
    filtered_comm,
    ~ select(.x, 
             scientific_name, Axis.1, Axis.2) %>% 
      column_to_rownames("scientific_name") %>% 
      as.matrix()
  )) %>% 
  mutate(hull_outside_species = map(
    filtered_comm,
    ~ nrow(.x)
  )) %>% 
  select(sample_ID, hull_outside_species) %>% 
  unnest(hull_outside_species) %>% 
  filter(hull_outside_species < 3)

hulls <- species_presence %>% 
  select(sample_ID, convex_hull) %>% 
  unnest(cols = c(convex_hull)) %>% 
  left_join(., comm_meta, by = "sample_ID")

spp_by_survey <- species_presence %>% 
  select(filtered_comm) %>% 
  unnest(cols = c(filtered_comm)) %>% 
  select(!presence) %>% 
  left_join(., comm_meta, by = "sample_ID")

during_space <- ggplot() +
  theme(panel.background = element_rect(fill = "lightgrey")) +
  # global species pool
  geom_polygon(data = ch,
               aes(x = x,
                   y = y),
               fill = "#FFFFFF") +
  geom_point(data = m,
             aes(x = Axis.1,
                 y = Axis.2),
             color = "grey",
             shape = 21,
             alpha = 0.5) +
  geom_polygon(data = hulls %>% 
                 filter(exp_dates == "during" & quality == "high"),
               aes(x = x,
                   y = y,
                   group = sample_ID),
               alpha = 0.025,
               fill = high_col
  ) +
  geom_polygon(data = hulls %>% 
                 filter(exp_dates == "during" & quality == "medium"),
               aes(x = x,
                   y = y,
                   group = sample_ID),
               alpha = 0.025,
               fill = medium_col
  ) +
  geom_polygon(data = hulls %>% 
                 filter(exp_dates == "during" & quality == "low"),
               aes(x = x,
                   y = y,
                   group = sample_ID),
               alpha = 0.03,
               fill = low_col
  ) +
  facet_grid(cols = vars(quality), vars(treatment)) +
  labs(title = "Removal period",
       x = "PCoA 1",
       y = "PCoA 2") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))

after_space <- ggplot() +
  theme(panel.background = element_rect(fill = "lightgrey")) +
  # global species pool
  geom_polygon(data = ch,
               aes(x = x,
                   y = y),
               fill = "#FFFFFF") +
  geom_point(data = m,
             aes(x = Axis.1,
                 y = Axis.2),
             color = "grey",
             shape = 21,
             alpha = 0.5) +
  geom_polygon(data = hulls %>% 
                 filter(exp_dates == "after" & quality == "high"),
               aes(x = x,
                   y = y,
                   group = sample_ID),
               alpha = 0.025,
               fill = high_col
  ) +
  geom_polygon(data = hulls %>% 
                 filter(exp_dates == "after" & quality == "medium"),
               aes(x = x,
                   y = y,
                   group = sample_ID),
               alpha = 0.025,
               fill = medium_col
  ) +
  geom_polygon(data = hulls %>% 
                 filter(exp_dates == "after" & quality == "low"),
               aes(x = x,
                   y = y,
                   group = sample_ID),
               alpha = 0.03,
               fill = low_col
  ) +
  facet_grid(cols = vars(quality), vars(treatment)) +
  labs(title = "Recovery period",
       x = "PCoA 1",
       y = "PCoA 2") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))

space_together <- during_space + after_space

ggsave(here::here("figures",
                  "trait-space",
                  paste0("exp-dates-comparison_reduced_", today(), ".jpg")),
       space_together,
       width = 24,
       height = 10,
       units = "cm",
       dpi = 200)


# ⟞ e. log response ratio -------------------------------------------------

log_response_ratios <- fd_metrics %>% 
  nest(.by = c(exp_dates, quality),
       data = everything()) %>% 
  mutate(spp_rich_lrr = map(
    data,
    ~ means_ratio(x = spp_rich ~ treatment,
                  adjust = TRUE,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_spp_rich = log_Means_ratio_adjusted,
             CI_low_spp_rich = CI_low,
             CI_high_spp_rich = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(fric_lrr = map(
    data,
    ~ means_ratio(x = fric ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_fric = log_Means_ratio_adjusted,
             CI_low_fric = CI_low,
             CI_high_fric = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(redund_lrr = map(
    data,
    ~ means_ratio(x = redund ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_redund = log_Means_ratio_adjusted,
             CI_low_redund = CI_low,
             CI_high_redund = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(raoq_lrr = map(
    data,
    ~ means_ratio(x = raoq ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_raoq = log_Means_ratio_adjusted,
             CI_low_raoq = CI_low,
             CI_high_raoq = CI_high) %>% 
      select(!CI)
  )) %>%   
  mutate(fdis_lrr = map(
    data,
    ~ means_ratio(x = fdis ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_fdis = log_Means_ratio_adjusted,
             CI_low_fdis = CI_low,
             CI_high_fdis = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(feve_lrr = map(
    data,
    ~ means_ratio(x = feve ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_feve = log_Means_ratio_adjusted,
             CI_low_feve = CI_low,
             CI_high_feve = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(npp_lrr = map(
    data,
    ~ means_ratio(x = total_npp ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_npp = log_Means_ratio_adjusted,
             CI_low_npp = CI_low,
             CI_high_npp = CI_high) %>% 
      select(!CI)
  )) %>%  
  # select(!data) %>% 
  mutate(aggregate = pmap(
    list(exp_dates, quality, 
         spp_rich_lrr, fric_lrr, redund_lrr, 
         raoq_lrr, fdis_lrr, feve_lrr, npp_lrr),
    bind_cols
  )) 

lrr_for_plots <- log_response_ratios %>% 
  select(aggregate) %>% 
  unnest(cols = c(aggregate)) %>% 
  rename(exp_dates = ...1,
         quality = ...2) %>% 
  pivot_longer(cols = !c(exp_dates, quality),
               names_to = "stat",
               values_to = "value") %>% 
  mutate(name = case_when(
    str_detect(stat, "spp_rich") ~ "spp_rich",
    str_detect(stat, "fric") ~ "fric",
    str_detect(stat, "redund") ~ "redund",
    str_detect(stat, "raoq") ~ "raoq",
    str_detect(stat, "fdis") ~ "fdis",
    str_detect(stat, "feve") ~ "feve",
    str_detect(stat, "npp") ~ "npp"
  )) %>% 
  mutate(stat = case_when(
    str_detect(stat, "lrr") ~ "lrr",
    str_detect(stat, "CI_low") ~ "CI_low",
    str_detect(stat, "CI_high") ~ "CI_high"
  )) %>% 
  pivot_wider(names_from = "stat",
              values_from = "value") %>% 
  mutate(name = fct_relevel(name, "spp_rich", "fric", "redund", "fdis", "raoq", "feve", "npp"))

lrr_npp <- lrr_for_plots %>% 
  filter(name %in% c("npp")) %>% 
  ggplot(aes(x = quality,
             y = lrr)) + 
  geom_hline(yintercept = 0, 
             linetype = 2,
             color = "grey") +
  geom_pointrange(aes(ymin = CI_low,
                      ymax = CI_high)) +
  labs(x = "Experimental period",
       y = "log response ratio \n (log continual/control)") +
  facet_grid(cols = vars(exp_dates))
  

lrr_func <- lrr_for_plots %>% 
  filter(name %in% c("spp_rich", "fric", "feve", "fdis")) %>% 
ggplot(aes(x = quality,
           y = lrr)) + 
  geom_hline(yintercept = 0, 
             linetype = 2,
             color = "grey") +
  geom_pointrange(aes(ymin = CI_low,
                      ymax = CI_high)) +
  labs(x = "Experimental period",
       y = "log response ratio (log continual/control)") +
  facet_grid(rows = vars(name),
             cols = vars(exp_dates),
             labeller = labeller(
               name = c("spp_rich" = "Species richness", 
                        "fric" = "Functional richness",
                        "feve" = "Functional evenness",
                        "fdis" = "Functional dispersion")))
# higher value = more in continual removal plot than control
# after: higher spp richness in removal plot
# after: higher functional richness in removal plot
# no difference in redundancy
# during: lower Rao Q

# % change = 100 * (e^LRR - 1)
# from https://jepusto.com/files/ABAI-2019-Log-response-ratios.pdf

ggsave(here::here("figures",
                  "log-response-ratio",
                  paste0("npp_site-quality_exp-dates_", today(), ".jpg")),
       lrr_npp,
       height = 6, 
       width = 15,
       units = "cm",
       dpi = 300)


ggsave(here::here("figures",
                  "log-response-ratio",
                  paste0("func_site-quality_exp-dates_", today(), ".jpg")),
       lrr_func,
       height = 18, 
       width = 15,
       units = "cm",
       dpi = 300)


log_response_ratios_reduced <- fd_metrics_reduced %>% 
  nest(.by = c(exp_dates, quality),
       data = everything()) %>% 
  mutate(spp_rich_lrr = map(
    data,
    ~ means_ratio(x = spp_rich ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_spp_rich = log_Means_ratio_adjusted,
             CI_low_spp_rich = CI_low,
             CI_high_spp_rich = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(fric_lrr = map(
    data,
    ~ means_ratio(x = fric ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_fric = log_Means_ratio_adjusted,
             CI_low_fric = CI_low,
             CI_high_fric = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(redund_lrr = map(
    data,
    ~ means_ratio(x = redund ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_redund = log_Means_ratio_adjusted,
             CI_low_redund = CI_low,
             CI_high_redund = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(raoq_lrr = map(
    data,
    ~ means_ratio(x = raoq ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_raoq = log_Means_ratio_adjusted,
             CI_low_raoq = CI_low,
             CI_high_raoq = CI_high) %>% 
      select(!CI)
  )) %>%   
  mutate(fdis_lrr = map(
    data,
    ~ means_ratio(x = fdis ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_fdis = log_Means_ratio_adjusted,
             CI_low_fdis = CI_low,
             CI_high_fdis = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(feve_lrr = map(
    data,
    ~ means_ratio(x = feve ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_feve = log_Means_ratio_adjusted,
             CI_low_feve = CI_low,
             CI_high_feve = CI_high) %>% 
      select(!CI)
  )) %>% 
  mutate(npp_lrr = map(
    data,
    ~ means_ratio(x = total_npp ~ treatment,
                  data = .x,
                  log = TRUE,
                  ci = 0.95) %>% 
      rename(lrr_npp = log_Means_ratio_adjusted,
             CI_low_npp = CI_low,
             CI_high_npp = CI_high) %>% 
      select(!CI)
  )) %>%  
  # select(!data) %>% 
  mutate(aggregate = pmap(
    list(exp_dates, quality, 
         spp_rich_lrr, fric_lrr, redund_lrr, 
         raoq_lrr, fdis_lrr, feve_lrr, npp_lrr),
    bind_cols
  )) 

lrr_for_plots_reduced <- log_response_ratios_reduced %>% 
  select(aggregate) %>% 
  unnest(cols = c(aggregate)) %>% 
  rename(exp_dates = ...1,
         quality = ...2) %>% 
  pivot_longer(cols = !c(exp_dates, quality),
               names_to = "stat",
               values_to = "value") %>% 
  mutate(name = case_when(
    str_detect(stat, "spp_rich") ~ "spp_rich",
    str_detect(stat, "fric") ~ "fric",
    str_detect(stat, "redund") ~ "redund",
    str_detect(stat, "raoq") ~ "raoq",
    str_detect(stat, "fdis") ~ "fdis",
    str_detect(stat, "feve") ~ "feve",
    str_detect(stat, "npp") ~ "npp"
  )) %>% 
  mutate(stat = case_when(
    str_detect(stat, "lrr") ~ "lrr",
    str_detect(stat, "CI_low") ~ "CI_low",
    str_detect(stat, "CI_high") ~ "CI_high"
  )) %>% 
  pivot_wider(names_from = "stat",
              values_from = "value") %>% 
  mutate(name = fct_relevel(name, "spp_rich", "fric", "redund", "fdis", "raoq", "feve", "npp"))

lrr_npp_reduced <- lrr_for_plots_reduced %>% 
  filter(name %in% c("npp")) %>% 
  ggplot(aes(x = quality,
             y = lrr)) + 
  geom_hline(yintercept = 0, 
             linetype = 2,
             color = "grey") +
  geom_pointrange(aes(ymin = CI_low,
                      ymax = CI_high)) +
  labs(x = "Experimental period",
       y = "log response ratio \n (log continual/control)") +
  facet_grid(cols = vars(exp_dates))


lrr_func_reduced <- lrr_for_plots_reduced %>% 
  filter(name %in% c("spp_rich", "fric", "feve")) %>% 
  ggplot(aes(x = quality,
             y = lrr)) + 
  geom_hline(yintercept = 0, 
             linetype = 2,
             color = "grey") +
  geom_pointrange(aes(ymin = CI_low,
                      ymax = CI_high)) +
  labs(x = "Experimental period",
       y = "log response ratio (log continual/control)") +
  facet_grid(rows = vars(name),
             cols = vars(exp_dates),
             labeller = labeller(
               name = c("spp_rich" = "Species richness", 
                        "fric" = "Functional richness",
                        "feve" = "Functional evenness")))
# higher value = more in continual removal plot than control
# after: higher spp richness in removal plot
# after: higher functional richness in removal plot
# no difference in redundancy
# during: lower Rao Q

# % change = 100 * (e^LRR - 1)
# from https://jepusto.com/files/ABAI-2019-Log-response-ratios.pdf

ggsave(here::here("figures",
                  "log-response-ratio",
                  paste0("npp_site-quality_exp-dates_", today(), ".jpg")),
       lrr_npp,
       height = 6, 
       width = 15,
       units = "cm",
       dpi = 300)


ggsave(here::here("figures",
                  "log-response-ratio",
                  paste0("func_site-quality_exp-dates_reduced_", today(), ".jpg")),
       lrr_func_reduced,
       height = 18, 
       width = 15,
       units = "cm",
       dpi = 300)









