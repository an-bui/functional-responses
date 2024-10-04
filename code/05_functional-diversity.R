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
    "calcification" ~ "N"
  ))

# reduced trait matrix
trait_matrix_reduced <- trait_matrix %>% 
  select(stipe, branching, blade_category, attachment, growth, calcification)

# reduced trait categories
algae_traits_cat_reduced <- algae_traits_cat %>% 
  filter(trait_name %in% colnames(trait_matrix_reduced))

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
trait_gower <- gowdis(trait_matrix)

# doing PCoA to get dimensions
trait_pcoa <- ape::pcoa(D = trait_gower)

# extracting proportion of inertia explained by each axis
proportion_inertia <- trait_pcoa$values %>% 
  rownames_to_column("axis") %>% 
  mutate(axis = as_factor(axis))

# plotting for visualization
prop_inertia_plot <- ggplot(data = proportion_inertia,
                            aes(x = axis,
                                y = Rel_corr_eig)) +
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

# messages:
# Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept. 
# FEVe: Could not be calculated for communities with <3 functionally singular species. 
# FDis: Equals 0 in communities with only one functionally singular species. 
# FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
# FRic: Dimensionality reduction was required. The last 16 PCoA axes (out of 18 in total) were removed. 
# FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) = 0.2849965 
# FDiv: Could not be calculated for communities with <3 functionally singular species. 

spp_ric <- algae_fd$nbsp %>% 
  enframe() %>% 
  left_join(., comm_meta_algae, by = c("name" = "sample_ID"))

fric <- algae_fd$FRic %>% 
  enframe() %>% 
  left_join(., comm_meta_algae, by = c("name" = "sample_ID"))

algae_div <- vegan::diversity(x = comm_mat_algae,
                              index = "simpson") %>% 
  enframe() %>% 
  rename(simpson = value)

raoq <- algae_fd$RaoQ %>% 
  enframe() %>% 
  rename(raoq = value) %>% 
  left_join(., comm_meta_algae, by = c("name" = "sample_ID")) %>% 
  left_join(., algae_div, by = "name") %>% 
  mutate(redund = simpson - raoq)

# ⟞ ⟞ ii. models ----------------------------------------------------------

rich_during <- glmmTMB(value ~ time_since_end*treatment + (1|site) + (1|year),
                       data = spp_ric %>% filter(exp_dates == "during"))

rich_after <- glmmTMB(value ~ time_since_end*treatment + (1|site) + (1|year),
                      data = spp_ric %>% filter(exp_dates == "after"))

redund_during <- glmmTMB(redund ~ time_since_end*treatment + (1|site) + (1|year),
                         data = raoq %>% filter(exp_dates == "during"))

redund_after <- glmmTMB(redund ~ time_since_end*treatment + (1|site) + (1|year),
                        data = raoq %>% filter(exp_dates == "after"))

fric_during <- glmmTMB(value ~ time_since_end*treatment + (1|site) + (1|year),
                       data = fric %>% filter(exp_dates == "during"))

fric_after <- glmmTMB(value ~ time_since_end*treatment + (1|site) + (1|year),
                      data = fric %>% filter(exp_dates == "after"))

plot(simulateResiduals(rich_during)) 
plot(simulateResiduals(rich_after))
plot(simulateResiduals(redund_during)) # bad
plot(simulateResiduals(redund_after)) 
plot(simulateResiduals(fric_during))
plot(simulateResiduals(fric_after))

summary(rich_during) # significant interaction between time and treatment
summary(rich_after) # significant interaction between time and treatment
summary(fric_during) # significant effects of time and treatment, but not interaction
summary(fric_after) # significant interaction between time and treatment
summary(redund_during) # significant effect of treatment only
summary(redund_after) # significant effect of time only

rich_pred_during <- ggpredict(rich_during,
                              terms = c("time_since_end", "treatment")) %>% 
  rename(time_since_end = x,
         treatment = group)

rich_pred_after <- ggpredict(rich_after,
                              terms = c("time_since_end", "treatment")) %>% 
  rename(time_since_end = x,
         treatment = group)

redund_pred_during <- ggpredict(redund_during,
                                terms = c("time_since_end", "treatment")) %>% 
  rename(time_since_end = x,
         treatment = group)

redund_pred_after <- ggpredict(redund_after,
                               terms = c("time_since_end", "treatment")) %>% 
  rename(time_since_end = x,
         treatment = group)

fric_pred_during <- ggpredict(fric_during,
                              terms = c("time_since_end", "treatment")) %>% 
  rename(time_since_end = x,
         treatment = group)

fric_pred_after <- ggpredict(fric_after,
                             terms = c("time_since_end", "treatment")) %>% 
  rename(time_since_end = x,
         treatment = group)

rich_pred_plot <- ggplot(mapping = aes(group = treatment,
                                       linetype = treatment,
                                       x = time_since_end)) +
  geom_point(data = spp_ric,
             aes(y = value,
                 color = treatment),
             alpha = 0.3,
             shape = 21) +
  geom_ribbon(data = rich_pred_during,
              aes(y = predicted,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  geom_ribbon(data = rich_pred_after,
              aes(y = predicted,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  geom_line(data = rich_pred_during,
            aes(y = predicted,
                color = treatment),
            linewidth = 1) +
  geom_line(data = rich_pred_after,
            aes(y = predicted,
                color = treatment),
            linewidth = 1) +
  model_preds_aesthetics +
  model_preds_theme() +
  labs(title = "Species richness")
rich_pred_plot

fric_pred_plot <- ggplot(mapping = aes(group = treatment,
                                       linetype = treatment,
                                       x = time_since_end)) +
  geom_point(data = fric,
             aes(y = value,
                 color = treatment),
             alpha = 0.3,
             shape = 21) +
  geom_ribbon(data = fric_pred_during,
              aes(y = predicted,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  geom_ribbon(data = fric_pred_after,
              aes(y = predicted,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  geom_line(data = fric_pred_during,
            aes(y = predicted,
                color = treatment),
            linewidth = 1) +
  geom_line(data = fric_pred_after,
            aes(y = predicted,
                color = treatment),
            linewidth = 1) +
  model_preds_aesthetics +
  model_preds_theme() +
  labs(title = "Functional richness")
fric_pred_plot

redund_pred_plot <- ggplot(mapping = aes(group = treatment,
                                       linetype = treatment,
                                       x = time_since_end)) +
  geom_point(data = raoq,
             aes(y = redund,
                 color = treatment),
             alpha = 0.3,
             shape = 21) +
  geom_ribbon(data = redund_pred_during,
              aes(y = predicted,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  geom_ribbon(data = redund_pred_after,
              aes(y = predicted,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  geom_line(data = redund_pred_during,
            aes(y = predicted,
                color = treatment),
            linewidth = 1) +
  geom_line(data = redund_pred_after,
            aes(y = predicted,
                color = treatment),
            linewidth = 1) +
  model_preds_aesthetics +
  model_preds_theme() +
  labs(title = "Functional redundancy")
redund_pred_plot

plots_together <- rich_pred_plot / fric_pred_plot / redund_pred_plot

ggsave(here::here("figures",
                  "model-predictions",
                  paste0("div-models_no-site-quality_",
                         today(),
                         ".jpg")),
       plots_together,
       height = 16,
       width = 12,
       units = "cm",
       dpi = 200)

# includes quality

rich_during <- glmmTMB(value ~ time_since_end*treatment*quality + (1|site) + (1|year),
               data = spp_ric %>% filter(exp_dates == "during"))

rich_after <- glmmTMB(value ~ time_since_end*treatment*quality + (1|site) + (1|year),
                       data = spp_ric %>% filter(exp_dates == "after"))

redund_during <- glmmTMB(redund ~ time_since_end*treatment*quality + (1|site) + (1|year),
                         data = raoq %>% filter(exp_dates == "during"))

redund_after <- glmmTMB(redund ~ time_since_end*treatment*quality + (1|site) + (1|year),
                         data = raoq %>% filter(exp_dates == "after"))

fric_during <- glmmTMB(value ~ time_since_end*treatment*quality + (1|site) + (1|year),
                       data = fric %>% filter(exp_dates == "during"))

fric_after <- glmmTMB(value ~ time_since_end*treatment*quality + (1|site) + (1|year),
                       data = fric %>% filter(exp_dates == "after"))

plot(simulateResiduals(rich_during))
plot(simulateResiduals(rich_after))
plot(simulateResiduals(redund_during)) # bad
plot(simulateResiduals(redund_after)) # ok but model didn't converge
plot(simulateResiduals(fric_during))
plot(simulateResiduals(fric_after))

summary(rich_during) # effect of treatment and habitat quality interaction, not of time
summary(rich_after) # effect of time, quality, and treatment interactions
summary(fric_during) # effect of habitat quality only
summary(fric_after) # effect of interaction between time and habitat quality
summary(redund_during) # effect of time and habitat quality interaction
summary(redund_after) # nothing!!!!

rich_pred_during <- ggpredict(rich_during,
                              terms = c("time_since_end", "treatment", "quality")) %>% 
  rename(time_since_end = x,
         treatment = group,
         quality = facet)

ggemmeans(rich_during,
          terms = c("treatment", "quality")) %>% plot()

rich_pred_after <- ggpredict(rich_after,
                             terms = c("time_since_end", "treatment", "quality")) %>% 
  rename(time_since_end = x,
         treatment = group,
         quality = facet)

ggemmeans(rich_after,
          terms = c("treatment", "quality")) %>% plot()

redund_pred_during <- ggpredict(redund_during,
                                terms = c("time_since_end", "treatment", "quality")) %>% 
  rename(time_since_end = x,
         treatment = group,
         quality = facet)

ggemmeans(redund_during,
          terms = c("treatment", "quality")) %>% plot()

redund_pred_after <- ggpredict(redund_after,
                               terms = c("time_since_end", "treatment", "quality")) %>% 
  rename(time_since_end = x,
         treatment = group,
         quality = facet)

ggemmeans(redund_after,
          terms = c("treatment", "quality")) %>% plot()

fric_pred_during <- ggpredict(fric_during,
                                terms = c("time_since_end", "treatment", "quality")) %>% 
  rename(time_since_end = x,
         treatment = group,
         quality = facet)

ggemmeans(fric_during,
          terms = c("treatment", "quality")) %>% plot()

fric_pred_after <- ggpredict(fric_after,
                               terms = c("time_since_end", "treatment", "quality")) %>% 
  rename(time_since_end = x,
         treatment = group,
         quality = facet)

ggemmeans(fric_after,
          terms = c("time_since_end", "quality")) %>% plot()

rich_time <- ggplot() +
  coord_cartesian(ylim = c(-0.01, 20)) +
  geom_point(data = spp_ric,
             aes(x = time_since_end,
                 y = value,
                 color = treatment),
             alpha = 0.2, 
             shape = 21) +
  geom_ribbon(data = rich_pred_during,
              aes(x = time_since_end,
                  ymin = conf.low,
                  ymax = conf.high,
                  group = treatment),
              alpha = 0.05) +
  geom_ribbon(data = rich_pred_after,
              aes(x = time_since_end,
                  ymin = conf.low,
                  ymax = conf.high,
                  group = treatment),
              alpha = 0.05) +
  geom_line(data = rich_pred_during,
            aes(x = time_since_end,
                y = predicted,
                group = treatment,
                color = treatment,
                linetype = treatment),
            linewidth = 1) +
  geom_line(data = rich_pred_after,
            aes(x = time_since_end,
                y = predicted,
                group = treatment,
                color = treatment,
                linetype = treatment),
            linewidth = 1) +
  scale_color_manual(values = c(control = control_col, continual = continual_col),
                     labels = c("Removal", "Reference")) +
  scale_linetype_manual(values = c(control = "22", continual = "solid"),
                        labels = c("Removal", "Reference")) +
  labs(title = "Species richness") +
  facet_wrap(~ quality)

redund_time <- ggplot() +
  geom_point(data = raoq,
             aes(x = time_since_end,
                 y = redund,
                 color = treatment),
             alpha = 0.2, 
             shape = 21) +
  geom_ribbon(data = redund_pred_during,
              aes(x = time_since_end,
                  ymin = conf.low,
                  ymax = conf.high,
                  group = treatment),
              alpha = 0.1) +
  # geom_ribbon(data = redund_pred_after,
  #             aes(x = time_since_end,
  #                 ymin = conf.low,
  #                 ymax = conf.high,
  #                 group = treatment),
  #             alpha = 0.05) +
  geom_line(data = redund_pred_during,
            aes(x = time_since_end,
                y = predicted,
                group = treatment,
                color = treatment,
                linetype = treatment),
            linewidth = 1) +
  geom_line(data = redund_pred_after,
            aes(x = time_since_end,
                y = predicted,
                group = treatment,
                color = treatment,
                linetype = treatment),
            linewidth = 1) +
  scale_color_manual(values = c("control" = control_col, "continual" = continual_col),
                     labels = c("Removal", "Reference")) +
  scale_linetype_manual(values = c("control" = "22", "continual" = "solid"),
                        labels = c("Removal", "Reference")) +
  labs(title = "Functional redundancy") +
  facet_wrap(~ quality)

fric_time <- ggplot() +
  coord_cartesian(ylim = c(-0.01, 0.32)) +
  geom_point(data = fric,
             aes(x = time_since_end,
                 y = value,
                 color = treatment),
             alpha = 0.2, 
             shape = 21) +
  geom_ribbon(data = fric_pred_during,
              aes(x = time_since_end,
                  ymin = conf.low,
                  ymax = conf.high,
                  group = treatment),
              alpha = 0.05) +
  geom_ribbon(data = fric_pred_after,
              aes(x = time_since_end,
                  ymin = conf.low,
                  ymax = conf.high,
                  group = treatment),
              alpha = 0.05) +
  geom_line(data = fric_pred_during,
            aes(x = time_since_end,
                y = predicted,
                group = treatment,
                color = treatment,
                linetype = treatment),
            linewidth = 1) +
  geom_line(data = fric_pred_after,
            aes(x = time_since_end,
                y = predicted,
                group = treatment,
                color = treatment,
                linetype = treatment),
            linewidth = 1) +
  scale_color_manual(values = c(control = control_col, continual = continual_col),
                     labels = c("Removal", "Reference")) +
  scale_linetype_manual(values = c(control = "22", continual = "solid"),
                        labels = c("Removal", "Reference")) +
  labs(title = "Functional richness") +
  facet_wrap(~ quality)
fric_time

plots_together <- rich_time / fric_time / redund_time

ggsave(filename = here::here(
  "figures",
  "model-predictions",
  paste0("div-models_", today(), ".jpg")),
  plots_together,
  height = 14,
  width = 14,
  units = "cm",
  dpi = 300)


ggplot(raoq %>% filter(exp_dates == "during"),
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
  facet_wrap(~treatment, ncol = 2) +
  labs(title = "Functional redundancy")

redund_plot <- ggplot(raoq %>% filter(treatment == "control"),
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
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Functional redundancy")

spric_plot <- ggplot(spp_ric %>% filter(treatment == "continual"),
                     aes(x = quality,
                         y = value)) +
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

fric_plot <- ggplot(fric %>% filter(treatment == "continual"),
                    aes(x = quality,
                        y = value)) +
  geom_violin() +
  geom_point(alpha = 0.05, 
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Functional richness (convex hull volume)")

cowplot::plot_grid(spric_plot, redund_plot, fric_plot, ncol = 3)

ggplot(raoq %>% filter(exp_dates == "during"),
       aes(x = time_since_end,
           y = value,
           color = quality)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(raoq %>% filter(exp_dates == "after"),
       aes(x = time_since_end,
           y = value,
           color = quality)) +
  geom_point() +
  geom_smooth(method = "lm")


# ⟞ b. `mFD` --------------------------------------------------------------

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


# computing distances
algae_trait_distance <- mFD::funct.dist(
  sp_tr         = trait_matrix,
  tr_cat        = algae_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = FALSE)

algae_quality <- mFD::quality.fspaces(
  sp_dist             = algae_trait_distance,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

# this gives the same axis values as trait_pcoa$vectors
sp_faxes_coord_algae <- algae_quality$"details_fspaces"$"sp_pc_coord"

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



alpha_fd_indices_algae <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_algae[ , c("PC1", "PC2")],
  asb_sp_w         = comm_mat_algae_reduced,
  ind_vect         = c("fdis", "feve", "fric", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_algae,
  plot_asb_nm = c("aque_continual_2011-01-11", "aque_control_2011-01-11"),
  ind_nm = c("fric")
)

plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_algae,
  plot_asb_nm              = c("aque_continual_2011-01-11", "aque_control_2011-01-11 "),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

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



library(adiv)

trait_gower <- gowdis(trait_matrix)

trait_div_adiv <- uniqueness(comm = comm_mat_algae_reduced,
           dis = trait_gower)

redund <- trait_div_adiv$red %>% 
  rownames_to_column("sample_ID") %>% 
  left_join(., comm_meta, by = "sample_ID") %>% 
  filter(treatment == "continual")
            
            
rich <- ggplot(redund,
       aes(x = quality,
           y = N)) +
  geom_point(alpha = 0.1,
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  labs(title = "Species richess") +
  facet_wrap(~exp_dates, ncol = 1)
  
redundancy <- ggplot(redund,
       aes(x = quality,
           y = Rstar)) +
  geom_point(alpha = 0.1,
             position = position_jitter(width = 0.1, seed = 666)) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  labs(title = "Functional redundancy") +
  facet_wrap(~exp_dates, ncol = 2)
  
  
cowplot::plot_grid(rich, redundancy, ncol = 2)


# ⟞ c. `fundiversity` -----------------------------------------------------

library(fundiversity)

fric_fundiversity <- fd_fric(traits = trait_pcoa$vectors[, 1:2], 
        sp_com = comm_mat_algae_reduced) %>% 
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



# ⟞ d. `betapart` ---------------------------------------------------------

library(betapart)

functional.beta.multi(x = comm_mat_algae_reduced_bin,
                      traits = trait_axes[, 1:2],
                      )










