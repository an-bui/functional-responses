##########################################################################-
# Community analyses
# last modified: 2024-04-24

# This is a script to use the clusters from `03_trait-clustering.R` to 
# attache the clusters to biomass data from LTE surveys.
##########################################################################-

##########################################################################-
# 1. source ---------------------------------------------------------------
##########################################################################-

# only need to do this once per session
source(here::here("code", "03_trait-clustering.R"))

##########################################################################-
# 2. site species ordination ----------------------------------------------
##########################################################################-

# Note: this is already in the LTE paper. Just recreating it here for reference

spp_nmds <- metaMDS(comm_mat_algae, distance = "altGower")

stressplot(spp_nmds)

adonis2(comm_mat_algae ~ time_since_end*treatment, 
        data = comm_meta_algae %>% mutate(time_since_end = as_factor(time_since_end)), 
        method = "altGower")

spp_gower_dist <- vegdist(comm_mat_algae, method = "altGower")

spp_disper_treatment <- betadisper(spp_gower_dist, comm_meta_algae$treatment)

anova(spp_disper_treatment)

spp_nmds_scores <- scores(spp_nmds, choices = c(1, 2), display = "sites") %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("sample_ID") %>% 
  left_join(comm_meta_algae, by = "sample_ID")

ggplot(data = spp_nmds_scores,
       aes(x = NMDS1,
           y = NMDS2,
           color = treatment,
           shape = site_full)) +
  geom_point(size = 2)

##########################################################################-
# 3. calculating per cluster biomass for each survey ----------------------
##########################################################################-

# match community data with clusters
cluster_match <- comm_df %>% 
  # only include algae
  filter(new_group == "algae") %>% 
  # join with clusters
  left_join(., pam_clusters_7, by = "sp_code")

# find the taxa that are not included in the clustering
cluster_missing <- cluster_match %>% 
  filter(is.na(scientific_name))

# function to calculate cluster biomass
# takes argument for the actual column name
calc_group_biomass <- function(grouping) {
  cluster_match %>% 
    filter(!is.na(cluster)) %>% 
    # group by sample_ID
    group_by(sample_ID, {{ grouping }}) %>% 
    # add up biomass within a cluster
    summarize(group_total = sum(dry_gm2, na.rm = TRUE)) %>% 
    # calculate total biomass and proportion of each cluster
    mutate(total = sum(group_total),
           prop = group_total/total,
           prop = case_when(
             prop == "NaN" ~ 0,
             TRUE ~ prop
           )) %>% 
    ungroup() %>% 
    separate_wider_delim(cols = sample_ID,
                         delim = "_",
                         names = c("site", "treatment", "date"),
                         cols_remove = FALSE) %>% 
    mutate(date = as_date(date)) %>% 
    select(!c(site, treatment, date)) %>% 
    left_join(., comm_meta, by = "sample_ID") %>% 
    filter(site == "aque" & date > "2010-04-26" |
             site == "napl" & date > "2010-04-27" |
             site == "mohk" & date > "2010-05-05" |
             site == "carp" & date > "2010-04-23") 
}

# calculate per group biomass in each survey
cluster_biomass <- calc_group_biomass(cluster)
gf_biomass <- calc_group_biomass(sd_growth_form)
ff_biomass <- calc_group_biomass(ll_func_form)

# function to plot timeseries
group_timeseries <- function(df, grouping, y) {
  df %>% 
    ggplot(aes(x = time_since_end,
               y = {{ y }},
               color = {{ grouping }},
               alpha = treatment)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 0,
               lty = 2) +
    scale_alpha_manual(values = c("control" = 0.4, "continual" = 1)) +
    facet_wrap(~site, ncol = 1, scales = "free_y") +
    theme(panel.grid = element_blank())
}

cluster_prop_timeseries <- group_timeseries(
  df = cluster_biomass, 
  grouping = cluster,
  y = prop
) +
  scale_color_manual(values = cluster_cols,
                     name = "Cluster")

cluster_prop_timeseries

ff_prop_timeseries <- group_timeseries(
  df = ff_biomass,
  grouping = ll_func_form,
  y = prop
)

ff_prop_timeseries

gf_prop_timeseries <- group_timeseries(
  df = gf_biomass,
  grouping = sd_growth_form,
  y = prop
)

gf_prop_timeseries

# ggsave(here("figures", "community-timeseries", paste0("cluster-prop-timeseries_", today(), ".jpg")),
#        cluster_prop_timeseries,
#        height = 12,
#        width = 8)

cluster_total_timeseries <- group_timeseries(
  df = cluster_biomass %>% 
    filter(sample_ID != "napl_control_2023-05-18"),
  grouping = cluster,
  y = group_total
) +
  scale_color_manual(values = cluster_cols,
                     name = "Cluster")

cluster_total_timeseries
  
ff_total_timeseries <- group_timeseries(
  df = ff_biomass %>% 
    filter(sample_ID != "napl_control_2023-05-18"),
  grouping = ll_func_form,
  y = group_total
) 

ff_total_timeseries

gf_total_timeseries <- group_timeseries(
  df = gf_biomass %>% 
    filter(sample_ID != "napl_control_2023-05-18"),
  grouping = sd_growth_form,
  y = group_total
) 

gf_total_timeseries

# ggsave(here("figures", "community-timeseries", paste0("cluster7-pam-total-timeseries_", today(), ".jpg")),
#        cluster_total_timeseries,
#        height = 12,
#        width = 8)

##########################################################################-
# 4. total biomass model --------------------------------------------------
##########################################################################-

## a. clusters ------------------------------------------------------------

# model 1: includes all predictors and random effects of site and year
cluster_total_biomass_during_model_1 <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment + (1|site) + (1|year),
  data = cluster_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1
  )

# model 2: includes all predictors and random effect of site
cluster_total_biomass_during_model_2 <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment + (1|site),
  data = cluster_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1
  )

# model 3: includes all predictors, no random effects
cluster_total_biomass_during_model_3 <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment,
  data = cluster_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 4: includes time since end and cluster as only predictors, no random effects
cluster_total_biomass_during_model_4 <- glmmTMB(
  group_total ~ time_since_end*cluster,
  data = cluster_biomass %>% filter(exp_dates == "after" & treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# best model is 1 with all predictors and random effect of site and year
MuMIn::AICc(cluster_total_biomass_during_model_1,
            cluster_total_biomass_during_model_2,
            cluster_total_biomass_during_model_3) %>% 
  arrange(AICc)

# looking at model residuals
simulateResiduals(cluster_total_biomass_during_model_1, plot = TRUE)

# model summary
summary(cluster_total_biomass_during_model_1)

cluster_total_model_summary <- tbl_regression(
  cluster_total_biomass_during_model_1,
  intercept = TRUE, 
  add_estimate_to_reference_rows = TRUE
  ) %>% 
  gtsummary::as_flex_table()

cluster_total_model_summary

# save table as doc
# save_as_docx(cluster_total_model_summary,
#              path = here("tables", "model-summaries", paste0("cluster-total-model-summary_", today(), ".docx")))

cluster_total_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment + (1|site),
  data = cluster_biomass %>% filter(exp_dates == "after"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

simulateResiduals(cluster_total_biomass_after_model, plot = TRUE)

testOutliers(cluster_total_biomass_after_model)

summary(cluster_total_biomass_after_model)

cluster_during_model_preds <- ggpredict(cluster_total_biomass_during_model_1, terms = c("time_since_end", "cluster", "treatment")) %>% 
  rename(time_since_end = x,
         cluster = group,
         treatment = facet)

cluster_after_model_preds <- ggpredict(cluster_total_biomass_after_model, terms = c("time_since_end", "cluster", "treatment")) %>% 
  rename(time_since_end = x,
         cluster = group,
         treatment = facet)

cluster_biomass_continual_plot <- ggplot() +
  # during removal 
  # 95% CI ribbons
  geom_ribbon(data = cluster_during_model_preds %>% filter(treatment == "continual"),
            aes(x = time_since_end,
                y = predicted, 
                ymin = conf.low,
                ymax = conf.high,
                # alpha = treatment,
                fill = cluster,
                group = cluster),
            alpha = 0.1) +
  # model predictions
  geom_line(data = cluster_during_model_preds %>% filter(treatment == "continual"),
             aes(x = time_since_end,
                 y = predicted,
                 color = cluster,
                 group = cluster),
            linewidth = 1) +
  # after removal 
  # 95% CI ribbon
  geom_ribbon(data = cluster_after_model_preds %>% filter(treatment == "continual"),
              aes(x = time_since_end,
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high,
                  # alpha = treatment,
                  fill = cluster,
                  group = cluster),
              alpha = 0.1) +
  # model prediction
  geom_line(data = cluster_after_model_preds %>% filter(treatment == "continual"),
            aes(x = time_since_end,
                y = predicted, 
                # alpha = treatment,
                color = cluster,
                group = cluster),
            linewidth = 1) +
  
  # actual data 
  geom_point(data = cluster_biomass %>% filter(treatment == "continual"),
             aes(x = time_since_end,
                 y = group_total,
                 color = cluster),
             alpha = 0.4) +
  
  # plot appearance 
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  geom_vline(xintercept = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 700)) +
  labs(x = "Time since end",
       y = "Cluster biomass (dry g/m\U00B2)") +
  theme(panel.grid = element_blank())

cluster_biomass_continual_plot

ggplot() +
  geom_ribbon(data = cluster_during_model_preds %>% filter(treatment == "control"),
              aes(x = time_since_end,
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high,
                  # alpha = treatment,
                  fill = cluster,
                  group = cluster),
              alpha = 0.1) +
  geom_line(data = cluster_during_model_preds %>% filter(treatment == "control"),
            aes(x = time_since_end,
                y = predicted, 
                # alpha = treatment,
                color = cluster,
                group = cluster)) +
  geom_ribbon(data = cluster_after_model_preds %>% filter(treatment == "control"),
              aes(x = time_since_end,
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high,
                  # alpha = treatment,
                  fill = cluster,
                  group = cluster),
              alpha = 0.1) +
  geom_line(data = cluster_after_model_preds %>% filter(treatment == "control"),
            aes(x = time_since_end,
                y = predicted, 
                # alpha = treatment,
                color = cluster,
                group = cluster)) +
  geom_point(data = cluster_biomass %>% filter(treatment == "control"),
             aes(x = time_since_end,
                 y = group_total,
                 color = cluster),
             alpha = 0.4) +
  geom_vline(xintercept = 0)

## b. growth forms --------------------------------------------------------

# model 1: all predictors and random effects of site and year
gf_total_biomass_during_model_1 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment + (1|site) + (1|year),
  data = gf_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 2: includes all predictors and random effect of site
gf_total_biomass_during_model_2 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment + (1|site),
  data = gf_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 3: includes all predictors, no random effects
gf_total_biomass_during_model_3 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment,
  data = gf_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 4: includes time since end and gf as only predictors, no random effects
# continual removal only, no treatment
gf_total_biomass_during_model_4 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form,
  data = gf_biomass %>% filter(exp_dates == "after" & treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# best model is 2 with all predictors and random effect of site only
# delta AIC is less than 1 so just going with both random effects
MuMIn::AICc(gf_total_biomass_during_model_1,
            gf_total_biomass_during_model_2,
            gf_total_biomass_during_model_3) %>% 
  arrange(AICc)

# residuals
simulateResiduals(gf_total_biomass_during_model_1,
                  plot = TRUE)
# larger residuals at higher values of model predictions?

gf_total_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment + (1|site),
  data = gf_biomass %>% filter(exp_dates == "after" & sample_ID != "napl_control_2023-05-18"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

gf_total_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*treatment + (1|site),
  data = gf_biomass %>% filter(exp_dates == "after" & sample_ID != "napl_control_2023-05-18"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

parameters::standard_error(gf_total_biomass_after_model)

simulateResiduals(gf_total_biomass_after_model, plot = TRUE)

summary(gf_total_biomass_after_model)

ggpredict(gf_total_biomass_during_model_1,
          terms = c("time_since_end", "sd_growth_form")) %>% 
  plot(show_data = TRUE)

gf_during_model_preds <- ggpredict(gf_total_biomass_during_model_1,
                                   terms = c("time_since_end", "sd_growth_form", "treatment")) %>% 
  rename(time_since_end = x,
         sd_growth_form = group,
         treatment = facet)

gf_after_model_preds <- ggpredict(gf_total_biomass_after_model,
                                   terms = c("time_since_end", "sd_growth_form", "treatment")) %>% 
  rename(time_since_end = x,
         sd_growth_form = group,
         treatment = facet)

gf_biomass_continual_plot <- ggplot() +
  geom_ribbon(data = gf_during_model_preds %>% filter(treatment == "continual"),
              aes(x = time_since_end,
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high,
                  # alpha = treatment,
                  fill = sd_growth_form,
                  group = sd_growth_form),
              alpha = 0.1) +
  geom_line(data = gf_during_model_preds %>% filter(treatment == "continual"),
            aes(x = time_since_end,
                y = predicted,
                color = sd_growth_form,
                group = sd_growth_form),
            linewidth = 1) +
  geom_ribbon(data = gf_after_model_preds %>% filter(treatment == "continual"),
              aes(x = time_since_end,
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high,
                  # alpha = treatment,
                  fill = sd_growth_form,
                  group = sd_growth_form),
              alpha = 0.1) +
  geom_line(data = gf_after_model_preds %>% filter(treatment == "continual"),
            aes(x = time_since_end,
                y = predicted, 
                # alpha = treatment,
                color = sd_growth_form,
                group = sd_growth_form),
            linewidth = 1)

gf_biomass_continual_plot
   
   +
  geom_point(data = cluster_biomass %>% filter(treatment == "continual"),
             aes(x = time_since_end,
                 y = group_total,
                 color = cluster),
             alpha = 0.4) +
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  geom_vline(xintercept = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 700)) +
  labs(x = "Time since end",
       y = "Cluster biomass (dry g/m\U00B2)") +
  theme(panel.grid = element_blank())

## c. functional forms ----------------------------------------------------

# model 1: all predictors and random effects of site and year
ff_total_biomass_during_model_1 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment + (1|site) + (1|year),
  data = ff_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 2: includes all predictors and random effect of site
ff_total_biomass_during_model_2 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment + (1|site),
  data = ff_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 3: includes all predictors, no random effects
ff_total_biomass_during_model_3 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment,
  data = ff_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 4: includes time since end and gf as only predictors, no random effects
# continual removal only, no treatment
ff_total_biomass_during_model_4 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form,
  data = ff_biomass %>% filter(exp_dates == "after" & treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# best model is 2 with all predictors and random effect of site only
# delta AIC is less than 1 so just going with both random effects
MuMIn::AICc(ff_total_biomass_during_model_1,
            ff_total_biomass_during_model_2,
            ff_total_biomass_during_model_3) %>% 
  arrange(AICc)

# residuals
simulateResiduals(ff_total_biomass_during_model_1,
                  plot = TRUE)
# larger residuals at higher values of model predictions?

ff_total_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment + (1|site),
  data = ff_biomass %>% filter(exp_dates == "after" & sample_ID != "napl_control_2023-05-18"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

simulateResiduals(ff_total_biomass_after_model,
                  plot = TRUE)

ff_during_model_preds <- ggpredict(ff_total_biomass_during_model_1, terms = c("time_since_end", "ll_func_form", "treatment")) %>% 
  rename(time_since_end = x,
         ll_func_form = group,
         treatment = facet)

ff_after_model_preds <- ggpredict(ff_total_biomass_after_model, terms = c("time_since_end", "ll_func_form", "treatment")) %>% 
  rename(time_since_end = x,
         ll_func_form = group,
         treatment = facet)

ff_biomass_continual_plot <- ggplot() +
  geom_ribbon(data = ff_during_model_preds %>% filter(treatment == "continual"),
              aes(x = time_since_end,
                  y = predicted,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = ll_func_form),
              alpha = 0.1) +
  geom_line(data = ff_during_model_preds %>% filter(treatment == "continual"),
            aes(x = time_since_end,
                y = predicted,
                color = ll_func_form),
            linewidth = 1) +
  geom_ribbon(data = ff_after_model_preds %>% filter(treatment == "continual"),
              aes(x = time_since_end,
                  y = predicted,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = ll_func_form),
              alpha = 0.1) +
  geom_line(data = ff_after_model_preds %>% filter(treatment == "continual"),
            aes(x = time_since_end,
                y = predicted,
                color = ll_func_form),
            linewidth = 1) +
  geom_point(data = ff_biomass %>% 
               # taking out super high points for visualization only (easier to see lines)
               filter(!(sample_ID %in% c("mohk_continual_2012-11-15", 
                                         "napl_control_2023-05-18"))),
             aes(x = time_since_end,
                 y = group_total,
                 color = ll_func_form),
             alpha = 0.1,
             shape = 21) +
  scale_color_manual(values = ff_cols) +
  scale_fill_manual(values = ff_cols) +
  labs(x = "Time since end (years)",
       y = "Biomass (dry g/m\U00B2)",
       fill = "Littler & Littler functional form",
       color = "Littler & Littler functional form") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6))

ff_biomass_continual_plot 

# ggsave(here("figures", "model-predictions", paste0("ll-func-form_total-biomass_", today(), ".jpg")),
#        ff_biomass_continual_plot,
#        width = 16,
#        height = 10,
#        units = "cm",
#        dpi = 300)

# coarsely branched above jointed calcareous during removal
# coarsely branched below jointed calcareous after removal


