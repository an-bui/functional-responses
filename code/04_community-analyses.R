##########################################################################-
# Community analyses
# last modified: 2024-06-14

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

# ⟞ a. LTE ----------------------------------------------------------------

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
calc_group_biomass <- function(type, grouping) {
  
  if(type == "LTE") {
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
  } else if(type == "benthics") {
    benthic_cluster_match %>% 
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
                           names = c("site", "year", "transect"),
                           cols_remove = FALSE) %>% 
      select(!c(site, year, transect)) %>% 
      left_join(., benthic_comm_meta, by = "sample_ID")  
  } else {
    warning("Your type or grouping might be wrong!")
    return(NA)
  }

}

# calculate per group biomass in each survey
cluster_biomass <- calc_group_biomass("LTE", cluster)
gf_biomass <- calc_group_biomass("LTE", sd_growth_form)
ff_biomass <- calc_group_biomass("LTE", ll_func_form)

check <- cluster_match %>% 
  filter(site == "aque" & date > "2010-04-26" |
           site == "napl" & date > "2010-04-27" |
           site == "mohk" & date > "2010-05-05" |
           site == "carp" & date > "2010-04-23") %>% 
  filter(ll_func_form == "thick_leathery") %>% 
  filter(!is.na(cluster)) %>% 
  filter(treatment == "continual") %>% 
  filter(time_since_end == 3) %>% 
  pull(dry_gm2) %>% 
  sum()
# this output should match up with total in cluster_thick_leathery for a given time_since_end

cluster_thick_leathery <- cluster_match %>% 
  filter(site == "aque" & date > "2010-04-26" |
           site == "napl" & date > "2010-04-27" |
           site == "mohk" & date > "2010-05-05" |
           site == "carp" & date > "2010-04-23") %>% 
  filter(ll_func_form == "thick_leathery") %>% 
  filter(!is.na(cluster)) %>% 
  filter(treatment == "continual") %>% 
  # group by sample_ID
  group_by(time_since_end, cluster) %>% 
  # add up biomass within a cluster
  summarize(group_total = sum(dry_gm2, na.rm = TRUE)) %>% 
  # calculate total biomass and proportion of each cluster
  mutate(total = sum(group_total),
         prop = group_total/total,
         prop = case_when(
           prop == "NaN" ~ 0,
           TRUE ~ prop
         )) %>% 
  ungroup()

cluster_coarsely_branched <- cluster_match %>% 
  filter(site == "aque" & date > "2010-04-26" |
           site == "napl" & date > "2010-04-27" |
           site == "mohk" & date > "2010-05-05" |
           site == "carp" & date > "2010-04-23") %>% 
  filter(ll_func_form == "coarsely_branched") %>% 
  filter(!is.na(cluster)) %>% 
  filter(treatment == "continual") %>% 
  # group by sample_ID
  group_by(time_since_end, cluster) %>% 
  # add up biomass within a cluster
  summarize(group_total = sum(dry_gm2, na.rm = TRUE)) %>% 
  # calculate total biomass and proportion of each cluster
  mutate(total = sum(group_total),
         prop = group_total/total,
         prop = case_when(
           prop == "NaN" ~ 0,
           TRUE ~ prop
         )) %>% 
  ungroup()

thick_leathery_prop <- ggplot(data = cluster_thick_leathery,
       aes(x = time_since_end,
           y = prop,
           fill = cluster)) +
  geom_area(color = "white") + 
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(cluster2, cluster3, cluster5, cluster6, cluster7)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Time since end (years)",
       y = "Proportion",
       title = "Proportion of thick leathery biomass by cluster",
       subtitle = "Sum cluster biomass across sites for a single sampling period")

thick_leathery_prop

# ggsave(here("figures", "community-timeseries",
#             paste0("cluster-thick-leathery-timeseries_", today(), ".jpg")),
#        thick_leathery_prop,
#        width = 8,
#        height = 4)

coarsely_branched_prop <- ggplot(data = cluster_coarsely_branched,
       aes(x = time_since_end,
           y = prop,
           fill = cluster)) +
  geom_area(color = "white") + 
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(cluster2, cluster3, cluster4, cluster6)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Time since end (years)",
       y = "Proportion",
       title = "Proportion of coarsely branched biomass by cluster",
       subtitle = "Sum cluster biomass across sites for a single sampling period")

coarsely_branched_prop

# ggsave(here("figures", "community-timeseries",
#             paste0("cluster-coarsely-branched-timeseries_", today(), ".jpg")),
#        coarsely_branched_prop,
#        width = 8,
#        height = 4)

# function to widen group biomass
widen_group_biomass <- function(df) {
  df %>% 
    select(-prop) %>% 
    pivot_wider(names_from = 2,
                values_from = 3) %>% 
    relocate(total, .after = last_col())
}

cluster_biomass_wide <- widen_group_biomass(cluster_biomass) %>% 
  # renames cluster columns with cluster_ prefix
  rename_with(.cols = `1`:`7`, ~ paste0("cluster_", .))
gf_biomass_wide <- widen_group_biomass(gf_biomass)
ff_biomass_wide <- widen_group_biomass(ff_biomass)

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
) +
  scale_color_manual(values = ff_cols,
                     name = "L&L functional form")

ff_prop_timeseries

gf_prop_timeseries <- group_timeseries(
  df = gf_biomass,
  grouping = sd_growth_form,
  y = prop
) +
  scale_color_manual(values = gf_cols,
                     name = "S&D growth form")

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

cluster_mean_summary <- cluster_biomass %>% 
  filter(sample_ID != "napl_control_2023-05-18") %>% 
  group_by(cluster, time_since_end) %>% 
  summarize(mean = mean(group_total, na.rm = TRUE))

cluster_thick_leathery_summary <- cluster_thick_leathery %>% 
  group_by(cluster, time_since_end) %>% 
  summarize(mean = mean(prop, na.rm = TRUE)) 

ggplot(data = cluster_thick_leathery_summary,
       aes(x = time_since_end,
           y = mean,
           fill = cluster)) +
  geom_area() +
  scale_fill_manual(values = c(cluster2, cluster3, cluster5, cluster6, cluster7))
  

ff_mean_summary <- ff_biomass %>% 
  filter(sample_ID != "napl_control_2023-05-18") %>% 
  group_by(ll_func_form, time_since_end) %>% 
  summarize(mean = mean(group_total, na.rm = TRUE))

ggplot(data = cluster_biomass %>% 
         filter(sample_ID != "napl_control_2023-05-18") %>% 
         filter(cluster %in% c(2, 3, 5, 6, 7)),
       aes(x = time_since_end,
           y = group_total)) +
  # geom_point(shape = 21,
  #            alpha = 0.3) +
  geom_line(data = cluster_mean_summary %>% 
              filter(cluster %in% c(2, 3, 5, 6, 7)),
            aes(x = time_since_end,
                y = mean,
                color = cluster),
            linewidth = 1,
            linetype = 2) +
  scale_color_manual(values = c(cluster2, cluster3, cluster5, cluster6, cluster7)) +
  geom_line(data = ff_mean_summary %>% 
              filter(ll_func_form == "thick_leathery"),
           aes(x = time_since_end,
               y = mean),
           color = thi_lea_col,
           linewidth = 1) +
  labs(x = "Time since end (years)",
       y = "Mean group biomass") +
  theme(panel.grid = element_blank())

ggplot(data = cluster_thick_leathery,
       aes(x = time_since_end,
           y = prop,
           fill = cluster)) +
  geom_area()

group_timeseries(
  df = cluster_biomass %>% 
    filter(sample_ID != "napl_control_2023-05-18"),
  grouping = cluster,
  y = group_total
)
  
ff_total_timeseries <- group_timeseries(
  df = ff_biomass %>% 
    filter(sample_ID != "napl_control_2023-05-18"),
  grouping = ll_func_form,
  y = group_total
) +
  scale_color_manual(values = ff_cols,
                     name = "L&L functional form")

ff_total_timeseries

gf_total_timeseries <- group_timeseries(
  df = gf_biomass %>% 
    filter(sample_ID != "napl_control_2023-05-18"),
  grouping = sd_growth_form,
  y = group_total
) +
  scale_color_manual(values = gf_cols,
                     name = "S&D growth form")

gf_total_timeseries

# ggsave(here("figures", "community-timeseries", paste0("cluster7-pam-total-timeseries_", today(), ".jpg")),
#        cluster_total_timeseries,
#        height = 12,
#        width = 8)

# ⟞ b. benthics -----------------------------------------------------------

# match community data with clusters
benthic_cluster_match <- benthic_comm_df %>% 
  # join with clusters
  left_join(., pam_clusters_7, by = "sp_code")

# find the taxa that are not included in the clustering
cluster_missing <- benthic_cluster_match %>% 
  filter(is.na(scientific_name))

# calculate per group biomass in each survey
benthic_cluster_biomass <- calc_group_biomass("benthics", cluster)
benthic_gf_biomass <- calc_group_biomass("LTE", sd_growth_form)
benthic_ff_biomass <- calc_group_biomass("LTE", ll_func_form)

benthic_check <- benthic_cluster_match %>% 
  filter(ll_func_form == "thick_leathery") %>% 
  filter(!is.na(cluster)) %>% 
  filter(year == 2010) %>% 
  pull(dry_gm2) %>% 
  sum()
# this output should match up with total in cluster_thick_leathery for a given time_since_end

benthic_cluster_thick_leathery <- benthic_cluster_match %>% 
  filter(ll_func_form == "thick_leathery") %>% 
  filter(!is.na(cluster)) %>% 
  # group by year and cluster
  group_by(year, cluster) %>% 
  # add up biomass within a cluster
  summarize(group_total = sum(dry_gm2, na.rm = TRUE)) %>% 
  # calculate total biomass and proportion of each cluster
  mutate(total = sum(group_total),
         prop = group_total/total,
         prop = case_when(
           prop == "NaN" ~ 0,
           TRUE ~ prop
         )) %>% 
  ungroup()

benthic_cluster_coarsely_branched <- benthic_cluster_match %>% 
  filter(ll_func_form == "coarsely_branched") %>% 
  filter(!is.na(cluster)) %>% 
  # group by year and cluster
  group_by(year, cluster) %>% 
  # add up biomass within a cluster
  summarize(group_total = sum(dry_gm2, na.rm = TRUE)) %>% 
  # calculate total biomass and proportion of each cluster
  mutate(total = sum(group_total),
         prop = group_total/total,
         prop = case_when(
           prop == "NaN" ~ 0,
           TRUE ~ prop
         )) %>% 
  ungroup()

benthic_thick_leathery_prop <- ggplot(
  data = benthic_cluster_thick_leathery,
  aes(x = year,
      y = prop,
      fill = cluster)) +
  geom_area(color = "white") + 
  scale_fill_manual(values = c(cluster2, cluster3, cluster5, cluster6, cluster7)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportion",
       title = "BENTHICS Proportion of thick leathery biomass by cluster",
       subtitle = "Sum cluster biomass across sites for a single sampling period")

benthic_thick_leathery_prop

# ggsave(here("figures", "community-timeseries",
#             paste0("benthic_cluster-thick-leathery-timeseries_", today(), ".jpg")),
#        benthic_thick_leathery_prop,
#        width = 8,
#        height = 4)

benthic_coarsely_branched_prop <- ggplot(
  data = benthic_cluster_coarsely_branched,
       aes(x = year,
           y = prop,
           fill = cluster)) +
  geom_area(color = "white") + 
  scale_fill_manual(values = c(cluster2, cluster3, cluster4, cluster6)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportion",
       title = "BENTHICS Proportion of coarsely branched biomass by cluster",
       subtitle = "Sum cluster biomass across sites for a single sampling period")

benthic_coarsely_branched_prop

# ggsave(here("figures", "community-timeseries",
#             paste0("benthic_cluster-coarsely-branched-timeseries_", today(), ".jpg")),
#        benthic_coarsely_branched_prop,
#        width = 8,
#        height = 4)

##########################################################################-
# 4. group biomass model --------------------------------------------------
##########################################################################-

## a. clusters ------------------------------------------------------------

# model 1: includes all predictors and random effects of site and year
cluster_group_biomass_during_model_1 <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment + (1|site) + (1|year),
  data = cluster_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1
  )

# model 2: includes all predictors and random effect of site
cluster_group_biomass_during_model_2 <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment + (1|site),
  data = cluster_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1
  )

# model 3: includes all predictors, no random effects
cluster_group_biomass_during_model_3 <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment,
  data = cluster_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 4: includes time since end and cluster as only predictors, no random effects
cluster_group_biomass_during_model_4 <- glmmTMB(
  group_total ~ time_since_end*cluster,
  data = cluster_biomass %>% filter(exp_dates == "after" & treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# best model is 1 with all predictors and random effect of site and year
model.sel(cluster_group_biomass_during_model_1,
          cluster_group_biomass_during_model_2,
          cluster_group_biomass_during_model_3)

# looking at model residuals
simulateResiduals(cluster_group_biomass_during_model_1, plot = TRUE)

# model summary
summary(cluster_group_biomass_during_model_1)

cluster_group_model_summary <- tbl_regression(
  cluster_group_biomass_during_model_1,
  intercept = TRUE, 
  add_estimate_to_reference_rows = TRUE
  ) %>% 
  gtsummary::as_flex_table()

cluster_group_model_summary

# save table as doc
# save_as_docx(cluster_group_model_summary,
#              path = here("tables", "model-summaries", paste0("cluster-total-model-summary_", today(), ".docx")))

cluster_group_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*cluster*treatment + (1|site),
  data = cluster_biomass %>% filter(exp_dates == "after"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

simulateResiduals(cluster_group_biomass_after_model, plot = TRUE)

testOutliers(cluster_group_biomass_after_model)

summary(cluster_group_biomass_after_model)

cluster_during_model_preds <- ggpredict(cluster_group_biomass_during_model_1, terms = c("time_since_end", "cluster", "treatment")) %>% 
  rename(time_since_end = x,
         cluster = group,
         treatment = facet)

cluster_after_model_preds <- ggpredict(cluster_group_biomass_after_model, terms = c("time_since_end", "cluster", "treatment")) %>% 
  rename(time_since_end = x,
         cluster = group,
         treatment = facet)

cluster_biomass_continual_plot <- ggplot() +
  # actual data 
  geom_point(data = cluster_biomass %>% 
               filter(treatment == "continual") %>% 
               filter(!(sample_ID %in% c("mohk_continual_2012-11-15",
                                         "napl_control_2023-05-18"))),
             aes(x = time_since_end,
                 y = group_total,
                 color = cluster),
             alpha = 0.2,
             shape = 21) +
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
  
  # plot appearance 
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  geom_vline(xintercept = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 420)) +
  labs(x = "Time since end",
       y = "Biomass (dry g/m\U00B2)") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6))

cluster_biomass_continual_plot
# cluster 7 is tall with stipes: eisenia, laminaria, pterygophora
# cluster 5 is stephanocystis, egregia, sargassum horneri/muticum
# cluster 2 is acrosorium, chondracanthus, callophyllis, cryptopleura, desmarestia, gloiocladia, mazzaella, osmundea, nienburgia, polyneura
# cluster 3 is cryptopleura, dictyota, gymnogrongrus, phycodrys, prionitis, rhodymenia, sarcodiotheca, stenogramma
# cluster 1 is amphiroa, bossiella, calliarthron, corallina, lithothrix

# ggsave(here("figures", "model-predictions", paste0("cluster_group-biomass_", today(), ".jpg")),
#        cluster_biomass_continual_plot,
#        width = 16,
#        height = 10,
#        units = "cm",
#        dpi = 300)

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
gf_group_biomass_during_model_1 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment + (1|site) + (1|year),
  data = gf_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 2: includes all predictors and random effect of site
gf_group_biomass_during_model_2 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment + (1|site),
  data = gf_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 3: includes all predictors, no random effects
gf_group_biomass_during_model_3 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment,
  data = gf_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 4: includes time since end and gf as only predictors, no random effects
# continual removal only, no treatment
gf_group_biomass_during_model_4 <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form,
  data = gf_biomass %>% filter(exp_dates == "after" & treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# best model is 2 with all predictors and random effect of site only
# delta AIC is less than 1 so just going with both random effects
model.sel(gf_group_biomass_during_model_1,
          gf_group_biomass_during_model_2,
          gf_group_biomass_during_model_3) 

# residuals
simulateResiduals(gf_group_biomass_during_model_1,
                  plot = TRUE)
# larger residuals at higher values of model predictions?

gf_group_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*sd_growth_form*treatment + (1|site),
  data = gf_biomass %>% filter(exp_dates == "after" & sample_ID != "napl_control_2023-05-18"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# gf_group_biomass_after_model <- glmmTMB(
#   group_total ~ time_since_end*treatment + (1|site),
#   data = gf_biomass %>% filter(exp_dates == "after" & sample_ID != "napl_control_2023-05-18"),
#   family = ziGamma(link = "log"),
#   ziformula = ~1)

parameters::standard_error(gf_group_biomass_after_model)

simulateResiduals(gf_group_biomass_after_model, plot = TRUE)

summary(gf_group_biomass_after_model)

ggpredict(gf_group_biomass_during_model_1,
          terms = c("time_since_end", "sd_growth_form")) %>% 
  plot(show_data = TRUE)

gf_during_model_preds <- ggpredict(gf_group_biomass_during_model_1,
                                   terms = c("time_since_end", "sd_growth_form", "treatment")) %>% 
  rename(time_since_end = x,
         sd_growth_form = group,
         treatment = facet)

gf_after_model_preds <- ggpredict(gf_group_biomass_after_model,
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
            linewidth = 1) +
  scale_color_manual(values = gf_cols,
                     name = "S&D growth form") +
  scale_fill_manual(values = gf_cols,
                     name = "S&D growth form")

gf_biomass_continual_plot

## c. functional forms ----------------------------------------------------

# model 1: all predictors and random effects of site and year
ff_group_biomass_during_model_1 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment + (1|site) + (1|year),
  data = ff_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 2: includes all predictors and random effect of site
ff_group_biomass_during_model_2 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment + (1|site),
  data = ff_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 3: includes all predictors, no random effects
ff_group_biomass_during_model_3 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment,
  data = ff_biomass %>% filter(exp_dates == "during"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# model 4: includes time since end and gf as only predictors, no random effects
# continual removal only, no treatment
ff_group_biomass_during_model_4 <- glmmTMB(
  group_total ~ time_since_end*ll_func_form,
  data = ff_biomass %>% filter(exp_dates == "after" & treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

# best model is 2 with all predictors and random effect of site only
# delta AIC is less than 1 so just going with both random effects
model.sel(ff_group_biomass_during_model_1,
          ff_group_biomass_during_model_2,
          ff_group_biomass_during_model_3)

# residuals
simulateResiduals(ff_group_biomass_during_model_1,
                  plot = TRUE)
# larger residuals at higher values of model predictions?

ff_group_biomass_after_model <- glmmTMB(
  group_total ~ time_since_end*ll_func_form*treatment + (1|site),
  data = ff_biomass %>% filter(exp_dates == "after" & sample_ID != "napl_control_2023-05-18"),
  family = ziGamma(link = "log"),
  ziformula = ~1)

simulateResiduals(ff_group_biomass_after_model,
                  plot = TRUE)

ff_during_model_preds <- ggpredict(ff_group_biomass_during_model_1, terms = c("time_since_end", "ll_func_form", "treatment")) %>% 
  rename(time_since_end = x,
         ll_func_form = group,
         treatment = facet)

ff_after_model_preds <- ggpredict(ff_group_biomass_after_model, terms = c("time_since_end", "ll_func_form", "treatment")) %>% 
  rename(time_since_end = x,
         ll_func_form = group,
         treatment = facet)

ff_biomass_continual_plot <- ggplot() +
  geom_point(data = ff_biomass %>% 
               # taking out super high points for visualization only (easier to see lines)
               filter(!(sample_ID %in% c("mohk_continual_2012-11-15", 
                                         "napl_control_2023-05-18"))),
             aes(x = time_since_end,
                 y = group_total,
                 color = ll_func_form),
             alpha = 0.2,
             shape = 21) +
  
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
  
  # plot appearance 
  scale_color_manual(values = ff_cols) +
  scale_fill_manual(values = ff_cols) +
  geom_vline(xintercept = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 420)) +
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


##########################################################################-
# 5. total biomass model --------------------------------------------------
##########################################################################-

## a. clusters ------------------------------------------------------------

ggplot(cluster_biomass_wide,
       aes(x = log(total))) +
  geom_histogram(bins = 15,
                 color = "black",
                 fill = "blue")

ggplot(cluster_biomass_wide,
       aes(x = total)) +
  geom_histogram(bins = 15,
                 color = "black",
                 fill = "blue")

data <- cluster_biomass_wide %>% 
  filter(exp_dates == "during",
         treatment == "continual") 
  # scaling cluster biomasses
  # mutate(across(18:24, ~ scale(.)[1:length(.), 1]))
  # mutate(total = scale(total))
  # mutate(across(.cols = cluster_1:cluster_2, ~ scale)) %>% 
  # rename_with(~str_remove_all(., fixed("[, 1]")))
  # filter(!(sample_ID %in% c("mohk_continual_2012-11-15",
  #                           "mohk_continual_2017-08-11",
  #                           "napl_continual_2015-05-19",
  #                           "napl_continual_2013-08-15", 
  #                           "mohk_continual_2014-11-18", 
  #                           "napl_continual_2014-05-21")))
  # filter(total > 0) %>% 
  # filter(!(sample_ID %in% outliers))

zeros <- data %>% 
  select(cluster_1:cluster_7) %>% 
  pivot_longer(cols = cluster_1:cluster_7,
               names_to = "cluster",
               values_to = "biomass") %>% 
  group_by(cluster) %>% 
  summarize(morethan = length(cluster[biomass>0]))

ggplot(data,
       aes(x = cluster_1,
           y = total)) +
  geom_point()

ggplot(data,
       aes(x = total)) +
  geom_histogram()

outliers <- c("carp_continual_2015-08-18", # 69
              "carp_continual_2016-02-22", # 71
              "carp_control_2011-01-25", # 80
              "carp_control_2011-10-17", # 82
              "carp_control_2013-08-16", # 86
              "carp_control_2016-02-22" # 93
              )

GGally::ggpairs(
  data = data,
  columns = c("total", "cluster_2", "cluster_3", "cluster_5", "cluster_7")
)

# model 1: includes all predictors and random effects of site and year
cluster_total_biomass_during_model_1 <- glmmTMB(
  total ~ cluster_2 + cluster_3 + cluster_5 + cluster_7 + (1|site) + (1|year),
  data = data
)

simulated_res <- simulateResiduals(cluster_total_biomass_during_model_1)

plot(simulated_res)

testZeroInflation(simulated_res)
testDispersion(simulated_res)

hist(simulated_res)

plotResiduals(simulated_res, data$time_since_end)

sim <- simulateResiduals(cluster_total_biomass_during_model_1)
which(residuals(sim) == 1 | residuals(sim) == 0)
# 56, 84, 90

summary(cluster_total_biomass_during_model_1)

clust2_plot <- ggpredict(cluster_total_biomass_during_model_1, terms = c("cluster_2")) %>% plot(show_data = TRUE) +
  labs(title = "cluster 2")

clust3_plot <- ggpredict(cluster_total_biomass_during_model_1, terms = c("cluster_3")) %>% plot(show_data = TRUE) +
  labs(title = "cluster 3")

clust5_plot <- ggpredict(cluster_total_biomass_during_model_1, terms = c("cluster_5")) %>% plot(show_data = TRUE) +
  labs(title = "cluster 5")

clust7_plot <- ggpredict(cluster_total_biomass_during_model_1, terms = c("cluster_7")) %>% plot(show_data = TRUE) +
  labs(title = "cluster 7")

(clust2_plot + clust3_plot) / (clust5_plot + clust7_plot)

performance::check_overdispersion(cluster_total_biomass_during_model_1)
outliers(cluster_total_biomass_during_model_1)


# observations 67 69 71 80 82 87 93 97

#
cluster_test <- glmmTMB(
  log(total) ~ time_since_end + cluster_1 + cluster_2 + cluster_3 + cluster_4 + cluster_5 + cluster_6 + cluster_7 + treatment + (1|site) + (1|year),
  data = cluster_biomass_wide %>% filter(exp_dates == "during"),
  family = inverse.gaussian()
)

cluster_test2 <- glmmTMB(
  total ~ time_since_end + cluster_1 + cluster_2 + cluster_3 + cluster_4 + cluster_5 + cluster_6 + cluster_7 + treatment + (1|site) + (1|year),
  data = cluster_biomass_wide %>% filter(exp_dates == "during"),
  family = tweedie(link = "log")
)

simulateResiduals(cluster_test2) %>% plot()

## b. growth forms --------------------------------------------------------

data <- gf_biomass_wide %>% 
  filter(exp_dates == "during") %>% 
  filter(total > 0) 
  

# model 1: includes all predictors and random effects of site and year
gf_total_biomass_during_model_1 <- glmmTMB(
  total ~ time_since_end + articulated_calcareous + corticated_foliose + corticated_macrophytes + crustose + filamentous_algae + foliose + leathery_macrophyte + treatment + (1|site) + (1|year),
  data = data,
  family = ziGamma(link = "log"),
  ziformula = ~1,
  na.action = na.pass
)

simulateResiduals(gf_total_biomass_during_model_1) %>% plot()

## c. functional forms ----------------------------------------------------

data <- ff_biomass_wide %>% 
  filter(exp_dates == "during") %>% 
  filter(total > 0)  

GGally::ggpairs(data = data,
                columns = c("total", "coarsely_branched", "jointed_calcareous", "thick_leathery"))

# testing this model: ignore
ff_total_biomass_during_model_1 <- glmmTMB(
  total ~ coarsely_branched + jointed_calcareous + thick_leathery + (1|site) + (1|year),
  data = data %>% filter(treatment == "continual")
)

# model 1: includes all predictors and random effects of site and year
ff_total_biomass_during_model_1 <- glmmTMB(
  total ~ time_since_end + coarsely_branched + jointed_calcareous + thick_leathery + (1|site) + (1|year),
  data = data %>% filter(treatment == "continual"),
  family = ziGamma(link = "log"),
  ziformula = ~ 1
)

simulateResiduals(ff_total_biomass_during_model_1) %>% plot()

cb_plot <- ggpredict(ff_total_biomass_during_model_1, terms = c("coarsely_branched")) %>% plot(show_data = TRUE) +
  labs(title = "coarsely branched")

jc_plot <- ggpredict(ff_total_biomass_during_model_1, terms = c("jointed_calcareous")) %>% plot(show_data = TRUE) +
  labs(title = "jointed calcareous")

tl_plot <- ggpredict(ff_total_biomass_during_model_1, terms = c("thick_leathery")) %>% plot(show_data = TRUE) +
  labs(title = "thick leathery")

cb_plot | jc_plot | tl_plot

model.sel(ff_total_biomass_during_model_1,
          cluster_total_biomass_during_model_1)

AICc(ff_total_biomass_during_model_1,
     cluster_total_biomass_during_model_1,
     gf_total_biomass_during_model_1) %>% 
  arrange(AICc)

