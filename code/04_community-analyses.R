##########################################################################-
# Community analyses
# last modified: 2024-04-23

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
  left_join(., groups, by = "sp_code")

# find the taxa that are not included in the clustering
cluster_missing <- cluster_match %>% 
  filter(is.na(scientific_name))
  
# calculate per cluster biomass in each survey
cluster_biomass <- cluster_match %>% 
  filter(!is.na(cluster)) %>% 
  # group by sample_ID
  group_by(sample_ID, cluster) %>% 
  # add up biomass within a cluster
  summarize(cluster_total = sum(dry_gm2, na.rm = TRUE)) %>% 
  # calculate total biomass and proportion of each cluster
  mutate(total = sum(cluster_total),
         prop = cluster_total/total,
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

cluster_prop_timeseries <- ggplot(cluster_biomass,
       aes(x = time_since_end,
           y = prop,
           color = cluster,
           alpha = treatment)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_color_manual(values = c("1" = chloro_col, "2" = ochro_col),
                     name = "Cluster") +
  scale_alpha_manual(values = c("control" = 0.4, continual = 1)) +
  facet_grid(rows = vars(site)) +
  theme(panel.grid = element_blank())

cluster_prop_timeseries

# ggsave(here("figures", "community-timeseries", paste0("cluster-prop-timeseries_", today(), ".jpg")),
#        cluster_prop_timeseries,
#        height = 12,
#        width = 8)

cluster_total_timeseries <- ggplot(cluster_biomass %>% 
         filter(sample_ID != "napl_control_2023-05-18"),
       aes(x = time_since_end,
           y = cluster_total,
           color = cluster,
           alpha = treatment)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_color_manual(values = c("1" = chloro_col, "2" = ochro_col),
                     name = "Cluster") +
  scale_alpha_manual(values = c("control" = 0.4, continual = 1)) +
  facet_wrap(~site, ncol = 1, scales = "free") +
  theme(panel.grid = element_blank())

cluster_total_timeseries

# ggsave(here("figures", "community-timeseries", paste0("cluster-total-timeseries_", today(), ".jpg")),
#        cluster_total_timeseries,
#        height = 12,
#        width = 8)


##########################################################################-
# 4. proportion biomass model ---------------------------------------------
##########################################################################-

cluster_biomass_during_model <- glmmTMB(prop ~ time_since_end*cluster*treatment,
                                 data = cluster_biomass %>% filter(exp_dates == "during"),
                                 family = nbinom2)

# very weird residuals
simulateResiduals(cluster_biomass_during_model, plot = TRUE)

summary(cluster_biomass_during_model)

cluster_biomass_after_model <- glmmTMB(prop ~ time_since_end*cluster*treatment,
                                        data = cluster_biomass %>% filter(exp_dates == "after"),
                                        family = nbinom2)

# very weird residuals
simulateResiduals(cluster_biomass_after_model, plot = TRUE)

summary(cluster_biomass_after_model)

cluster_total_biomass_during_model <- glmmTMB(cluster_total ~ time_since_end*cluster*treatment,
                                              data = cluster_biomass %>% filter(exp_dates == "during"),
                                              family = ziGamma(link = "log"),
                                              ziformula = ~1)

simulateResiduals(cluster_total_biomass_during_model, plot = TRUE)

summary(cluster_total_biomass_during_model)

