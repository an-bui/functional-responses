##########################################################################-
# Community analyses
# last modified: 2024-04-09

# This is a script to get the LTE community data and trait data together
# into a site by trait matrix, then do PERMANOVA to determine effect of
# treatment and time period on trait composition at each site.
##########################################################################-

##########################################################################-
# 1. source ---------------------------------------------------------------
##########################################################################-

# only need to do this once per session
source(here::here("code", "02_data-cleaning.R"))

##########################################################################-
# 2. Gower dissimilarity --------------------------------------------------
##########################################################################-

# creating gower dissimilarity matrix
trait_gower <- gowdis(trait_matrix)
trait_gower_daisy <- daisy(trait_matrix,
                           metric = "gower") # outputs are the same

##########################################################################-
# 2. trait clustering -----------------------------------------------------
##########################################################################-

# Darling et al: gower dist, then wards clustering
# create clusters using Ward clustering using gower distance
trait_clust <- hclust(d = trait_gower_daisy, 
                      method = "ward.D2") 

# determine the best number of clusters
trait_groups <- NbClust(diss = trait_gower_daisy, 
                        distance = NULL, 
                        min.nc = 2,
                        max.nc = 45, 
                        method = "ward.D2", 
                        index = "silhouette")

trait_groups
trait_groups$All.index
plot(trait_groups$All.index)
# 9 groups, though not sure if the "max" number should be different

# vlsualizing the histogram
fviz_dend(trait_clust, k = 35)

# extract clusters
groups <- cutree(trait_clust, k = 9) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  rename(cluster = value)

##########################################################################-
# 2. getting species x trait matrix using Gower dissimilarity -------------
##########################################################################-


# doing PCoA to get dimensions
trait_pcoa <- wcmdscale(d = trait_gower)
trait_nmds <- metaMDS(trait_gower_daisy)

# extracting scores
trait_pcoa_scores <- scores(trait_pcoa, choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  left_join(., algae_taxa, by = "scientific_name")
trait_nmds_scores <- scores(trait_nmds, choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  left_join(., algae_taxa, by = "scientific_name")

# plotting PCoA axes
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

# plotting NMDS axes
trait_nmds_plot <- ggplot(trait_nmds_scores,
                          aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = taxon_phylum, shape = taxon_phylum),
             size = 3,
             alpha = 0.9) +
  scale_color_manual(values = c("Chlorophyta" = chloro_col,
                                "Ochrophyta" = ochro_col,
                                "Rhodophyta" = rhodo_col)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, 0.55)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.45, 0.55)) +
  theme(
    legend.title = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid = element_blank()
  )

trait_nmds_plot

# ggsave(filename = here("figures",
#                        "trait-ordination",
#                        paste0("gower-traits_", today(), ".jpg")),
#        plot = trait_pcoa_plot,
#        dpi = 300,
#        width = 6,
#        height = 6)

# ggsave(filename = here("figures",
#                        "trait-ordination",
#                        paste0("gower-traits_nmds_", today(), ".jpg")),
#        plot = trait_nmds_plot,
#        dpi = 300,
#        width = 8,
#        height = 6)

##########################################################################-
# 3. site species ordination ----------------------------------------------
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
# 5. creating site x trait matrix -----------------------------------------
##########################################################################-

# creating species x trait data
spp_trait_data <- scores(trait_nmds, choices = c(1, 2), display = "sites") %>% 
  as_tibble(rownames = NA) %>% 
  # mutate(Dim1_new = Dim1 + 1,
  #        Dim2_new = Dim2 + 2) %>% 
  # select(scientific_name, Dim1_new, Dim2_new) %>% 
  # column_to_rownames("scientific_name") %>% 
  as.matrix()

# monitoring data
comm_mat_algae_matrix <- comm_mat_algae %>% 
  # putting the columns in the right order
  select(rownames(spp_trait_data)) %>% 
  as.matrix()

# double checking that species are in the right order
rownames(spp_trait_data) == colnames(comm_mat_algae_matrix)

# site x trait data
site_by_trait <- comm_mat_algae_matrix %*% spp_trait_data

##########################################################################-
# 6. doing an NMDS --------------------------------------------------------
##########################################################################-

site_trait_nmds <- metaMDS(comm = site_by_trait)

site_trait_scores <- scores(site_trait_nmds, choices = c(1, 2), display = "sites") %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("sample_ID") %>% 
  left_join(., comm_meta_algae, by = "sample_ID") %>% 
  # taking out weird points for now
  filter(!(sample_ID %in% c("napl_control_2023-05-18", "mohk_continual_2012-11-15")))

ggplot(site_trait_scores,
       aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = treatment))

# ggsave(filename = here("figures",
#                        "trait-ordination",
#                        paste0("site-trait_nmds_", today(), ".jpg")),
#        dpi = 300,
#        width = 8,
#        height = 6)



