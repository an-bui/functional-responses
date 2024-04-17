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
trait_groups <- NbClust(trait_gower_daisy,
                        min.nc = 2,
                        max.nc = 45, 
                        method = "ward.D2", 
                        index = "gap")

trait_groups

# new function
my_fviz_nbclust(trait_groups)

# vlsualizing the dendrogram
trait_dendro <- fviz_dend(trait_clust, 
                          k = 2,
                          k_colors = c(chloro_col, ochro_col)) +
  theme(text = element_text(size = 12),
        axis.title = element_blank()) +
  coord_flip()

trait_clusters <- cutree(trait_clust, 
                         k = 2) %>% 
  enframe(name = "label",
          value = "cluster") %>% 
  mutate(cluster = fct_relevel(as_factor(cluster), "1", "2"))

trait_dendro <- dendro_data(trait_clust, 
                            type = "rectangle")

trait_dendro_segments <- segment(trait_dendro) %>% 
  select(x, xend, y, yend) %>% 
  rownames_to_column("segment_num") %>% 
  mutate(cluster = case_when(
    x > 17.5 ~ 2,
    x < 17.5 ~ 1,
  )) %>% 
  mutate(cluster = case_when(
    segment_num == 3 ~ 2,
    segment_num == 1 ~ 1,
    TRUE ~ cluster
  ),
  cluster = fct_relevel(as_factor(cluster), "1", "2"))

trait_dendro_labels <- label(trait_dendro) %>% 
  left_join(., trait_clusters, by = "label")


trait_dendro_plot <- ggplot() + 
  geom_segment(data = trait_dendro_segments, 
               aes(x = x, 
                   y = y, 
                   xend = xend, 
                   yend = yend,
                   color = cluster),
               linewidth = 1) +
  geom_text(data = trait_dendro_labels,
            aes(x = x, 
                y = y, 
                label = label,
                hjust = 0,
                color = cluster),
            size = 12) +
  coord_flip() +
  scale_color_manual(values = c("1" = chloro_col, "2" = ochro_col)) +
  scale_y_reverse(expand = c(0.2, 0), 
                  labels = scales::wrap_format(10)) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


# ggsave(here("figures", "trait-ordination", paste0("dendrogram_", today(), ".jpg")),
#        trait_dendro_plot,
#        height = 20,
#        width = 30)

# extract clusters
groups <- cutree(trait_clust, k = 2) %>% 
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

trait_pcoa_scores_only <- scores(trait_pcoa, choices = c(1, 2, 3)) %>% 
  as_tibble(rownames = NA)
trait_nmds_scores <- scores(trait_nmds, 
                            choices = c(1, 2),
                            tidy = TRUE) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  left_join(., algae_taxa, by = "scientific_name") %>% 
  left_join(., trait_clusters, by = c("scientific_name" = "label"))
trait_nmds_scores_only <- scores(trait_nmds, 
                                 choices = c(1, 2)) %>% 
  as_tibble(rownames = NA)
  

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
  geom_point(aes(color = cluster, shape = taxon_phylum),
             size = 3,
             alpha = 0.9) +
  scale_color_manual(values = c("1" = chloro_col, "2" = ochro_col),
                     name = "Cluster") +
  scale_shape_manual(values = c(16, 15, 17),
                     name = "Phylum") +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, 0.55)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.45, 0.55)) +
  theme(
    # legend.title = element_blank(),
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



