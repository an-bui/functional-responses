##########################################################################-
# Trait clustering
# last modified: 2024-04-23

# This is a script to cluster categorical and continuous traits based on 
# Gower dissimilarity. It depends on `02_data-cleaning.R`, which
# organizes the trait data into a matrix. It feeds into `04_community_
# analyses.R`, which takes the clusters and attaches them to the biomass
# data from the LTE surveys. 
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
# 3. trait clustering -----------------------------------------------------
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

trait_groups_pam <- cluster::pam(x = trait_gower_daisy,
                                 k = 2,
                                 metric = "euclidean")

trait_groups_pam

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

trait_dendro_plot

# ggsave(here("figures", "trait-clustering", paste0("dendrogram_", today(), ".jpg")),
#        trait_dendro_plot,
#        height = 20,
#        width = 30)

# extract clusters
groups <- cutree(trait_clust, k = 2) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("scientific_name") %>% 
  rename(cluster = value) %>% 
  mutate(cluster = as_factor(cluster)) %>% 
  left_join(., algae_taxa, by = c("scientific_name"))

##########################################################################-
# 4. visualizing Gower dissimilarity with clusters ------------------------
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

