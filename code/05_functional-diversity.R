# this is junk!!!

zero_spp <- c("ANPA", "FR", "EA", "IR", "NEO", "SELO")

trait_matrix_reduced <- trait_matrix %>% 
  select(stipe, midrib, branching, blades, calcification)

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

algae_traits_cat_reduced <- algae_traits_cat %>% 
  filter(trait_name %in% colnames(trait_matrix_reduced))

comm_mat_bin <- comm_mat_algae %>% 
  as.data.frame() %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) %>% 
  as.matrix()
  

algae_fd <- dbFD(x = trait_matrix,
                 a = comm_mat_algae,
                 corr = "none")

spp_ric <- algae_fd$nbsp %>% 
  enframe() %>% 
  left_join(., comm_meta_algae, by = c("name" = "sample_ID"))

fric <- algae_fd$FRic %>% 
  enframe() %>% 
  left_join(., comm_meta_algae, by = c("name" = "sample_ID"))

algae_div <- vegan::diversity(x = comm_mat_algae,
                              index = "shannon") %>% 
  enframe() %>% 
  rename(shannon = value)

raoq <- algae_fd$RaoQ %>% 
  enframe() %>% 
  left_join(., comm_meta_algae, by = c("name" = "sample_ID")) %>% 
  left_join(., algae_div, by = "name") %>% 
  mutate(redund = shannon - value)

redund_plot <- ggplot(raoq,
                      aes(x = quality,
                          y = redund)) +
  geom_violin() +
  geom_point(alpha = 0.05) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Functional redundancy")

spric_plot <- ggplot(spp_ric,
                     aes(x = quality,
                         y = value)) +
  geom_violin() +
  geom_point(alpha = 0.05) +
  stat_summary(geom = "pointrange",
               fun.data = "mean_cl_boot",
               color = "red") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~exp_dates, ncol = 1) +
  labs(title = "Species richness")

fric_plot <- ggplot(fric,
                    aes(x = quality,
                        y = value)) +
  geom_violin() +
  geom_point(alpha = 0.05) +
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

sp_faxes_coord_algae <- algae_quality$"details_fspaces"$"sp_pc_coord"

alpha_fd_indices_algae <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_algae[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = comm_mat_algae,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

algae_fe <- sp.to.fe(
  sp_tr = trait_matrix_reduced,
  tr_cat = algae_traits_cat_reduced,
  fe_nm_type = "fe_rank",
  check_input = TRUE
)

# functional entity names
algae_fe$"fe_nm"

# how each species is distributed into a functional entity
algae_fe$"sp_fe"

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

fred <- ggplot(data = df,
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

frich <- ggplot(data = df,
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





