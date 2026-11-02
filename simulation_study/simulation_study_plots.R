# Simulation study plots
# setup workspace --------------------------------------------------------- ####
library(here)
# load custom functions
source(here("functions/utility_functions.R"))
sourceDir(here("functions"))
# Load necessary packages
load_lib(dplyr, tidyr, stringr,          # data
         mvnfast, ambient,               # simulate example track
         ggplot2, scales, patchwork)     # plots

# settings ---------------------------------------------------------------- ####
sim_result_dir <- here("simulation_study/outputs")

# custom functions
# standardise output file names: out("Tmax_estimates.png") etc.
out <- function(suffix) file.path(results_dir, paste0("sim_study_", suffix))
# format method title
method_title <- function(method) {
  case_when(method == "bbis"  ~ "BBIS",
            method == "euler" ~ "Euler",
            method == "ssf"   ~ "SSF")
}
# plotmath labels per parameter (parsed by label_parsed)
coef_label <- function(par) {
  case_when(par == "beta1"   ~ "beta[1]",
            par == "beta2"   ~ "beta[2]",
            par == "beta3"   ~ "beta[3]",
            par == "gammasq" ~ "gamma^2")
}

# truth per parameter (beta, speed are set by the simulation script)
speed <- 5              # speed parameter for Langevin model
beta  <- c(4, 2, -0.1)  # covariate coefficients
true_pars <- data.frame(par  = c("beta1", "beta2", "beta3", "gammasq"),
                        true = c(beta[1], beta[2], beta[3], speed))

# long format + truth + titles; one row per sim x method x factor x parameter
to_long <- function(df) {
  df %>%
    pivot_longer(c(beta1, beta2, beta3, gammasq),
                 names_to = "par", values_to = "estimate") %>%
    left_join(true_pars, by = "par") %>%
    mutate(method = method_title(method),
           par    = coef_label(par))
}

# define plot specs ------------------------------------------------------- ####
results_dir <- here("simulation_study/figures")
dodge       <- 0.4
make_path(results_dir)

# Sim 1: varying delta_t, fixed number of observations -------------------- ####
# as Tmax grows the fix spacing grows with it (same n_obs); x-axis is Tmax in days
# change titles
pp_s1 <- read.csv(file.path(sim_result_dir, "varying_thin_estimates.csv")) %>%
  to_long() %>%
  mutate(x = floor(Tmax / 24))

# generate plot
p_s1 <- ggplot(pp_s1, aes(factor(x), estimate, colour = method, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +  # zero line
  geom_boxplot(outlier.size = 0.3, alpha = 0.4,
               outlier.stroke = 0.3,
               linewidth = 0.3,
               position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept = true), linetype = "dashed", colour = "#C43733") +   # truth
  scale_colour_viridis_d(end = 0.9) +
  scale_fill_viridis_d(end = 0.9) +
  facet_wrap(~ par, scales = "free_y", labeller = label_parsed) +
  labs(x = expression(T[max] ~ "(days)"), y = "Estimate",
       colour = "Method", fill = "Method") +
  theme_bw()
# p_s1

# save
ggsave(out("sim_vary_Tmax.png"), p_s1, width = 7, height = 4); browseURL(out("sim_vary_Tmax.png"))

# Sim 2: varying delta_t, fixed maximum time ------------------------------ ####
# Tmax fixed, so coarser delta means fewer fixes; x-axis is delta in hours
# change titles
pp_s2 <- read.csv(file.path(sim_result_dir, "varying_thin_estimates_fixed_Tmax.csv")) %>%
  to_long()

# generate plot
p_s2 <- ggplot(pp_s2, aes(factor(delta), estimate, colour = method, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +  # zero line
  geom_boxplot(outlier.size = 0.3, alpha = 0.4,
               outlier.stroke = 0.3,
               linewidth = 0.3,
               position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept = true), linetype = "dashed", colour = "#C43733") +   # truth
  scale_colour_viridis_d(end = 0.9) +
  scale_fill_viridis_d(end = 0.9) +
  facet_wrap(~ par, scales = "free_y", labeller = label_parsed) +
  labs(x = expression(Delta[t] ~ "(hours)"), y = "Estimate",
       colour = "Method", fill = "Method") +
  theme_bw()
# p_s2

# save
ggsave(out("sim_vary_delta_fixedTmax.png"), p_s2, width = 7, height = 4); browseURL(out("sim_vary_delta_fixedTmax.png"))

# Sim 3: varying number of bridges (M) ------------------------------------ ####
# TODO: not implemented in the original base-R script; data in
# simulation_study/varying_M_estimates.Rda if adding later

# example simulated track figure ------------------------------------------ ####
# simulate single track
set.seed(123)
## track pars 
speed <- 5              # speed parameter for Langevin model
dt    <- 1/60           # temporal resolution of simulated tracks (fraction of 1 hour)
beta  <- c(4, 2, -0.1)  # covariate coefficients
loc0  <- c(0, 0)        # starting location of tracks

## default estimation pars
n_sim  <- 2    # number of simulations per simulation study
ncores <- parallel::detectCores() - 2 # number of cores used in parallel computations
thin   <- 100      # thinning
N      <- 1/(dt*15) - 1 # default nodes
M      <- 100       # default number of bridges
n_obs  <- 2*365*24 / dt # #hours / dt

## covariate pars
res  <- 1  # resolution of covariates 
ncov <- 2  # number of covariates
scal = 100 #kilometers
ext  <- c(-1, 1, -1, 1)*scal  # extent of study area
perlin_f <- 0.05  # Perlin noise frequency

# simulate covariates with Perlin noise 
covlist <- list()
xgrid  <- seq(ext[1], ext[2], by = res)
ygrid  <- seq(ext[3], ext[4], by = res)
grid   <- expand.grid(x = xgrid, y = ygrid)   # cell centres, shared by all layers
# two Perlin-noise covariates
for (i in 1:ncov) {
  vals <- 3 * noise_perlin(c(length(xgrid), length(ygrid)), frequency = perlin_f)
  covlist[[i]] <- list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}
# squared distance to map centre as a third covariate
dist2 <- (grid$x^2 + grid$y^2) / 100
covlist[[3]] <- list(x = xgrid, y = ygrid,
                     z = matrix(dist2, length(xgrid), length(ygrid)))
# cov_names <- c("Covariate 1", "Covariate 2", "Squared distance to centre")
cov_names <- c("Covariates", "", "")

# simulate one high-res track + stationary UD (load if exisits) ----------- ####
example_track_path <- out("example_track.csv")
# import if exists
if (file.exists(example_track_path)) {
  X_full <- read.csv(example_track_path)
} else {  # simulate
  X_full <- data.frame()
  for (i in seq(n_sim)) {
    # high-res track (dt = 1 min); first row is loc0
    track <- simLMM(dt, speed, covlist, beta, loc0, n_obs) %>% 
      as.data.frame() %>% transmute(x = V1, y = V2, ID = i)
    X_full <- rbind(X_full, track) # combine with previous
  }
  # save
  write.csv(X_full, out("example_track.csv"), row.names = FALSE)
}

# stationary UD on the grid: pi(s) proportional to exp(sum_k beta_k cov_k(s))
ud_z  <- exp(beta[1] * covlist[[1]]$z +
             beta[2] * covlist[[2]]$z +
             beta[3] * covlist[[3]]$z)
ud_df <- data.frame(grid, ud = as.vector(log(ud_z / sum(ud_z))))

# thin the track to the experiment's resolutions 
# "1 min" (= full track) through 24 h; thinning factor is interval / dt
int_lab <- c("1 min", "30 min", "6 h", "24 h") 
int_hrs <-  c(1/60, 0.5, 6, 24)
# int_lab <- c("1 min", "30 min", "1 h", "2 h", "6 h", "12 h", "24 h")
# int_hrs <- c(dt, 0.5, 1, 2, 6, 12, 24)
tracks_full <- Map(function(h, lab) {
  xt <- thinTrack(X_full, round(h / dt))
  data.frame(x = xt[, 1], y = xt[, 2], interval = lab, ID = xt$ID)
}, int_hrs[1], int_lab[1]) %>%
  bind_rows() %>%
  mutate(interval = factor(interval, levels = int_lab))  # fine -> coarse

# extract segment of track
delay <- 5*24  # (hours)
days <- 7.1

nrow_start <- floor(days*24/dt)  # number of rows to get 7 days 
tracks_start <- Map(function(h, lab) {
  xt <- thinTrack(X_full[1:nrow_start + delay/dt,], round(h / dt))
  data.frame(x = xt[, 1], y = xt[, 2], interval = lab)
}, int_hrs, int_lab) %>%
  bind_rows() %>%
  mutate(interval = factor(interval, levels = int_lab))  # fine -> coarse

# crop to sample area
## track ext
ext_x <- round(range(tracks_start$x)) + c(-2, 2)
ext_y <- round(range(tracks_start$y)) + c(-2, 2)
## square around centre
centre <- c(mean(ext_x), mean(ext_y))
ext_rad <- round(max(diff(ext_x), diff(ext_y))/2)
ext_x <- c(centre[1] + c(-ext_rad, ext_rad))
ext_y <- c(centre[2] + c(-ext_rad, ext_rad))
## crop
ud_crop_df <- filter(ud_df, x >= ext_x[1], x <= ext_x[2], 
                     y >= ext_y[1], y <= ext_y[2])
## define inset
ext_path <- data.frame(x = ext_x[c(1,1,2,2,1)],
                       y = ext_y[c(1,2,2,1,1)])

# plot specs 
# white (full track) ramping to dark blue (coarsest thinning)
# track_cols <- colorRampPalette(c("white", "#08306B"))(length(int_lab))
# track_cols <- colorRampPalette(c("#EAEEF9", "#A2B5E4", "#466DCD", "#042197"))(length(int_lab))
# track_cols <- colorRampPalette(c("#7D6F5E", "#AB6235", "#CA5425", "#823614"))(length(int_lab))
# track_cols <- colorRampPalette(c("#A29A88", "#85934A", "#29582D", "#122703"))(length(int_lab))
# track_cols <- colorRampPalette(c("#7D6F5E", "#847EB2", "#4B4C7E", "#132A4F"))(length(int_lab))
track_cols <- colorRampPalette(c("#A29A88", "#54A790", "#427B89", "#1B3C42"))(length(int_lab))
# raster_col <- scale_fill_viridis_c(direction = -1, name = "Log UD")
# raster_col <- scale_fill_continuous(palette = "blues", name = "Log UD")
raster_col <- scale_fill_gradient(low = "#FBFBF3", high = "#B1A493",
                                  limits = c(-33, -5), name = "Log UD")

# left panel: the three covariates stacked, each on its own colour scale
cov_panel <- function(i) {
  ggplot(data.frame(grid, z = as.vector(covlist[[i]]$z))) +
    geom_raster(aes(x, y, fill = z)) +
    # scale_fill_viridis_c(name = NULL) +
    scale_fill_gradient(low = "#FBFBF3", high = "#B1A493", name = NULL) + 
    coord_equal(expand = FALSE) +
    # coord_cartesian(expand = FALSE) +
    # labs(subtitle = cov_names[i]) +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(c(0,1,1,1)[i], 0,
                               c(1,1,1,0)[i], 0),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none")
          # legend.key.width  = unit(0.3, "lines"),
          # legend.key.height = unit(0.5, "lines"))
}

# Panel 1: raw covs 
p_covs <- cov_panel(1) / cov_panel(2) / cov_panel(3)

# Panel 2: UD background with the track drawn at fine resolution
p_track <- ggplot(filter(tracks_full, # %>% slice(1:1000), 
                         ID == 1), aes(x, y)) +
  geom_raster(data = ud_df, aes(x, y, fill = ud), inherit.aes = FALSE) +
  # scale_fill_viridis_c(direction = -1, name = "Log UD") +
  raster_col + 
  # path
  # geom_path(colour = "white") +
  geom_segment(aes(xend = lead(x), yend = lead(y)),
               alpha = 0.1, colour = "#54A790") +
  # inset
  geom_path(data = ext_path, colour = "black") +  # inset
  # scale_colour_manual(values = track_cols[c(1, length(int_lab))], name = expression(Delta[t])) +
  # scale_alpha_discrete(range = c(1, 0.5), name = expression(Delta[t])) +
  # scale_linewidth_discrete(range = c(0.5, 1.5), name = expression(Delta[t])) +
  coord_equal(expand = FALSE) +
  theme_void(base_size = 8) + 
  theme(plot.margin = margin(0,4,0,4),
        axis.title = element_blank(), 
        legend.position = "none")

# Panel 3: UD background with the track drawn at each resolution
p_track_seg <- ggplot(tracks_start, aes(x, y)) +
  geom_raster(data = ud_crop_df, aes(x, y, fill = ud), inherit.aes = FALSE) +
  scale_fill_gradient(low = "#FBFBF3", high = "#B1A493", name = "Log UD", 
                      guide  = "none") + 
  geom_path(aes(colour = interval, linewidth = interval,
                alpha = interval,
                group = interval)) +
  scale_colour_manual(values = track_cols, name = expression(Delta[t])) +
  scale_alpha_manual(values = c(0.7, 0.8, 0.7, 0.6), name = expression(Delta[t])) +
  scale_linewidth_discrete(range = c(0.4, 1.6), name = expression(Delta[t])) +
  # geom_point(data = tibble(x = 0, y = 0), size = 2.5) +
  coord_equal(expand = FALSE) +
  theme_void(base_size = 8) + 
  theme(plot.margin = margin(0,4,0,4),
        axis.title = element_blank(),
        panel.grid = element_blank())

# combine panels
p_example <- (wrap_elements(p_covs) | p_track | p_track_seg) +
  plot_layout(widths = c(0.76, 2, 1.96)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10))  # same tag size on wrapped a) as b), c)

# save
ggsave(out("example_track.png"), 
       p_example, width = 8, height = 3.25); browseURL(out("example_track.png"))

# TODO: if the right panel is too busy, add a cropped inset of a sub-region
# as a 3rd panel via patchwork::inset_element()
