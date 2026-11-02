# reindeer simulation plots
# setup workspace --------------------------------------------------------- ####
library(here)
# load custom functions
source(here("functions/utility_functions.R"))
sourceDir(here("functions"))                             # ud_from_beta(), eval_ud(), jack_se(), ...
sourceDir(here("case_studies/reindeer/R"))
# Load necessary packages
load_lib(terra, dplyr, tidyr, lubridate, stringr,  # data
         ggplot2, patchwork, scales)               # plots

# settings ---------------------------------------------------------------- ####
results_dir <- here("simulation_study/outputs/reindeer_simulation")
rd_results  <- here("case_studies/reindeer/results")   # reindeer fits for comparison
rd_exp      <- "thinning_fixed_hour"                   # matching reindeer experiment

# custom functions
# standardise output file names: out("estimates.png") etc.
out <- function(suffix) file.path(results_dir, paste0("reindeer_simulation_", suffix))
# format method title
method_title <- function(method) {
  case_when(method == "bbis"  ~ "BBIS",
            method == "euler" ~ "Euler",
            method == "ssf"   ~ "SSF")
}
# format covariate title
coef_title <- function(coef) {
  coef <- str_replace_all(coef, "_", " ") %>% str_to_sentence() %>%
    str_replace("Resevoirs", "Reservoirs")  # fix typo in source column name
  ifelse(str_detect(coef, "2$"),  # quadratics
         str_replace(coef, "2$", "²"), coef)
}

# read results ------------------------------------------------------------ ####
# estimates of all replicates (fitted on the server): one row per
# sim x method x sampling interval
est <- list.files(results_dir, pattern = "^reindeer_landscape_estimates_sim\\d+\\.csv$",
                  full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  # optim code 0 = converged (bbis); euler is closed form, stored as 1
  mutate(converged = method == "euler" | convergence == 0)
true_pars <- read.csv(file.path(results_dir, "true_pars.csv"))

# quick checks: replicates, convergence and compute time per method x interval
est %>%
  group_by(method, delta) %>%
  summarise(n_sim = n(), converged = mean(converged), n_obs = mean(n_obs),
            median_mins = median(time) / 60, .groups = "drop")

# long format: one row per sim x method x interval x parameter, with the truth
est_long <- est %>%
  pivot_longer(c(starts_with("beta_"), gammasq), names_to = "par", values_to = "estimate") %>%
  left_join(true_pars, by = "par") %>%
  mutate(par = sub("^beta_", "", par),
         bias = estimate - true)

# reindeer comparison estimates: the same change from 3 h in the real data (BBIS);
# the fixed-hour thinning run if it has been summarised, else the first thinning
# run (M = 200; its 6-24 h levels used the old rolling thinning, ~21 h at "24 h")
rd_file <- file.path(rd_results, paste0(rd_exp, "_beta_summary.csv"))
rd <- if (file.exists(rd_file)) {
  read.csv(rd_file) %>%
    filter(method == "bbis")
} else {
  read.csv(file.path(rd_results, "thinning_beta_summary.csv")) %>%
    filter(method == "bbis", lvl_M == 200) %>%
    rename(interval = lvl_interval)
}
rd_shift <- rd %>%
  group_by(coef) %>%
  mutate(d_3h = mean - mean[interval == min(interval)]) %>%
  ungroup() %>%
  filter(interval != min(interval)) %>%
  select(par = coef, delta = interval, d_3h)

# habitat betas only (drop dist2_centre so it matches the reindeer vector)
hab <- setdiff(sub("^beta_", "", grep("^beta_", true_pars$par, value = TRUE)), "dist2_centre")

# define plot specs ------------------------------------------------------- ####
dodge <- 0.3  # horizontal spread so methods don't overlap
sim_var <- unique(est_long$delta)  # sampling intervals (h)
# x-axis factor levels and labels: sim_var (hours) -> "3 h" etc
x_lvls <- sim_var
x_labs <- paste(sim_var, "h")

# plot output directory
results_dir <- here("simulation_study/figures")

# plot: estimates vs sampling interval ------------------------------------ ####
# every converged replicate as a box; dashed line = true value
# change titles + drop extreme outliers (per facet, > 6 * MAD from median)
cut <- est_long %>%
  filter(converged, par %!in% c("dist2_centre", "gammasq")) %>%
  mutate(par = coef_title(par), method = method_title(method)) %>%
  group_by(par) %>%
  mutate(med = median(estimate), lim = 6 * mad(estimate)) %>%
  ungroup()
pp_est <- filter(cut, abs(estimate - med) <= lim)
# tally the trimmed points per facet x delta x method, above vs below
excl <- cut %>%
  filter(abs(estimate - med) > lim) %>%
  mutate(side = ifelse(estimate > med, "hi", "lo")) %>%
  count(par, delta, method, side)

# generate plot
p_est <- ggplot(pp_est, aes(factor(delta, x_lvls, x_labs), estimate,
                            colour = method, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +               # zero line
  geom_boxplot(outlier.size = 0.3, alpha = 0.4,
               outlier.stroke = 0.3, 
               box.linewidth = 0.3, median.linewidth = 0.4,
               position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept = true), linetype = "dashed", colour = "#C43733") +   # truth
  geom_text(data = filter(excl, side == "hi"),
            aes(factor(delta, x_lvls, x_labs), Inf,
                label = paste0("\u25B2", n), colour = method, group = method),
            inherit.aes = FALSE, position = position_dodge(width = 0.8),
            vjust = 1.1, size = 2.2, show.legend = FALSE) +
  geom_text(data = filter(excl, side == "lo"),
            aes(factor(delta, x_lvls, x_labs), -Inf,
                label = paste0("\u25BC", n), colour = method, group = method),
            inherit.aes = FALSE, position = position_dodge(width = 0.8),
            vjust = -0.1, size = 2.2, show.legend = FALSE) +
  scale_colour_viridis_d(end = 0.9) +
  scale_fill_viridis_d(end = 0.9) +
  facet_wrap(~ par, scales = "free_y") +
  labs(x = "Sampling interval", y = "Estimate",
       colour = "Method", fill = "Method") +
  theme_bw()
# p_est

# save
ggsave(out("beta_est.png"), p_est, width = 7, height = 5); browseURL(out("beta_est.png"))
# 
# # plot: simulated tracks over a covariate --------------------------------- ####
# ids     <- c(1:4)             # simulated animals to plot (1..n_animals)
# cov     <- "mountain_veg"     # covariate drawn under the tracks
# id_cols <- c(`1` = "white",   `2` = "#E41A1C",
#              `3` = "#FF7F00", `4` = "#F781BF")   # colour per animal
# 
# # selected animals at each sampling schedule
# X_sel <- X_full[X_full$ID %in% ids, ]
# tracks_plot <- bind_rows(
#   mutate(X_fine[X_fine$ID %in% ids, ], schedule = "simulated (10 min)"),
#   mutate(thin_track(X_sel, 3),  schedule = "3 h"),
#   mutate(thin_track(X_sel, 12), schedule = "12 h"),
#   mutate(thin_track(X_sel, 24), schedule = "24 h")) %>%
#   mutate(schedule = factor(schedule, c("simulated (10 min)", "3 h", "12 h", "24 h")),
#          ID = factor(ID))
# 
# # fixes per animal at each schedule, for the subtitle
# n_fix <- table(droplevels(tracks_plot$schedule[tracks_plot$schedule != "simulated (10 min)"])) /
#   length(ids)
# 
# # covariate cropped to the tracks (+ 3 km), as a data.frame for ggplot
# bb     <- ext(range(tracks_plot$x) + c(-3, 3), range(tracks_plot$y) + c(-3, 3))
# cov_df <- as.data.frame(crop(covs[[cov]], bb), xy = TRUE) %>%
#   setNames(c("x", "y", "value"))
# 
# # generate plot
# p_tracks <- ggplot(tracks_plot, aes(x, y, colour = ID, group = ID)) +
#   # covariate, greyscale so the tracks stand out
#   geom_raster(data = cov_df, aes(x, y, fill = value), inherit.aes = FALSE) +
#   # scale_fill_gradient(low = "black", high = "grey90", name = cov) +
#   # tracks: path for all, points for the sampled fixes only
#   geom_path(linewidth = 0.2, alpha = 0.8) +
#   geom_point(data = filter(tracks_plot, schedule != "simulated (10 min)"), size = 0.4) +
#   # start locations
#   geom_point(data = slice(group_by(tracks_plot, schedule, ID), 1),
#              shape = 17, size = 2) +
#   scale_colour_manual(values = id_cols, name = "Simulated ID") +
#   facet_wrap(~ schedule, nrow = 1) +
#   coord_equal(expand = FALSE) +
#   labs(x = "Easting (km)", y = "Northing (km)",
#        subtitle = sprintf("replicate %d, %g-day tracks; fixes per animal: %s", ik, track_days,
#                           paste(names(n_fix), round(n_fix), sep = " = ", collapse = ", "))) +
#   theme_bw()
# # p_tracks
# 
# # save
# ggsave(out("tracks.png"), p_tracks, width = 11, height = 4); browseURL(out("tracks.png"))
