# reindeer simulation plots
# setup workspace --------------------------------------------------------- ####
library(here)
# load custom functions
source(here("functions/utility_functions.R"))
sourceDir(here("functions"))                             # ud_from_beta(), eval_ud(), jack_se(), ...
sourceDir(here("case_studies/reindeer/R/"))  
# Load necessary packages
load_lib(terra, dplyr, tidyr, lubridate, stringr,  # data
         ggplot2, patchwork, scales, grid)  # plots

# settings ---------------------------------------------------------------- ####
exp_name       <- "thinning_fixed_hour"   # as in reindeer_fit.R
fits_dir       <- here("case_studies/reindeer/fitted_estimates/experiments", exp_name)
score_dir      <- file.path(fits_dir, "scores")  # UD scores cached here between runs
metrics        <- c("logdens", "kstar", "logdens_kstar",  # magnitude-sensitive (main)
                    "auc", "cbi", "tss")#,                  # rank-only (supplementary)
# "UDcor")                              # ud overlap
ref_interval   <- 3                       # native GPS interval: reference for drift
anchor_legacy  <- 12                      # as in reindeer_fit.R (fits without anchor_h)
converged_only <- TRUE                    # drop non-converged fits from summaries
make_path(results_dir)
make_path(score_dir)

# custom functions
# standardise output file names
out <- function(suffix) file.path(results_dir, paste0(exp_name, "_", suffix))
# format method title 
method_title <- function(method) {
  case_when(method == "bbis"  ~ "BBIS",
            method == "euler" ~ "Euler",
            method == "ssf"   ~ "SSF")
}
# format covariate title 
coef_title <- function(coef) {
  coef <- str_replace_all(coef, "_", " ") %>% str_to_sentence()
  ifelse(str_detect(coef, "2$"),  # quadratics
         str_replace(coef, "2$", "\u00b2"), coef)
}

# load data --------------------------------------------------------------- ####
covs_raw    <- rast(here("case_studies/reindeer/data/reindeer_covariates.tiff"))
covs_smooth <- rast(here("case_studies/reindeer/data/smoothed/reindeer_covariates_smooth_0.5km.tif")) |>
  setNames(names(covs_raw))
cov_fields  <- list(`0` = covs_raw, `0.5` = covs_smooth)  # as in reindeer_fit.R
tracks <- load_reindeer_tracks(covs_raw)

# import sim result summaries --------------------------------------------- ####
results_dir  <- here("case_studies/reindeer/results")
fits         <- read.csv(out("fits.csv"))
betas        <- read.csv(out("betas.csv"))
beta_sm      <- read.csv(out("beta_summary.csv"))
beta_phase   <- read.csv(out("beta_phase.csv"))
beta_drift   <- read.csv(out("beta_drift.csv"))
beta_shift   <- read.csv(out("beta_shift.csv"))
scores       <- read.csv(out("scores.csv"))
score_sm     <- read.csv(out("score_summary.csv"))
score_rel    <- read.csv(out("score_vs_bbis.csv"))
diag_sm      <- read.csv(out("diagnostic_summary.csv"))
# import UD rasters
ud_stack     <- rast(out("mean_UDs.tif"))

# define plot specs ------------------------------------------------------- ####
dodge   <- 0.3  # horizontal spread so methods don't overlap
sub     <- "phase-averaged; mean +/- 2 jackknife SE over leave-one-out folds"
x_lvls  <- sort(unique(fits$interval))
x_labs  <- paste(x_lvls, "h")

# plot output directory
results_dir    <- here("case_studies/reindeer/figures")

# plot: Estimate vs interval -------------------------------- ####
# error bars are ± 2 jackknife SE over leave-one-out folds (phase-averaged)
# change titles
beta_sm <- beta_sm %>%
  mutate(coef = coef_title(coef), method = method_title(method)) 

# generate plot
p_betas <- ggplot(beta_sm, aes(factor(interval, x_lvls, x_labs), mean,
                               colour = method, group = method)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +   # zero line
  geom_line(position = position_dodge(width = dodge)) +
  geom_pointrange(aes(ymin = mean - 2 * se_jack, ymax = mean + 2 * se_jack),
                  position = position_dodge(width = dodge), size = 0.1) +
  scale_colour_viridis_d(end = 0.9) +
  facet_wrap(~ coef, scales = "free_y") +
  labs(x = "Fix interval", y = "Estimate", 
       colour = "Method") +
  theme_bw()
# p_betas

# save
ggsave(out("betas.png"), p_betas, width = 7, height = 5); #browseURL(out("betas.png"))

# plot: time-of-day effect ------------------------------------------------ ####
# betas per clock phase of the coarsest schedule (one fix a day, at that UTC hour)
# error bars are ± 2 jackknife SE over leave-one-out folds

# change titles
beta_phase <- beta_phase %>%
  mutate(coef = coef_title(coef), method = method_title(method)) 

p_phase <- ggplot(filter(beta_phase, interval == max(x_lvls)),
                  aes(anchor_h, mean, colour = method, group = method)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_line(position = position_dodge(width = 1)) +
  geom_pointrange(aes(ymin = mean - 2 * se_jack, ymax = mean + 2 * se_jack),
                  position = position_dodge(width = 1), size = 0.1) +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  scale_colour_viridis_d(end = 0.9) +
  facet_wrap(~ coef, scales = "free_y") +
  labs(x = "Start hour (UTC)", y = "Estimate", colour = "Method") +
  theme_bw()
p_phase

# save
ggsave(out("beta_phase.png"), p_phase, width = 7, height = 5); #browseURL(out("beta_phase.png"))

# plot: out-of-sample UD performance vs interval -------------------------- ####
# mean score +/- 2 SE over withheld animals (phase-averaged)
# mean(log(UD)) 
# (mean log(predicted UD) across locations of left-out animal); 
# mean ± 2*se calculated across 12 folds and phases: mean/se(group_by(score_fold, method, interval, metric))

pp <- filter(score_sm, metric == "logdens") %>% 
  mutate(metric = "Log density",
         method = method_title(method)) 

p_perf <- ggplot(pp, aes(factor(interval, x_lvls, x_labs), mean,
                         colour = method, group = method)) +
  geom_line(position = position_dodge(width = dodge)) +
  geom_pointrange(aes(ymin = mean - 2 * se, ymax = mean + 2 * se),
                  position = position_dodge(width = dodge), size = 0.1) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = "Fix interval", y = "Mean out-of-sample\nlog-density", colour = "Method") +
  theme_bw()
# p_perf

# save
ggsave(out("performance.png"), p_perf, width = 3.5, height = 2); #browseURL(out("performance.png"))

# rank-only metrics (blind to the size of beta), for the supplement
# mean +/- 2 SE over withheld animals (phase-averaged)
pp_rank <- filter(score_sm, metric %in% 
                    c("UDCor_BA", "UDCor_D", "UDCor_spearman", 
                      "auc", "cbi", "tss")) %>% 
  mutate(method = method_title(method),
         metric = str_to_upper(metric)) 

p_perf_rank <- ggplot(pp_rank,
                      aes(factor(interval, x_lvls, x_labs), mean,
                          colour = method, group = method)) +
  geom_line(position = position_dodge(width = dodge)) +
  geom_pointrange(aes(ymin = mean - 2 * se, ymax = mean + 2 * se),
                  position = position_dodge(width = dodge), size = 0.1) +
  scale_colour_viridis_d(end = 0.9) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Fix interval", y = "Metric value", colour = "Method") +
  theme_bw()
# p_perf_rank

# save
ggsave(out("performance_rank.png"), p_perf_rank, width = 5.5, height = 2); #browseURL(out("performance_rank.png"))

# plot: mean UD maps ------------------------------------------------------ ####

methods <- c("bbis", "euler", "ssf")
# Rasters: row 1 = 3 h, row 2 = 24 h - 3 h
r3h <- lapply(methods, function(m) ud_stack[[paste(m, "| 3 h")]])
dUD <- lapply(methods, function(m) {
  ud_stack[[paste(m, "| 24 h")]] - ud_stack[[paste(m, "| 3 h")]]
})

# Raster -> long data frame (x/y are needed by geom_raster; axes are hidden)
to_df <- function(r) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "value"
  df
}

df_row1 <- lapply(r3h, to_df)
df_row2 <- lapply(dUD, to_df)

# Shared scale bars
fill_ud <- scale_fill_gradientn(
  colours  = hcl.colors(100, "viridis"),
  limits   = c(0, 9e-6),
  breaks   = seq(0, 9e-6, by = 3e-6),
  labels   = label_scientific(digits = 1),
  name     = NULL,
  na.value = "transparent",
  guide    = guide_colourbar(barheight = unit(2, "in"), barwidth = unit(0.25, "in"))
)

fill_diff <- scale_fill_gradientn(
  colours  = hcl.colors(101, "Blue-Red 3"),
  limits   = c(-3.5e-6, 3.5e-6),
  breaks   = seq(-3e-6, 3e-6, by = 1.5e-6),
  labels   = label_scientific(digits = 1),
  name     = NULL,
  na.value = "transparent",
  guide    = guide_colourbar(barheight = unit(2, "in"), barwidth = unit(0.25, "in"))
)

# One panel: raster only, no axes, optional title
panel <- function(df, sc, title = NULL) {
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    sc +
    coord_equal(expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          plot.margin = margin(2, 1, 2, 2))
}

p_row1 <- Map(function(df, m) panel(df, fill_ud, method_title(m)), df_row1, methods)
p_row2 <- lapply(df_row2, function(df) panel(df, fill_diff))

# Each row collects its shared legend onto the right
row1 <- wrap_plots(p_row1, nrow = 1, guides = "collect") & theme(legend.position = "right")
row2 <- wrap_plots(p_row2, nrow = 1, guides = "collect") & theme(legend.position = "right")

# Rotated row titles on the left
row_lab <- function(lbl) {
  wrap_elements(full = textGrob(lbl, rot = 90, gp = gpar(fontsize = 10)))
}

top <- (row_lab(expression("3 h UD")) | row1) + plot_layout(widths = c(0.01, 1))
bot <- (row_lab(expression(Delta * "(UD): 24 h - 3 h")) | row2) + plot_layout(widths = c(0.01, 1))
fig <- top / bot

ggsave(out("dUD_panels_gg.png"), plot = fig,
       width = 6, height = 4.5, dpi = 150, bg = "white") #; browseURL(out("dUD_panels_gg.png"))


