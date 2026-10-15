# Simulation study plots
# setup workspace --------------------------------------------------------- ####
library(here)
# load custom functions
source(here("functions/utility_functions.R"))
sourceDir(here("functions"))
# Load necessary packages
load_lib(dplyr, tidyr, stringr,     # data
         ggplot2, scales)           # plots

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
