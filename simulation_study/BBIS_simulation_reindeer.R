# Simulation study: reindeer landscape
# simulation replication of reindeer case study
# simulate reindeer movement over landscape using parameters from case study, 
# using random points from empirical tracks as starting locations
# 100 simulations X 6 individuals X 180 days, at 1 minute fix rate
# fit at 3 h, 6 h, 12 h, and 24 h
# thin tracks to different sampling intervals, fit Euler, BBIS, and SSF to data 
# 
# simulate and fit on Sigma2 servers using:
# Rscript simulation_study/BBIS_simulation_reindeer.R [sim]
#   [sim] one replicate (e.g. SLURM_ARRAY_TASK_ID); omit to run all n_sim

message(sprintf("%s: %s", Sys.time(), "reindeer landscape simulation study"))
# prep workspace ---------------------------------------------------------- ####
library(here)
# custom functions
source(here("functions/utility_functions.R"))  # custom general perpose functions
sourceDir("functions")  # custom function to load all functions in folder
source(here("case_studies/reindeer/R/reindeer_data.R"))  # load_reindeer_tracks()
load_lib(mvnfast, parallel, terra, dplyr,
         Rcpp, amt, survival)
cpp_path <- here("functions/compute_lik_grad_full.cpp")
Rcpp::sourceCpp(cpp_path)

output_path <- "simulation_study/outputs/reindeer_simulation"
make_path(output_path)
# Define parameters up front ---------------------------------------------- ####
## track pars
dt          <- 1/60  # temporal resolution of simulated tracks (fraction of 1 hour)
beta_centre <- -0.5  # coefficient of squared distance to range centre (see covariates)
n_animals   <- 6     # simulated animals per replicate
track_days  <- 180    # duration of each simulated track (days)
# speed and habitat betas are taken from the reindeer fit (below)

## default estimation pars (same as case_studies/reindeer/reindeer_fit.R)
n_sim   <- 100  # number of simulations
ncores  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                 unset = as.character(parallel::detectCores() - 2)))  # cores for BBIS
M       <- 200  # number of bridges
dt_max  <- 0.5  # max hours between BBIS nodes
max_retries <- 4  # BBIS refits with fresh bridges if not converged
sim_var <- c(3, 6, 12, 24)  # sampling intervals (h); 3 h = reindeer GPS schedule
oneimpact_lib <- "/cluster/projects/nn11017k/BioDivMapping/R"  # oneimpact on Saga
if (!dir.exists(oneimpact_lib)) oneimpact_lib <- NULL

# iSSF
ssf_n_control <- 10     # random steps per used step
ssf_tol_frac  <- 1 / 6  # resampling tolerance, fraction of the interval

# reindeer covariates ----------------------------------------------------- ####
# the 0.5 km smoothed field used by the reindeer BBIS & Euler fits
covs <- rast(here("case_studies/reindeer/data/smoothed/reindeer_covariates_smooth_0.5km.tif"))
names(covs) <- names(rast(here("case_studies/reindeer/data/reindeer_covariates.tiff")))

# reindeer tracks: simulated animals start at random observed fixes
tracks <- load_reindeer_tracks(covs, lib = oneimpact_lib)

# Include squared distance to centre of the reindeer range as covariate (as in
# BBIS_simulations.R). The landscape is only 63 x 88 km and, at the reindeer
# speed, an unconfined animal drifts ~30 km per axis over a year and walks off
# the raster. With beta_centre = -0.5 the range sd is ~10 km at equilibrium.
xc <- mean(tracks$x)
yc <- mean(tracks$y)
dist2 <- ((init(covs[[1]], "x") - xc)^2 + (init(covs[[1]], "y") - yc)^2) / 100
covs  <- c(covs, setNames(dist2, "dist2_centre"))
covlist <- spatRast_to_covlist(covs)

# reindeer movement parameters -------------------------------------------- ####
# BBIS at the native 3 h, M = 200, mean over the leave-one-out folds of the
# first thinning experiment (the 3 h level keeps every fix, so it is unaffected
# by the change to fixed-clock-hour thinning)
fit_3h <- read.csv(here("case_studies/reindeer/results/thinning_betas.csv")) %>%
  filter(method == "bbis", lvl_interval == 3, lvl_M == 200, converged) %>%
  group_by(coef) %>%
  summarise(beta = mean(beta))
speed <- read.csv(here("case_studies/reindeer/results/thinning_diagnostics.csv")) %>%
  filter(method == "bbis", lvl_interval == 3, lvl_M == 200, converged) %>%
  pull(sigma) %>%
  mean()  # gamma2 (km^2 per hour)
beta <- c(setNames(fit_3h$beta, fit_3h$coef)[names(covlist)[-length(covlist)]],
          dist2_centre = beta_centre)  # covariate coefficients, in covlist order

# save true values alongside the estimates
write.csv(data.frame(par = c(paste0("beta_", names(beta)), "gammasq"),
                     true = c(beta, speed)),
          file.path(output_path, "true_pars.csv"), row.names = FALSE)

# define result data.frame column names
col_names <- c(
  "sim",
  "method",  # (euler/bbis)
  "dt",
  "Tmax",   # total hours tracked (all animals)
  "delta",  # sampling interval (h)
  "n_obs",
  "N",      # approx. nodes per step (delta / dt_max)
  "M",
  "convergence",
  "iterations",
  "time",   # compute time (s)
  paste0("beta_", names(covlist)),
  "gammasq"
)

# Sim: reindeer landscape, varying delta_t (fixed duration) ---------------- ####
print("reindeer landscape, varying delta_t, fixed duration")
args <- commandArgs(trailingOnly = TRUE)
sims <- if (length(args)) as.integer(args[1]) else seq_len(n_sim)
message(sprintf("Starting %d replicate(s): %s", length(sims), paste(sims, collapse = ", ")))
for (ik in sims) {
  set.seed(123 + ik)
  sim_results <- data.frame()  # refresh result target
  
  # start locations (and clock times) drawn from the observed reindeer fixes
  starts <- tracks[sample(nrow(tracks), n_animals), c("x", "y", "time")]
  n_obs_sim <- round(track_days * 24 / dt) + 1
  
  # simulate one track per animal
  X_full <- data.frame()
  X_fine <- data.frame()  # 10-min copy of the simulated tracks, for plotting
  message(sprintf("Simulating %d animals for replicate %d", n_animals, ik))
  for (a in seq_len(n_animals)) {
    # re-simulate the (rare) track that leaves the raster (NA gradient)
    repeat {
      X <- simLMM(dt, speed, covlist, beta, c(starts$x[a], starts$y[a]), n_obs_sim)
      if (!anyNA(X)) break
    }
    X <- data.frame(x = X[, 1], y = X[, 2],
                    time = starts$time[a] + round((seq_len(n_obs_sim) - 1) * dt * 3600),
                    ID = a)
    X_fine <- rbind(X_fine, X[seq(1, nrow(X), 10), ])  # every 10th step
    # reindeer GPS schedule: one fix every 3 h on the clock
    X_full <- rbind(X_full, thin_track(X, 3))
  }
  Tmax <- n_animals * track_days * 24
  print(nrow(X_full))
  
  for (jk in seq_along(sim_var)) {
    message(sprintf("Fitting replicate %d, delta = %g h", ik, sim_var[jk]))
    # set up simulation parameters
    delta <- sim_var[jk]
    N_sim <- delta / dt_max
    
    # thinning track: fixed clock hours, as for the reindeer
    X_thin <- thin_track(X_full, delta)
    
    # step duration = the thinned interval 
    # estimate with ssf
    ssf <- fit_ssf(X_thin, covs, resample_rate = delta,
                   tolerance = delta * ssf_tol_frac,
                   n_control = ssf_n_control, seed = ik)
    sim_results <- data.frame(ik, "ssf",   # sim & method
                              dt, Tmax,  # sim conditions
                              delta, nrow(X_thin), NA, NA,   # fit conditions
                              ssf$convergence, NA,     # convergence, iterations
                              ssf$time_secs,  # compute time
                              matrix(c(ssf$beta, NA), nrow = 1)) %>%
      setNames(col_names) %>%
      rbind(sim_results, .)
    
    # estimate with euler
    euler <- fit_euler(X_thin, covlist)
    ## extract & store euler outputs
    sim_results <- data.frame(ik, "euler",   # sim & method
                              dt, Tmax,  # sim conditions
                              delta, nrow(X_thin), NA, NA,   # fit conditions
                              1, NA,     # convergence, iterations
                              euler$time_secs,  # compute time
                              matrix(c(euler$beta, euler$gamma2), nrow = 1)) %>%
      setNames(col_names) %>%
      rbind(sim_results, .)
    
    # estimate with bbis (Euler gamma2 seeds init_sigma, as in reindeer_fit.R)
    bbis <- fit_langevin_bbis_retry(X_thin, covlist, M = M, dt_max = dt_max,
                                    dt_units = "hours", seed = 1000 * ik + delta,
                                    ncores = ncores, init_sigma = euler$gamma2,
                                    max_retries = max_retries, diagnostics = FALSE)
    
    ## extract & store bbis outputs
    sim_results <- data.frame(ik, "bbis",   # sim, method
                              dt, Tmax,  # sim conditions
                              delta, nrow(X_thin), N_sim, M,  # fit conditions
                              bbis$fit$convergence,  # convergence
                              as.numeric((bbis$fit$counts)[1]),  # iterations
                              as.numeric(bbis$fit$time, units = "secs"),  # compute time
                              matrix(bbis$fit$par, nrow = 1)) %>%   # estimates
      setNames(col_names) %>%
      rbind(sim_results, .)
  }
  write.csv(sim_results,
            file = file.path(output_path, sprintf("reindeer_landscape_estimates_sim%03d.csv", ik)),
            row.names = FALSE)
}
message(sprintf("Finished %d replicate(s): %s", length(sims), paste(sims, collapse = ", ")))
