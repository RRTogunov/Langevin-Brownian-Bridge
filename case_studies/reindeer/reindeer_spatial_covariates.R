# Reindeer case study - prep spatial covariates

# prep workspace ---------------------------------------------------------- ####
source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, terra, sf, oneimpact)

# create raster data ------------------------------------------------------ ####
# based on Panzacchi et al 2016. 
# load
load("case_studies/reindeer/data/env_raster_brick.rda")
covs <- rast(rb)
rm(rb)
# create smoothed layers of categorical cov 
foo <- as.factor(covs$newNORUT_resamp)
# replace na with cliff
foo <- ifel(is.na(foo), 3, foo)
# define level namse
levels(foo)[[1]] <- 
  data.frame(ID = 0:9, 
             newNORUT_resamp = 
               c("NA",             # LC0: NA
                 "forest",    # LC1: forest veg.	Proportion of Land Cover a class: forest vegetation
                 "bog",            # LC2: bog	Prop. of LC: bog
                 "cliff",          # LC3: cliff
                 "mountain_rock",  # LC4: mountain not edible veg.	Prop. of LC: other mountain vegetation
                 "mountain_veg",      # LC5: mountain edible veg.	Prop. of LC: mountain veg. edible for reindeer in spring
                 "fields",         # LC6: fields	Prop. of LC: fields
                 "urban",    # LC7: urban areas	Prop. of LC: urban areas
                 "lakes",  # LC8: natural lakes	Prop. of LC: non dammed lakes
                 "resevoirs"))    # LC9: reservoirs
covs$newNORUT_resamp <- foo

# create smoothed versions of categorical covs 
## retain only covs in Panzacchi et al 2016. final model
levels <- c(2, 4, 5, 8, 9)
cov_lvl <- levels(covs$newNORUT_resamp)[[1]]
## smooth 
for (i in levels) {
  out <- potential_GPU(covs$newNORUT_resamp == i, shape = "gaus", 
                       alphas = dist_to_alpha(500, shape = "gaus"))
  # add smoothed layer to covs & assign names
  covs <- c(covs, out) |> 
    setNames(c(names(covs), cov_lvl[with(cov_lvl, ID == i), 2]))
}

# other changes ----------------------------------------------------------- ####
# square slope as in Panzacchi et al 2016. 
covs$Slope_25 <- covs$Slope_25^2
# rescale solar
covs$SolarRadiat <- covs$SolarRadiat/100000

# filter covs to include in model (based on top model in Panzacchi et al 2016.)
covs <- covs |> 
  tidyterra::select(Slope_25, Ker_Trail_1k, ker1_WRoPubl, SolarRadiat,
                    bog, mountain_rock, mountain_veg, lakes,
                    resevoirs) |> 
  setNames(c("slope2", "trail_density", "road_density", "solar_radiation",
           "bog", "mountain_rock", "mountain_veg", 
           "lakes", "resevoirs"))

# # standardise covariates
# covs <- covs |> standardise_raster() 

# Change the resolution and extent from m to km for numerical stability
crs_km <- gsub("units=m", "units=km", crs(covs, proj = TRUE))
r <- rast(nrows = nrow(covs), ncols = ncol(covs),
          ext   = as.vector(ext(covs)) / 1000,
          crs   = crs_km)                 # define template raster

# reproject covs to km
covs <- project(covs, r)              # transform raster

# save -------------------------------------------------------------------- ####
writeRaster(covs, "case_studies/reindeer/data/reindeer_covariates.tiff", 
            overwrite = TRUE)

# # archived version using own covariates ----------------------------------- ####
# # 0. define paths to the raw covariates ----------------------------------- ####
# bdm <- "C:/Users/rtogunov/Research/BioDivMapping/data/temp"
# cov_paths <- c(
#   summer_precipitation = file.path(bdm, "met/summer_precipitation_150fe9d5b80b1fa5363f64b5474d3772.tiff"),
#   summer_temperature   = file.path(bdm, "met/summer_temperature_150fe9d5b80b1fa5363f64b5474d3772.tiff"),
#   elevation            = file.path(bdm, "geonorge/elevation_EPSG3045_X-454505_1512595_Y6054905_9348605.tiff"),
#   slope                = file.path(bdm, "geonorge/slope_EPSG3045_X-454505_1512595_Y6054905_9348605.tiff"),
#   density_roads        = file.path(bdm, "geonorge/density_roads_af60dbdc89d49e380c883b1b5b579cc1.tiff"),
#   forest_volume        = file.path(bdm, "SkogRover/forest_volume.tiff")
# )
# 
# # The raster whose CRS & (cropped) extent define the common target grid.
# target_cov <- "density_roads"
# target_res <- 100  # metres
# 
# # 1. define focal region (buffered extent of the reindeer data) ----------- ####
# data("reindeer")
# reindeer <- sf::st_as_sf(reindeer, coords = c("x", "y"), crs = 25833)
# reindeer_v <- terra::vect(reindeer)
# 
# e   <- ext(reindeer_v)
# xr  <- e[2] - e[1]          # x range
# yr  <- e[4] - e[3]          # y range
# buf <- as.numeric(max(xr, yr) * 0.2)   # buffer = 5% of the larger range
# e_focal <- ext(e[1] - buf, e[2] + buf, e[3] - buf, e[4] + buf)
# 
# # extent -> polygon, densified (vertex ~ every 1 km) for faithful reprojection
# focal_poly <- as.polygons(e_focal, crs = crs(reindeer_v)) |>
#   densify(interval = 1000)
# 
# # 2. Crop each covariate in their own CRS --------------------------------- ####
# # For each covariate: project the focal polygon into the covariate's CRS, take
# # its extent (a rectangle), and crop the raster to that extent.
# cropped <- lapply(cov_paths, function(p) {
#   r <- rast(p)
#   r <- crop(r, ext(project(focal_poly, r)))
#   ifel(is.na(r), 0, r)
# })
# names(cropped) <- names(cov_paths)
# 
# # 3. Project covs to common raster (crs, res, ext) ------------------------ ####
# target_crs <- crs(cropped[[target_cov]])
# template   <- rast(ext(cropped[[target_cov]]),
#                    resolution = target_res,
#                    crs        = target_crs)
# 
# # approximate a raster's ground resolution in metres (handles lon/lat sources)
# res_in_m <- function(r) {
#   if (is.lonlat(r)) {
#     yc <- mean(as.vector(ext(r))[3:4])
#     c(res(r)[1] * 111320 * cos(yc * pi / 180), res(r)[2] * 111320)
#   } else {
#     res(r)
#   }
# }
# 
# common <- lapply(names(cropped), function(nm) {
#   r      <- cropped[[nm]]
#   srcres <- min(res_in_m(r))
#   # finer output (coarse source)  -> "near"
#   # coarser output (fine source)  -> "average" (mean of sub-cells)
#   method <- if (srcres < target_res * 0.999) "average" else "near"
#   message(sprintf("Projecting %-20s (source ~%4.0f m) with method '%s'",
#                   nm, srcres, method))
#   project(r, template, method = method)
# })
# covs <- rast(common) |> setNames(names(cropped))
# 
# # mask to the focal region so standardisation uses focal-region cells only
# covs <- mask(covs, project(focal_poly, target_crs))
# 
# # 4. Standardise continuous covariates (z-score) -------------------------- ####
# sds  <- global(covs, "sd", na.rm = TRUE)[, 1]
# drop <- names(covs)[is.na(sds) | sds == 0]
# if (length(drop)) {
#   warning(sprintf("Dropping zero-variance covariate(s): %s",
#                   paste(drop, collapse = ", ")))
#   covs <- covs[[setdiff(names(covs), drop)]]
# }
# 
# covs_std <- scale(covs)   # z-score each layer (mean 0, unit variance)
# 
# # 5. Save ----------------------------------------------------------------- ####
# make_path(here("case_studies/reindeer"))
# writeRaster(covs_std,
#             here("case_studies/reindeer/reindeer_covariates.tif"),
#             overwrite = TRUE)
# writeVector(project(focal_poly, target_crs),
#             here("case_studies/reindeer/focal_region.gpkg"),
#             overwrite = TRUE)
# 
# message("Saved reindeer_covariates.tif with layers: ",
#         paste(names(covs_std), collapse = ", "))
