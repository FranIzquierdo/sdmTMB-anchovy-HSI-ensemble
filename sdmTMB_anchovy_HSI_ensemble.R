#==========================================================
# sdmTMB_anchovy_HSI_ensemble.R
# Engraulis encrasicolus — Anchovy HSI
# Author: Francisco Izquierdo
# Contributor: Gabriella Lo Cicero
# Last edit: 21 May 2026
#==========================================================

# Anchovy (Engraulis encrasicolus) HSI — western Mediterranean.
# Presence-only data fitted with a Bernoulli spatiotemporal GLM (sdmTMB).
# No environmental covariates: spatial (omega_s) and spatiotemporal (epsilon_st)
# random fields capture all structure. 
# Epsilon_st is the quantity of interest, exported as a raster stack for a later
#  climate velocity analysis.
# A sensitivity analysis across PA methods, ratios, and temporal structures
# is performed over 10 seeds (120 runs) to identify the best configuration.
# Finally, best selected (PA + ratio + model structure) convergent seeds models
# are combined into an ensemble to reduce sensitivity to pseudo-absence placement.

# OUTLINE:

#   Data              — load, project (UTM 33N km), subset to 2014-2025
#   Mesh              — INLA triangulation mesh (RangeGuess = 130 km)
#   Mask              — inner marine domain for PA sampling and prediction
#   Bathymetry        — NOAA bathy download; reference shelf mask (450 m)
#   Pseudo-absences   — Buffer and RWPAS methods, 3 ratios, 10 seeds
#   Model selection   — 120 sdmTMB fits; stability diagnostics; best config
#   Ensemble setup    — select seeds_ok; define paths and shared objects
#   Model fitting     — refit ensemble (anisotropy = TRUE); cache per seed
#   Anisotropy AIC    — confirm anisotropy improvement across all 7 seeds
#   DHARMa            — residual diagnostics per sub-model (KS, spatial AC)
#   Individual CV     — spatial blockCV per sub-model (k=5, 200 km); flag weak models
#   Ensemble CV       — average fold predictions → ensemble AUC, Boyce
#   Ensemble preds    — full-domain predictions for all years; HSI, omega, epsilon
#   Parameter summary — extract and tabulate sdmTMB parameters across seeds
#   Epsilon export    — epsilon_st raster stacks (UTM + WGS84, .tif + .zip)

# OUTPUT FOLDERS:

#   input/dataset/                              raw data + df_sub_anchovy.rds
#   input/exploratory/                          coverage maps, bathymetry plots
#   input/mesh/                                 range guess + mesh plot
#   input/masks/                                mask raster for predictions
#   input/pseudo_absences/                      PA method plots

#   output/model_selection/                     one subfolder per run + tabla_comp.csv
#   output/model_selection/pseudo_absences/buffer/
#   output/model_selection/pseudo_absences/rwpas/
#   output/model_selection/stability/           Test 1/2/3 plots + summary
#   output/ensemble/models/                     fitted models (anisotropy = TRUE)
#   output/ensemble/cv/                         CV objects, metrics, ROC, Boyce
#   output/ensemble/predictions/                ensemble HSI, omega, epsilon
#   output/ensemble/parameters/                 parameter table + plot
#   output/ensemble/dharma/                     DHARMa QQ plots + summary
#   output/ensemble/epsilon_export/             epsilon raster stacks (UTM + WGS84 zips)

# Start here -------------------------------------------------------------------

## Packages --------------------------------------------------------------------

library(sp)
library(sf)
library(INLA)
library(dplyr)
library(marmap)
library(terra)
library(visreg)
library(Hmisc)
library(sdmTMB)
library(ggplot2)
library(DHARMa)
library(corrplot)
library(sdmTMBextra)
library(rnaturalearth)
library(openxlsx)
library(flexsdm)
library(pROC)
library(future)
library(patchwork)
library(GGally)
library(tidyterra)
library(ggridges)
library(blockCV)

## Directories -----------------------------------------------------------------

path_data        <- "input/dataset"
path_exploratory <- "input/exploratory"
path_mesh        <- "input/mesh"
path_masks       <- "input/masks"
path_pa          <- "input/pseudo_absences"

path_out         <- "output"
path_model_sel   <- "output/model_selection"
path_pa_buffer   <- "output/model_selection/pseudo_absences/buffer"
path_pa_rwpas    <- "output/model_selection/pseudo_absences/rwpas"
path_stab        <- "output/model_selection/stability"
path_best        <- "output/best_model"

dirs_to_create <- c(
  path_data, path_exploratory, path_mesh, path_masks, path_pa,
  path_out, path_pa_buffer, path_pa_rwpas,
  path_model_sel, path_stab
)
lapply(dirs_to_create, function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE))

## CRS -------------------------------------------------------------------------

# All model coordinates in UTM (km) Zone 33N 
# Equal-distance projection keeps range/mesh/PA distance parameters in km

crs_km    <- "+proj=utm +zone=33 +units=km +ellps=WGS84 +datum=WGS84 +no_defs"
crs_m_str <- "EPSG:32633"
CRS.new   <- CRS("+proj=utm +zone=33 +units=km +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")


# Data -------------------------------------------------------------------------

## Load raw data ---------------------------------------------------------------

# ./input/dataset/Occ_EE_IAS.xlsx contains all anchovy (Engraulis encrasicolus) survey records

df_full_anchovy      <- read.xlsx(file.path(path_data, "Occ_EE_IAS.xlsx"))
df_full_anchovy$Year <- as.integer(sub(".*?(\\d{4})$", "\\1", df_full_anchovy$Campagna))
table(df_full_anchovy$Year)

coast <- ne_countries(scale = "medium", returnclass = "sf")

cat("=== Unique surveys ===\n")
print(table(df_full_anchovy$Campagna))

df_full_anchovy$Survey <- trimws(gsub("\\d", "", df_full_anchovy$Campagna))
cat("\n=== Detected surveys ===\n")
print(table(df_full_anchovy$Survey))

## Map presences by year -------------------------------------------------------

p_by_year <- ggplot() +
  geom_sf(data = coast, fill = "grey80", color = "grey60", linewidth = 0.2) +
  geom_point(data = df_full_anchovy,
             aes(x = Lon_Media, y = Lat_Media, color = Survey),
             size = 1.2, alpha = 0.8) +
  facet_wrap(~Year, nrow = 3) +
  coord_sf(xlim = range(df_full_anchovy$Lon_Media) + c(-0.5, 0.5),
           ylim = range(df_full_anchovy$Lat_Media) + c(-0.5, 0.5)) +
  scale_color_manual(values = c("#e0ab5b", "#225ea8"), name = "Survey") +
  labs(title    = "Spatial coverage by year",
       subtitle = "Engraulis encrasicolus — each panel = one sampling year",
       x = "Longitude", y = "Latitude") +
  theme_light() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

print(p_by_year)
ggsave(file.path(path_exploratory, "df_full_coverage_by_year.png"),
       p_by_year, width = 14, height = 10, dpi = 900)

## Map latitudinal distribution by year ----------------------------------------

p_lat <- ggplot(df_full_anchovy, aes(x = factor(Year), y = Lat_Media, fill = Survey)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  labs(title = "Latitudinal distribution by year", x = "Year", y = "Latitude (N)") +
  scale_fill_manual(values = c("#e0ab5b", "#225ea8"), name = "Survey") +
  theme_light() + theme(legend.position = "bottom")

print(p_lat)
ggsave(file.path(path_exploratory, "coverage_latitude_by_year.png"),
       p_lat, width = 10, height = 5, dpi = 900)

## Subset to modelling years ---------------------------------------------------

# Keep only years with comparable sampling coverage (2014-2025)
# Earlier years have heterogeneous effort and are excluded

years_sel      <- 2014:2025
df_sub_anchovy <- subset(df_full_anchovy, Year %in% years_sel)
df_sub_anchovy$pr_ab <- 1

# Faceted map — one panel per modelling year
p_sub_year <- ggplot() +
  geom_sf(data = coast, fill = "grey80", color = "grey60", linewidth = 0.2) +
  geom_point(data = df_sub_anchovy,
             aes(x = Lon_Media, y = Lat_Media, color = Survey),
             size = 1.5, alpha = 0.9) +
  facet_wrap(~Year, nrow = 3) +
  coord_sf(xlim = range(df_sub_anchovy$Lon_Media) + c(-0.5, 0.5),
           ylim = range(df_sub_anchovy$Lat_Media) + c(-0.5, 0.5)) +
  scale_color_manual(values = c("#e0ab5b", "#225ea8"), name = "Survey") +
  labs(title    = "Modelling dataset — spatial coverage by year",
       subtitle = paste0("Engraulis encrasicolus | ", min(years_sel), "–", max(years_sel),
                         " | n = ", nrow(df_sub_anchovy), " presences"),
       x = "Longitude", y = "Latitude") +
  theme_light() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

print(p_sub_year)
ggsave(file.path(path_exploratory, "subset_coverage_by_year.png"),
       p_sub_year, width = 14, height = 10, dpi = 900)

# General overview map — all years pooled
p_overview <- ggplot() +
  geom_sf(data = coast, fill = "grey80", color = "grey60", linewidth = 0.3) +
  geom_point(data = df_sub_anchovy,
             aes(x = Lon_Media, y = Lat_Media, color = Survey),
             size = 1.8, alpha = 0.7) +
  coord_sf(xlim = range(df_sub_anchovy$Lon_Media) + c(-0.5, 0.5),
           ylim = range(df_sub_anchovy$Lat_Media) + c(-0.5, 0.5)) +
  scale_color_manual(values = c("#e0ab5b", "#225ea8"), name = "Survey") +
  labs(title    = "Study area — Engraulis encrasicolus",
       subtitle = paste0("Western Mediterranean | ", min(years_sel), "–", max(years_sel),
                         " | n = ", nrow(df_sub_anchovy), " presences across ",
                         length(years_sel), " years"),
       x = "Longitude", y = "Latitude") +
  theme_light() +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5))

print(p_overview)
ggsave(file.path(path_exploratory, "study_area_overview.png"),
       p_overview, width = 8, height = 8, dpi = 900)

## UTM projection (km) ---------------------------------------------------------

# Project lon/lat to UTM Zone 33N in KILOMETRES (not metres)

d <- data.frame(lon = df_sub_anchovy$Lon_Media, lat = df_sub_anchovy$Lat_Media)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
dutm      <- spTransform(d, CRS.new)
coords_km <- as.matrix(coordinates(dutm))
cat("coords_km range - X:", round(range(coords_km[, 1])), "| Y:", round(range(coords_km[, 2])), "\n")

df_sub_anchovy$X <- coords_km[, 1]
df_sub_anchovy$Y <- coords_km[, 2]
df_sub_anchovy   <- na.omit(df_sub_anchovy)

## Exploratory pairs plot ------------------------------------------------------

ppair <- ggpairs(df_sub_anchovy[, c("pr_ab", "X", "Y", "Lat_Media", "Lon_Media")],
                 upper = list(continuous = wrap("cor", size = 3)),
                 lower = list(continuous = wrap("points", size = 0.5, alpha = 0.5))) +
  theme_light()
print(ppair)
ggsave(file.path(path_exploratory, "ggpairs.png"), ppair, width = 10, height = 5, dpi = 900)

saveRDS(df_sub_anchovy, file.path(path_data, "df_sub_anchovy.rds"))


# Mesh -------------------------------------------------------------------------

## Range guess from inter-site distances ---------------------------------------

# The spatial range is the distance beyond which two locations 
# can be considered uncorrelated

# Before fitting the model, we estimate a minimal plausible range value 
# (rangue guess) following Zuur et al. (2017)

D       <- dist(coords_km)
df_dist <- data.frame(dist = as.vector(D))

p1 <- ggplot(df_dist, aes(x = dist)) +
  geom_histogram(bins = 30, fill = "#a1dab4", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 130, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  scale_x_continuous(breaks = seq(0, max(df_dist$dist), by = 100)) +
  theme_light() +
  labs(title = "Distance distribution between sites",
       x = "Distance between points (km)", y = "Frequency") +
  theme(panel.grid.minor = element_blank())

p2 <- ggplot(df_dist, aes(x = dist)) +
  stat_ecdf(geom = "step", color = "#225ea8", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, max(df_dist$dist), by = 100)) +
  theme_light() +
  labs(title = "Cumulative proportion of distances",
       x = "Distance between points (km)", y = "Cumulative proportion")

p_final <- (p1 | p2) +
  plot_annotation(title    = "Guess spatial range analysis",
                  subtitle = "Euclidean distances between sampling points (km)",
                  theme    = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
print(p_final)
ggsave(file.path(path_mesh, "guess_range_analysis.png"),
       p_final, width = 12, height = 5, dpi = 900)

# First peak ~130 km = within-survey range → used as RangeGuess
RangeGuess <- 130              # prior spatial range (km)
MaxEdge    <- RangeGuess / 3  #  max triangle edge (~2.5 triangles per range)
cutoff_val <- MaxEdge / 2       # min distance between mesh nodes

## INLA mesh (reference period 2014–2020) --------------------------------------

# Mesh built on a fixed reference period so convex hull geometry is stable
# and independent of which years are added later. All models share this mesh
cat("\nPresences per year:\n"); print(table(df_sub_anchovy$Year))

mesh_years  <- 2014:2020
df_mesh_ref <- df_sub_anchovy[df_sub_anchovy$Year %in% mesh_years, ]
coords_mesh <- as.matrix(df_mesh_ref[, c("X", "Y")])
cat("Years used for mesh:", paste(range(mesh_years), collapse = "-"),
    "(n =", nrow(df_mesh_ref), "presences)\n")

ConvHull  <- inla.nonconvex.hull(coords_mesh, convex = -0.05)
mesh_inla <- inla.mesh.2d(boundary = ConvHull,
                           max.edge = c(1, 6) * MaxEdge,
                           cutoff   = cutoff_val,
                           offset   = c(MaxEdge, 1.6 * MaxEdge))
cat("Mesh nodes:", mesh_inla$n, "\n")
saveRDS(mesh_inla, file.path(path_mesh, "mesh_inla.rds"))

plot(mesh_inla, axes = TRUE, xlab = "X (km)", ylab = "Y (km)",
     main = paste0("INLA mesh — ", mesh_inla$n, " nodes | RangeGuess = ", RangeGuess, " km"))
points(coords_km[, 1], coords_km[, 2], pch = 21, bg = "red", col = "darkred", cex = 0.3)

coast_utm_km <- ne_countries(scale = "medium", returnclass = "sf") %>% st_transform(crs_km)

png(file.path(path_mesh, "mesh_anchovy.png"), width = 8, height = 8, units = "in", res = 900)
plot(mesh_inla, axes = TRUE, xlab = "X (km)", ylab = "Y (km)",
     main = paste0("INLA mesh (", paste(range(mesh_years), collapse = "-"), ") — ",
                   mesh_inla$n, " nodes | MaxEdge = ", MaxEdge, " km"))
plot(st_geometry(coast_utm_km), add = TRUE, col = "grey66", border = "white")
points(coords_mesh[, 1], coords_mesh[, 2], pch = 21, bg = "red", col = "darkred", cex = 0.4)
dev.off()

saveRDS(mesh_inla, file.path(path_mesh, "mesh_inla.rds"))

## Spatial coverage check ------------------------------------------------------

## https://github.com/FranIzquierdo/sdm-spatial-coverage/blob/main/mesh_coverage.R

# Assigns each observation to its mesh triangle and summarises how many points
# fall per triangle. This is a visual check on whether the mesh resolution is
# appropriate for the available data

tri_idx    <- mesh_inla$graph$tv
node_xy    <- mesh_inla$loc[, 1:2]
n_tri      <- nrow(tri_idx)

tri_polys <- lapply(seq_len(n_tri), function(i) {
  pts <- node_xy[tri_idx[i, ], ]
  st_polygon(list(rbind(pts, pts[1, ])))
})
tri_sf <- st_sf(triangle_id = seq_len(n_tri),
                geometry    = st_sfc(tri_polys, crs = crs_km))

pres_sf_km   <- st_as_sf(df_sub_anchovy, coords = c("X", "Y"), crs = crs_km)
hits         <- st_intersects(tri_sf, pres_sf_km)
tri_sf$n_pts <- lengths(hits)

pts     <- tri_sf$n_pts
pts_occ <- pts[pts > 0]
n_occ   <- sum(pts > 0)
ratio   <- sum(pts) / max(n_occ, 1)
p95     <- as.numeric(quantile(pts_occ, 0.95))
pct_occ <- n_occ / n_tri

cat("\n── Mesh coverage ────────────────────────────────────────\n")
cat("  Triangles total / occupied:", n_tri, "/", n_occ,
    "(", round(pct_occ * 100, 1), "%)\n")
cat("  Ratio obs/occupied triangle:", round(ratio, 2), "\n")
cat("  P95 pts per triangle       :", round(p95, 2), "\n")
cat("─────────────────────────────────────────────────────────\n")

# Histogram + CDF + spatial map
df_hist <- data.frame(n_pts = pts_occ) %>%
  mutate(category = cut(n_pts, breaks = c(-Inf, 3, 8, Inf),
                        labels = c("Low (1-3)", "Moderate (4-8)", "High (>8)")))

p_cov_hist <- ggplot(df_hist, aes(x = n_pts, fill = category)) +
  geom_histogram(binwidth = 1, color = "white", boundary = 0.5) +
  scale_fill_manual(values = c("Low (1-3)"      = "#a1dab4",
                               "Moderate (4-8)" = "#41b6c4",
                               "High (>8)"      = "#225ea8"),
                    name = NULL) +
  scale_x_continuous(breaks = seq(0, max(pts_occ), by = 1)) +
  labs(title    = "Sampling points per triangle",
       subtitle = paste0("Occupied: ", n_occ, " / ", n_tri,
                         " (", round(pct_occ * 100, 1), "%)  |  Nodes: ", mesh_inla$n),
       x = "Points per triangle", y = "Count (triangles)") +
  theme_light() + theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_cov_cdf <- ggplot(data.frame(n_pts = pts_occ), aes(x = n_pts)) +
  stat_ecdf(geom = "step", color = "aquamarine3", linewidth = 1) +
  geom_vline(xintercept = p95, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(0, max(pts_occ), by = 1)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title    = "Cumulative distribution",
       subtitle = paste0("P95 = ", round(p95, 2), " pts  |  ratio = ", round(ratio, 2)),
       x = "Points per triangle", y = "Cumulative proportion") +
  theme_light() + theme(panel.grid.minor = element_blank())

bbox_tri <- st_bbox(tri_sf)
p_cov_map <- ggplot() +
  geom_sf(data = tri_sf, aes(fill = n_pts), color = "grey75", linewidth = 0.15) +
  geom_sf(data = coast_utm_km, fill = "grey85", color = "grey60", linewidth = 0.2) +
  geom_point(data = data.frame(X = df_sub_anchovy$X, Y = df_sub_anchovy$Y),
             aes(x = X, y = Y), color = "firebrick", size = 0.6, alpha = 0.8) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "N points") +
  coord_sf(xlim = c(bbox_tri["xmin"], bbox_tri["xmax"]),
           ylim = c(bbox_tri["ymin"], bbox_tri["ymax"]), datum = NULL) +
  labs(title    = "Spatial distribution",
       subtitle = paste0("MaxEdge = ", round(MaxEdge, 1), " km  |  ",
                         mesh_inla$n, " nodes"),
       x = "X (km)", y = "Y (km)") +
  theme_light() + theme(panel.grid = element_blank())

p_cov_final <- p_cov_hist + p_cov_cdf + p_cov_map +
  plot_layout(design = "AACC\nBBCC") +
  plot_annotation(
    title    = paste0("Mesh coverage "),
    subtitle = "Anchovy presences (2014-2020) over INLA mesh triangles",
    theme = theme(plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5)))

print(p_cov_final)
ggsave(file.path(path_mesh, "mesh_coverage.png"),
       p_cov_final, width = 16, height = 11, dpi = 900)

# Note: coverage computed over the full mesh (inner + outer buffer); the inner
# domain only is used as the prediction area (coverage there will be higher)

# Mask -------------------------------------------------------------------------

## Inner mesh marine domain ----------------------------------------------------

# Two domains are extracted from the INLA mesh after removing land:
#   outer_sf     — exterior convex hull of the survey domain
#   inner_mar_sf — accessible interior domain (our model prediction area)

# The inner mesh domain is the correct space for generating pseudo-absences: 
# points in outer mesh would be "not sampled" = not true absences

bnd_idx      <- mesh_inla$segm$bnd$idx[, 1]
bnd_coords_m <- mesh_inla$loc[bnd_idx, 1:2] * 1000
outer_sf     <- st_polygon(list(rbind(bnd_coords_m, bnd_coords_m[1, ]))) %>%
  st_sfc(crs = 32633) %>% st_sf() %>% st_make_valid()

bbox_wgs84 <- st_bbox(st_transform(outer_sf, 4326))
land_crop  <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>% st_crop(bbox_wgs84) %>%
  st_union() %>% st_make_valid() %>% st_transform(32633)

area_mar_sf <- st_difference(outer_sf, land_crop) %>% st_make_valid()
st_write(area_mar_sf, file.path(path_masks, "area_mar_dominio.gpkg"), delete_dsn = TRUE)

int_idx      <- mesh_inla$segm$int$idx[, 1]
int_coords_m <- mesh_inla$loc[int_idx, 1:2] * 1000
inner_sf     <- st_polygon(list(rbind(int_coords_m, int_coords_m[1, ]))) %>%
  st_sfc(crs = 32633) %>% st_sf() %>% st_make_valid()
inner_mar_sf <- st_difference(inner_sf, land_crop) %>% st_make_valid()

v_inner          <- terra::vect(inner_mar_sf)
mask_inner_final <- terra::rasterize(v_inner,
  terra::rast(terra::ext(v_inner), res = 2000, crs = crs_m_str), field = 1)
terra::writeRaster(mask_inner_final, file.path(path_masks, "mask_inner_mesh.tif"), overwrite = TRUE)

p_mask_pretty <- ggplot() +
  geom_spatraster(data = mask_inner_final, alpha = 0.4) +
  scale_fill_gradient(low = "#a1dab4", high = "#a1dab4", na.value = "transparent", guide = "none") +
  geom_sf(data = land_crop,   fill = "#e0e0e0", color = "grey70", linewidth = 0.1) +
  geom_sf(data = area_mar_sf, fill = NA,        color = "#225ea8", linewidth = 0.4) +
  theme_light() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
        plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5)) +
  coord_sf(crs = st_crs(crs_km), datum = st_crs(crs_km), expand = FALSE) +
  labs(title    = "Domain mask",
       subtitle = "Marine area for pseudo-absence sampling and HSI projection",
       x = "Longitude (UTM 33N)", y = "Latitude")

print(p_mask_pretty)
ggsave(file.path(path_masks, "mask_inner_mesh_pro.png"),
       p_mask_pretty, width = 8, height = 9, dpi = 900)


# Bathymetry -------------------------------------------------------------------

# NOAA data (2-min resolution). Downloaded once and saved to disk
# Used for visualisation and the -450 m shelf mask (reference only;
# predictions use the inner mesh domain)

coast_wgs84 <- ne_countries(scale = "medium", returnclass = "sf")

bathy_raw    <- getNOAA.bathy(lon1 = 5, lon2 = 18, lat1 = 34, lat2 = 47,
                              resolution = 2, keep = FALSE)
bathy_marine <- rast(as.xyz(bathy_raw), type = "xyz", crs = "EPSG:4326")
bathy_marine[bathy_marine > 0] <- NA

bathy_df <- as.data.frame(bathy_marine, xy = TRUE)
colnames(bathy_df) <- c("lon", "lat", "depth")

p_bathy_full <- ggplot() +
  geom_raster(data = bathy_df, aes(x = lon, y = lat, fill = depth)) +
  geom_sf(data = coast_wgs84, fill = "grey70", color = "white", linewidth = 0.2) +
  geom_point(data = df_sub_anchovy, aes(x = Lon_Media, y = Lat_Media),
             color = "red", size = 1.2, alpha = 0.7) +
  scale_fill_viridis_c(name = "Depth (m)", option = "mako", direction = -1) +
  coord_sf(xlim = c(5, 18), ylim = c(34, 47)) +
  labs(title = "Bathymetry (full marine) + anchovy presences",
       x = "Longitude", y = "Latitude") + theme_light()
ggsave(file.path(path_exploratory, "map_bathymetry_presences.png"), p_bathy_full, width = 8, height = 8, dpi = 900)

inner_wgs84 <- st_transform(inner_mar_sf, 4326)
bathy_450   <- bathy_marine; bathy_450[bathy_450 < -450] <- NA
poly_450    <- as.polygons(bathy_450 > -Inf) %>% st_as_sf() %>% st_union() %>% st_make_valid()
mask_450_sf <- st_intersection(poly_450, inner_wgs84) %>%
  st_difference(st_union(st_make_valid(coast_wgs84))) %>% st_make_valid()

p_bathy_450 <- ggplot() +
  geom_sf(data = mask_450_sf, fill = "#a1dab4", alpha = 0.6, color = "#a1dab4", linewidth = 0.3) +
  geom_sf(data = coast_wgs84, fill = "grey70",  color = "white",  linewidth = 0.2) +
  geom_point(data = df_sub_anchovy, aes(x = Lon_Media, y = Lat_Media),
             color = "#225ea8", size = 1.2, alpha = 0.8) +
  coord_sf(xlim = c(6.5, 16.5), ylim = c(35.5, 45.5)) +
  labs(title    = "Study area",
       subtitle = "Presences (2014-2025) | Engraulis encrasicolus | <450m depth",
       x = "Longitude", y = "Latitude") + theme_light()
print(p_bathy_450)
ggsave(file.path(path_exploratory, "study_area_450m.png"), p_bathy_450, width = 8, height = 8, dpi = 900)

mask_450_utm        <- st_transform(mask_450_sf, 32633) %>% st_make_valid()
v_450               <- terra::vect(mask_450_utm)
mask_bathy_450_rast <- terra::rasterize(v_450,
  terra::rast(terra::ext(v_450), res = 2000, crs = crs_m_str), field = 1)
terra::writeRaster(mask_bathy_450_rast, file.path(path_masks, "mask_bathy_450.tif"), overwrite = TRUE)


# Pseudo-absences (2 methods x 3 ratios x 10 seeds) ----------------------------

# Method 1 — Buffer : hard 20 km exclusion zone; PAs drawn uniformly outside.

# Method 2 — RWPAS  : PAs drawn by weighted sampling from the full domain,
#                     using a bias layer as sampling probability:
#                       dist < RangeGuess : p = p_min (0.20) — low probability near presences
#                       dist > RangeGuess : p increases sigmoidally from 0.20 → 1.0
#                     This reflects the ecological reasoning that within the
#                     species' range, a tow may still return an absence 2/10 times


## Pseudo absence methods  settings --------------------------------------------

ratios_pa <- c(1, 3, 5)          # 1:1, 1:3, 1:5 pres:abs balance
seeds_pa  <- c(101, 202, 303, 404, 505, 606, 707, 808, 909, 1010) # 10

## Shared helper objects -------------------------------------------------------

# Buffer method only:
#   buffer_km / buffer_width_m — hard exclusion zone around presences

# RWPAS method only:
#   range_matern_m     — sigmoid midpoint (= RangeGuess in metres)
#   suavizado_m        — sigmoid slope (60 km)
#   p_min              — bias-layer floor probability in [0, RangeGuess] zone
#   candidatos_sf_full — full candidate grid (all domain cells) for PA sampling

# Shared by both methods:
#   df_pres_m          — presences in metres (st_distance works in metres)
#   m_to_km            — convert PA coords from metres back to km
#   costa_km, mask_df_km, xlim_km, ylim_km — plotting helpers

# Buffer parameters
# buffer_km = 20: approximately the minimum inter-trawl distance
buffer_km      <- 20
buffer_width_m <- buffer_km * 1000

# RWPAS parameters
range_matern_m <- RangeGuess * 1000   # m — sigmoid midpoint (= RangeGuess km)
suavizado_m    <- 60 * 1000           # 60 km — sigmoid slope
p_min          <- 0.20                # bias-layer floor probability in [0, RangeGuess] zone

# Shared
df_pres_m  <- df_sub_anchovy %>% mutate(X = X * 1000, Y = Y * 1000)
m_to_km    <- function(df) df %>% mutate(X = X / 1000, Y = Y / 1000)
costa_km   <- st_transform(land_crop, crs_km)
mask_df_km <- as.data.frame(mask_inner_final, xy = TRUE) %>%
  mutate(x = x / 1000, y = y / 1000)
xlim_km    <- range(df_sub_anchovy$X) + c(-30, 30)
ylim_km    <- range(df_sub_anchovy$Y) + c(-30, 30)

# candidatos_sf_full: full candidate pool for RWPAS pseudo-absence sampling.
# Every cell of the inner marine domain raster is converted
# to an SF point, these are all the locations where a pseudo-absence could
# potentially be placed. Computed once here because the domain does not change
# across years, ratios, or seeds; recomputing it inside the loop would repeat
# the same expensive raster-to-SF conversion ~30 times for no reason.
candidatos_sf_full <- as.data.frame(mask_inner_final, xy = TRUE, na.rm = TRUE) %>%
  st_as_sf(coords = c("x", "y"), crs = 32633)
cat("Candidate grid:", nrow(candidatos_sf_full), "cells\n")

## ─────────────────────────────────────────────────────────────────────────────
## Method 1: Buffer (20 km exclusion)
## ─────────────────────────────────────────────────────────────────────────────

# For each combination of ratio and seed, pseudo-absences are sampled randomly 
# from the inner mesh domain, excluding any location within 20 km of a presence 

# One representative plot is saved per ratio (first seed only)

for (ratio_pa in ratios_pa) {
  for (seed_i in seeds_pa) {
    fname <- paste0("df_anchovy_buffer_ratio1_", ratio_pa, "_seed_", seed_i, ".rds")
    if (file.exists(file.path(path_pa_buffer, fname))) {
      cat("  Skip (exists):", fname, "\n"); next
    }
    set.seed(seed_i)
    cat("\nBUFFER | Ratio 1:", ratio_pa, "| Seed:", seed_i, "\n")

    lista_buffer <- list()
    for (a in sort(unique(df_pres_m$Year))) {
      pres_year <- df_pres_m %>% filter(Year == a)
      pa <- sample_pseudoabs(
        data   = pres_year,
        x      = "X", y = "Y",
        n      = nrow(pres_year) * ratio_pa,
        method = c("geo_const", width = buffer_width_m),
        rlayer = mask_inner_final
      )
      pa$Year  <- a
      pa$pr_ab <- 0
      lista_buffer[[as.character(a)]] <- pa
    }

    df_out <- bind_rows(
      df_sub_anchovy %>% select(X, Y, Year, pr_ab),
      bind_rows(lista_buffer) %>% m_to_km() %>% select(X, Y, Year, pr_ab)
    )

    fname <- paste0("df_anchovy_buffer_ratio1_", ratio_pa, "_seed_", seed_i, ".rds")
    saveRDS(df_out, file.path(path_pa_buffer, fname))
    cat("  Saved:", fname, "| n:", nrow(df_out), "\n")

    # Representative map — one per ratio (first seed only, avoid 15 near-identical plots)
    if (seed_i == seeds_pa[1]) {
      p_rep <- ggplot() +
        geom_raster(data = mask_df_km, aes(x = x, y = y), fill = "aliceblue") +
        geom_sf(data = costa_km, fill = "grey90", color = "grey70") +
        geom_point(data = filter(df_out, pr_ab == 0), aes(x = X, y = Y),
                   color = "#225ea8", size = 0.8, alpha = 0.4) +
        geom_point(data = filter(df_out, pr_ab == 1), aes(x = X, y = Y),
                   color = "#e0ab5b", size = 1.2) +
        facet_wrap(~Year, nrow = 3) +
        coord_sf(xlim = xlim_km, ylim = ylim_km, crs = st_crs(crs_km), datum = st_crs(crs_km)) +
        theme_light() +
        theme(strip.text = element_text(face = "bold", size = 8)) +
        labs(title    = paste0("Buffer (", buffer_km, " km) | Ratio 1:", ratio_pa, " | Seed ", seed_i),
             subtitle = "Orange = presences | Blue = pseudo-absences",
             x = "UTM 33N (km)", y = "")
      ggsave(file.path(path_pa, paste0("buffer_ratio1_", ratio_pa, "_seed_", seed_i, ".png")),
             p_rep, width = 14, height = 10, dpi = 900)
    }
  }
}

## ─────────────────────────────────────────────────────────────────────────────
## Method 2: RWPAS (weighted sampling by bias layer)
## ─────────────────────────────────────────────────────────────────────────────

# For each combination of ratio and seed, pseudo-absences are sampled from the 
# inner mesh domain using a sigmoid probability weight based on distance to the 
# nearest presence. 

# One representative plot is saved per ratio (first seed only)

# Bias layer drives sampling probability:
#   dist < RangeGuess : p = p_min (0.20) — low probability near presences
#   dist = RangeGuess : sigmoid midpoint
#   dist > RangeGuess : p increases sigmoidally → 1.0 far from presences
# Applied per year on the full inner domain grid.

for (ratio_pa in ratios_pa) {
  for (seed_i in seeds_pa) {
    fname <- paste0("df_anchovy_rwpas_ratio1_", ratio_pa, "_seed_", seed_i, ".rds")
    if (file.exists(file.path(path_pa_rwpas, fname))) {
      cat("  Skip (exists):", fname, "\n"); next
    }
    set.seed(seed_i)
    cat("\nRWPAS | Ratio 1:", ratio_pa, "| Seed:", seed_i, "\n")

    lista_rwpas <- list()
    for (y in sort(unique(df_pres_m$Year))) {
      pres_year <- df_pres_m %>% filter(Year == y)
      pres_sf   <- st_as_sf(pres_year, coords = c("X", "Y"), crs = 32633)

      # Bias layer — used as sampling probability (p_min within RangeGuess, sigmoid beyond)
      min_dist    <- st_distance(candidatos_sf_full, pres_sf) %>% apply(1, min) %>% as.numeric()
      sigmoid_far <- p_min + (1 - p_min) * (2 * plogis((min_dist - range_matern_m) / suavizado_m) - 1)
      bias_layer  <- ifelse(min_dist < range_matern_m, p_min, sigmoid_far)

      # Weighted sampling — prob proportional to bias_layer (0.2 near presences, sigmoid beyond RangeGuess)
      idx      <- sample(seq_len(nrow(candidatos_sf_full)),
                         size    = nrow(pres_year) * ratio_pa,
                         replace = FALSE,
                         prob    = bias_layer)
      coords_f <- st_coordinates(candidatos_sf_full[idx, ])

      lista_rwpas[[as.character(y)]] <- data.frame(
        X     = coords_f[, 1], Y     = coords_f[, 2],
        Year  = y,             pr_ab = 0,
        bias  = bias_layer[idx])   # stored for optional use as covariate
    }

    df_out <- bind_rows(
      df_sub_anchovy %>% select(X, Y, Year, pr_ab) %>% mutate(bias = NA_real_),
      bind_rows(lista_rwpas) %>% m_to_km() %>% select(X, Y, Year, pr_ab, bias)
    )

    fname <- paste0("df_anchovy_rwpas_ratio1_", ratio_pa, "_seed_", seed_i, ".rds")
    saveRDS(df_out, file.path(path_pa_rwpas, fname))
    cat("  Saved:", fname, "| n:", nrow(df_out), "\n")

    if (seed_i == seeds_pa[1]) {
      p_rep <- ggplot() +
        geom_raster(data = mask_df_km, aes(x = x, y = y), fill = "aliceblue") +
        geom_sf(data = costa_km, fill = "grey90", color = "grey70") +
        geom_point(data = filter(df_out, pr_ab == 0), aes(x = X, y = Y),
                   color = "#225ea8", size = 0.8, alpha = 0.4) +
        geom_point(data = filter(df_out, pr_ab == 1), aes(x = X, y = Y),
                   color = "#e0ab5b", size = 1.2) +
        facet_wrap(~Year, nrow = 3) +
        coord_sf(xlim = xlim_km, ylim = ylim_km, crs = st_crs(crs_km), datum = st_crs(crs_km)) +
        theme_light() +
        theme(strip.text = element_text(face = "bold", size = 8)) +
        labs(title    = paste0("RWPAS | Ratio 1:", ratio_pa, " | Seed ", seed_i),
             subtitle = "Orange = presences | Blue = pseudo-absences (weighted by bias layer)",
             x = "UTM 33N (km)", y = "")
      ggsave(file.path(path_pa, paste0("rwpas_ratio1_", ratio_pa, "_seed_", seed_i, ".png")),
             p_rep, width = 14, height = 10, dpi = 900)
    }
  }
}

## PA sampling probability comparison ------------------------------------------

# Buffer : p = 0 within 20 km buffer, p = 1 (uniform) beyond

# RWPAS  : continuous sigmoid from p_min (0.20, at dist=0) → 1.0, centred at RangeGuess.
#          PAs are drawn by WEIGHTED sampling (prob = bias_layer), so cells
#          closer to presences are sampled with lower probability

dist_seq     <- seq(0, 800, by = 1) * 1000
sig_far     <- p_min + (1 - p_min) * (2 * plogis((dist_seq - range_matern_m) / suavizado_m) - 1)
p_rwpas_new <- ifelse(dist_seq < range_matern_m, p_min, sig_far)

df_curves <- data.frame(
  dist_km = rep(dist_seq / 1000, 2),
  prob    = c(ifelse(dist_seq < buffer_width_m, 0, 1), p_rwpas_new),
  Method  = rep(c("Buffer (hard exclusion)", "RWPAS (bias layer)"),
                each = length(dist_seq))
)

p_pa_curves <- ggplot(df_curves, aes(x = dist_km, y = prob, colour = Method)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = RangeGuess, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  annotate("text", x = RangeGuess + 5, y = 0.08, hjust = 0,
           label = paste0("RangeGuess = ", RangeGuess, " km"), color = "grey40", size = 3.2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_colour_manual(values = c("#a1dab4", "#225ea8")) +
  labs(title    = "PA bias-layer probability vs distance from nearest presence",
       subtitle = paste0("RWPAS: flat ", p_min, " [0, ", RangeGuess,
                         " km], sigmoid 0.5→1 beyond ", RangeGuess,
                         " km (random sampling)"),
       x = "Distance to nearest presence (km)", y = "Bias layer probability",
       colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
p_pa_curves
ggsave(file.path(path_pa, "pa_sampling_probability_curves.png"),
       p_pa_curves, width = 8, height = 5, dpi = 900)


# Model selection -------------------------------------------------------------

## Selection loop --------------------------------------------------------------

## Overview --------------------------------------------------------------------

# Maximum: 2 PA methods x 3 ratios x 10 seeds x 2 models = 120 runs.

# Results are saved incrementally after every run to tabla_comp.csv
# If the script is interrupted, it resumes from where it left off (run_id check)

# ═══════════════════════════════════════════════════════════════════════════════
# USER SELECTORS — edit this block before running the loop
# ═══════════════════════════════════════════════════════════════════════════════

## Option A — trial run (quick range check: 1 method x 1 ratio x 1 seed x 1 model)
# methods_to_run <- "rwpas"
# ratios_to_run  <- 1
# seeds_to_run   <- 101
# models_to_run  <- "ar1"

## Option B — full run (all methods, ratios, seeds, models)
methods_to_run <- c("rwpas", "buffer")
ratios_to_run  <- c(1, 3, 5)
seeds_to_run   <- c(101, 202, 303, 404, 505, 606, 707, 808, 909, 1010)
models_to_run  <- c("rw", "ar1")

## Anisotropy ------------------------------------------------------------------

# Fixed to FALSE for the selection loop — faster and more stable
# Will be set to TRUE when refitting the best model structure in the final selection
run_anisotropy <- FALSE

# ═══════════════════════════════════════════════════════════════════════════════

## Model formula and family ----------------------------------------------------

# Intercept-only: no environmental covariates. The spatial (omega) and
# spatiotemporal (epsilon) random fields capture all structure. This is
# intentional: epsilon_st is the output fed into the climate velocity model

formula_base <- pr_ab ~ 1
family_base  <- binomial(link = "logit")

## Priors ----------------------------------------------------------------------

# Penalised-complexity priors on the Matérn covariance for both spatial (omega)
# and spatiotemporal (epsilon) fields

# range_gt = 50 km: P(range < 50) = 0.05 — anchovy operates at ~100 km scale.
# sigma_lt = 5: P(sigma > 5) = 0.05 — effectively uninformative on the logit
# scale; just prevents the optimiser from wandering to absurd values

# No hard constraint on sigma_E so epsilon_st remains free to express its
# full magnitude (the quantity we want to obtain)

mi_control <- sdmTMBcontrol(
  priors = sdmTMBpriors(
    matern_s  = pc_matern(range_gt = 50, sigma_lt = 3),
    matern_st = pc_matern(range_gt = 50, sigma_lt = 3)
  ),
  parallel     = 10,
  newton_loops = 1
)

dataset_paths <- list(
  buffer = path_pa_buffer,
  rwpas  = path_pa_rwpas
)[methods_to_run]

# Resume-safe: load existing results if available
tabla_comp_file <- file.path(path_model_sel, "tabla_comp.csv")
if (file.exists(tabla_comp_file)) {
  tabla_comp <- read.csv(tabla_comp_file, stringsAsFactors = FALSE)
  tabla_comp$seed <- as.character(tabla_comp$seed)
  cat("Resuming — existing rows:", nrow(tabla_comp), "\n")
} else {
  tabla_comp <- data.frame()
}

safe_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x[1]))
}

## Model selection loop --------------------------------------------------------

# Iterates over all PA datasets (method × ratio × st structure x seed), 
# fits one sdmTMB model per dataset × temporal structure, extracts parameters 
# and appends each run to tabla_comp.csv. Already-completed runs are skipped by
# run_id, safe to interrupt and resume

cat("\n================ MODEL SELECTION LOOP ================\n")
cat("Methods       :", paste(names(dataset_paths), collapse = ", "), "\n")
cat("Ratios        :", paste(as.character(ratios_pa),  collapse = ", "), "\n")
cat("Seeds         :", paste(as.character(seeds_pa),   collapse = ", "), "\n")
cat("Models        :", paste(models_to_run,  collapse = ", "), "\n")
cat("Anisotropy    :", run_anisotropy, "\n")
cat("Expected runs :", length(dataset_paths) * length(ratios_pa) *
                       length(seeds_pa) * length(models_to_run), "\n")
cat("======================================================\n")
gc()

for (metodo in names(dataset_paths)) {

  archivos_rds <- list.files(dataset_paths[[metodo]], pattern = "\\.rds$", full.names = TRUE)
  cat("\nMETHOD:", metodo, "| Datasets found:", length(archivos_rds), "\n")

  for (file_i in archivos_rds) {
    fname_i <- basename(file_i)
    ratio_i <- sub(".*ratio1_([0-9]+).*",     "\\1", fname_i)
    seed_i  <- sub(".*seed_([0-9]+)\\.rds$",  "\\1", fname_i)

    if (!(as.numeric(ratio_i) %in% ratios_to_run)) next
    if (!(as.numeric(seed_i)  %in% seeds_to_run))  next

    for (modelo in models_to_run) {
      run_id <- paste(metodo,
                      paste0("ratio1_", ratio_i),
                      paste0("seed_",  seed_i),
                      modelo, sep = "_")

      # Skip already-completed runs (resume safety)
      if (nrow(tabla_comp) > 0 && run_id %in% tabla_comp$run_id) {
        cat("  SKIP:", run_id, "\n")
        next
      }

      cat("\n>>> ", run_id, "\n")
      df_mod   <- readRDS(file_i)

      ## Link dataset INLA mesh ------------------------------------------------
      
      # make_mesh() must be called per dataset even though all share mesh_inla
      
      # It builds the projection matrix A that maps the dataset's observation
      # locations onto the fixed mesh nodes. A changes with the dataset (different
      # PA locations), but the mesh geometry (nodes, triangles) stays identical
      
      mesh_mod <- make_mesh(df_mod, c("X", "Y"), mesh = mesh_inla)

      t0    <- Sys.time()
      mod_i <- tryCatch(
        sdmTMB(
          data           = df_mod,
          formula        = formula_base,
          mesh           = mesh_mod,
          family         = family_base,
          spatial        = "on",
          time           = "Year",
          spatiotemporal = modelo,
          anisotropy     = run_anisotropy,
          silent         = TRUE,
          control        = mi_control
        ),
        error = function(e) { cat("  FAILED:", e$message, "\n"); NULL }
      )
      runtime_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)

      # Build row (failed model gets NAs but is still recorded)
      if (is.null(mod_i)) {
        fila_i <- data.frame(
          run_id = run_id, method = metodo, ratio = paste0("1:", ratio_i),
          seed = seed_i, model = modelo, runtime_min = runtime_min,
          converged = FALSE, AIC = NA_real_, sanity_ok = NA_real_,
          hessian_ok = NA_real_, se_ok = NA_real_, sigma_ok = NA_real_,
          sigma_O = NA_real_, sigma_E = NA_real_, rho = NA_real_,
          range = NA_real_
        )
      } else {
        path_run_i <- file.path(path_model_sel, run_id)
        dir.create(path_run_i, recursive = TRUE, showWarnings = FALSE)
        saveRDS(mod_i, file.path(path_run_i, "model.rds"))

        ## Extract pars. and run sanity  ---------------------------------------
        
        # tidy("ran_pars") returns: range, sigma_O, sigma_E, rho (if AR1).
        # sanity() checks: hessian invertibility, SE magnitudes, sigma bounds.

        pars_df <- tryCatch(sdmTMB::tidy(mod_i, "ran_pars"), error = function(e) NULL)
        s_i     <- tryCatch(sanity(mod_i, silent = TRUE),
                            error = function(e) list(all_ok = FALSE))

        get_par <- function(nm) {
          if (is.null(pars_df)) return(NA_real_)
          v <- pars_df$estimate[pars_df$term == nm]
          if (length(v) == 0) NA_real_ else round(as.numeric(v[1]), 4)
        }

        fila_i <- data.frame(
          run_id      = run_id,
          method      = metodo,
          ratio       = paste0("1:", ratio_i),
          seed        = seed_i,
          model       = modelo,
          runtime_min = runtime_min,
          converged   = TRUE,
          AIC         = tryCatch(round(AIC(mod_i), 2), error = function(e) NA_real_),
          sanity_ok   = safe_scalar(s_i$all_ok),
          hessian_ok  = safe_scalar(s_i$hessian_ok),
          se_ok       = safe_scalar(s_i$se_magnitude_ok),
          sigma_ok    = safe_scalar(s_i$sigmas_ok),
          sigma_O     = get_par("sigma_O"),
          sigma_E     = get_par("sigma_E"),
          rho         = get_par("rho"),
          range       = get_par("range")
        )
        cat("  OK | AIC:", fila_i$AIC, "| sanity:", s_i$all_ok,
            "| range:", fila_i$range, "| time:", runtime_min, "min\n")
      }

      tabla_comp <- dplyr::bind_rows(tabla_comp, fila_i)
      saveRDS(tabla_comp, file.path(path_model_sel, "tabla_comp.rds"))
      write.csv(tabla_comp, tabla_comp_file, row.names = FALSE)
      gc()
    }
  }
}

cat("\n================ LOOP FINISHED ================\n")
cat("Total rows in tabla_comp:", nrow(tabla_comp), "\n")
print(tabla_comp)


## Stability diagnostics -------------------------------------------------------

# Visualises convergence rate, parameter distributions, and AIC comparisons
# across the 120 runs to identify the best-performing configuration

## Load tabla_comp -------------------------------------------------------------

# Load results table from disk if not already in memory (e.g. after a fresh session).
if (!exists("tabla_comp") || nrow(tabla_comp) == 0) {
  tabla_comp <- read.csv(file.path(path_model_sel, "tabla_comp.csv"),
                         stringsAsFactors = FALSE)
}
tc <- tabla_comp %>% filter(converged == TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 1: Sanity convergence rate by method × ratio × model
## Which PA ratio produce more convergent fits?
## ─────────────────────────────────────────────────────────────────────────────

sanity_pct <- tc %>%
  group_by(method, ratio, model) %>%
  summarise(pct_sanity = mean(sanity_ok == 1, na.rm = TRUE) * 100,
            .groups = "drop") %>%
  mutate(ratio = factor(ratio, levels = c("1:1", "1:3", "1:5")),
         model = toupper(model))

p_sanity <- ggplot(sanity_pct,
                   aes(x = ratio, y = pct_sanity, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = paste0(round(pct_sanity), "%")),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2) +
  facet_wrap(~model) +
  scale_fill_manual(values = c("buffer" = "#225ea8", "rwpas" = "#a1dab4"),
                    name = "PA method") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Sanity check pass rate by PA method, ratio, and temporal model",
       x = "Presence:absence ratio", y = "% sanity OK") +
  theme_light() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
p_sanity

ggsave(file.path(path_stab, "plot1_sanity_rate.png"),
       p_sanity, width = 10, height = 5, dpi = 900)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 2: Parameter distribution — σ_ε, σ_ω, range (ratio 1:3, sanity OK)
## Presence-absence method selection
## ─────────────────────────────────────────────────────────────────────────────

tc_13 <- tc %>%
  filter(ratio == "1:3", sanity_ok == 1) %>%
  mutate(model = toupper(model))

# Compute outlier bounds per parameter
bounds <- tc_13 %>%
  summarise(
    sigE_lo  = mean(sigma_E, na.rm = TRUE) - 2 * sd(sigma_E, na.rm = TRUE),
    sigE_hi  = mean(sigma_E, na.rm = TRUE) + 2 * sd(sigma_E, na.rm = TRUE),
    sigO_lo  = mean(sigma_O, na.rm = TRUE) - 2 * sd(sigma_O, na.rm = TRUE),
    sigO_hi  = mean(sigma_O, na.rm = TRUE) + 2 * sd(sigma_O, na.rm = TRUE),
    range_lo = mean(range,   na.rm = TRUE) - 2 * sd(range,   na.rm = TRUE),
    range_hi = mean(range,   na.rm = TRUE) + 2 * sd(range,   na.rm = TRUE)
  )

tc_13 <- tc_13 %>%
  mutate(
    out_sigE  = sigma_E < bounds$sigE_lo  | sigma_E > bounds$sigE_hi,
    out_sigO  = sigma_O < bounds$sigO_lo  | sigma_O > bounds$sigO_hi,
    out_range = range   < bounds$range_lo | range   > bounds$range_hi,
    outlier   = out_sigE | out_sigO | out_range
  )

make_strip <- function(data, yvar, ylab, lo, hi) {
  ggplot(data, aes(x = method, y = .data[[yvar]], colour = outlier, shape = model)) +
    geom_hline(yintercept = c(lo, hi), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = c("FALSE" = "#225ea8", "TRUE" = "#e0ab5b"),
                        labels = c("Normal", "Outlier (±2 SD)"),
                        name = NULL) +
    scale_shape_manual(values = c("RW" = 16, "AR1" = 17), name = "Model") +
    labs(x = NULL, y = ylab) +
    theme_light() +
    theme(legend.position = "bottom")
}

p2a <- make_strip(tc_13, "sigma_E", expression(sigma[epsilon]),
                  bounds$sigE_lo, bounds$sigE_hi)
p2b <- make_strip(tc_13, "sigma_O", expression(sigma[omega]),
                  bounds$sigO_lo, bounds$sigO_hi)
p2c <- make_strip(tc_13, "range",   "Spatial range (km)",
                  bounds$range_lo, bounds$range_hi)

p_params <- p2a + p2b + p2c +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
p_params

ggsave(file.path(path_stab, "plot2_param_distribution.png"),
       p_params, width = 12, height = 5, dpi = 900)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 3a: AIC comparison — RW vs AR1, RWPAS, ratio 1:3, per seed
## AIC is only valid within the same dataset (same seed/method/ratio)
## Here seeds share the same PA draws → AIC comparisons are valid
## ─────────────────────────────────────────────────────────────────────────────

tc_aic <- tc %>%
  filter(method == "rwpas", ratio == "1:3", model %in% c("rw", "ar1"),
         sanity_ok == 1, !is.na(AIC)) %>%
  mutate(model = toupper(model)) %>%
  group_by(seed) %>%
  filter(n_distinct(model) == 2) %>%   # only seeds where both models converged
  ungroup()

p_aic <- ggplot(tc_aic, aes(x = model, y = AIC, group = seed)) +
  geom_line(colour = "grey60", linewidth = 0.5) +
  geom_point(aes(colour = model), size = 3.5) +
  scale_colour_manual(values = c("RW" = "#225ea8", "AR1" = "#a1dab4"),
                      name = "Temporal model") +
  labs(title    = "AIC per seed — RW vs AR1 (RWPAS, ratio 1:3, sanity OK)",
       subtitle = "Lines connect the same seed; lower AIC = better fit",
       x = NULL, y = "AIC") +
  theme_light() +
  theme(legend.position = "bottom")

p_aic
ggsave(file.path(path_stab, "plot3a_aic_rw_vs_ar1.png"),
       p_aic, width = 6, height = 5, dpi = 900)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 3b: AR1 structure: rho — RWPAS models at ratio 1:3,  boundary failures
## Shows what proportion of AR1 models hit the ±1 boundary
## ─────────────────────────────────────────────────────────────────────────────

tc_ar1 <- tc %>%
  filter(method == "rwpas", ratio == "1:3", model == "ar1", !is.na(rho)) %>%
  mutate(rho_status = case_when(
    abs(rho) > 0.9  ~ "boundary (|ρ| > 0.9)",
    TRUE            ~ "acceptable"
  ))

n_boundary <- sum(abs(tc_ar1$rho) > 0.9)
n_total    <- nrow(tc_ar1)

p_rho <- ggplot(tc_ar1, aes(x = factor(seed), y = rho, colour = rho_status)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =  0.9, ymax =  1,
           fill = "#e0ab5b", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1,   ymax = -0.9,
           fill = "#e0ab5b", alpha = 0.18) +
  geom_hline(yintercept = c(-0.9, 0.9), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_point(size = 4, alpha = 0.9) +
  annotate("text", x = Inf, y =  0.955,
           label = "≈ RW (ρ → +1)",
           hjust = 1.1, size = 3, colour = "grey35") +
  annotate("text", x = Inf, y = -0.955,
           label = "oscillating (ρ → −1)",
           hjust = 1.1, size = 3, colour = "grey35") +
  scale_colour_manual(
    values = c("boundary (|ρ| > 0.9)" = "#e0ab5b", "acceptable" = "#225ea8"),
    name   = NULL) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(title    = expression("AR1 autocorrelation ("*rho*") — RWPAS, ratio 1:3"),
       subtitle = paste0(n_boundary, " / ", n_total,
                         " seeds in boundary zone (|ρ| > 0.9)"),
       x = "Seed", y = expression(rho)) +
  theme_light() +
  theme(legend.position = "bottom")
p_rho

ggsave(file.path(path_stab, "plot3b_rho_distribution.png"),
       p_rho, width = 7, height = 5, dpi = 900)

## Best configuration — USER SETS THIS ----------------------------------------

# After inspecting tabla_comp.csv and the stability test plots (output/model_
# selection/stability/), set the three variables below. seeds_ok is derived
# automatically: all seeds that converged AND passed sanity() for this
# configuration are included in the ensemble.

# WHY AN ENSEMBLE?
#   Pseudo-absences are not real absences: they are randomly sampled points
#   that define the "available habitat" for the model. The specific PA set
#   affects model parameters and predictions. 

# Fitting the same model structure once per seed and averaging predictions 
# reduces this sensitivity. Note that we average PREDICTIONS, not parameters

selected_method <- "rwpas"   # "rwpas" | "buffer"
selected_ratio  <- 3         # 1 | 3 | 5
selected_model  <- "rw"      # "rw" | "ar1"

# Load objects created by the model selection section from
# disk in case the script is re-entered here from a new R session

if (!exists("tabla_comp") || nrow(tabla_comp) == 0)
  tabla_comp <- read.csv(file.path(path_model_sel, "tabla_comp.csv"),
                         stringsAsFactors = FALSE)
if (!exists("mesh_inla"))
  mesh_inla <- readRDS(file.path(path_mesh, "mesh_inla.rds"))
if (!exists("mask_inner_final"))
  mask_inner_final <- terra::rast(file.path(path_masks, "mask_inner_mesh.tif"))
if (!exists("formula_base"))
  formula_base <- pr_ab ~ 1
if (!exists("family_base"))
  family_base <- binomial(link = "logit")
if (!exists("mi_control"))
  mi_control <- sdmTMBcontrol(
    priors = sdmTMBpriors(
      matern_s  = pc_matern(range_gt = 50, sigma_lt = 3),
      matern_st = pc_matern(range_gt = 50, sigma_lt = 3)
    ),
    parallel     = 10,
    newton_loops = 1
  )

# Derive seeds_ok from tabla_comp (converged + sanity OK for this configuration)
seeds_ok <- tabla_comp |>
  filter(
    method    == selected_method,
    ratio     == paste0("1:", selected_ratio),
    model     == selected_model,
    converged == TRUE,
    sanity_ok == 1
  ) |>
  pull(seed) |>
  as.numeric() |>
  sort()

cat("\nSelected configuration:", selected_method,
    "| ratio 1:", selected_ratio, "| model:", selected_model, "\n")
cat("Seeds passing sanity (", length(seeds_ok), "):",
    paste(seeds_ok, collapse = ", "), "\n")

if (length(seeds_ok) < 3)
  warning("Fewer than 3 seeds passed sanity — inspect tabla_comp before proceeding.")

# Ensemble model --------------------------------------------------------------

## Ensemble setup --------------------------------------------------------------

# Selected_method, selected_ratio, selected_model and seeds_ok are set in the
# transition block above

# ref_seed: used to define spatial CV blocks (Step 3a) and align presence
#   predictions across models (Step 3d). Set to the first seed after sorting
#   seeds_ok (i.e. the seed with the lowest index among those that passed sanity)

# Crs_m: UTM Zone 33N in metres (EPSG:32633). blockCV requires a metric CRS
#   in metres. Different from crs_km (same projection, units = km)

# Formula_base, family_base, mi_control, crs_km are inherited from the model
# selection 

ref_seed <- seeds_ok[1]
crs_m    <- "EPSG:32633"

# path_pa_sel: PA dataset folder for the selected method.
#   Uses path_pa_rwpas / path_pa_buffer defined at the top of the script
path_pa_sel <- switch(selected_method, rwpas = path_pa_rwpas, buffer = path_pa_buffer)

## Ensemble paths --------------------------------------------------------------

path_ens    <- "output/ensemble"
path_mods   <- file.path(path_ens, "models")
path_cv     <- file.path(path_ens, "cv")
path_preds  <- file.path(path_ens, "predictions")
path_params <- file.path(path_ens, "parameters")
path_dh     <- file.path(path_ens, "dharma")
path_eps    <- file.path(path_ens, "epsilon_export")

for (p in c(path_ens, path_mods, path_cv, path_preds, path_params, path_dh, path_eps))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)

## Shared objects --------------------------------------------------------------

# mesh_inla and mask_inner_final are already in memory from the Mesh and Masks
# sections. coast is re-projected to crs_km for map plotting (overrides the
# unprojected version loaded in the Data section)

coast <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(crs_km)

# Prediction grid: built from mask_inner_final (already in memory)
# Same domain used throughout the model selection 
pred_grid_base <- as.data.frame(mask_inner_final, xy = TRUE, na.rm = TRUE) |>
  transmute(X = x / 1000, Y = y / 1000)

xlim_km <- range(pred_grid_base$X) + c(-5,  5)
ylim_km <- range(pred_grid_base$Y) + c(-5, 10)

# Reference dataset: used to define CV blocks and extract years_sel for the ensemble

# Presences are identical across all seeds; only pseudo-absences differ
ds_ref    <- paste0("df_anchovy_", selected_method, "_ratio1_", selected_ratio,
                    "_seed_", ref_seed, ".rds")
df_ref    <- readRDS(file.path(path_pa_sel, ds_ref))
years_sel <- sort(unique(df_ref$Year))

cat("Seeds:", paste(seeds_ok, collapse = ", "), "\n")
cat("Years:", paste(years_sel, collapse = ", "), "\n")

# WHY REFIT?
#   The model selection loop ran with anisotropy = FALSE to speed up the grid
#   search. Once the best configuration was identified (RWPAS, 1:3, RW),
#   anisotropy was tested separately (check_anisotropy.R) and confirmed:
#   delta AIC = -18 strongly favours anisotropy = TRUE. All ensemble models
#   are therefore refit with anisotropy = TRUE.

# CACHING:
#   Each fitted model is saved to output/ensemble/models/seed_XXX/model.rds.
#   On subsequent runs (or if the script is interrupted), models are loaded
#   from disk instead of refitted
#   seed 707 may already exist from check_anisotropy.R or check_blockcv.R.

# SANITY CHECK:
#   sdmTMB::sanity() checks for: gradient > 0.01 (convergence), NaN parameters,
#   Hessian positive definiteness, and large SE. A model that fails sanity
#   should be inspected before including it in the ensemble.

## Model fitting ---------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("Model fitting\n")
cat(strrep("=", 60), "\n")

# This loop fits one sdmTMB model per seed using the selected configuration 
# (RWPAS, 1:3, RW, anisotropy = TRUE)

# If the model already exists on disk it is loaded instead of refitted. 
# Each fitted model and its dataset are stored in `models_list` 
# and `datasets_list` for use in subsequent sections (DHARMa, CV, predictions)

datasets_list <- list()
models_list   <- list()

for (seed in seeds_ok) {

  seed_label <- as.character(seed)
  seed_dir   <- file.path(path_mods, paste0("seed_", seed))
  if (!dir.exists(seed_dir)) dir.create(seed_dir, recursive = TRUE)

  ds_fname <- paste0("df_anchovy_", selected_method, "_ratio1_", selected_ratio,
                     "_seed_", seed, ".rds")
  df <- readRDS(file.path(path_pa_sel, ds_fname))
  datasets_list[[seed_label]] <- df   # store for later use in CV and DHARMa

  mod_path <- file.path(seed_dir, "model.rds")
  if (file.exists(mod_path)) {
    mod <- readRDS(mod_path)
    cat("Loaded  seed", seed, "\n")
  } else {
    cat("Fitting seed", seed, "...\n")
    mesh_s <- make_mesh(df, c("X", "Y"), mesh = mesh_inla)
    mod <- sdmTMB(
      data           = df,
      formula        = formula_base,
      mesh           = mesh_s,
      family         = family_base,
      spatial        = "on",          # omega_s: time-invariant spatial field
      time           = "Year",
      spatiotemporal = selected_model, # epsilon_st: random walk over years
      anisotropy     = TRUE,           # directional spatial range
      silent         = TRUE,
      control        = mi_control
    )
    saveRDS(mod, mod_path)
    cat("  Saved:", mod_path, "\n")
  }

  s <- sanity(mod)
  cat("  Sanity seed", seed, ":", ifelse(s$all_ok, "OK", "FAIL"), "\n")
  models_list[[seed_label]] <- mod
}

# PURPOSE: confirm that anisotropy = TRUE improves fit consistently across ALL
# 7 seeds, not just seed 707 (which was tested in check_anisotropy.R).

# INTERPRETATION:
#   delta_AIC = AIC(aniso=TRUE) - AIC(aniso=FALSE)
#   delta_AIC < -5 for all seeds → anisotropy consistently improves fit
#   delta_AIC > 0 for any seed  → that model did not benefit; investigate

# AIC penalises model complexity (anisotropy adds 2 parameters: range_max and
# an angle)

## Anisotropy AIC check --------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("Anisotropy AIC check\n")
cat(strrep("=", 60), "\n")

aic_rows <- lapply(seeds_ok, function(seed) {
  run_id  <- paste(selected_method,
                   paste0("ratio1_", selected_ratio),
                   paste0("seed_",   seed),
                   selected_model, sep = "_")
  mod_no  <- readRDS(file.path(path_model_sel, run_id, "model.rds"))
  mod_yes <- models_list[[as.character(seed)]]
  data.frame(
    seed         = seed,
    AIC_no_aniso = round(AIC(mod_no),  2),
    AIC_aniso    = round(AIC(mod_yes), 2),
    delta_AIC    = round(AIC(mod_yes) - AIC(mod_no), 2)
  )
})

aic_tab <- do.call(rbind, aic_rows)
cat("\n")
print(aic_tab)
cat("\nMean delta AIC:", round(mean(aic_tab$delta_AIC), 2), "\n")
cat("All seeds improve with anisotropy:",
    all(aic_tab$delta_AIC < -2), "\n")

write.csv(aic_tab, file.path(path_ens, "aic_comparison.csv"), row.names = FALSE)
cat("Saved: aic_comparison.csv\n")

# Runs on each sub-model independently. Two diagnostics:
#   KS uniformity (testUniformity): scaled residuals should be Uniform[0,1].
#     p > 0.05 → no evidence of misfit

#   Spatial autocorrelation (testSpatialAutocorrelation): Moran's I on residuals
#     from the middle survey year (one year only to avoid duplicate coordinates)
#     p > 0.05 → spatial field has captured the residual structure

## DHARMa diagnostics ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("DHARMa diagnostics\n")
cat(strrep("=", 60), "\n")

dharma_summary <- list()

for (seed in seeds_ok) {
  seed_label <- as.character(seed)
  cat("DHARMa: seed", seed, "...\n")

  mod    <- models_list[[seed_label]]
  df     <- datasets_list[[seed_label]]
  dh_dir <- file.path(path_dh, paste0("seed_", seed))
  if (!dir.exists(dh_dir)) dir.create(dh_dir, recursive = TRUE)

  set.seed(1)
  sim <- tryCatch(
    simulate(mod, nsim = 500, type = "mle-mvn"),
    error = function(e) { cat("  simulate failed:", conditionMessage(e), "\n"); NULL }
  )

  if (is.null(sim)) {
    dharma_summary[[seed_label]] <- data.frame(seed = seed, status = "FAILED")
    next
  }

  fitted_prob <- predict(mod, newdata = df)$est |> plogis()

  dh_obj <- DHARMa::createDHARMa(
    simulatedResponse       = sim,
    observedResponse        = df$pr_ab,
    fittedPredictedResponse = fitted_prob
  )

  png(file.path(dh_dir, "dharma_qq_resid.png"), width = 1000, height = 500)
  plot(dh_obj, main = paste0("DHARMa — seed ", seed))
  dev.off()

  png(file.path(dh_dir, "dharma_tests.png"), width = 1100, height = 700, res = 120)
  DHARMa::testResiduals(dh_obj)
  dev.off()

  yr_mid    <- years_sel[ceiling(length(years_sel) / 2)]
  idx_yr    <- which(df$Year == yr_mid)
  dh_yr     <- DHARMa::recalculateResiduals(dh_obj, sel = idx_yr)
  spat_test <- tryCatch(
    DHARMa::testSpatialAutocorrelation(dh_yr,
                                       x = df$X[idx_yr], y = df$Y[idx_yr],
                                       plot = FALSE),
    error = function(e) NULL
  )

  dharma_summary[[seed_label]] <- data.frame(
    seed         = seed,
    status       = "OK",
    KS_p         = round(DHARMa::testUniformity(dh_obj, plot = FALSE)$p.value, 4),
    dispersion_p = round(DHARMa::testDispersion(dh_obj, plot = FALSE)$p.value, 4),
    spatial_p    = if (!is.null(spat_test)) round(spat_test$p.value, 4) else NA
  )

  cat("  KS p =",         dharma_summary[[seed_label]]$KS_p,
      "| Dispersion p =", dharma_summary[[seed_label]]$dispersion_p,
      "| Spatial p =",    dharma_summary[[seed_label]]$spatial_p, "\n")
}

dh_tab <- do.call(rbind, dharma_summary)
cat("\n=== DHARMa summary ===\n")
print(dh_tab)
write.csv(dh_tab, file.path(path_dh, "dharma_summary.csv"), row.names = FALSE)

## Independent CV --------------------------------------------------------------

# Spatial block CV for all sub-models

cat("\n", strrep("=", 60), "\n")
cat("Individual CV\n")
cat(strrep("=", 60), "\n")
cat("Spatial block CV (k = 5, 200 km) per sub-model independently.\n",
    "Ensemble CV reuses these fold predictions — no additional fits needed.\n\n")

# Skip both Individual CV and Ensemble CV if output files already exist
# They share in-memory fold predictions so must run or be skipped together
# Delete either CSV to force a rerun of both
path_icv_out <- file.path(path_cv, "individual_cv_metrics.csv")
path_ecv_out <- file.path(path_cv, "ensemble_cv_metrics.csv")

if (file.exists(path_icv_out) && file.exists(path_ecv_out)) {
  cat("CV already computed — skipping Individual CV and Ensemble CV.\n",
      "Delete", path_icv_out, "or", path_ecv_out, "to rerun.\n")
} else {

dir.create(path_cv, showWarnings = FALSE, recursive = TRUE)

# ── Block design ──────────────────────────────────────────────────────────────

pres_unique <- df_ref |> dplyr::filter(pr_ab == 1) |> dplyr::distinct(X, Y)
n_pres      <- nrow(pres_unique)

pres_sf_m <- pres_unique |>
  sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(crs_km)) |>
  sf::st_transform(crs_m)

set.seed(42)
cv_obj <- blockCV::cv_spatial(
  x         = pres_sf_m,
  size      = 200000,
  k         = 5,
  selection = "random",
  iteration = 100,
  seed      = 42,
  progress  = FALSE,
  report    = FALSE
)

fold_df <- pres_unique |> dplyr::mutate(fold = cv_obj$folds_ids)
k_vals  <- sort(unique(fold_df$fold))

cat("Block design: k =", length(k_vals), "| size = 200 km | n_pres =", n_pres, "\n")
cat("Presences per fold:\n"); print(table(fold_df$fold))

# Assign fold to any (X, Y) data frame via nearest-presence lookup
# Each point inherits the fold of its spatially nearest presence, placing
# it in the same geographic block
assign_folds_nn <- function(xy_df, fold_df) {
  pts  <- xy_df  |> sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(crs_km))
  pres <- fold_df |> sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(crs_km))
  fold_df$fold[sf::st_nearest_feature(pts, pres)]
}

# Block map — hexagonal block polygons + presences + pseudo-absences
# Transform blockCV hexagons and points to WGS84 for geographic display
blocks_wgs84 <- cv_obj$blocks |> sf::st_transform(4326)

# Pool all unique PA from all 7 seeds for the background overlay
all_pa_plot <- dplyr::bind_rows(lapply(seeds_ok, function(s) {
  datasets_list[[as.character(s)]] |>
    dplyr::filter(pr_ab == 0) |>
    dplyr::select(X, Y)
})) |>
  dplyr::distinct(X, Y) |>
  sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(crs_km)) |>
  sf::st_transform(4326)

pres_plot_wgs84 <- pres_unique |>
  sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(crs_km)) |>
  sf::st_transform(4326)

coast_wgs84_cv <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

fold_pal <- c("1" = "#f8766d", "2" = "#619cff", "3" = "#00ba38",
              "4" = "#9a4be1", "5" = "#f0a30a")

# Combine points with type column for unified legend
pts_plot <- dplyr::bind_rows(
  pres_plot_wgs84  |> dplyr::mutate(type = "Presence"),
  all_pa_plot      |> dplyr::mutate(type = "Pseudo-absence")
)

p_blocks <- ggplot() +
  geom_sf(data = coast_wgs84_cv, fill = "white", color = "grey70", linewidth = 0.2) +
  geom_sf(data = blocks_wgs84,
          aes(fill = factor(folds)), alpha = 0.45, color = "white", linewidth = 0.15) +
  geom_sf(data = dplyr::filter(pts_plot, type == "Pseudo-absence"),
          aes(shape = type), colour = "grey30", size = 0.7, alpha = 0.5) +
  geom_sf(data = dplyr::filter(pts_plot, type == "Presence"),
          aes(shape = type), colour = "black", size = 1.4, alpha = 0.9) +
  scale_fill_manual(values = fold_pal, name = "Fold") +
  scale_shape_manual(values = c("Presence" = 16, "Pseudo-absence" = 1), name = NULL) +
  coord_sf(xlim = c(7.5, 17.5), ylim = c(33.5, 46.5)) +
  labs(title    = "CV spatial blocks — k = 5, size = 200 km",
       subtitle = "Same blocks applied to all 7 ensemble models") +
  guides(
    fill  = guide_legend(order = 1, override.aes = list(alpha = 0.6, size = 5)),
    shape = guide_legend(order = 2, override.aes = list(size  = 2.5))
  ) +
  theme_light() +
  theme(legend.position = "right")

ggsave(file.path(path_cv, "cv_blocks.png"), p_blocks,
       width = 8, height = 11, dpi = 300)
cat("Saved: cv_blocks.png\n")

# Save fold summary (presences + background per fold) for Quarto table
fold_summary <- data.frame(
  fold       = k_vals,
  n_pres     = as.integer(table(fold_df$fold)[as.character(k_vals)]),
  n_bg_union = sapply(pa_test_by_fold, nrow)
)
fold_summary <- rbind(
  fold_summary,
  data.frame(fold = 0L, n_pres = sum(fold_summary$n_pres),
             n_bg_union = sum(fold_summary$n_bg_union))
)
write.csv(fold_summary, file.path(path_cv, "cv_fold_summary.csv"), row.names = FALSE)

# ── Test background per fold ──────────────────────────────────────────────────

# Pool all unique PA locations from all 7 seeds, assign to folds, split by fold
# Each fold's test background = union of PA from the 7 seeds in that block,
# deduplicated. Same fold model predicts both test presences and test background

all_pa_unique <- dplyr::bind_rows(lapply(seeds_ok, function(seed) {
  datasets_list[[as.character(seed)]] |>
    dplyr::filter(pr_ab == 0) |>
    dplyr::select(X, Y)
})) |>
  dplyr::distinct(X, Y)

all_pa_unique <- all_pa_unique |>
  dplyr::mutate(fold = assign_folds_nn(all_pa_unique, fold_df))

# Use ALL unique PA per fold, no subsampling
# AUC is rank-based (ratio-invariant); Boyce benefits from more background
# points (smoother P/E curve). The 1:3 ratio is a training design choice,
# not a constraint for evaluation metrics
pa_test_by_fold <- lapply(k_vals, function(k) {
  all_pa_unique |> dplyr::filter(fold == k) |> dplyr::select(X, Y)
})
names(pa_test_by_fold) <- as.character(k_vals)

cat("Background locations per fold (union of all 7 seeds, all unique — no subsampling):\n")
cat("  Presences per fold: "); print(table(fold_df$fold))
cat("  Background per fold:\n")
print(sapply(pa_test_by_fold, nrow))

# ── Combined loop: held-out presences + held-out background ───────────────────

  # For each seed × fold:
#   1. Train on all data NOT in fold k (presences + PA)
#   2. Predict test presences × all years → average over years
#   3. Predict test background × all years → average over years
# Both predictions from the same fold model → consistent scale

cat("\nFitting", length(seeds_ok), "models ×", length(k_vals), "folds =",
    length(seeds_ok) * length(k_vals), "total fits...\n")

pres_pred_by_fold <- lapply(k_vals, function(k) {
  matrix(NA_real_, nrow = sum(fold_df$fold == k), ncol = length(seeds_ok))
})
bg_pred_by_fold <- lapply(k_vals, function(k) {
  matrix(NA_real_, nrow = nrow(pa_test_by_fold[[as.character(k)]]),
         ncol = length(seeds_ok))
})

for (si in seq_along(seeds_ok)) {
  seed <- seeds_ok[si]
  df_s <- datasets_list[[as.character(seed)]] |>
    dplyr::mutate(fold = assign_folds_nn(dplyr::select(datasets_list[[as.character(seed)]], X, Y), fold_df))

  for (ki in seq_along(k_vals)) {
    k <- k_vals[ki]

    df_train     <- df_s |> dplyr::filter(fold != k) |> dplyr::select(-fold)
    test_pres_xy <- fold_df |> dplyr::filter(fold == k) |> dplyr::select(X, Y)
    test_bg_xy   <- pa_test_by_fold[[as.character(k)]]

    if (nrow(df_train) < 10 || nrow(test_pres_xy) == 0) next

    mesh_fold <- sdmTMB::make_mesh(df_train, c("X", "Y"), mesh = mesh_inla)

    mod_k <- tryCatch(
      sdmTMB::sdmTMB(
        data           = df_train,
        formula        = pr_ab ~ 1,
        mesh           = mesh_fold,
        family         = binomial(link = "logit"),
        spatial        = "on",
        time           = "Year",
        spatiotemporal = "rw",
        anisotropy     = TRUE,
        silent         = TRUE,
        control        = sdmTMB::sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1)
      ),
      error = function(e) NULL
    )
    if (is.null(mod_k)) { cat("  seed", seed, "fold", k, ": SKIP\n"); next }

    # Predict test presences
    new_pres <- sdmTMB::replicate_df(test_pres_xy, "Year", time_values = years_sel)
    hsi_pres <- new_pres |>
      dplyr::mutate(HSI = plogis(predict(mod_k, newdata = new_pres)$est)) |>
      dplyr::group_by(X, Y) |>
      dplyr::summarise(HSI = mean(HSI), .groups = "drop")
    pres_pred_by_fold[[ki]][, si] <- dplyr::left_join(
      test_pres_xy, hsi_pres, by = c("X", "Y"))$HSI

    # Predict test background (same fold model → same scale)
    if (nrow(test_bg_xy) > 0) {
      new_bg  <- sdmTMB::replicate_df(test_bg_xy, "Year", time_values = years_sel)
      hsi_bg  <- new_bg |>
        dplyr::mutate(HSI = plogis(predict(mod_k, newdata = new_bg)$est)) |>
        dplyr::group_by(X, Y) |>
        dplyr::summarise(HSI = mean(HSI), .groups = "drop")
      bg_pred_by_fold[[ki]][, si] <- dplyr::left_join(
        test_bg_xy, hsi_bg, by = c("X", "Y"))$HSI
    }

    cat("  seed", seed, "| fold", k, "OK\n")
  }
}

# ── Individual model metrics ──────────────────────────────────────────────────

# For each seed: pool held-out predictions across all folds → AUC, Boyce
# No new fitting needed — reuses pres_pred_by_fold / bg_pred_by_fold from above

cat("\n--- Computing individual CV metrics per sub-model ---\n")

indiv_cv <- lapply(seq_along(seeds_ok), function(si) {
  seed <- seeds_ok[si]

  pres_s <- do.call(c, lapply(pres_pred_by_fold, function(m) m[, si]))
  bg_s   <- do.call(c, lapply(bg_pred_by_fold,   function(m) m[, si]))

  pres_v <- !is.na(pres_s)
  bg_v   <- !is.na(bg_s)

  if (sum(pres_v) < 5 || sum(bg_v) < 5) {
    cat("  seed", seed, ": insufficient valid predictions — SKIP\n")
    return(data.frame(seed = seed, status = "SKIP",
                      AUC = NA, Boyce = NA, flag = "SKIP"))
  }

  obs_s  <- c(rep(1L, sum(pres_v)), rep(0L, sum(bg_v)))
  pred_s <- c(pres_s[pres_v], bg_s[bg_v])

  auc_s <- tryCatch({
    round(as.numeric(pROC::auc(pROC::roc(obs_s, pred_s, quiet = TRUE))), 4)
  }, error = function(e) NA_real_)

  boyce_s <- tryCatch({
    round(ecospat::ecospat.boyce(
      fit      = c(bg_s[bg_v], pres_s[pres_v]),
      obs      = pres_s[pres_v],
      nclass   = 0, window.w = "default", res = 100, PEplot = FALSE
    )$cor, 4)
  }, error = function(e) NA_real_)

  flag <- ifelse(!is.na(auc_s) & !is.na(boyce_s) & auc_s >= 0.7 & boyce_s >= 0.5,
                 "OK", "REVIEW")

  cat("  seed", seed, "| AUC =", auc_s, "| Boyce =", boyce_s, "|", flag, "\n")

  data.frame(seed = seed, status = "OK",
             AUC = auc_s, Boyce = boyce_s, flag = flag)
})

indiv_cv_tab <- do.call(rbind, indiv_cv)
cat("\n=== Individual model CV summary ===\n")
print(indiv_cv_tab)
write.csv(indiv_cv_tab, file.path(path_cv, "individual_cv_metrics.csv"), row.names = FALSE)

indiv_cv_long <- indiv_cv_tab |>
  dplyr::filter(!is.na(AUC)) |>
  tidyr::pivot_longer(cols = c(AUC, Boyce), names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    seed   = factor(seed),
    metric = factor(metric, levels = c("AUC", "Boyce")),
    flag   = rep(indiv_cv_tab$flag[!is.na(indiv_cv_tab$AUC)], each = 2)
  )

threshold_df <- data.frame(
  metric    = factor(c("AUC", "Boyce"), levels = c("AUC", "Boyce")),
  threshold = c(0.7, 0.5)
)

p_indiv <- ggplot(indiv_cv_long, aes(x = seed, y = value, fill = flag)) +
  geom_col(width = 0.65, colour = "white") +
  geom_hline(data = threshold_df, aes(yintercept = threshold),
             linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  scale_fill_manual(values = c("OK" = "steelblue", "REVIEW" = "firebrick"),
                    name = "Status") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title    = "Individual sub-model CV performance (blockCV k = 5, 200 km)",
       subtitle = "Dashed line = minimum threshold (AUC ≥ 0.7, Boyce ≥ 0.5)",
       x = "Seed (sub-model)", y = "Value") +
  theme_light() +
  theme(strip.text = element_text(face = "bold"))
print(p_indiv)
ggsave(file.path(path_cv, "individual_cv_barplot.png"), p_indiv,
       width = 8, height = 4, dpi = 300)
cat("Saved: individual_cv_metrics.csv, individual_cv_barplot.png\n")

## Ensemble CV -----------------------------------------------------------------

  cat("\n", strrep("=", 60), "\n")
cat("Ensemble CV\n")
cat(strrep("=", 60), "\n")
cat("Reuses fold predictions from Individual CV — no new fits.\n",
    "Sub-model predictions are averaged per location before computing AUC and Boyce.\n\n")

# Average over seeds, collect across folds
pres_ens_prob <- do.call(c, lapply(pres_pred_by_fold, function(m) rowMeans(m, na.rm = TRUE)))
bg_ens_prob   <- do.call(c, lapply(bg_pred_by_fold,   function(m) rowMeans(m, na.rm = TRUE)))

pres_valid   <- !is.na(pres_ens_prob)
bg_valid     <- !is.na(bg_ens_prob)
n_pres_valid <- sum(pres_valid)
n_bg_valid   <- sum(bg_valid)

cat("\nHeld-out presence HSI : n =", n_pres_valid, "| NAs =", sum(!pres_valid), "\n")
cat("Held-out background HSI: n =", n_bg_valid,   "| NAs =", sum(!bg_valid),   "\n")

# ── AUC + Boyce ───────────────────────────────────────────────────────────────

obs_vec  <- c(rep(1L, n_pres_valid), rep(0L, n_bg_valid))
pred_vec <- c(pres_ens_prob[pres_valid], bg_ens_prob[bg_valid])
roc_obj  <- pROC::roc(obs_vec, pred_vec, quiet = TRUE)
auc_val  <- round(as.numeric(pROC::auc(roc_obj)), 4)

boyce <- ecospat::ecospat.boyce(
  fit      = c(bg_ens_prob[bg_valid], pres_ens_prob[pres_valid]),
  obs      = pres_ens_prob[pres_valid],
  nclass   = 0,
  window.w = "default",
  res      = 100,
  PEplot   = FALSE
)
boyce_val <- round(boyce$cor, 4)

cat("\n=== Ensemble CV metrics ===\n")
cat("  Models        :", length(seeds_ok), "\n")
cat("  Folds (k)     : 5\n")
cat("  Block size    : 200 km\n")
cat("  Presences     :", n_pres_valid, "\n")
cat("  Background    :", n_bg_valid, "(RWPAS test fold union)\n")
cat("  AUC           :", auc_val, "\n")
cat("  Boyce index   :", boyce_val, "\n")

cv_metrics <- data.frame(
  n_models  = length(seeds_ok),
  k_folds   = 5L,
  block_km  = 200,
  n_pres    = n_pres_valid,
  n_bg      = n_bg_valid,
  bg_method = "rwpas_pool7seeds_allunique",
  AUC       = auc_val,
  Boyce     = boyce_val
)
write.csv(cv_metrics, file.path(path_cv, "ensemble_cv_metrics.csv"), row.names = FALSE)

# ROC plot
roc_df <- data.frame(fpr = 1 - roc_obj$specificities,
                     tpr = roc_obj$sensitivities)
p_roc  <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.7, y = 0.12,
           label = paste0("AUC = ", auc_val), size = 5, colour = "steelblue",
           fontface = "bold") +
  labs(title    = "Ensemble ROC — blockCV (k = 5, 200 km)",
       subtitle = paste0(length(seeds_ok), " models | RWPAS background pooled (all unique)"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_light()
print(p_roc)
ggsave(file.path(path_cv, "cv_roc.png"), p_roc, width = 5, height = 5, dpi = 300)

# Boyce plot
png(file.path(path_cv, "cv_boyce.png"), width = 600, height = 500)
ecospat::ecospat.boyce(
  fit      = c(bg_ens_prob, pres_ens_prob[pres_valid]),
  obs      = pres_ens_prob[pres_valid],
  nclass   = 0,
  window.w = "default",
  res      = 100,
  PEplot   = TRUE
)
title(main = paste0("Boyce = ", boyce_val,
                    "  |  blockCV k=5 | ", length(seeds_ok), " models | RWPAS pooled"))
dev.off()

cat("Saved: cv_blocks.png, cv_roc.png, cv_boyce.png, ensemble_cv_metrics.csv\n")

} # end combined CV skip (sections 4 + 5)

# Each of the 7 sub-models predicts on the full prediction grid for all years.
# predict() returns: est (linear predictor), omega_s (spatial field),
# epsilon_st (spatiotemporal field), and other components.
#
# WHAT EACH FIELD MEANS:
#   est        → total linear predictor: intercept + omega_s + epsilon_st
#   plogis(est) → predicted HSI on [0,1] scale
#   omega_s    → time-invariant spatial random field: persistent habitat features
#                (e.g. bathymetry, coastal structure) not captured by predictors
#   epsilon_st → spatiotemporal random field: year-to-year deviations from
#                the mean spatial pattern (e.g. interannual environmental variability)

# ENSEMBLE STATISTICS (per grid cell × year):
#   HSI        → mean plogis(est) across 7 models: ensemble HSI map
#   HSI_sd     → SD of plogis(est) across 7 models: inter-model uncertainty
#                (high SD = models disagree on suitability at that cell)
#   omega_s    → mean omega_s: ensemble spatial random field
#   epsilon_st → mean epsilon_st: ensemble spatiotemporal field

# These fields are used downstream for climate velocity analysis:
#   omega_s + epsilon_st decompose the HSI into its stable vs dynamic components.

## Ensemble predictions --------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("Ensemble predictions\n")
cat(strrep("=", 60), "\n")

pred_grid_years <- sdmTMB::replicate_df(pred_grid_base, "Year",
                                         time_values = years_sel)

cat("Predicting on full grid (", nrow(pred_grid_base), "cells ×",
    length(years_sel), "years)...\n")

preds_all <- lapply(seeds_ok, function(seed) {
  cat("  seed", seed, "\n")
  p       <- predict(models_list[[as.character(seed)]], newdata = pred_grid_years)
  p$seed  <- seed
  p
}) |> bind_rows()

# Ensemble statistics per cell × year
preds_ens_yr <- preds_all |>
  group_by(X, Y, Year) |>
  summarise(
    HSI        = mean(plogis(est)),   # ensemble mean HSI
    HSI_sd     = sd(plogis(est)),     # inter-model uncertainty
    omega_s    = mean(omega_s),       # ensemble spatial field
    epsilon_st = mean(epsilon_st),    # ensemble spatiotemporal field
    .groups = "drop"
  )

# Time-averaged: mean over all years per cell
preds_ens_mean <- preds_ens_yr |>
  group_by(X, Y) |>
  summarise(
    HSI_mean    = mean(HSI),        # long-term mean suitability
    HSI_sd_mean = mean(HSI_sd),     # mean inter-model uncertainty
    omega_s     = mean(omega_s),    # static spatial pattern
    .groups = "drop"
  )

saveRDS(preds_ens_yr,   file.path(path_preds, "predictions_ensemble_by_year.rds"))
saveRDS(preds_ens_mean, file.path(path_preds, "predictions_ensemble_mean.rds"))
cat("Predictions saved\n")

## Maps ------------------------------------------------------------------------

ylim_km_tight <- range(pred_grid_base$Y) + c(-5, 5)

map_base <- list(
  geom_sf(data = coast, fill = "grey80", colour = "grey60", linewidth = 0.2),
  coord_sf(xlim = xlim_km, ylim = ylim_km_tight,
           crs = st_crs(crs_km), datum = st_crs(crs_km)),
  theme_light(),
  theme(legend.position = "right",
        plot.margin = margin(2, 2, 2, 2, "pt"))
)

# Mean HSI: long-term average suitability across the study period
p_hsi_mean <- ggplot() +
  geom_tile(data = preds_ens_mean, aes(x = X, y = Y, fill = HSI_mean)) +
  map_base +
  scale_fill_steps2(name = "HSI", low = "#deebf7", mid = "white", high = "#e0ab5b",
                    midpoint = 0.5, n.breaks = 7, limits = c(0, 1), na.value = NA) +
  labs(title    = "Ensemble mean HSI",
       subtitle = paste0(min(years_sel), "–", max(years_sel),
                         " | ", length(seeds_ok), " models"))
print(p_hsi_mean)
ggsave(file.path(path_preds, "hsi_mean.png"), p_hsi_mean,
       width = 6, height = 8, dpi = 300)

# HSI by year: interannual variability in suitability
p_hsi_yr <- ggplot() +
  geom_tile(data = preds_ens_yr, aes(x = X, y = Y, fill = HSI)) +
  map_base +
  scale_fill_steps2(name = "HSI", low = "#deebf7", mid = "white", high = "#e0ab5b",
                    midpoint = 0.5, n.breaks = 7, limits = c(0, 1), na.value = NA) +
  facet_wrap(~Year, nrow = 3) +
  labs(title = "Ensemble HSI by year",
       x = "UTM 33N (km)", y = "") +
  theme(strip.text = element_text(face = "bold", size = 8))
print(p_hsi_yr)
ggsave(file.path(path_preds, "hsi_by_year.png"), p_hsi_yr,
       width = 10, height = 12, dpi = 300)

# Uncertainty: where do the 7 models disagree most?
# High SD = PA configuration strongly influences predictions at that location.
p_uncert <- ggplot() +
  geom_tile(data = preds_ens_mean, aes(x = X, y = Y, fill = HSI_sd_mean)) +
  map_base +
  scale_fill_steps(name = "SD",
                   low = "white", high = "#b03a00",
                   n.breaks = 7,
                   limits = c(0, max(preds_ens_mean$HSI_sd_mean, na.rm = TRUE)),
                   na.value = NA) +
  labs(title    = "Ensemble uncertainty (inter-model SD)",
       subtitle = paste0(length(seeds_ok), " models"))
print(p_uncert)
ggsave(file.path(path_preds, "hsi_uncertainty.png"), p_uncert,
       width = 6, height = 8, dpi = 300)

# Omega: time-invariant spatial random field — persistent habitat structure
# Orange = persistently above average (persistent hotspot)
# Blue   = persistently below average (persistent cold spot)
p_omega <- ggplot() +
  geom_tile(data = preds_ens_mean, aes(x = X, y = Y, fill = omega_s)) +
  map_base +
  scale_fill_steps2(low = "#225ea8", mid = "white", high = "#e0ab5b",
                    midpoint = 0, n.breaks = 7, name = expression(omega[s]),
                    na.value = NA) +
  labs(title = "Ensemble mean spatial random field (omega)")
print(p_omega)
ggsave(file.path(path_preds, "omega_mean.png"), p_omega,
       width = 6, height = 8, dpi = 300)

# Epsilon by year: year-to-year deviations from the mean spatial pattern
# Orange = above-average suitability that year; Blue = below-average
# This field captures interannual environmental variability (used in climate velocity)
p_eps <- ggplot() +
  geom_tile(data = preds_ens_yr, aes(x = X, y = Y, fill = epsilon_st)) +
  map_base +
  scale_fill_steps2(low = "#225ea8", mid = "white", high = "#e0ab5b",
                    midpoint = 0, n.breaks = 7, name = expression(epsilon[s*","*t]),
                    na.value = NA) +
  facet_wrap(~Year, nrow = 3) +
  labs(title = "Ensemble spatiotemporal random field (epsilon) by year",
       x = "UTM 33N (km)", y = "") +
  theme(strip.text = element_text(face = "bold", size = 8))
print(p_eps)
ggsave(file.path(path_preds, "epsilon_by_year.png"), p_eps,
       width = 10, height = 12, dpi = 300)

cat("All prediction maps saved to", path_preds, "\n")


# PURPOSE: assess consistency of parameter estimates across the 7 models.
#   If parameters are stable (low SD across seeds) → the ensemble is robust
#   to pseudo-absence choice. High SD in any parameter → that parameter is
#   sensitive to the specific PA locations → interpret with caution

# PARAMETERS:
#   intercept → overall mean log-odds of presence (before spatial fields)
#   range_max → maximum spatial range (km), direction of strongest autocorrelation
#   range_min → minimum spatial range (km), perpendicular direction
#   sigma_O   → SD of the spatial field (omega_s): magnitude of persistent habitat structure
#   sigma_E   → SD of the spatiotemporal field (epsilon_st): magnitude of interannual variation
#   phi       → overdispersion (NA for binomial — only applies to continuous families)

# NOTE: parameters are not averaged across models. Each model was fitted to a
# different dataset (different PA) so their likelihoods are not comparable.
# The mean ± SD row summarises variability, not a meaningful "ensemble parameter"

## Parameter summary -----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("Parameter summary\n")
cat(strrep("=", 60), "\n")

param_rows <- lapply(seeds_ok, function(seed) {
  mod    <- models_list[[as.character(seed)]]
  rp     <- sdmTMB::tidy(mod, "ran_pars")
  fp     <- sdmTMB::tidy(mod, "fixed")
  aic_v  <- round(AIC(mod), 2)

  get_est <- function(df, nm) {
    v <- df$estimate[df$term == nm]
    if (length(v) == 0) NA_real_ else round(v, 4)
  }

  data.frame(
    seed      = seed,
    intercept = get_est(fp,  "(Intercept)"),
    range_max = get_est(rp,  "range_max"),
    range_min = get_est(rp,  "range_min"),
    sigma_O   = get_est(rp,  "sigma_O"),
    sigma_E   = get_est(rp,  "sigma_E"),
    phi       = get_est(rp,  "phi"),   # will be NA for binomial
    AIC       = aic_v
  )
})

param_tab <- do.call(rbind, param_rows)

# Geometric mean range: sqrt(range_max × range_min) — summarises the effective
# isotropic range when anisotropy is present. More informative than either axis alone
param_tab <- param_tab |>
  mutate(range_mean = round(sqrt(range_max * range_min), 4))

# Summary row: mean ± SD across the 7 seeds
num_cols <- setdiff(names(param_tab), "seed")
summ_row <- as.data.frame(t(sapply(num_cols, function(col) {
  vals <- param_tab[[col]]
  paste0(round(mean(vals, na.rm = TRUE), 4),
         " ± ", round(sd(vals, na.rm = TRUE), 4))
})))
summ_row <- cbind(data.frame(seed = "mean ± SD"), summ_row)
names(summ_row) <- names(param_tab)

param_out <- rbind(
  param_tab |> mutate(across(all_of(num_cols), as.character)),
  summ_row
)

cat("\n")
print(param_out)
write.csv(param_out, file.path(path_params, "parameter_table.csv"), row.names = FALSE)
cat("Saved: parameter_table.csv\n")

# Parameter plot: one point per seed per parameter, dashed line = ensemble mean
# range_mean (geometric mean of range_max and range_min) summarises the
# anisotropic range as a single value alongside the two axes
plot_terms <- c("range_mean", "sigma_O", "sigma_E")

p_params <- param_tab |>
  tidyr::pivot_longer(all_of(plot_terms), names_to = "term", values_to = "estimate") |>
  ggplot(aes(x = factor(seed), y = estimate, colour = term)) +
  geom_point(size = 2.5) +
  geom_hline(
    data = param_tab |>
      tidyr::pivot_longer(all_of(plot_terms), names_to = "term", values_to = "estimate") |>
      group_by(term) |> summarise(m = mean(estimate), .groups = "drop"),
    aes(yintercept = m, colour = term),
    linetype = "dashed", linewidth = 0.6
  ) +
  facet_wrap(~term, scales = "free_y") +
  labs(title = "Parameter estimates across ensemble models",
       x = "Seed", y = "Estimate") +
  theme_light() +
  theme(legend.position = "none")
print(p_params)
ggsave(file.path(path_params, "parameter_plot.png"), p_params,
       width = 10, height = 5, dpi = 300)

# PURPOSE: export the ensemble spatiotemporal random field (epsilon_st) as a
#   multi-layer GeoTIFF (one band per year) in two coordinate systems:
#   UTM Zone 33N km (native model CRS) and WGS84 (degrees). These stacks are
#   the primary input for the climate velocity model

# epsilon_st captures interannual deviations from the mean spatial pattern:
#   positive → above-average suitability that year
#   negative → below-average suitability
# The temporal structure of epsilon_st drives the velocity estimates

# OUTPUTS (output/ensemble/epsilon_export/):
#   epsilon_utm.tif    — multi-band GeoTIFF, UTM 33N km, one band per year
#   epsilon_utm.zip    — compressed archive
#   epsilon_wgs84.tif  — multi-band GeoTIFF, WGS84, one band per year
#   epsilon_wgs84.zip  — compressed archive

## Epsilon export --------------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("Epsilon export\n")
cat(strrep("=", 60), "\n")

years_sorted <- sort(years_sel)

# Build one raster layer per year by rasterizing ensemble epsilon_st onto the
# mask template (same extent and resolution as mask_inner_final).
cat("Building epsilon raster stack (", length(years_sorted), "layers)...\n")

eps_layers <- lapply(years_sorted, function(yr) {
  df_yr <- preds_ens_yr |>
    filter(Year == yr) |>
    select(X, Y, epsilon_st)
  pts <- terra::vect(df_yr, geom = c("X", "Y"), crs = crs_km)
  r   <- terra::rasterize(pts, mask_inner_final, field = "epsilon_st", fun = mean)
  names(r) <- paste0("Y", yr)
  r
})

eps_stack_utm <- terra::rast(eps_layers)
cat("Stack dimensions:", dim(eps_stack_utm), "\n")

# UTM export
tif_utm <- file.path(path_eps, "epsilon_utm.tif")
terra::writeRaster(eps_stack_utm, tif_utm, overwrite = TRUE)
zip(file.path(path_eps, "epsilon_utm.zip"), files = tif_utm, flags = "-j")
cat("Saved: epsilon_utm.tif + .zip\n")

# WGS84 export
eps_stack_wgs84 <- terra::project(eps_stack_utm, "EPSG:4326")
tif_wgs84 <- file.path(path_eps, "epsilon_wgs84.tif")
terra::writeRaster(eps_stack_wgs84, tif_wgs84, overwrite = TRUE)
zip(file.path(path_eps, "epsilon_wgs84.zip"), files = tif_wgs84, flags = "-j")
cat("Saved: epsilon_wgs84.tif + .zip\n")

cat("Epsilon export complete:", length(years_sorted), "layers | years:",
    paste(range(years_sorted), collapse = "–"), "\n")
