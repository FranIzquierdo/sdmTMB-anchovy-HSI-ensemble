# sdmTMB anchovy HSI ensemble

Spatiotemporal Habitat Suitability Index (HSI) model for anchovy (*Engraulis encrasicolus*) in the western Mediterranean using **sdmTMB**. 

## Objective

Model the HSI of anchovy from presence-only acoustic survey data (2014–2025) using a Bernoulli spatiotemporal GLM with spatial (ω_s) and spatiotemporal (ε_{s,t}) random fields (and no environmental covariates).

The key output is the spatiotemporal random field ε_{s,t}, which captures interannual distributional anomalies and is used as input for a complementary climate velocity analysis.

A systematic sensitivity analysis is run across 2 pseudo-absence (PA) methods (buffer and rwpas), 3 presence:absence ratios (1:1, 1:3, 1:5), 2 spatiotemporal structures (RW and AR1) and 10 random seeds (120 model runs). 

Finally, the best PA configuration, ratio and spatiotemporal structure for 10 seeds is combined into an ensemble of 7 models (average predictions, not parameters) to reduce sensitivity to pseudo-absence placement/generation and quantify inter-model uncertainty.

## Dependencies

R packages: `sdmTMB`, `INLA`, `blockCV`, `DHARMa`, `pROC`, `ecospat`, `flexsdm`, `terra`, `sf`, `ggplot2`, `patchwork`, `rnaturalearth`, `openxlsx`

## Folder structure

```
input/
  dataset/            raw data (Occ_EE_IAS.xlsx) + df_sub_anchovy.rds
  exploratory/        spatial coverage maps, bathymetry plots
  mesh/               range guess analysis, INLA mesh plots, mesh coverage
  masks/              inner domain mask (GeoTIFF + GPKG)
  pseudo_absences/    representative PA plots per method and ratio

output/
  model_selection/
    tabla_comp.csv    results table for all 120 runs (AIC, sanity, parameters, runtime)
    stability/        diagnostic plots (sanity rate, parameter distribution, AIC/rho)
    pseudo_absences/
      buffer/         PA datasets — Buffer method
      rwpas/          PA datasets — RWPAS method

  ensemble/
    models/           fitted models with anisotropy = TRUE (one .rds per seed)
    aic_comparison.csv  AIC: anisotropy TRUE vs FALSE
    parameters/       parameter table and plot across seeds
    dharma/           DHARMa diagnostic plots and summary table
    cv/               spatial blocks map, fold summary, sub-model CV metrics, ensemble CV metrics, ROC curve, Boyce plot
    predictions/      ensemble HSI, HSI_sd, omega, epsilon (by year and mean)
    epsilon_export/   epsilon_st raster stacks in UTM and WGS84 (.tif + .zip)
```

## What this script does

1. Projects and filters presence data (UTM Zone 33N, 2014–2025)
2. Builds an INLA triangulation mesh (RangeGuess = 130 km)
3. Generates pseudo-absences: Buffer and RWPAS methods, 3 ratios, 10 seeds
4. Fits 120 sdmTMB models; selects best configuration via stability diagnostics
5. Refits the selected configuration (RWPAS, 1:3, RW) as an ensemble of 7 sub-models with anisotropy
6. Confirms anisotropy improves fit via AIC comparison across all 7 seeds
7. Runs DHARMa residual diagnostics on each sub-model
8. Evaluates each sub-model via spatial block CV (k = 5, 200 km): AUC + Boyce 
9. Evaluates the ensemble with the same spatial block CV (35 fits): AUC + Boyce
10. Generates ensemble predictions: mean HSI, inter-model uncertainty, spatial and spatiotemporal random fields
11. Exports ε_{s,t} raster stacks (UTM and WGS84) for downstream climate velocity analysis
