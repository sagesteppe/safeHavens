## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "",
  eval=FALSE
)

## ----load libraries, message =F, warning = F----------------------------------
# library(safeHavens)
# library(sf) ## vector operations
# library(terra) ## raster operations
# library(geodata) ## environmental variables
# library(dplyr) ## general data handling
# library(tidyr) ## general data handling
# library(ggplot2) ## plotting
# library(patchwork) ## multiplots
# library(pROC)
# library(DHARMa)
# library(tidybayes)

## ----Download species occurrence data, echo = FALSE---------------------------
# cols = c('decimalLatitude', 'decimalLongitude', 'dateIdentified', 'species', 'acceptedScientificName', 'datasetName',
#   'coordinateUncertaintyInMeters', 'basisOfRecord', 'institutionCode', 'catalogNumber')
# 
# ## download species data using scientificName, can use keys and lookup tables for automating many taxa.
# hemi <- rgbif::occ_search(scientificName = "Helianthella microcephala")
# hemi <- hemi[['data']][,cols]  |>
#   drop_na(decimalLatitude, decimalLongitude) |> # any missing coords need dropped.
#   distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) |> # no dupes can be present
#   st_as_sf(coords = c( 'decimalLongitude', 'decimalLatitude'), crs = 4326, remove = F)

## ----prep basemap, echo = F, warning=F----------------------------------------
# western_states <- spData::us_states |> ## for making a quick basemap.
#   dplyr::filter(NAME %in%
#     c('Utah', 'Arizona', 'Colorado', 'New Mexico', 'Wyoming', 'Nevada', 'Idaho', 'California')) |>
#   dplyr::select(NAME, geometry) |>
#   st_transform(4326)
# 
# bb <- st_bbox(
#   c(
#     xmin = -114.5,
#     xmax = -106,
#     ymax = 42,
#     ymin = 33.5
#     ),
#     crs = st_crs(4326)
#     )
# 
# western_states <- st_crop(western_states, bb)
# 
# theme_navy <- function(base_size = 14, base_family = "sans") {
#   ggplot2::theme_grey(base_size = base_size, base_family = base_family) +
#   ggplot2::theme(
#     plot.title        = ggplot2::element_text(face = "bold", size = rel(1.2),
#                           hjust = 0.5, colour = "#dce3f0"),
#     text              = ggplot2::element_text(colour = "#dce3f0"),
#     panel.background  = ggplot2::element_rect(fill = "#0d1220", colour = NA),
#     plot.background   = ggplot2::element_rect(fill = "#0a0e1a", colour = NA),
#     panel.border      = ggplot2::element_rect(fill = NA, colour = NA),
#     axis.title        = ggplot2::element_text(face = "bold", size = rel(1),
#                           colour = "#dce3f0"),
#     axis.text         = ggplot2::element_text(colour = "#8a99b8"),
#     axis.line.x       = ggplot2::element_line(colour = "#3a4a6b"),
#     axis.line.y       = ggplot2::element_line(colour = "#3a4a6b"),
#     axis.ticks        = ggplot2::element_line(colour = "#3a4a6b"),
#     panel.grid.major  = ggplot2::element_line(colour = "#1c2540"),
#     panel.grid.minor  = ggplot2::element_blank(),
#     legend.background = ggplot2::element_rect(fill = "#0a0e1a", colour = NA),
#     legend.key        = ggplot2::element_rect(fill = "#0a0e1a", colour = NA),
#     legend.position   = "bottom",
#     legend.direction  = "horizontal",
#     strip.background  = ggplot2::element_rect(fill = "#1c2540", colour = "#1c2540"),
#     strip.text        = ggplot2::element_text(face = "bold", colour = "#dce3f0")
#   )
# }
# 
# bmap <- ggplot() +
#     geom_sf(data = western_states, fill = NA, lwd = 1.25) +
#     geom_sf(data = hemi, col = 'white') +
#     theme_minimal() +
#     coord_sf(
#         xlim = c(bb[['xmin']], bb[['xmax']]),
#         ylim = c(bb[['ymin']], bb[['ymax']])
#         )  +
#     theme_navy() +
#     theme(
#       axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank()
#       )
# 
# bmap +
#     labs(title = 'Helianthella microcephala\noccurrence records')

## ----remove cols and western states, echo = F---------------------------------
# rm(cols, western_states)

## ----download data, echo = FALSE----------------------------------------------
# # Download WorldClim bioclim at ~10 km
# bio_current <- worldclim_global(var="bioc", res=2.5)
# 
# bio_future <- cmip6_world(
#   model = "CNRM-CM6-1", ## modelling method
#   ssp   = "245", ## "Middle of the Road" scenario
#   time  = "2041-2060", # time period
#   var   = "bioc", # just use the bioclim variables
#   res   = 2.5
# )
# 
# # Crop to domain
# bbox <- ext(bb)
# 
# bio_current <- crop(bio_current, bbox)
# bio_future <- crop(bio_future, bbox)

## ----standardize raster layer names, echo = FALSE-----------------------------
# simplify_names <- function(x){
#     paste0('bio_', sprintf("%02d", as.numeric(gsub('\\D+','', names(x)))))
# }
# 
# names(bio_current) <- gsub('^.*5m_', '', names(bio_current))
# names(bio_future) <- gsub('^.*2060_', '', names(bio_future))
# 
# names(bio_current) <- simplify_names(bio_current)
# names(bio_future) <- simplify_names(bio_future)
# 
# hemi <- select(hemi)
# 
# rm(simplify_names)

## ----Create Species Distribution Model, message=F, warning=F, eval=F----------
# sdModel <- bayesianSDM(
#   x = hemi,
#   pca_predictors = FALSE,
#   predictors = bio_current,
#   quantile_v = 0.025,
#   planar_proj = 5070,
#   resample = TRUE,
#   feature_selection = 'ffs',
#   min_ffs_var = 3,
#   fact = 3,
#   silent = 2 # print fewer messages during model fit.
#   )

## ----save model, eval = F, echo = F-------------------------------------------
# setwd('~/Documents/assoRted/safeHavens/vignettes')
# sdModel$RasterPredictions <- terra::wrap(sdModel$RasterPredictions)
# sdModel$Predictors <- NULL ## need to save some memory
# sdModel$RasterPredictions_sd <- terra::wrap(sdModel$RasterPredictions_sd)
# sdModel$AOA <- terra::wrap(sdModel$AOA)
# saveRDS(sdModel, file.path('..', 'inst', 'extdata', 'BayesSDM.rds'))

## ----load pre-computed sdmodel------------------------------------------------
# sdModel <- readRDS(
#   file.path(system.file(package="safeHavens"), 'extdata',  'BayesSDM.rds')
#   )
# 
# coefs <- rownames(as.data.frame(brms::fixef(sdModel$Model, summary = TRUE)))
# sdModel$Predictors <- bio_current[[ coefs[grep('bio_', coefs)]]]
# 
# # add the above back on, had to save memory
# sdModel$RasterPredictions <- terra::unwrap(sdModel$RasterPredictions)
# sdModel$RasterPredictions_sd <- terra::unwrap(sdModel$RasterPredictions_sd)
# sdModel$AOA <- terra::unwrap(sdModel$AOA)

## ----Explore SDM output - Different alpha-------------------------------------
# sdm_td <- sdModel$TrainData
# 
# par(mfrow = c(1, 2))
# plot(sdModel$RasterPredictions, main = 'mean prediction')
# points(vect(sdm_td[sdm_td$occurrence==1,]),col = 'black', cex = 0.5)
# points(vect(sdm_td[sdm_td$occurrence==0,]),col = 'red', cex = 0.5)
# 
# plot(sdModel$RasterPredictions_sd, main = 'sd prediction')
# points(vect(sdm_td[sdm_td$occurrence==1,]),col = 'black', cex = 0.5)
# points(vect(sdm_td[sdm_td$occurrence==0,]),col = 'red', cex = 0.5)
# par(mfrow = c(1, 1))

## ----general moel overview----------------------------------------------------
# sdModel$Model

## ----plot aoa-----------------------------------------------------------------
# plot(sdModel$AOA)

## ----check diagnostics--------------------------------------------------------
# sdModel$Diagnostics

## ----Explore SDM output - Confusion Matrix------------------------------------
# knitr::kable(sdModel$ConfusionMatrix$byClass)

## ----Threshold the SDM output-------------------------------------------------
# threshold_rasts <- PostProcessSDM(
#   rast_cont = sdModel$RasterPredictions,
#   test = sdModel$TestData,
#   train = sdModel$TrainData,
#   planar_proj = 5070,
#   thresh_metric = 'equal_sens_spec',
#   quant_amt = 0.2
#   )
# knitr::kable(threshold_rasts$Threshold)

## ----Compare Threshold Results------------------------------------------------
# terra::plot(threshold_rasts$FinalRasters)

## ----Rescale Predictor Variables----------------------------------------------
# rr <- RescaleRasters_bayes(
#   model = sdModel$Model,
#   predictors = sdModel$Predictors,
#   training_data = sdModel$TrainData,
#   pred_mat = sdModel$PredictMatrix,
#   include_uncertainty = TRUE
#   )
# 
# terra::plot(rr$RescaledPredictors)

## ----plot coef uncertainty layer----------------------------------------------
# plot(subset(rr$RescaledPredictors, 'coef_uncertainty'))

## ----Show beta Coefficients from model----------------------------------------
# knitr::kable(
#   rr$BetaCoefficients[
#     order(abs(rr$BetaCoefficients$Estimate),decreasing = T),
#     ])

## ----posterior clustering-----------------------------------------------------
# current_cluster <- PosteriorCluster(
#   model = sdModel$Model,
#   predictors = subset(rr$RescaledPredictors, 'coef_uncertainty', negate = TRUE),
#   f_rasts = threshold_rasts$FinalRasters,
#   lyr = 'Supplemented',
#   pred_mat = sdModel$PredictMatrix,
#   training_data = sdModel$Train,
#   planar_proj = 5070,
#   coord_wt = 0.5,
#   n_pts = 5000
# )

## ----plot stability raster with sample points---------------------------------
# plot(current_cluster$SummaryRaster$stability, main = 'Cluster Stability')
# points(current_cluster$SamplePoints, cex = 0.3)

## ----show posterior clusters, warning=F---------------------------------------
# bmap +
#   geom_sf(data = current_cluster$Geometry, aes(fill = factor(ID))) +
#   labs(fill = 'Cluster', title = 'Consensus Clusters')

## ----view top 3 ranks side by side--------------------------------------------
# par(mfrow = c(1, 3))
# plot(as.factor(current_cluster$SummaryRaster$consensus), main = 'Consensus')
# plot(current_cluster$RankRaster$rank2_final, main = '2nd most')
# plot(current_cluster$RankRaster$rank3_final, main = '3rd most')
# par(mfrow = c(1,1))

## ----Save SDM results, eval = F-----------------------------------------------
# bp <- '~/Documents/assoRted/StrategizingGermplasmCollections'
# 
# writeSDMresults(
#   cv_model = sdModel$CVStructure,
#   model = sdModel$Model,
#   cm = sdModel$ConfusionMatrix,
#   coef_tab = rr$BetaCoefficients,
#   f_rasts = threshold_rasts$FinalRasters,
#   thresh = threshold_rasts$Threshold,
#   file.path(bp, 'results', 'SDM'), 'hemi_test')
# 
# # we can see that the files were placed here using this.
# list.files( file.path(bp, 'results', 'SDM'), recursive = TRUE )

## ----projected clusters into future time point, message = F-------------------
# future_rescaled <- projectClustersBayes(
#   bSDM_object = sdModel,
#   posterior_clusters = current_cluster,
#   future_predictors = bio_future,
#   current_predictors = bio_current,
#   threshold_rasts = threshold_rasts,
#   planar_proj = "epsg:5070",
# )
# 
# nCl = seq(max(future_rescaled$Geometry$ID))
# full_pal = c(
#   '#a6cee3','#1f78b4','#b2df8a','#33a02c',
#   '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
#   '#cab2d6','#6a3d9a','#ffff99','#b15928'
#   )
# cus_pal <-full_pal[nCl]
# names(cus_pal) <- nCl
# 
# bmap +
#   geom_sf(
#     data = current_cluster$Geometry,
#     aes(fill = factor(ID),
#     col = NA)
#     ) +
#   scale_fill_manual(values = cus_pal) +
#   theme(legend.position = 'none', NULL) +
#   labs(title = 'Current') +
# bmap +
#   geom_sf(
#     data = future_rescaled$Geometry,
#     aes(fill = factor(ID)),
#     col = NA
#     ) +
#   scale_fill_manual(values = cus_pal) +
#   labs(title = 'Future', fill = 'Cluster')
# 

## ----clean, echo=F------------------------------------------------------------
# rm(nCl, cus_pal)

## ----Clean up SDM variables, warning=FALSE, echo = FALSE----------------------
# rm(rr, files, bio_future, bio_current, bmap, bbox, bb, current_cluster, threshold_rasts)

## ----dark theme, echo = F-----------------------------------------------------
# update_geom_defaults("line", list(colour = "#dce3f0"))

## ----pp checks 1, message = F-------------------------------------------------
# brms::pp_check(sdModel$Model, ndraws = 10) +
#   labs(title = 'stuff') +
#   theme_navy() +
#   labs(
#     title    = "Posterior Predictive Check - Density",
#     subtitle = "Simulated datasets from the model (light) overlaid on observed data (dark)",
#     caption  = "Lines closely following the observed distribution suggest the model captures the data generating process well."
#   )

## ----pp checks 2, message = F-------------------------------------------------
# brms::pp_check(sdModel$Model, type = "bars", ndraws = 10) +
#   labs(title = 'more') +
#   theme_navy() +
#   labs(
#     title    = "Posterior Predictive Check - Class Balance",
#     subtitle = "Predicted counts of 0s and 1s vs observed across posterior draws",
#     caption  = "Bars heights within the uncertainty range (point and line) indicates the model reproduces the observed presence/absence ratio."
#   )

## ----pp checks 3, message = F-------------------------------------------------
# brms::pp_check(sdModel$Model, type = "stat", stat = "mean") +
#   labs(title = 'more') +
#   theme_navy()  +
#   labs(
#     title    = "Posterior Predictive Check - Mean",
#     subtitle = "Distribution of predicted means across draws vs the observed mean (vertical line)",
#     caption  = "The observed mean should fall near the center of the predicted distribution. Models outside of this area are likely biased."
#   )

## ----posterior predicted probabilities----------------------------------------
# loo_probs <- tibble(
#   pred_prob = brms::loo_epred(sdModel$Model, moment_match = TRUE)[, 1],
#   observed  = factor(sdModel$TrainData$occurrence)
# )
# # Get posterior predicted probabilities
# pred_probs <- sdModel$TrainData |>
#   add_epred_draws(sdModel$Model, ndraws = 100) |>
#   group_by(.row, .draw) |>
#   summarise(prob = mean(.epred), .groups = "drop")

## ----plot loo predicted probabilities-----------------------------------------
# loo_probs |>
#   mutate(bin = cut(pred_prob, breaks = 10)) |>
#   group_by(bin) |>
#   summarise(mean_pred = mean(pred_prob),
#             mean_obs  = mean(as.numeric(as.character(observed)))) |>
#   ggplot(aes(x = mean_pred, y = mean_obs)) +
#   geom_point(col = '#dce3f0') +
#   geom_line(col = '#dce3f0') +
#   geom_abline(linetype = "dashed", col = '#dce3f0') +
#   labs(
#     x = "Mean LOO predicted probability",
#     y = "Observed frequency",
#     title    = "LOO Calibration",
#     subtitle = "Leave-one-out predicted probabilities vs observed frequencies, binned across the probability range",
#     caption  = "Points close to the diagonal indicate well-calibrated predictions. Deviation suggests the model is over or underconfident in certain probability ranges."
#     ) +
#   theme_navy()

## ----auc plot, message = F, warning = F---------------------------------------
# pred_probs |>
#   group_by(.draw) |>
#   summarise(
#     auc = as.numeric(
#       auc(
#         roc(sdModel$TrainData$occurrence[.row],
#         prob,
#         quiet = TRUE))
#         )
#   ) |>
#   ggplot(aes(x = auc)) +
#   stat_halfeye(col = '#8B2635', fill = '#E0E2DB') +
#   labs(
#     x = "AUC",
#     title    = "Posterior Distribution of AUC",
#     subtitle = "Uncertainty in discrimination ability estimated from posterior predicted probabilities on training data",
#     caption  = "AUC closer to 1 indicates better separation between presences and absences.\nThis is estimated on training data and will be optimistic use the matrix data from function."
# ) +
#   theme_navy()

## ----loo roc plot-------------------------------------------------------------
# update_geom_defaults("line",   list(colour = "#dce3f0"))
# roc(loo_probs$observed, loo_probs$pred_prob) |>
#   ggroc() +
#   geom_abline(slope = 1, intercept = 1, linetype = "dashed", col = "#dce3f0") +
#   labs(
#     title    = "LOO ROC Curve",
#     subtitle = "Receiver operating characteristic curve using leave-one-out predicted probabilities",
#     caption  = "The curve shows the tradeoff between true positive and false positive rates across all classification thresholds.\nCurves closer to the top-left corner indicate strong out-of-sample discrimination."
#   ) +
#   theme_navy()

## ----loo predicted probs by class---------------------------------------------
# loo_probs |>
#   ggplot(aes(x = pred_prob, fill = observed)) +
#   stat_halfeye(alpha = 0.6) +
#   scale_fill_manual(values = c('#8B2635', '#dce3f0')) +
#   labs(x = "LOO predicted probability", fill = "Observed") +
#   theme_navy() +
#   labs(
#     title    = "Predicted Probability by Class",
#     subtitle = "LOO predicted probability distributions for pseudo-absences (0) and observed presences (1)",
#     caption  = "Less overlap between distributions indicates the model better separates suitable from unsuitable habitat.\nA high-performing model shows two distinct, separated peaks."
#     )

## ----qqplot from dharma-------------------------------------------------------
# sims <- brms::posterior_predict(sdModel$Model, ndraws = 500)
# 
# dharma_res <- createDHARMa(
#   simulatedResponse = t(sims),
#   observedResponse  = sdModel$TrainData$occurrence,
#   fittedPredictedResponse = apply(sims, 2, mean)
# )
# 
# plot(dharma_res)# QQ + residuals vs fitted

## ----posterior distribution per predictor-------------------------------------
# spr_draw <- sdModel$Model |>
#   spread_draws(`b_.*`, regex = TRUE) |>
#   select(starts_with('b_bio'))
# 
# spr_draw |>
#   pivot_longer(starts_with("b_"), names_to = "param", values_to = "value") |>
#   ggplot(aes(x = value)) +
#   stat_dots(alpha = 0.7, normalize = "panels") +
#   facet_wrap(~ param) +
#   theme(legend.position = "none") +
#   theme_navy() +
#   labs(
#     title    = "Posterior Distribution per Predictor",
#     subtitle = "Each panel shows the full posterior sample of a predictor's coefficient on its own scale",
#     caption  = "Wider distributions indicate greater uncertainty in a variable's effect.\nDistributions overlapping zero suggest weaker or less certain contributions to the model."
#   )

## ----standardised posterior effects plot, warning = F-------------------------
# spr_draw |>
#   mutate(across(everything(), scale)) |>
#   pivot_longer(
#     everything(),
#     names_to = "param",
#     names_transform = list(param = \(x) sub("^b_", "", x)),
#     values_to = "value"
#   )|>
#   ggplot(aes(x = value, y = reorder(param, value))) +
#   stat_halfeye(col = '#8B2635', fill = '#9191E9') +
#   xlim(-5, 5) + ## tighten this up for visualizing on vignette html page
#   geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, col = '#E0E2DB') +
#   labs(
#     x = "Posterior estimate",
#     y = NULL,
#     title    = "Standardised Posterior Effects",
#     subtitle = "Predictor coefficients scaled to unit variance, allowing direct comparison of effect sizes",
#     caption  = "Variables further from zero have stronger effects on predicted occurrence.\nTighter distributions indicate more certainty — wide distributions crossing zero warrant caution."
#     ) +
#   theme_navy()

## ----marginal effects plot, warnings = F--------------------------------------
# plots <- purrr::map(
#   gsub("^b_", '', colnames(spr_draw) ), ~
#   marginaleffects::plot_predictions(sdModel$Model, condition = .x) +
#     labs(title = .x) +
#     xlab(NULL) +
#     ylim(0,1) +
#     theme_navy() +
#     theme(legend.position = 'none') +
#     ggplot2::update_geom_defaults("ribbon", list(fill = "#8b2635")) +
#     ggplot2::update_geom_defaults("line", list(colour = "white"))
# )
# 
# patchwork::wrap_plots(plots, ncol = 2) +
#   plot_annotation(
#     title    = "Marginal Effects per Predictor",
#     subtitle = "Predicted occurrence probability across each variable's range, averaged over all other predictors",
#     caption  = "The shape of each curve reveals whether a variable has a linear, threshold, or unimodal relationship with occurrence.\nWider ribbons indicate greater uncertainty, often at the extremes of a variable's range."
# ) +
#   plot_annotation(theme = theme(
#     plot.background = element_rect(fill = "#0a0e1a", colour = NA),
#     plot.title      = element_text(face = "bold", size = rel(1.2),colour = "#dce3f0"),
#     plot.subtitle   = element_text(colour = "#dce3f0", size = rel(1.1)),
#     plot.caption   = element_text(colour = "#dce3f0", size = rel(1))
#   ))
# 
# update_geom_defaults("line", list(colour = "black"))

## ----eval=F-------------------------------------------------------------------
# rm(theme_navy, spr_draw, sdModel, plots, sims, dharma_res, loo_probs,pred_probs)

