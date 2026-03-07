# Bayesian Approaches

## Species Distribution Modelling

Please refer to all frequentist methods before viewing this vignette.
The rationale is roughly similar, but this approach provides better
estimates of uncertainty for the processing.

### Background

This vignette details the steps required to create a bayesian Species
Distribution Model (SDM) using the functions provided in `safeHavens`,
which form a wrapper around brms::brms(). The SDM is then post-processed
to create a binary raster map of suitable and unsuitable habitat, that
is used to rescale the environmental predictor variables according to
the parameter posteriors for the model. These parameter posteriors are
then used in a clustering algorithm to partition the species range into
environmentally distinct regions for germplasm sampling.

The goal of most SDM’s is to create a model that accurately predicts
where a species will be located in environmental space and can then be
predicted in geographic space. The goal of these models is to understand
the degree and direction to which various environmental features
correlate with a species observed range.

### About

Here we use bayesian hierarchical models to estimate coefficient
uncertainty to use for clustering climatically similar areas via an SDM
approach.

#### prep data

The below steps are all present at the start of the ‘Predictive
Provenance’ vignette. We will ‘hide’ some of the processing steps to
keep this vignette entry relatively short.

``` r
library(safeHavens)
library(terra)
library(geodata)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
```

![](BayesianApproaches_files/figure-html/prep%20basemap-1.png)

### fit the model

Here we fit as SDM using `bayesianSDM`. For arguments it requires
occurrence data `x`, a raster stack of `predictors`, a `quantile_v` that
is an offset used to create pseudo-absence data, and a `planar_proj`
used to calculate distances between occurrence data and possible
pseudo-absence points.

``` r
sdModel <- bayesianSDM( 
  x = hemi,
  pca_predictors = FALSE,
  predictors = bio_current,
  quantile_v = 0.025,
  planar_proj = 5070, 
  resample = TRUE,
  feature_selection = 'ffs',
  min_ffs_var = 5, 
  fact = 3, 
  silent = 2 # print fewer messages during model fit.
  )
```

``` r
#setwd('~/Documents/assoRted/safeHavens/vignettes')
#sdModel$RasterPredictions <- terra::wrap(sdModel$RasterPredictions)
#sdModel$Predictors <- terra::wrap(sdModel$Predictors)
#sdModel$RasterPredictions_sd <- terra::wrap(sdModel$RasterPredictions_sd)
#sdModel$AOA <- terra::wrap(sdModel$AOA)
#saveRDS(sdModel, file.path('..', 'inst', 'extdata', 'BayesSDM.rds'))
```

For the sake of this vignette, we will load a run of the above model.

``` r
sdModel <- readRDS(
  file.path(system.file(package="safeHavens"), 'extdata',  'BayesSDM.rds')
  )

sdModel$RasterPredictions <- terra::unwrap(sdModel$RasterPredictions)
sdModel$Predictors <- terra::unwrap(sdModel$Predictors)
sdModel$RasterPredictions_sd <- terra::unwrap(sdModel$RasterPredictions_sd)
sdModel$AOA <- terra::unwrap(sdModel$AOA)
```

Under the hood this function uses `brms`, which dispatches to Stan, to
fit a bayesian hierarchical model (GLMM), where the distances between
occurrence records are the fixed effect, useful to reduce the effects of
spatial autocorrelation to increase interpret-ability of the models
parameters.

#### explore the output

Plot the predictions

``` r
sdm_td <- sdModel$TrainData

par(mfrow = c(1, 2))
plot(sdModel$RasterPredictions, main = 'mean prediction')
points(vect(sdm_td[sdm_td$occurrence==1,]),col = 'black', cex = 0.5)
points(vect(sdm_td[sdm_td$occurrence==0,]),col = 'red', cex = 0.5)

plot(sdModel$RasterPredictions_sd, main = 'sd prediction')
points(vect(sdm_td[sdm_td$occurrence==1,]),col = 'black', cex = 0.5)
points(vect(sdm_td[sdm_td$occurrence==0,]),col = 'red', cex = 0.5)
```

![](BayesianApproaches_files/figure-html/Explore%20SDM%20output%20-%20Different%20alpha-1.png)

The brms:: model summary can be accessed from the sdModel list.

``` r
sdModel$Model
Warning: There were 32 divergent transitions after warmup. Increasing
adapt_delta above 0.99 may help. See
http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
 Family: bernoulli 
  Links: mu = logit 
Formula: occurrence ~ bio_02 + bio_03 + bio_04 + bio_09 + bio_14 + s(gp_x, gp_y, bs = "gp", k = 50) 
   Data: sf::st_drop_geometry(train) (Number of observations: 415) 
  Draws: 4 chains, each with iter = 5000; warmup = 2000; thin = 1;
         total post-warmup draws = 12000

Smoothing Spline Hyperparameters:
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sds(sgp_xgp_y_1)  1076.12    420.52   479.08  2109.27 1.00     2349     4826

Regression Coefficients:
            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept       5.15      7.15    -3.30    21.86 1.00     2420     5747
bio_02         -0.03      0.11    -0.32     0.09 1.00     5916     5747
bio_03         -0.06      0.10    -0.31     0.03 1.00     3753     6509
bio_04         -0.00      0.00    -0.02     0.00 1.00     2778     6755
bio_09          0.00      0.02    -0.03     0.05 1.00     9355     9702
bio_14         -0.05      0.06    -0.19     0.01 1.00     2793     7626
sgp_xgp_y_1    -0.00      0.25    -0.28     0.28 1.00    11289     7949
sgp_xgp_y_2    -0.01      0.26    -0.32     0.25 1.00    12436     8781

Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

rHat \<1.01 ideal, \<1.05 acceptable, both ESS should be greater than
400 increase iterations to achieve this. Divergence of 0 is ideal, any
divergences are flags for a possibly misspecified model.

Diagnostics can be accessed from a list, and the stan website has some
information on them here as well as warnings (as may be depicted above)
<https://mc-stan.org/misc/warnings/>

rHat \<1.01 ideal, \<1.05 acceptable, both ESS should be greater than
400 increase iterations to achieve this. Divergence of 0 is ideal, any
divergences are flags for a possibly misspecified model.

Similar to `elasticSDM` we also have area of applicabality and other
surfaces from CAST

``` r
plot(sdModel$AOA)
```

![](BayesianApproaches_files/figure-html/plot%20aoa-1.png)

``` r
sdModel$Diagnostics
$max_Rhat
[1] 1.004849

$min_BulkESS
[1] 2131.685

$min_TailESS
[1] 4453.781

$n_divergent
[1] 32
```

To view how well the model is able to predict suitable habitat of the
hold out sample the confusion matrix can be accessed from the test and
train splits. While balanced accuracy for these is generally low, as
these models are not intended for prediction, so much as understanding
the relationships between environmental variables driving the species
distribution, the results show the model was able to learn some
important characteristics of the species.

``` r
sdModel$ConfusionMatrix
$positive
[1] "1"

$table
          Reference
Prediction  0  1
         0 71 11
         1  7 13

$overall
      Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
    0.82352941     0.47959184     0.73551513     0.89192685     0.76470588 
AccuracyPValue  McnemarPValue 
    0.09683511     0.47950012 

$byClass
         Sensitivity          Specificity       Pos Pred Value 
           0.5416667            0.9102564            0.6500000 
      Neg Pred Value            Precision               Recall 
           0.8658537            0.6500000            0.5416667 
                  F1           Prevalence       Detection Rate 
           0.5909091            0.2352941            0.1274510 
Detection Prevalence    Balanced Accuracy 
           0.1960784            0.7259615 

$mode
[1] "sens_spec"

$dots
list()

attr(,"class")
[1] "confusionMatrix"
```

Here we can see how our selected model works on predicting the state of
test data. Note these accuracy results are slightly higher than those
from the CV folds. This is not a bug, CV folds are testing on their own
holdouts, while this is a brand new set of holdouts. The main reason the
confusion matrix results are likely to be higher is due to spatial
auto-correlation which are not address in a typical ‘random split’ of
test and train data. In order to minimize the effects of
spatial-autocorrelation on our model we use `CAST` under the hood which
allows for spatially informed cross validation.

Consider the output from CVStructure to be a bit more realistic.

#### binarize the output

The bayesian sdm is post-processed using the same function as is used by
the `elasticSDM` workflow.

``` r
threshold_rasts <- PostProcessSDM(
  rast_cont = sdModel$RasterPredictions, 
  test = sdModel$TestData,
  train = sdModel$TrainData,
  planar_proj = 5070,
  thresh_metric = 'equal_sens_spec', 
  quant_amt = 0.2
  )
knitr::kable(threshold_rasts$Threshold)
```

|            |    kappa | spec_sens | no_omission | prevalence | equal_sens_spec | sensitivity |
|:-----------|---------:|----------:|------------:|-----------:|----------------:|------------:|
| thresholds | 0.220913 |  0.220913 |   0.0510473 |  0.2444197 |       0.2866278 |   0.2444197 |

We can compare the results of applying this function side by side using
the output from the function.

``` r
terra::plot(threshold_rasts$FinalRasters)
```

![](BayesianApproaches_files/figure-html/Compare%20Threshold%20Results-1.png)

#### rescale predictor variables

The `bayesianSDM` requires a different function from `elasticSDM` for
creating weighted surfaces, with the suffix `RescaleRasters_bayes`.

``` r
rr <- RescaleRasters_bayes(
  model = sdModel$Model,
  predictors = sdModel$Predictors, 
  training_data = sdModel$TrainData, 
  pred_mat = sdModel$PredictMatrix,
  include_uncertainty = TRUE
  )

terra::plot(rr$RescaledPredictors)
```

![](BayesianApproaches_files/figure-html/Rescale%20Predictor%20Variables-1.png)

An optionally returned layer (via include_uncertainty = TRUE), depicts
areas with greater uncertainty towards the prediction.

``` r
plot(subset(rr$RescaledPredictors, 'coef_uncertainty'))
```

![](BayesianApproaches_files/figure-html/plot%20coef%20uncertainty%20layer-1.png)

We can also show the variables that contributed the most to the sdm
models predictions, along with their signs.

``` r
knitr::kable(
  rr$BetaCoefficients[ 
    order(abs(rr$BetaCoefficients$Estimate),decreasing = T), 
    ])
```

|     | Variable |   Estimate | Est.Error |       Q2.5 |     Q97.5 |
|:----|:---------|-----------:|----------:|-----------:|----------:|
| 2   | bio_03   | -0.0597135 | 0.1015725 | -0.3129111 | 0.0253125 |
| 5   | bio_14   | -0.0481234 | 0.0588451 | -0.1918073 | 0.0109816 |
| 1   | bio_02   | -0.0287647 | 0.1124815 | -0.3225050 | 0.0859241 |
| 3   | bio_04   | -0.0047051 | 0.0048937 | -0.0158306 | 0.0014900 |
| 4   | bio_09   |  0.0041970 | 0.0163085 | -0.0253077 | 0.0466395 |

To delineate putative seed collection zones we use the
`PosteriorCluster` function. It will use the raster surfaces of mean
beta coefficient importance

``` r
obby <- PosteriorCluster(
  model = sdModel$Model,
  predictors = subset(rr$RescaledPredictors, 'coef_uncertainty', negate = TRUE),
  f_rasts = threshold_rasts$FinalRasters,
  lyr = 'Supplemented',
  pred_mat = sdModel$PredictMatrix,
  training_data = sdModel$Train,
  planar_proj = 5070
)
Drawing 100 posterior beta samples ...
Sampling fixed point set from prediction surface ...
  500 of 500 requested points have complete predictor data.
Clustering over 100 posterior draws ...
  Draw 25 / 100
  Draw 50 / 100
  Draw 75 / 100
  Draw 100 / 100
Computing co-occurrence matrix ...
Deriving consensus clustering ...
Computing top-3 cluster assignments per point ...
Projecting consensus clusters to raster ...
Loading required package: lattice
```

Using the Bayesian posteriors and clustering allows us to determine
areas that have higher and lower uncertainty for clustering assignments.
Values with a higher proportion had greater colocation scores than areas
with lower proportions.

``` r
par(mfrow=c(1,1))
plot(obby$SummaryRaster$stability)
points(obby$SamplePoints, cex = 0.5)
```

![](BayesianApproaches_files/figure-html/plot%20stability%20raster%20with%20sample%20points-1.png)

As with the other sampling functions, spatial polygons are returned for
use with other applications.

``` r
bmap +
  geom_sf(data = obby$Geometry, aes(fill = factor(ID))) + 
  labs(fill = 'Cluster')
 [1m [22mCoordinate system already present.
 [36mℹ [39m Adding new coordinate system, which will replace the existing one.
```

![](BayesianApproaches_files/figure-html/show%20posterior%20clusters-1.png)

Interested users can see what the second most and third most stable
clustering scenarios would look like via the following ranks. Do note
that the renumbering of clusters from NW to SE (left to right, like many
language groups read books) has not occurred on these rasters yet, so
the numbers differ from the consensus surface.

``` r
par(mfrow = c(1, 3))
plot(obby$SummaryRaster$consensus)
plot(obby$RankRaster$rank2_final)
plot(obby$RankRaster$rank3_final)
```

![](BayesianApproaches_files/figure-html/view%20top%203%20ranks%20side%20by%20side-1.png)

#### save results

writeSDMresults can be used to save most of the core SDM output from the
bayesian work flow as well.

``` r
bp <- '~/Documents/assoRted/StrategizingGermplasmCollections'

writeSDMresults(
  cv_model = sdModel$CVStructure, 
  model = sdModel$Model, 
  cm = sdModel$ConfusionMatrix, 
  coef_tab = rr$BetaCoefficients, 
  f_rasts = threshold_rasts$FinalRasters,
  thresh = threshold_rasts$Threshold,
  file.path(bp, 'results', 'SDM'), 'hemi_test')

# we can see that the files were placed here using this. 
list.files( file.path(bp, 'results', 'SDM'), recursive = TRUE )
```

#### wrapping up

While the bayesian hierarchical GLM does not produce excellent
predictions of suitable habitat, it has several advantages. 1) It can
incorporate spatial smoothers (splines) from mcgv to reduce the effects
of spatial-autocorrelation on model parameter estimates. 2) It provides
uncertainty around the contribution of each selected predictor
(i.e. many predictors can be weighed differently to get similar results)
3)
