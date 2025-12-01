# Rare Species Sampling Schema

## prepare data

Load the required packages.

``` r
library(safeHavens)
library(ggplot2)
library(patchwork)
set.seed(99)
```

Here we will use the Bradypus data included in the `dismo` package
again.

``` r
x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
```

And we will create the same base map used in `GettingStarted`.

    #> Warning: attribute variables are assumed to be spatially constant throughout
    #> all geometries

While all other functions in the package handle `sf` objects directly,
this function will actually just use a simple data frame of sites, to
simplfy handing data off to C++ for optimization routines.

The input to the `maximizeDispersion` function is a list with two
elements: a distance matrix, and a data frame of site locations and
attributes. The data frame must contain the following columns.

``` r
n_sites <- nrow(x) 
df <- data.frame(
  site_id = seq_len(n_sites),
  required = FALSE,
  coord_uncertainty = 0, 
  lon = sf::st_coordinates(x)[,1], 
  lat = sf::st_coordinates(x)[,2]
)

head(df)
#>   site_id required coord_uncertainty      lon      lat
#> 1       1    FALSE                 0 -65.4000 -10.3833
#> 2       2    FALSE                 0 -65.3833 -10.3833
#> 3       3    FALSE                 0 -65.1333 -16.8000
#> 4       4    FALSE                 0 -63.6667 -17.4500
#> 5       5    FALSE                 0 -63.8500 -17.4000
#> 6       6    FALSE                 0 -64.4167 -16.0000
```

The second required element, the distance matrix, can be calculated with
the `greatCircleDistance` function in the package. Please use this
rather than st_distance from `sf` for consistency, as the units differ
slightly. If you want to use
[`sf::st_distance`](https://r-spatial.github.io/sf/reference/geos_measures.html),
make sure to convert the units to match the scale of the
`greatCircleDistance` function, otherwise the results will be incorrect.

``` r
dist_mat <- sapply(1:nrow(df), function(i) {
   greatCircleDistance(
     df$lat[i], df$lon[i],
     df$lat, df$lon
   )
 })
```

The optimization routine requires at least one ‘required’ site to be
specified. Here we will select the site closest to the geographic center
of all sites as the required site.

Normally this can refer to existing accessions, or administrative units,
or preserves which are helping to implement the germplasm collection,
and are fortunate enough to already have some samples or at least
guaranteed access.

``` r
dists2c <- greatCircleDistance(
  median(df$lat), 
  median(df$lon), 
  df$lat, 
  df$lon
)
df[order(dists2c)[1],'required'] <- TRUE
```

This function not only bootstraps sites to simulate the true
distribution distribution of the species, but it also bootstraps
coordinate uncertainty for each site. Here we will randomly assign 20%
of the sites to have coordinate uncertainty between 1 km and 40 km. Note
that he argument is always in meters.

``` r
uncertain_sites <- sample(
  setdiff(seq_len(n_sites), 
  which(df$required)), 
  size = round(n_sites*0.2, 0)
  )
df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 1000, 40000) # meters
```

## Run dispersion maximization based only on geographic distances

The input to the function is the distance matrix, and the site data.

``` r
test_data <- list(
  distances = dist_mat,
  sites = df
  )

str(test_data)
#> List of 2
#>  $ distances: num [1:116, 1:116] 0 1.83 714.09 807.71 797.93 ...
#>  $ sites    :'data.frame':   116 obs. of  5 variables:
#>   ..$ site_id          : int [1:116] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..$ required         : logi [1:116] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ coord_uncertainty: num [1:116] 0 0 0 0 0 ...
#>   ..$ lon              : num [1:116] -65.4 -65.4 -65.1 -63.7 -63.9 ...
#>   ..$ lat              : num [1:116] -10.4 -10.4 -16.8 -17.4 -17.4 ...

rm(x, n_sites, uncertain_sites, dists2c)
```

The funtion `maximizeDispersion` has several parameters to control the
optimization routine.

``` r
st <- system.time( {
    geo_res <- maximizeDispersion(  ## reduce some parameters for faster run. 
      input_data = test_data,
      n_sites = 10,
      lambda_var = 0.05,
      n_bootstrap = 500,
      objective = "sum",
      n_local_search_iter = 50,
      n_restarts = 2
    )
  }
)
#> Sites: 116 | Seeds: 1 | Requested: 10 | Coord. Uncertain: 23 | BS Replicates: 500
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
```

The function operates relatively quick with few bootstraps and few
sites, but will take considerably longer with more complex scenarios. We
recommened using at least 999 bootstraps for real world applications.

``` r
st
#>    user  system elapsed 
#>  25.106   0.005  25.113
rm(st)
```

### return output structure

Various elements are returned in the output list.

``` r
str(geo_res)
#> List of 5
#>  $ input_data     :'data.frame': 116 obs. of  9 variables:
#>   ..$ site_id          : int [1:116] 47 7 14 16 27 29 61 73 75 76 ...
#>   ..$ required         : logi [1:116] TRUE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ coord_uncertainty: num [1:116] 0 0 5301 0 0 ...
#>   ..$ lon              : num [1:116] -74.3 -63.2 -40.1 -49.5 -49.5 ...
#>   ..$ lat              : num [1:116] 4.58 -17.8 -19.33 -1 -2.25 ...
#>   ..$ cooccur_strength : num [1:116] 4501 4500 4500 4500 4500 ...
#>   ..$ is_seed          : logi [1:116] TRUE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ selected         : logi [1:116] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   ..$ sample_rank      : int [1:116] 1 2 2 2 2 2 2 2 2 2 ...
#>  $ selected_sites : int [1:10] 7 14 16 27 29 47 61 73 75 76
#>  $ stability_score: num 1
#>  $ stability      :'data.frame': 116 obs. of  3 variables:
#>   ..$ site_id         : int [1:116] 47 7 14 16 27 29 61 73 75 76 ...
#>   ..$ cooccur_strength: num [1:116] 4501 4500 4500 4500 4500 ...
#>   ..$ is_seed         : logi [1:116] TRUE FALSE FALSE FALSE FALSE FALSE ...
#>  $ settings       :'data.frame': 1 obs. of  6 variables:
#>   ..$ n_sites     : num 10
#>   ..$ n_bootstrap : num 500
#>   ..$ objective   : chr "sum"
#>   ..$ lambda      : num 0.05
#>   ..$ dropout_prob: num 0.15
#>   ..$ n_uncertain : int 23
```

The stability score shows how often the most frquently selected network
of sites was selected from the bootstrapped runs.

``` r
head(geo_res$stability_score)
#> [1] 1
```

The stability data frame shows how often each site was selected across
all bootstrap runs.

``` r
head(geo_res$stability)
#>    site_id cooccur_strength is_seed
#> 47      47             4501    TRUE
#> 7        7             4500   FALSE
#> 14      14             4500   FALSE
#> 16      16             4500   FALSE
#> 27      27             4500   FALSE
#> 29      29             4500   FALSE
```

Many users may find the combindation of their input data with a few
columns, to be all they need to carry on after the results.

``` r
head(geo_res$input_data)
#>    site_id required coord_uncertainty      lon      lat cooccur_strength
#> 47      47     TRUE             0.000 -74.3000   4.5833             4501
#> 7        7    FALSE             0.000 -63.1667 -17.8000             4500
#> 14      14    FALSE          5301.306 -40.0667 -19.3333             4500
#> 16      16    FALSE             0.000 -49.5000  -1.0000             4500
#> 27      27    FALSE             0.000 -49.5000  -2.2500             4500
#> 29      29    FALSE             0.000 -46.7333 -23.4500             4500
#>    is_seed selected sample_rank
#> 47    TRUE     TRUE           1
#> 7    FALSE     TRUE           2
#> 14   FALSE     TRUE           2
#> 16   FALSE     TRUE           2
#> 27   FALSE     TRUE           2
#> 29   FALSE     TRUE           2
```

Run parameters are saved in the settings element.

``` r
head(geo_res$settings)
#>   n_sites n_bootstrap objective lambda dropout_prob n_uncertain
#> 1      10         500       sum   0.05         0.15          23
```

### visualize the selection results

``` r
map + 
  geom_point(data = geo_res$input_data, 
  aes(
    x = lon, 
    y = lat, 
    shape = required, 
    size = cooccur_strength,
    color = selected
    )
  ) +
 # ggrepel::geom_label_repel(aes(label = site_id), size = 4) + 
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites; Geographic Distances')
```

![](RareSpecies_files/figure-html/first%20selection%20of%20target%20sites-1.png)

``` r
map + 
  geom_point(data = geo_res$input_data, 
    aes(
      x = lon, 
      y = lat, 
      shape = required, 
      size = -sample_rank,
      color = sample_rank
      )
    ) +
 # ggrepel::geom_label_repel(aes(label = sample_rank), size = 4) +
  theme_minimal()   
```

![](RareSpecies_files/figure-html/priority%20ranking%20plot-1.png)

## run dispersion maximization with geographic and environmental distances

### extract prep environmental distances

``` r
files <- list.files(
  path = file.path(system.file(package="dismo"), 'ex'), 
  pattern = 'grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the independent variables
rm(files)
```

For our environmental distances, we will use a PCA transformation of the
environmental variables. We will simply scrape 100 random points from
the raster layers to calculate the PCA. Then predict the PCA raster
layers across the entire study area. We will take the first two layers,
and calculate environmental distances based on these two layers.

``` r
pts <- terra::spatSample(predictors, 100, na.rm = TRUE)
pts <- pts[, names(pts)!='biome' ] # remove categorical variable for distance calc

pca_results <- stats::prcomp(pts, scale = TRUE)
round(pca_results$sdev^2 / sum(pca_results$sdev^2), 2) # variance explained
#> [1] 0.55 0.29 0.09 0.05 0.02 0.00 0.00 0.00
pca_raster <- terra::predict(predictors, pca_results)

terra::plot(terra::subset(pca_raster, c(1:2))) # prediction of the pca onto a new raster
```

![](RareSpecies_files/figure-html/run%20with%20geographic%20and%20environmental%20distances-1.png)

``` r
rm(pts, predictors, pca_results)
```

we keep the first two PCA layers for environmental distance calculation.
More layers will increase dimenstionality, and may lead to less useful
results. Note that it’s fine to use a euclidean distance calculation for
these, as the values are truly in the position of the plot.

``` r
env_values <- terra::extract(pca_raster, 
  sf::st_coordinates(
    sf::st_as_sf(
      df, 
      coords = c('lon', 'lat'), 
      crs = 4326
    )
  )
)[,1:2]
plot(env_values, main = 'environmental distance of points from first two PCA axis')
```

![](RareSpecies_files/figure-html/extract%20environmental%20values%20and%20calculate%20distance%20matrix-1.png)

``` r

env_dist_mat <- as.matrix(
    dist(env_values)
  )

rm(pca_raster)
```

``` r
test_data <- list(
  distances = simplify2array(list(dist_mat, env_dist_mat)),
  sites = df
  )

st <- system.time( 
  {
    env_res <- maximizeDispersion(  ## reduce some parameters for faster run. 
      input_data = test_data,
      n_sites = 10,
      lambda_var = 0.1,
      weight_1 = 0.1, 
      weight_2 = 0.9,
      n_bootstrap = 99, # considerably slower function!!!
      objective = "sum",
      n_local_search_iter = 50,
      n_restarts = 2
    )
  }
)
#> Sites: 116 | Seeds: 1 | Requested: 10 | Coord. Uncertain: 23 | BS Replicates: 99
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

rm(dist_mat, env_dist_mat)
```

Adding the second distance matrix into the mix slows things down
considerable.

Note that when this function has jittered points, not only is the point
jittered randomly for the geographic component, the environmental
distance matrix is also jittered. The jitter for the latter matrix come
from fitting a linear model where environmental distance is fit as a
function of geographic distance. The conditional mean is then used as a
noise to alter the environmental distance to/from each record. However,
this is unrelated to the slowdown.

``` r
st
#>    user  system elapsed 
#>  81.608   0.004  81.622
rm(st)
```

``` r
head(env_res$stability_score)
#> [1] 0.8787879
```

``` r
map + 
  geom_point(data = env_res$input_data, 
    aes(
      x = lon, 
      y = lat, 
      shape = required, 
      size = cooccur_strength,
      color = selected
      )
    ) +
 # ggrepel::geom_label_repel(aes(label = site_id), size = 4) + 
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites; Environmental')
```

![](RareSpecies_files/figure-html/unnamed-chunk-9-1.png)

## alternative methods for required central points

In the example above we use a point at the median geographic center of
the populations.

We can also identify the population which is *most* near the highest
density of populations. Intuitively, this would be suggested as a
population with a very high genetic diversity.

``` r
dens <- with(df, MASS::kde2d(lon, lat, n = 200))
max_idx <- which(dens$z == max(dens$z), arr.ind = TRUE)[1,]
max_point <- c(dens$x[max_idx[1]], dens$y[max_idx[2]])

pops_centre <- sweep(df[c('lon', 'lat')], 2, max_point, "-")
pop_centered_id <- which.min(rowSums(abs(pops_centre^2)))

rm(dens, max_idx, max_point, pops_centre)
```

Likewise we can identify the *population* which is most near the
‘center’ of the environmental variable space.

``` r
env_centered <- sweep(env_values, 2, sapply(env_values, median), "-")
env_centered_id <- which.min(rowSums(abs(env_centered^2)))

rm(env_values)
```

All three of these options can be combined to feed in three center
points for required sampling locations, providing a robust set of ‘core’
diversity for the species. Personally I would consider the ‘pop
centered’ population to be the most important required site to center a
design off of.

Here we will just showcase their positions

``` r
# geographic centroid was pt 47
centers <- df[ c(env_centered_id, pop_centered_id, 47), ] 
centers$type <- c('Environmental', 'Population', 'Geographic')

map +
  geom_point(
    data = df, 
    aes(x = lon, y = lat)
    ) + 
  geom_point(
    data = centers,  
    aes(x = lon, y = lat),
    col = '#FF1493', size = 4
    ) + 
  ggrepel::geom_label_repel(
    data = centers, 
    aes(label = type, x = lon, y = lat)
    ) + 
  theme_minimal() + 
  labs(title = 'Possbilities for required centers')
```

![](RareSpecies_files/figure-html/plot%20alternative%20centroids-1.png)

``` r

rm(env_centered_id, env_centered, pop_centered_id)
```

Using these three centers will give results as below

``` r
test_data$sites$required[centers$site_id] <- TRUE

env_res <- maximizeDispersion(  ## reduce some parameters for faster run. 
  input_data = test_data,
  n_sites = 10,
  lambda_var = 0.1,
  weight_1 = 0.3, 
  weight_2 = 0.7,
  n_bootstrap = 99, 
  objective = "sum",
  n_local_search_iter = 50,
  n_restarts = 2
)
#> Sites: 116 | Seeds: 3 | Requested: 10 | Coord. Uncertain: 23 | BS Replicates: 99
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

map + 
  geom_point(data = env_res$input_data, 
    aes(
      x = lon, 
      y = lat, 
      shape = required, 
      size = cooccur_strength,
      color = selected
      )
    ) +
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites with 3 required sites')
```

![](RareSpecies_files/figure-html/run%20multi%20required%20sites%20and%20plot-1.png)

``` r

rm(centers)
```

## show differences in lambda parameters

`maximizeDispersion` optimizes on two parameters. The first if seeking
to maximize the distance between selected populations. If this optimizer
runs alone ($\lambda$ = 0.0), it results in *only* drawing a set of
populations along the margin of the species. By introducting a
regularization term ($\lambda$), a balance of optimizing distance, with
also decreasing the variance between points is achieved.

Theoretically maintaining a low $\lambda$ helps maximize the geographic
coverage of populations at a small sample size. However, we find that
the best selection of populations are likely to come from a two-stage
approach. One where a user specifies the: - 1) required sites, - 2) runs
the algorithm to indentify ‘fringing’ sites, which may have more allelic
diversity, - 3) subjectively placing a couple sites across the range as
they see fit.

However, I believe that grid-based approaches may still be the most
suitable for most practicioners situations.

``` r
lambda_viewer <- function(lambda_val){

  env_res <- maximizeDispersion(   
    input_data = test_data,
    n_sites = 10,
    lambda_var = lambda_val,
    n_bootstrap = 50, 
    objective = "sum",
    n_local_search_iter = 50,
    n_restarts = 2, 
    verbose = FALSE
  )

  m <- map + 
    geom_point(data = env_res$input_data, 
      aes(
        x = lon, 
        y = lat, 
        shape = required, 
        size = cooccur_strength,
        color = selected
        )
      ) +
    labs(title = paste('lambda: ', lambda_val), x = NULL, y = NULL)

}

lambda_view <- lapply(c(0, 0.025, 0.05, 0.075, 0.1, 0.15), lambda_viewer)

lambda_view[[1]] + lambda_view[[2]] + lambda_view[[3]] + 
  lambda_view[[4]]  +  lambda_view[[5]] + lambda_view[[6]] + 
  plot_layout(ncol = 3)
```

![](RareSpecies_files/figure-html/unnamed-chunk-10-1.png)

``` r

rm(lambda_view, lambda_viewer, map)
```

## closing thoughts

This approach is highly experimental. Small permutations in lambda may
lead to large differences in the nature of results. Further, results are
often ‘ring like’, a consequence of both maximizing distance between
sites while minimizing variance.
