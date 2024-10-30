setwd('~/Documents/assoRted/StrategizingGermplasmCollections')
source('./scripts/createPCNM_fitModel.R')
source('./scripts/WriteSDMresults.R')
source('./scripts/postProcessSDM.R')
source('./scripts/trainKNN.R')
source('./scripts/EnvironmentalBasedSample.R')
source('./scripts/RescaleRasters.R')

x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- dplyr::distinct(x, .keep_all = )

files <- list.files(
  path = file.path(system.file(package="dismo"), 'ex'), 
  pattern = 'grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the independent variables

# Step 0 define spatial extent of study area. 

pts_plan <-   sf::st_transform(
  sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326), 
  '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')

bb <- sf::st_bbox(pts_plan)
buff_dist <- as.numeric( # here we get the mean distance of the XY distances of the bb
  ((bb[3] - bb[1]) + (bb[4] - bb[2])) / 2 
) / 2 # the mean distance * 0.25 is how much we will enlarge the area of analysis. 

bb1 <- sf::st_union(pts_plan) |>
  sf::st_buffer(buff_dist) |>
  terra::vect() |>
  terra::project(terra::crs(predictors)) |>
  terra::ext()

p1 <- terra::mask(predictors, bb1)

rm(pts_plan, bb, buff_dist, bb1)
# Step 1 Select Background points - let's use SDM package `envidist` for this

pa <- sdm::background(x = p1, n = nrow(x), sp = x, method = 'eDist') |>
  dplyr::select(lon = x,  lat = y)

pa$occurrence <- 0 ; x$occurrence <- 1
x <- dplyr::bind_rows(x, pa) |> # combine the presence and pseudoabsence points
  sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)  |>
  dplyr::mutate(occurrence = factor(occurrence))

brady.df <- data.frame(Species = 'Species', data.frame(sf::st_coordinates(x)))

dists <- sf::st_distance(x[ sf::st_nearest_feature(x), ], x, by_element = TRUE)
thinD <- as.numeric(quantile(dists, c(0.05)) / 1000) # ARGUMENT TO FN @PARAM 

thinned <- spThin::thin(
  loc.data = brady.df, thin.par = thinD,
  spec.col = 'Species',
  lat.col = 'Y', long.col = 'X', reps = 100, 
  locs.thinned.list.return = TRUE, 
  write.files = FALSE, 
  write.log.file = FALSE)

thinned <- data.frame(thinned[ which.max(unlist(lapply(thinned, nrow)))]) |>
  sf::st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326)

x <- x[lengths(sf::st_intersects(x, thinned))>0,]

rm(thinned, brady.df, pa, thinD, dists, files)

# Step 1.3 - Extract data to points for modelling
x <- terra::extract(predictors, x, bind = TRUE) |>
  sf::st_as_sf() 

# Step 1.2 - create a data split for testing the residuals of the glmnet model
# It's not ideal to do a simple split of these data, because spatial autocorrelation
# will mean that our results could be overly optimistic. 

index <- unlist(caret::createDataPartition(x$occurrence, p=0.8)) # @ ARGUMENT TO FN @PARAM
train <- x[index,]
test <- x[-index,]

rm(index)

# Develop CV folds for modelling
indices_knndm <- CAST::knndm(train, predictors, k=5)

# Recursive feature elimination using CAST developed folds
train_dat <- sf::st_drop_geometry(train[, -which(names(train) %in% c("occurrence"))])
ctrl <- caret::rfeControl(
  method = "LG0CV",
  repeats = 5,
  number = 10,
  functions = caret::lrFuncs,
  index = indices_knndm$indx_train,
  verbose = FALSE)

lmProfile <- caret::rfe(
  method = 'glmnet',
  sizes = c(3:ncol(train_dat)), 
  x = train_dat,
  y = sf::st_drop_geometry(train)$occurrence,
  rfeControl = ctrl)

# Step 4. Model fitting using CAST developed folds
train1 <- dplyr::mutate(
  train, # requires a YES OR NO or T/F just not anything numeric alike. 
  occurrence = dplyr::if_else(occurrence==1, 'YES', 'NO'))

cv_model <- train(
  x = sf::st_drop_geometry(train1[,predictors(lmProfile)]), 
  sf::st_drop_geometry(train)$occurrence, 
  method = "glmnet", 
  family = 'binomial', 
  index = indices_knndm$indx_train) 

sub <- train_dat[,predictors(lmProfile)]

rm(train1, train_dat)

# now fit the model just using glmnet::glmnet in order that we can get the 
# type of response for type='prob' rather than log odds or labelled classes
# which we need to work with terra::predict. 
mod <- glmnet::glmnet(
  x = sub, 
  sf::st_drop_geometry(train)$occurrence, 
  family = 'binomial', 
  keep = TRUE,
  lambda = cv_model$bestTune$lambda, alpha = cv_model$bestTune$alpha
)

rm(cv_model)

obs <- createPCNM_fitModel(
    x = train, 
    planar_proj = 
      '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')

mod <- obs$mod; cv_model <- obs$cv_model; pcnm <- obs$pcnm
pred_mat <- obs$pred_mat

predictors <- c(predictors, pcnm)

# get the variables to extract from the rasters for creating a matrix for 
# predictions, glmnet predict is kind of wonky and needs exact matrix dimensions. 

vars <- rownames(coef(mod)); vars <- vars[2:length(vars)]

# now we need just the COORDINATES FOR TEST and will extract the data from
# this set of predictors to them... 
predict_mat <- predictors[[vars]]
predict_mat <- as.matrix(
  terra::extract(predict_mat, test, ID = FALSE) 
)

cm <- caret::confusionMatrix(
  data = as.factor(predict(mod, newx = predict_mat, type = 'class')), 
  reference = test$occurrence,
  positive="1")

# rm(predict_mat)
## Predict our model onto a gridded surface (raster) ## This will allow for downstream
# use with the rest of the safeHavens workflow. 
preds <- predictors[[vars]]
predfun <- function(model, data, ...){
  predict(model, newx=as.matrix(data), type = 'response')
}

rast_cont <- terra::predict(preds, model = mod, fun=predfun, na.rm=TRUE)

rm(lmProfile, predfun)

ob <- postProcessSDM(rast_cont, thresh_metric = 'sensitivity', quant_amt = 0.25)
f_rasts <- ob$f_rasts
thresh <- ob$thresh

########### CREATE A COPY OF THE RASTER PREDICTORS WHERE WE HAVE 
# STANDARDIZED EACH VARIABLE - SO IT IS EQUIVALENT TO THE INPUT TO THE GLMNET
# FUNCTION, AND THEN MULTIPLIED IT BY IT'S BETA COEFFICIENT FROM THE FIT MODEL

rr <- RescaleRasters(model = mod, predictors = preds, training_data = train)

pred_rescale <- rr$rescaled_predictors
coef_tab <- rr$coefficient_table

# write out the results of the SDM process. 
writeSDMresults(
  file.path( 'results', 'SDM'), 'Bradypus_test')

# perform the clustering 
EnvironmentalBasedSample(
  pred_rescale = pred_rescale, 
  path = file.path( 'results', 'SDM'),
  taxon = 'Bradypus_test', 
  f_rasts = f_rasts, n = 20, fixedClusters = TRUE, n_pts = 500, 
  planar_proj = 
    '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs',
  buffer_d = 3, prop_split = 0.8)

