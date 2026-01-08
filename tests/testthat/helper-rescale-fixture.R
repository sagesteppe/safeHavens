rescale_fixture <- function(seed = 1) {

  set.seed(seed)

  r1 <- terra::rast(nrows = 5, ncols = 10, vals = runif(50))
  r2 <- terra::rast(nrows = 5, ncols = 10, vals = runif(50))
  predictors <- c(r1, r2)
  names(predictors) <- c("bio1", "bio2")

  training_data <- data.frame(
    occurrence = sample(c(0, 1), size = 30, replace = TRUE)
  )

  pred_mat <- data.frame(
    bio1 = rnorm(30),
    bio2 = rnorm(30)
  )

  model <- glmnet::glmnet(
    as.matrix(pred_mat),
    training_data$occurrence,
    family = "binomial",
    alpha = 0.5
  )

  list(
    predictors = predictors,
    training_data = training_data,
    pred_mat = pred_mat,
    model = model
  )
}

