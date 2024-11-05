#' partition data and train a simple KNN model 
#' 
#' Simply use this to partition your test data and quickly train a simple model
#' @param x the weighted matrix including a class ID in the column 'ID'
#' @param split prop for the data partitions. 
#' @keywords internal
#' @export
trainKNN <- function(x, split_prop){
  
  # partition data to ensure no leakage. 
  index <- unlist(caret::createDataPartition(x$ID, p = split_prop)) 
  train <- x[index,]
  test <- x[-index,]
  
  # next we will use a split which ensures that a few members of each class 
  # are represented in each fold
  trainControl <- caret::trainControl(
    method="repeatedcv", number=10, repeats=5)
  
  fit.knn <- caret::train(ID ~ ., data=train, method="knn",
                          trControl = trainControl, metric = 'Accuracy')
  
  # run the model and calculate results on the test data. 
  predicted <- predict(fit.knn, newdata = test)
  cm <- confusionMatrix(predicted, test$ID)
  
  return(
    list(
      fit.knn = fit.knn, 
      confusionMatrix = cm
      )
    )
  
}
