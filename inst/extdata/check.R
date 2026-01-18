library(terra)
setwd('~/Documents/assoRted/safeHavens/inst/extdata')

list.files()

dat <- readRDS('sdModel.rds')


names(dat)
plot(dat$RasterPredictions)
