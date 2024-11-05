# safeHavens <img src="inst/man/figures/logo.png" align="right" height="138" />
An R package for applying standardized file naming, directory structure, mapping, and metadata conventions to empirical Seed Transfer Zones (eSTZs). 

## Overview
The goal of this package is to provide germplasm curators with easily referrable spatial data sets to help prioritze field collection efforts.

## Installation
`safeHavens` is available only on github. 
It can be installed using `remotes` or `devtools` like so:
```
install.packages('devtools')
devtools::install_github('sagesteppe/safeHavens')

install.packages('remotes') # remotes is very similar and  a good alternativ for this use case.
remotes::install_github('sagesteppe/safeHavens')
```

Once installed it can be attached for use like any other package from github or CRAN

```
library(safeHavens)
```

## Description
This package helps germplasm curators communicate areas of interest to collection teams for them to collect new material for accession. It provides seven different sampling approaches for curators to choose from for each individual taxon they hope to process. 

## Usage  
`safeHavens` has only seven user facing functions for generating the sampling schemes.  

|        Function           |              Description               | Comp.| Envi.|
|---------------------------|----------------------------------------|------|------|
| `GridBasedSample`         | Creates and merges *n* grids over area |  L   |   L  |
| `PointBasedSample`        | Creates points to make pieces over area|  L   |   L  |
| `EqualAreaSample`         | Breaks area into similar size pieces   |  L   |   L  |
| `OpportunisticSample`     | Using PBS with existing records        |  L   |   L  |
| `IBDBasedSample`          | Breaks species range into clusters     |  H   |   M  |
| `EcoregionBasedSample`    | Using existing ecoregions to sample    |  L   |   H  |
| `EnvironmentalBasedSample`| Uses correlations from SDM to sample   |  H   |   H  |


The species distribution modelling section has a couple functions which are essential for achieving the `EnvironmentalBasedSample` design, these are: `elasticSDM`, `PostProcessSDM`, `RescaleRasters` and `writeSDMresults`.
