# safeHavens

The goal of this package is to provide germplasm curators with easily
referable spatial data sets to help prioritize field collection efforts.

## Overview

It provides functionality for seven sampling schemes which various
curators are interested in, many of which are likely to outperform
others for certain species or areas. The package also creates species
distribution models, but with the goal of germplasm sampling, rather
than predicting ranges at fine resolutions, or making inference; if you
are interested in this functionality R has several dozen other packages
which are tailored for these purposes.

## Description

This package helps germplasm curators communicate areas of interest to
collection teams for them to collect new material for accession. It
provides seven different sampling approaches for curators to choose from
for each individual taxon they hope to process.

## Installation

`safeHavens` is available only on github. It can be installed using
`remotes` or `devtools` like so:

``` R
install.packages('remotes') 
remotes::install_github('sagesteppe/safeHavens')
#devtools::install_github('sagesteppe/safeHavens') # if needed
```

Once installed it can be attached for use like any other package from
github or CRAN

``` R
library(safeHavens)
```

## Usage

`safeHavens` has only seven user facing functions for generating the
sampling schemes.

### Available Sampling Schemes

The following table shows the eight sampling approaches available in
safeHavens, with their computational complexity (Comp.) and
environmental data requirements (Envi.): L = Low, M = Medium, H = High.

| Function                   | Description                             | Comp. | Envi. |
|----------------------------|-----------------------------------------|-------|-------|
| `PointBasedSample`         | Creates points to make pieces over area | L     | L     |
| `EqualAreaSample`          | Breaks area into similar size pieces    | L     | L     |
| `OpportunisticSample`      | Using PBS with existing records         | L     | L     |
| `KMedoidsBasedSample`      | Use ecoregions or STSz for sample       | L     | M     |
| `IBDBasedSample`           | Breaks species range into clusters      | H     | M     |
| `PolygonBasedSample`       | Using existing ecoregions to sample     | L     | H     |
| `EnvironmentalBasedSample` | Uses correlations from SDM to sample    | H     | H     |

The species distribution modelling section has a couple functions which
are essential for achieving the `EnvironmentalBasedSample` design, these
are: `elasticSDM`, `PostProcessSDM`, `RescaleRasters` and
`writeSDMresults`.

``` r
remotes::install_github('sagesteppe/safeHavens') 
# devtools::install_github('sagesteppe/safeHavens')
```

An overview of the functionality in the package is below.

``` mermaid
%%{init: {'theme':'dark', 'themeVariables': {
    'primaryTextColor':'#ffffff',
    'lineColor':'#ffffff',
    'lineWidth':'24px',
    'fontSize':'18px',
    'fontFamily':'Arial'}}}%%
flowchart LR
A[/Species Occurrence Data/]
B[/Environmental Covariates/]

A --> D(PointBasedSample)
A --> E(EqualAreaSample)
A --> F(OpportunisticSample)
A --> G(IBDBasedSample)
A --> J[elasticSDM]
A --> K[bayesianSDM]
A --> L(KMedoidsBasedSample)
B -.-> L
B --> J
B --> K
D --> N[PrioritizeSample]
E --> N
F --> N
G --> N
A --> U[populationResistance]
B --> H[buildResistanceSurface]
H --> U
U --> V[IBRSurface]
V --> M[PolygonBasedSample]
M --> N
J --> O{postProcessSDM}
K --> O
O --> P{RescaleRaster}
O --> S{RescaleRasterBayes}
P --> Q(EnvironmentalBasedSample)
Q --> R[PredictiveProvenance]
Q --> M
R --> M
S --> T[PosteriorCluster]
T --> R
classDef geoColor fill:#d95f02,color:#FFFFFF
classDef polyColor fill:#66a61e,color:#FFFFFF
classDef envColor fill:#1b9e77,color:#FFFFFF
classDef dataColor fill:#a6761d,color:#FFFFFF
classDef decisionColor fill:#7570b3,color:#FFFFFF
classDef rareColor fill:#e6ab02,color:#FFFFFF
classDef ibrColor fill:#e7298a,color:#FFFFFF
class D,E,F,G geoColor
class H,M polyColor
class J,K,O,P,S,T,U,Q,R envColor
class A,B dataColor
class N decisionColor
class L rareColor
class H,U,V ibrColor
```

## Acknowledgements

Edzer Pebesma, and Krzysztof Dyba for help with seamless installations
on Ubuntu, and getting the pkgdown website running.
