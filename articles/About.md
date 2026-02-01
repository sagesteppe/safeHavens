# About

## about `safeHavens`

This package helps germplasm curators identify and communicate areas of
interest to collection teams for targeting new accessions. For common
species, it provides multiple sampling approaches for curators to choose
from for each individual taxon they aim to collect. It also provides an
*additional* sampling approach for rare species. The package allows for
integration into existing workflows, that is reading in accession data
from databases or Excel, comparing them to occurrence data, and writing
out results into both tabular (e.g.. CSV, a database, Excel, etc.) and
simple spatial data formats (e.g., geopackages, shapefiles, etc.), that
can be shared with collection teams and work on handheld electronic
devices (e.g., Qfield (free, preferred), ESRI products (paid, legacy)).
It allows germplasm curators to have sophisticated spatial capabilities
without the need to be GIS experts or paid software.

### ‘common’ species sampling approaches

Each approach is based on fundamental genetics and ecological theory,
and loosely aligns with the concepts in landscape genetics/genomics. In
practice, all functions attempt to capture the geographic and
environmental variation across a species range. To this end, they evolve
from simply partitioning a species range into `n` equal area, or
regularly spaced areas, to modeling isolation by distance, resistance,
and environment, or sampling from existing spatial data such as
ecoregions or seed transfer zones. The geographic distance functions
largely rely on Sewall Wright’s concept of Isolation by Distance (1943),
the now empirically supported idea that populations that are more
geographically proximate will be more genetically similar than
populations that are further apart. All approaches in the package are *a
priori*, with the exception of the `EnvironmentalBasedSample` function,
which relies on a species distribution model (SDM) to inform sampling
locations and is based on the concept of Isolation by Environment (Wang
& Bradburd 2014); however this still does not have field valdiation.
None of these sampling approaches are meant to supplant genetic
diversity-based sampling, but we recognize that having these data is
rare, and even gathering them would often miss the opportunity to make
opportunistic seed collections along the way.

The methods have various trade-offs in terms of computational and
environmental complexity, although all are designed to be run on a
standard desktop or laptop computer, and processing up to a couple
thousand of species overnight, with the faster methods, should be
possible. The table below presents the currently implemented sampling
scheme and the user-facing function associated with them.

| Function                   | Description                             | Comp. | Envi. |
|----------------------------|-----------------------------------------|-------|-------|
| `PointBasedSample`         | Creates points to make grids over area  | L     | L     |
| `EqualAreaSample`          | Breaks area into similar size pieces    | L     | L     |
| `OpportunisticSample`      | Using PBS with existing records         | L     | L     |
| `IBDBasedSample`           | Breaks species range into clusters      | H     | L     |
| `IBRSurface`               | Breaks species range into clusters      | H     | M     |
| `PolygonBasedSample`       | Using existing ecoregions or STZs       | L     | H     |
| `EnvironmentalBasedSample` | Uses correlations from SDM to sample    | H     | H     |
| `KMedoidsBasedSample`      | For rare species with known occurrences | M     | M     |

### geographic distance

The first two functions, `PointBasedSample`, and `EqualAreaSample`, are
flavors of the same process, where we partition the species range into
geographic clusters of similar sizes.  
The `OpportunisticSample` method is a special case of
`PointBasedSample`, where existing collection records are used to help
guide the sampling scheme, in other words additional samples are
designed around the existing samples. `OpportunisticSample` uses the
same underlying logic as `PointBasedSample`, but first ‘removes’ areas
around existing collection records, and then fills in the remaining gaps
with new sampling locations The `PointBasedSample` and `EqualAreaSample`
functions would be used when *establishing* a new collection strategy
for a species, while `OpportunisticSample` would be used when
*augmenting* an existing collection strategy.

`IBDBasedSample` also relies on geographic distance, but in lieu of
using the *continuity* of geographic space as its primary method, it
focuses on the *discontinuity* of space and uses distance matrices and
clustering to identify portions of the species range that are closer to
each other than to others. This method is slightly more computationally
expensive than the previous geographic distance methods.

### using regions

The `PolygonBasedSample` may be the most commonly encountered method in
North America, and in various formats, is driving two major germplasm
banking projects in both the Midwest and Southeastern United States, as
well as at a high level, composing the way that numerous native seed
collection coordinators are structured in the West. This method uses
environmental variation as an implicit guide to target populations for
seed collections; that is, the different ecoregions serve as
stratification agents. In broad strokes, the general thinking is that
these regions represent continuous transitions in the environment faced
by the species, and populations across these ranges will be differently
adapted to these environments. It can be used with either ecoregion or
seed transfer zone based data. However, it relies on existing spatial
data products, which may or may not be relevant to the ecology of the
species.

### resistance distance

The `IBRSurface` workstream allows users to parameterize a general cost
surface for a landscape. The least-cost paths between occupied portions
of the species range will then calculated and clustered. Spatial data
from the end of this process can be passed into `PolygonBasedSample` if
more than one collection per zone are desired.

### environmental distance

`EnvironmentalBasedSample`, is both the most computationally expensive
and the most environmentally explicit. This function fits a Species
Distribution Model (SDM), generated via a generalized linear model
(glmnet), supported by this package, to cluster populations based on
environmental variables related to their observed distributions and the
spatial configuration and distance between them. On paper, this draws
together all aspects of the above functions; however no empirical
testing of this approach has been implemented.

Spatial data from the end of this process can be passed into
`PolygonBasedSample` if more than one collection per zone are desired.

## rare species sampling approach

All of the above methods are designed for common species, that is
species with more than 50 or 100 unique occurrence records. Each method
returns a polygon geometry that provides guidelines for the general
regions where the species should be sampled. However, for rare species,
where the occurrence of individual populations are tracked, and which
are generally collected along maternal lines, requiring considerably
more field effort to gather seed, a supplemental approach exists that
suggests the priority in which individual populations (‘points’) can be
sampled.

The `KMedoidsBasedSample` method utilizes either geographic distances
alone or geographic distances in conjunction with a distance matrix
developed from a few key environmental variables to suggest which
populations should be prioritized for seed collections. This method is
based on the concept of maximizing the dispersion of selected points in
geographic and environmental spaces to best capture the overall
variation present in the species range.

### installation

`safeHavens` can be installed directly from github, using either
devtools or remotes. Note that the software relies on a few other
packages that may require additional system dependencies, such as
`rgeos` and `rgdal`. It also requires a working installation of `GDAL`,
`PROJ`, and `GEOS` on your system. This can occasionally be finicky to
install, but once you have them effectively all open-source geospatial
softwares are available to you.

`safeHavens` can be installed from GitHub with `remotes` or `devtools`.
Some users have issues with devtools on Windows, but remotes tends to
work well.

``` r
remotes::install_github('sagesteppe/safeHavens') 

# install.packages('devtools') 
# devtools::install_github('sagesteppe/safeHavens')
```
