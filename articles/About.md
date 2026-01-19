# About

## about `safeHavens`

This package helps germplasm curators communicate areas of interest to
collection teams for targeting new accessions. For common species, it
provides seven different sampling approaches for curators to choose from
for each individual taxon they hope to process. It also provides an
*additional* sampling approach for rare species. The package allows for
easy integration into existing workflows, that is reading in accession
data from databases or Excel, comparing them to occurrence data, and
then writing out run results into both tabular and simple spatial data
formats, which can be shared with collection teams and work on handheld
electronic devices (e.g. Qfield). It allows germplasm curators to have
sophisticated spatial abilities without the need to be GIS experts or
have access to any fee-based software.

#### ‘common’ species sampling approaches

Each approach is based on fundamental ecological theory. In practice,
they all attempt to capture the geographic and environmental variation
across a species range and largely rely on Sewall Wright’s concept of
Isolation by Distance (1943), the idea that populations that are more
geographically proximate will be more genetically similar than
populations that are further apart. Hence, they are all a priori
approaches, with the exception of the `EnvironmentalBasedSample`
function, which relies on a species distribution model (SDM) to inform
sampling locations and is based on the concept of Isolation by
Environment (Wang & Bradburd 2014). While none of these are meant to
supplant genetic diversity-based sampling, we recognize that having
these data is rare, and even gathering them would often miss the
opportunity to make opportunistic seed collections along the way.

The methods have various trade-offs in terms of computational and
environmental complexity, although all are designed to be run on a
standard desktop or laptop computer, and processing many thousands of
species overnight is possible. The table below presents the currently
implemented sampling scheme and the user-facing function associated with
them.

| Function                   | Description                             | Comp. | Envi. |
|----------------------------|-----------------------------------------|-------|-------|
| `PointBasedSample`         | Creates points to make grids over area  | L     | L     |
| `EqualAreaSample`          | Breaks area into similar size pieces    | L     | L     |
| `OpportunisticSample`      | Using PBS with existing records         | L     | L     |
| `KMedoidsBasedSample`      | For rare species with known occurrences | M     | L     |
| `IBDBasedSample`           | Breaks species range into clusters      | H     | M     |
| `PolygonBasedSample`       | Using existing ecoregions or STZs       | L     | H     |
| `EnvironmentalBasedSample` | Uses correlations from SDM to sample    | H     | H     |

The first two functions, `PointBasedSample`, and `EqualAreaSample`, are
flavors of the same process, where we try to partition the species range
into geographic chunks (polygons) of similar sizes.  
Each requires minimal computational power but features essentially no
environmental information of context beyond the species current range
restriction. The `OpportunisticSample` method is a special case of the
first `PointBasedSample`, where existing collection records are used to
help guide the sampling scheme. `OpportunisticSample` uses the same
underlying logic as `PointBasedSample`, but first ‘removes’ areas around
existing collection records, and then fills in the remaining gaps with
new sampling locations The first `PointBasedSample`, and
`EqualAreaSample`, functions would be used when *establishing* a new
collection strategy for a species, while `OpportunisticSample` would be
used when trying to *augment* an existing collection strategy.

`IBDBasedSample`, is largely in a class of its own; in lieu of using the
*continuity* of geographic space as its primary method, it focuses on
the discontinuity of space and uses distance matrices and clustering to
determine which patches of range are closer to each other than to other
patches. This method is more computationally expensive than the previous
four, and while it incorporates some environmental information, it is
not explicit in the way that the final two methods are.

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

`EnvironmentalBasedSample`, is both the most computationally expensive
and the most environmentally explicit. This function fits a Species
Distribution Model (SDM), generated via a generalized linear model
(glmnet), supported by this package, to cluster populations based on
environmental variables related to their observed distributions and the
spatial configuration and distance between them. On paper, this draws
together all aspects of the above functions; however no empirical
testing of this approach has been implemented.

#### rare species sampling approach

All of the above methods are designed for common species, that is
species with more than 50 or 100 unique occurrence records. Each method
returns a polygon geometry that provides guidelines for the general
regions where the species should be sampled. However, for rare species,
where the occurrence of individual populations is well tracked, and
which are generally collected along maternal lines, requiring
considerably more field effort to gather seed, a supplemental approach
exists that suggests the priority in which individual populations
(‘points’) can be sampled.

The `KMedoidsBasedSample` method utilizes either geographic distances
alone or geographic distances in conjunction with a distance matrix
developed from a few key environmental variables to suggest which
populations should be prioritized for seed collections. This method is
based on the concept of maximizing the dispersion of selected points in
geographic and environmental spaces to best capture the overall
variation present in the species range.

#### installation

`safeHavens` can be installed directly from github, using either
devtools or remotes. Note that the software relies on a few other
packages that may require additional system dependencies, such as
`rgeos` and `rgdal`. It also requires a working installation of `GDAL`,
`PROJ`, and `GEOS` on your system, as well as Rcpp.

Generally, these installations go off without a hitch but may require
additional attention from some users. The trade off is that these are
all *free* tools with enormous power and are quite commonly used in
ecological and geographical modelling.

`safeHavens` can be installed from GitHub with `remotes` or `devtools`.
Some users have issues with devtools on windows, but remotes tends to
work well.

``` r
remotes::install_github('sagesteppe/safeHavens') 

# install.packages('devtools') 
# devtools::install_github('sagesteppe/safeHavens')
```
