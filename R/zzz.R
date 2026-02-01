utils::globalVariables(c(
  # data for nse
  ".data",

  # common spatial / tidy columns
  "X",
  "Y",
  "ID",
  "geometry",
  "cell",
  "to_cluster",

  # allocation / assignment flags
  "Assigned",
  "NEWID",
  "Supplemented",

  # clustering / modeling objects
  "kmeans_centers",
  "sdModel",
  "occurrence",
  "s0",
  "pop_res_graphs",

  # polygon / area summaries
  "poly_area",
  "polygon_ct",
  "total_area_m2",

  # PCNM / SDM helpers
  "predictors",

  # NbClust / clustering parameters
  "min.nc",
  "max.nc",

  # ranking / prioritization
  "Level",

  # misc NSE / modeling symbols
  "x",
  "y",

  # dplyr placeholder (rare but legit)
  ".",

  # random
  'NewAssigned',
  'epsilon'
))
