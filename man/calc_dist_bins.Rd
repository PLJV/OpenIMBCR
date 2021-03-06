% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attributing.R
\name{calc_dist_bins}
\alias{calc_dist_bins}
\title{accepts a Spatial*DataFrame object with a distance observation
field name. Will attempt to arbitrarily bin the distances into
a number of distance bin intervals (specified by breaks). Will
reclass raw distances to bin identifier (e.g., distance class 3).}
\usage{
calc_dist_bins(df = NULL, p = 0.95, breaks = 10,
  force_shoulder = F)
}
\description{
accepts a Spatial*DataFrame object with a distance observation
field name. Will attempt to arbitrarily bin the distances into
a number of distance bin intervals (specified by breaks). Will
reclass raw distances to bin identifier (e.g., distance class 3).
}
