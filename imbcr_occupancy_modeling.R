#
# Single-season occupancy modeling work for 2016 data
#
# Author: Kyle Taylor
#

stripCommonName <- function(x) tolower(gsub(x,pattern=" |'|-",replacement=""))

spp <-
c(
  'Grasshopper sparrow'
  'Lark bunting',
  'Swainson’s Hawk',
  'Killdeer',
  'Upland sandpiper',
  'Ring-necked Pheasant',
  'Wild-turkey',
  'Mallard',
  'Redhead',
  'Great blue heron'
)

#
# Estimate naive occupancy for each species
#

#
# Fit a basic model that estimates probability of detection using IMBCR metadata
# across 1-km transects and evaluate its accuracy using an 80/20 k-fold cross-validation
#

#
# Fit a single-season occupancy model that allows for heterogeneity in detection probability
# across transects, using Andy Royle's (2012) model specification.
#
