#
# Workflow for Density (birds/km2) estimation for Ring-necked Pheasant using 2016 IMBCR data
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(OpenIMBCR)
require(unmarked)
require(rgdal)
require(raster)

# read-in our IMBCR source data
s <- imbcrTableToShapefile(recursiveFindFile(name="RawData_PLJV_IMBCR_20161201.csv", root="/home/ktaylora/Incoming"))
# determine our state covariates
vars <- list.files("/global_workspace/ring_necked_pheasant_imbcr_models/raster",pattern="tif$",full.name=F)
  vars <- gsub(vars,pattern="2016_|.tif",replacement="")

r <- list.files("/global_workspace/ring_necked_pheasant_imbcr_models/raster",pattern="tif$",full.name=T)
  r <- raster::stack(r)
    names(r) <- vars

# parse our dataset for RNEP records as an unmarked data.frame (distance)
umdf <- OpenIMBCR:::buildUnmarkedDistanceDf(r=r,s=s,spp="Ring-necked Pheasant")

#
# unmarked distance model fitting (~detection~abundance)
#

# find a decent null model (accounting for nusance detection parameters, if possible)
#m <- distsamp(~1~1,umf,keyfun="hazard",output="density",unitsOut="kmsq")
m <- distsamp(~doy+starttime~1,umdf,keyfun="hazard",output="density",unitsOut="kmsq")

#
# Define a random walk algorithm that seeks to minimize AIC against a randomized subset of
# predictors, recording along a gradient as it goes. Don't just try to select covariates
# based on minimum AIC. Actually consider the biological significance of covariates across
# the full walk of models.
#

# define the model space of all possible combinations of predictors
# parallelize our runs across nCores processors (defined at top)
vars_11x11 <- vars[grepl(vars,pattern="11x11")]
minimum_11x11 <- OpenIMBCR:::randomWalk_dAIC(vars=vars_11x11, umdf=umdf,
                           umFunction=unmarked::distsamp)
m_11x11_optimal <-
  distsamp(as.formula(as.character(minimum_11x11$formula[2])),umdf,keyfun="hazard",output="density",unitsOut="kmsq")


vars_33x33 <- vars[grepl(vars,pattern="33x33")]
minimum_33x33 <- OpenIMBCR:::randomWalk_dAIC(vars=vars_33x33, umdf=umdf,
                          umFunction=unmarked::distsamp)
m_33x33_optimal <-
  distsamp(as.formula(as.character(minimum_33x33$formula[2])),umdf,keyfun="hazard",output="density",unitsOut="kmsq")

vars_107x107 <- vars[grepl(vars,pattern="107x107")]
minimum_107x107 <- OpenIMBCR:::randomWalk_dAIC(vars=vars_107x107, umdf=umdf,
                          umFunction=unmarked::distsamp)
m_107x10_optimal <-
  distsamp(as.formula(as.character(minimum_107x107$formula[2])),umdf,keyfun="hazard",output="density",unitsOut="kmsq")

model_11x11_first <-
distsamp(as.formula("~doy~ crp_11x11+hay_11x11+hay_alfalfa_11x11+pasture_11x11+row_crop_11x11+small_grains_11x11+topographic_roughness_11x11")
         ,umdf,keyfun="hazard",output="density",unitsOut="kmsq")

model_11x11_second <-
distsamp(as.formula("~doy~ crp_11x11+hay_11x11+pasture_11x11+row_crop_11x11+small_grains_11x11+topographic_roughness_11x11")
         ,umdf,keyfun="hazard",output="density",unitsOut="kmsq")
# look over the random walk AICs and model, discuss with friends
# optimal <- as.character(minimum$formula[which(minimum$AIC == min(minimum$AIC))[1]])
m_final <- distsamp(as.formula("~doy~crp_11x11+row_crop_33x33+topographic_roughness_11x11+small_grains_33x33+hay_alfalfa_33x33+pasture_11x11")
                    ,umdf,keyfun="hazard",output="density",unitsOut="kmsq")
# m_final <- distsamp(as.formula(optimal),umdf,keyfun="hazard",output="density",unitsOut="kmsq")
# finish-up
#write.csv(minimum, "riph_models_selected.csv", rownames=F)
