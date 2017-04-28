#
# Workflow for Density (birds/km2) estimation for Ring-necked Pheasant using 2016 IMBCR data
#
# We are going to take a really large variable pool and slowly step-down variable selection
# using AIC to select variables for inclusion in a final model.
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(OpenIMBCR)
require(unmarked)
require(rgdal)
require(raster)


formulaToCovariates <- function(formula){
  vars <- as.character(formula)
  vars <- vars[length(vars)]
  vars <- unlist(strsplit(vars,split="[~]"))
  vars <- gsub(unlist(strsplit(vars[length(vars)], split="[+]")), pattern=" ", replacement="")
  return(vars)
}

# read-in our IMBCR source data
s <- imbcrTableToShapefile(recursiveFindFile(name="RawData_PLJV_IMBCR_20161201.csv", root="/home/ktaylora/Incoming"))
# append spatial covariates to our data.frame
long_lat <- data.frame(spTransform(s,CRS(projection("+init=epsg:4326")))@coords)
  names(long_lat) <- c("lon","lat")
    s@data <- cbind(s@data,long_lat)

# determine our state covariates
vars <- list.files("/global_workspace/ring_necked_pheasant_imbcr_models/raster",pattern="tif$",full.name=F)
  vars <- gsub(vars,pattern="2016_|.tif",replacement="")

r <- list.files("/global_workspace/ring_necked_pheasant_imbcr_models/raster",pattern="tif$",full.name=T)
  r <- raster::stack(r)
    names(r) <- vars

# parse our dataset for RNEP records as an unmarked data.frame (distance)
umdf <- OpenIMBCR:::buildUnmarkedDistanceDf(r=r,s=s,spp="Ring-necked Pheasant",vars=c("doy","starttime","lon","lat"))

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
  vars_11x11 <- append(vars_11x11,"crp_age") # let's test crp time-series here

minimum_11x11 <- OpenIMBCR:::randomWalk_dAIC(vars=vars_11x11, umdf=umdf,
                           umFunction=unmarked::distsamp, step=10)
m_11x11_optimal <-
  distsamp(as.formula(as.character(minimum_11x11$formula[nrow(minimum_11x11)])),umdf,keyfun="hazard",output="density",unitsOut="kmsq")


vars_33x33 <- vars[grepl(vars,pattern="33x33")]
minimum_33x33 <- OpenIMBCR:::randomWalk_dAIC(vars=vars_33x33, umdf=umdf,
                          umFunction=unmarked::distsamp, step=10)
m_33x33_optimal <-
  distsamp(as.formula(as.character(minimum_33x33$formula[nrow(minimum_33x33)])),umdf,keyfun="hazard",output="density",unitsOut="kmsq")

vars_107x107 <- vars[grepl(vars,pattern="107x107")]
minimum_107x107 <- OpenIMBCR:::randomWalk_dAIC(vars=vars_107x107, umdf=umdf,
                          umFunction=unmarked::distsamp, step=10)
m_107x107_optimal <-
  distsamp(as.formula(as.character(minimum_107x107$formula[nrow(minimum_107x107)])),umdf,keyfun="hazard",output="density",unitsOut="kmsq")

# model_11x11_first <-
# distsamp(as.formula("~doy~ crp_11x11+hay_11x11+hay_alfalfa_11x11+pasture_11x11+row_crop_11x11+small_grains_11x11+topographic_roughness_11x11")
#          ,umdf,keyfun="hazard",output="density",unitsOut="kmsq")
#
# model_11x11_second <-
# distsamp(as.formula("~doy~ crp_11x11+hay_11x11+pasture_11x11+row_crop_11x11+small_grains_11x11+topographic_roughness_11x11")
#          ,umdf,keyfun="hazard",output="density",unitsOut="kmsq")

# look over the random walk AICs and model, discuss with friends
# optimal <- as.character(minimum$formula[which(minimum$AIC == min(minimum$AIC))[1]])
final_vars <- c("crp_age","row_crop_11x11","small_grains_11x11",
                "topographic_roughness_11x11","row_crop_33x33","crp_11x11",
                "crp_33x33","small_grains_33x33","hay_107x107","hay_alfalfa_107x107"
                )

cat(" -- fitting a first-iteration of our 'final' variable series\n")
minimum_final <- OpenIMBCR:::randomWalk_dAIC(vars=final_vars, umdf=umdf,
                          umFunction=unmarked::distsamp, step=50)

final_vars <- formulaToCovariates(minimum_final$formula)

# keep our crp covariates in for our last round of testing
if (sum(grepl(final_vars,pattern="crp_age")) == 0){
  warning("we need our CRP covariate; manually adding for our last round of variable selection")
  minimum_final[,1] <- as.character(minimum_final[,1])
  minimum_final[nrow(minimum_final), 1] <-
    paste(as.character(minimum_final$formula[nrow(minimum_final)]),
                                  "+crp_age",sep="")
  final_vars <- formulaToCovariates(minimum_final$formula)
}

m_final <- distsamp(as.formula(as.character(minimum_final$formula[nrow(minimum_final)])),
                    umdf,keyfun="hazard",output="density",unitsOut="kmsq")
print(m_final)

cat(" -- intercept test:")
keep <- as.vector(OpenIMBCR:::bartuszevige_intercept_test(m_final))
  keep <- keep[2:length(keep)] # lose the intercept
    final_vars <- final_vars[keep]

cat(final_vars,"\n")
cat(" -- fitting a final model\n")

minimum_final <- OpenIMBCR:::randomWalk_dAIC(vars=final_vars, umdf=umdf,
                          umFunction=unmarked::distsamp, step=5)

# re-fit by dropping any lurking variables that failed the intercept test
m_final <- distsamp(as.formula(as.character(minimum_final$formula[nrow(minimum_final)])),
                    umdf,keyfun="hazard",output="density",unitsOut="kmsq")

# tag some non-linear spatial covariates onto our model for mapping purposes
# this will be overfit and this model won't be used for gleaning information about
# habitat relationships, that's m_final's job. The spatial model's purpose is just
# to "show where birds are at in 2016" with as little error as possible.
spatial_model <- as.character(minimum_final[nrow(minimum_final), 1])
  spatial_model <- paste(spatial_model,"+poly(lon,3)+poly(lat,3)",sep="")

m_final_for_mapping <-
  distsamp(as.formula(spatial_model),
           umdf,keyfun="hazard",output="density",unitsOut="kmsq")
# finish-up
write.csv(minimum_final, "riph_models_selected.csv", row.names=F)
save.image(paste("riph_final_model_",paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="_"),".rdata",sep=""))
