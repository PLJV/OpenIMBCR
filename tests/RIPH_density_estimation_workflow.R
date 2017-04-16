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
umdf <- buildUnmarkedDistanceDf(r=r,s=s,spp="Ring-necked Pheasant")

#
# unmarked distance model fitting (~detection~abundance)
#

# find a decent null model (accounting for nusance detection parameters, if possible)
#m <- distsamp(~1~1,umf,keyfun="hazard",output="density",unitsOut="kmsq")
#m <- distsamp(~doy+starttime~1,umf,keyfun="hazard",output="density",unitsOut="kmsq")

#
# Define a random walk algorithm that seeks to minimize AIC against a randomized subset of
# predictors, recording along a gradient as it goes. Don't just try to select covariates
# based on minimum AIC. Actually consider the biological significance of covariates across
# the full walk of models.
#

# define the model space of all possible combinations of predictors
# parallelize our runs across nCores processors (defined at top)
minimum <- randomWalk_dAIC(vars=vars, umdf=umdf,
                           umFunction=unmarked::distsamp)

# step <- 1000
# total_runs <- 1:nrow(models)
# cat(" -- starting a random walk:\n")
# # assign a null model AIC to beat (below)
# minimum <- data.frame(formula="~doy~1",AIC=m@AIC) # begin with our null (intercept) model
# # iterate over total_runs and try and minimize AIC as you go
# while(length(total_runs)>1){
#   # randomly sample total_runs that the cluster will consider for this run
#   focal_runs <- sample(total_runs,
#                        replace=F,
#                        size=ifelse(length(total_runs) > step, step, length(total_runs))
#                        )
#   # build models for this run
#   runs <- lapply(as.list(models[focal_runs,1]),FUN=as.formula)
#     runs <- parLapply(cl=cl, runs, fun=distsamp, data=umf,keyfun="hazard", output="density", unitsOut="kmsq")
#       runs <- unlist(lapply(runs,FUN=function(x){x@AIC}))
#   # if we beat the running lowest AIC, append it to the random walk table
#   if(runs[which(runs == min(runs))[1]] < min(minimum$AIC)){
#     minimum <- rbind( minimum,
#                       data.frame( formula=models[focal_runs[which(runs == min(runs))[1]],1],
#                                       AIC=runs[which(runs == min(runs))[1]] )
#                     )
#   }
#   total_runs <- total_runs[!(total_runs %in% focal_runs)]
#   cat(paste("[jobs remaining:",length(total_runs),"]",sep=""));
# }; cat("\n");

# look over the random walk AICs and model, discuss with friends
optimal <- as.character(minimum$formula[which(minimum$AIC == min(minimum$AIC))[1]])
m_final <- distsamp(as.formula(optimal),umdf,keyfun="hazard",output="density",unitsOut="kmsq")
# finish-up
write.csv(minimum, "riph_models_selected.csv", rownames=F)
save.image(gsub(paste("RIPH_density_estimation_workspace_",paste(strsplit(as.character(date()),split=" ")[[1]][c(2,3,5)],collapse="_"),".rdata",collapse=""),
                pattern=" ",
                replacement="")
          )
