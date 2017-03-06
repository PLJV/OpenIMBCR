#
# Workflow for Density (birds/km2) estimation for Ring-necked Pheasant using 2016 IMBCR data
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(OpenIMBCR)
require(unmarked)
require(rgdal)
require(raster)
require(parallel)

nCores <- parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)

# define our distance breaks through exploratory analysis
d <- c(0,100,200,300,400,500,600,700,800)

# read-in our IMBCR source data
s <- imbcrTableToShapefile(recursiveFindFile(name="RawData_PLJV_IMBCR_20161201.csv", root="/home/ktaylora/Incoming"))
# determine our state covariates
n <- list.files("/global_workspace/ring_necked_pheasant_imbcr_models/raster",pattern="tif$",full.name=F)
  n <- gsub(n,pattern="2016_|.tif",replacement="")

r <- list.files("/global_workspace/ring_necked_pheasant_imbcr_models/raster",pattern="tif$",full.name=T)
  r <- raster::stack(r)
    names(r) <- n
# merge the raster stack containing our state covariates with our IMBCR source table
s <- suppressWarnings(raster::extract(r,s,sp=T))
  s$doy <- as.numeric(strftime(as.POSIXct(as.Date(as.character(s$date), "%m/%d/%Y")),format="%j")) # convert date -> doy
    s@data <- s@data[,!grepl(names(s@data),pattern="FID")]

# kludging to select the covariates we will aggregate and use at the site level
n <- c(names(r),"doy","starttime") # append our covariates on detection

# parse our dataset for RNEP records
t <- s[s$common.name == "Ring-necked Pheasant",]@data
# build a distance table
distances <- data.frame(distance=t$radialdistance,transect=t$transectnum)
y <- unmarked::formatDistData(distances, distCol="distance",transectNameCol="transect",dist.breaks=d)
# build a target matrix
stateCovariates <- matrix(NA,ncol=length(n),nrow=length(levels(t$transectnum))) # e.g., 300 total transects x n state covariates
  rownames(stateCovariates) <- levels(t$transectnum)
# aggregate by field
for(i in 1:length(n)){
  stateCovariates[,i] <- aggregate(s@data[,n[i]], by=list(Category=s@data$transectnum), FUN=mean, na.rm=T)[,2]
}
# specify covariate names
colnames(stateCovariates) <- n
  stateCovariates <- data.frame(stateCovariates) # unmarked expects this as a data.frame

# format our training data as umf
umf <- unmarkedFrameDS(y=as.matrix(y), siteCovs=stateCovariates, survey="point", dist.breaks=d, unitsIn="m")

#
# unmarked distance model fitting (~detection~abundance)
#

# find a decent null model (accounting for nusance detection parameters, if possible)
#m <- distsamp(~1~1,umf,keyfun="hazard",output="density",unitsOut="kmsq")
#m <- distsamp(~doy+starttime~1,umf,keyfun="hazard",output="density",unitsOut="kmsq")
m <- distsamp(~doy~1,umf,keyfun="hazard",output="density",unitsOut="kmsq") # use AIC to determine an optimal detection function
# clean-up our state variables so we don't use our det covs estimate density
n <- n[!grepl(n,pattern="time|doy")]
# calculate : combinations w/o repetition (n!/(r!(n-r)!)
m_len <- vector(); for(i in 1:length(n)) { m_len <- append(m_len,dim(combn(n,m=i))[2]) }
 m_len <- sum(m_len)

# free-up some RAM
rm(r,s,t)
gc()

#
# check for a range of AIC support across all model combinations -- use best judgement on "optimal" variables selected here,
# because this approach is really hackish and not at all canonical. This will produce a table of optimal models randomly
# selected in blocks of 500. The smallest model from each block that beats our null (intercept) model (above) is included
# in this table
#

# minimize AIC across m_len possible models
models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
k <- 1
cat(" -- building a table of all possible models:")
for(i in 1:length(n)){
 combinations <- combn(n,m=i)
 for(j in 1:ncol(combinations)){
   models[k,1] <- paste("~doy~",paste(combinations[,j],collapse="+"),collapse="")
   k <- k+1
 }; cat(".");
}; cat("\n");

# parallelize our runs across nCores processors (defined at top)
total_runs <- 1:nrow(models)
focal_runs <- sample(total_runs,replace=F,size=500) # let's do our runs in blocks of 500
# prime the pump
runs <- lapply(as.list(models[focal_runs,1]),FUN=as.formula)
  runs <- parLapply(cl=cl, runs, fun=distsamp, data=umf,keyfun="hazard", output="density", unitsOut="kmsq")
    runs <- unlist(lapply(runs,FUN=function(x){x@AIC}))
# assign an AIC to beat (below)
minimum <- models[focal_runs[which(runs == min(runs))[1]],]
minimum$AIC <- runs[which(runs == min(runs))[1]]
# trim our total runs
total_runs <- total_runs[!(total_runs %in% focal_runs)]
# iterate over total_runs and try and minimize AIC as you go
while(length(total_runs)>1){
  focal_runs <- sample(total_runs,replace=F,size=100)
  runs <- lapply(as.list(models[focal_runs,1]),FUN=as.formula)
    runs <- parLapply(cl=cl, runs, fun=distsamp, data=umf,keyfun="hazard", output="density", unitsOut="kmsq")
      runs <- unlist(lapply(runs,FUN=function(x){x@AIC}))
  if(runs[which(runs == min(runs))[1]] < minimum$AIC){
    minimum <- rbind(minimum,
                     data.frame(formula=models[focal_runs[which(runs == min(runs))[1]],1],AIC=runs[which(runs == min(runs))[1]]))
  }
  total_runs <- total_runs[!(total_runs %in% focal_runs)]
  cat(paste("[jobs remaining:",length(total_runs),"]",sep=""));
}; cat("\n");

# take a peak at our models table

# fit our "final" optimized model
optimal <- models$formula[which(models$AIC == min(models$AIC))]
m_final <- distsamp(as.formula(optimal),umf,keyfun="hazard",output="density",unitsOut="kmsq")
# finish-up
save.image(paste("RIPH_density_estimation_workspace_",paste(strsplit(as.character(date()),split=" ")[[1]][c(2,3,5)],collapse="_"),".rdata",collapse=""))
