#
# exploratory Modeling Workflow for Density (birds/km2) Estimation for a Given Species of Interest
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(OpenIMBCR)
require(unmarked)
require(rgdal)
require(raster)
require(parallel)

system("clear")

# RUNTIME arguments

four_letter_code <- NULL
 spp_common_name <- NULL
          nCores <- parallel::detectCores()-1
              cl <- parallel::makeCluster(nCores)

cat(" -- species density estimation workflow\n")

# process any arguments passed by user

argv <- commandArgs(trailingOnly=T)

for(i in 1:(length(argv)-1)){
  if( grepl(tolower(argv[i]),pattern="-s|--spp") ){
    spp_common_name <- gsub(tolower(argv[i+1]),pattern="'",replacement="")
    cat(paste(" -- focal species:",spp_common_name,"\n"))
  } else if( grepl(tolower(argv[i]),pattern="-c|--code") ) {
    four_letter_code <- toupper(argv[i+1])
  }
}

if( grepl(class(four_letter_code),pattern="NULL") | is.na(four_letter_code) ){
  code <- unlist(strsplit(spp_common_name,split=" "))
  if(length(code) == 1){
    four_letter_code <- toupper(substr(code,1,4))
  } else if(length(code) == 2){
    four_letter_code <- paste(toupper(substr(code[1],1,2)), toupper(substr(code[2],1,2)), collapse="")
  }
  cat(" -- guessing species four-letter code as (use -c to force specificaiton):", four_letter_code, "\n")
}

cat(" -- parsing explantory (raster) variables\n")
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

# parse our dataset for focal species records
t <- s[OpenIMBCR:::stripCommonName(s$common.name) == OpenIMBCR:::stripCommonName(spp_common_name),]@data

if(nrow(t)==0) stop("couldn't find any records for species: ",spp_common_name)

# build a distance table
cat(" -- building a distance table and fitting a null (intercept) model\n")
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
  stateCovariates <- data.frame(stateCovariates) # unmarked expects this ahtos a data.frame

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

#
# Define a random walk algorithm that seeks to minimize AIC against a randomized subset of
# predictors, recording along a gradient as it goes. Don't just try to select covariates
# based on minimum AIC. Actually consider the biological significance of covariates across
# the full walk of models.
#

# define the model space of all possible combinations of predictors
models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
k <- 1
cat(" -- deriving model combinations:")
for(i in 1:length(n)){
 combinations <- combn(n,m=i)
 for(j in 1:ncol(combinations)){
   models[k,1] <- paste("~doy~",paste(combinations[,j],collapse="+"),collapse="")
   k <- k+1
 }; cat(".");
}; cat("\n");
# parallelize our runs across nCores processors (defined at top)
step <- 1000
total_runs <- 1:nrow(models)
cat(" -- starting a random walk:\n")
# assign a null model AIC to beat (below)
minimum <- data.frame(formula="~doy~1",AIC=m@AIC) # begin with our null (intercept) model
# iterate over total_runs and try and minimize AIC as you go
while(length(total_runs)>1){
  # randomly sample total_runs that the cluster will consider for this run
  focal_runs <- sample(total_runs,
                       replace=F,
                       size=ifelse(length(total_runs) > step, step, length(total_runs))
                       )
  # build models for this run
  runs <- lapply(as.list(models[focal_runs,1]),FUN=as.formula)
    runs <- parLapply(cl=cl, runs, fun=distsamp, data=umf,keyfun="hazard", output="density", unitsOut="kmsq")
      runs <- unlist(lapply(runs,FUN=function(x){x@AIC}))
  # if we beat the running lowest AIC, append it to the random walk table
  if(runs[which(runs == min(runs))[1]] < min(minimum$AIC)){
    minimum <- rbind( minimum,
                      data.frame( formula=models[focal_runs[which(runs == min(runs))[1]],1],
                                      AIC=runs[which(runs == min(runs))[1]] )
                    )
  }
  total_runs <- total_runs[!(total_runs %in% focal_runs)]
  cat(paste("[jobs remaining:",length(total_runs),"]",sep=""));
}; cat("\n");
# look over the random walk AICs and model, discuss with friends
optimal <- as.character(minimum$formula[which(minimum$AIC == min(minimum$AIC))[1]])
m_final <- distsamp(as.formula(optimal),umf,keyfun="hazard",output="density",unitsOut="kmsq")
# finish-up
cat(" -- cleaning-up and writing to disk\n")
parallel::stopCluster(cl)
write.csv(minimum, paste(four_letter_code,"_models_selected.csv", sep=""), row.names=F)
save.image(gsub(paste(four_letter_code,"_density_estimation_workspace_",paste(strsplit(as.character(date()),split=" ")[[1]][c(2,3,5)],collapse="_"),".rdata",collapse=""),
                pattern=" ",
                replacement="")
          )
