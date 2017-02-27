#
# Density (birds/km2) estimation for Ring-necked Pheasant using the 2016 IMBCR data
#

require(OpenIMBCR)
require(unmarked)
require(rgdal)
require(raster)

# define our distance breaks through exploratory analysis
d=c(0,100,200,300,400,500,600,700,800)

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
  s <- suppressWarnings(raster::extract(raster::raster("/global_workspace/aquifer_depletion/elevation.tif"),s,sp=T))
    s$doy <- as.numeric(strftime(as.POSIXct(as.Date(as.character(s$date), "%m/%d/%Y")),format="%j")) # convert date -> doy
      colnames(s@data) <- gsub(names(s@data),pattern="dem_",replacement="elev_")

# kludging to select the covariates we will aggregate and use at the site level
n <- names(s@data)
  n <- n[(length(n)-(nlayers(r)+1)) : length(n) ]
    n <- append(n,c("starttime")) # append our covariates on detection

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
m <- distsamp(~doy+timeofday~1,umf,keyfun="hazard",output="density",unitsOut="kmsq")
m <- distsamp(~doy~1,umf,keyfun="hazard",output="density",unitsOut="kmsq") # use AIC to determine an optimal detection function
# clean-up our state variables so we don't use our det covs estimate density
n <- n[!grepl(n,pattern="time|doy")]
# calculate : combinations w/o repetition (n!/(r!(n-r)!)
m_len <- vector(); for(i in 1:length(n)) { m_len <- append(m_len,dim(combn(n,m=i))[2]) }
 m_len <- sum(m_len)
# minimize AIC across m_len possible models
models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
k <- 1
cat(" -- optimizing AIC:\n")
for(i in 1:length(n)){
 combinations <- combn(n,m=i)
 cat(paste(" -- n=", i, " ways:", collapse=""))
 for(j in 1:ncol(combinations)){
   f <- as.formula(paste("~doy~",paste(combinations[,j],collapse="+"),collapse=""))
   m_2 <- distsamp(f,umf,keyfun="hazard",output="density",unitsOut="kmsq")
   models[k,1] <- paste("~doy~",paste(combinations[,j],collapse="+"),collapse="")
   models[k,2] <- m_2@AIC
   k <- k+1
   cat(paste("[",round(k/nrow(models),2),"%]",sep=""))
 }; cat("\n")
}
# check for a range of AIC support in our combinations
if(sum(abs(models$AIC-min(models$AIC)) < 5)>1){
  cat(" -- warning: support for more than 1 optimal model\n")
}

optimal <- models$formula[which(models$AIC == min(models$AIC))]
m_final <- distsamp(as.formula("~1~tree_107x107+wetland_107x107"),umf,keyfun="hazard",output="density",unitsOut="kmsq")
 m_2 <- distsamp(as.formula(paste("~1~",paste(n,collapse="+"),collapse="")),umf,keyfun="hazard",output="density",unitsOut="kmsq")

save.image(paste("RIPH_density_estimation_workspace_",paste(strsplit(as.character(date()),split=" ")[[1]][c(2,3,5)],collapse="_"),".rdata",collapse=""))
