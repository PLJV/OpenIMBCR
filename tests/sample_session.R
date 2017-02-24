#
# MAIN
#

require(rgdal)
require(raster)
require(landscapeAnalysis)
require(OpenIMBCR)

spp <-
c(
  "Grasshopper sparrow",
  "Lark bunting",
  "Swainson’s Hawk",
  "Cassin's Sparrow",
  "McCown's Longspur",
  "Long-billed curlew",
  "Killdeer",
  "Upland sandpiper",
  "Ring-necked Pheasant",
  "Northern bobwhite",
  "Wild-turkey",
  "Great blue heron"
)

# data frame -> spatial points data frame
s <- OpenIMBCR::imbcrTableToShapefile(filename=OpenIMBCR:::recursiveFindFile(name="RawData_PLJV_IMBCR_20161201.csv",root="/home/ktaylora/Incoming")[1])
# number of unique 1-km2 transects?
M <- length(as.character(unique(s$transectnum))) # number of transects

#
# parse our count observations and metadata into something we can use for est. p-detection
# and occupancy
#

transect_habitat_covs <- OpenIMBCR:::parseHabitatMetadataByTransect(s) # transect-level habitat covariates
        detectionHist <- OpenIMBCR:::parseStationCountsAsOccupancy(OpenIMBCR:::parseStationLevelMetadata(s,spp=spp[1])) # covariates for probability of detection

# testing input : do not run
# pairs(transect_habitat_covs,cex=0.7,pch=15,log=T)
# abs(cor(transect_habitat_covs[,2:ncol(transect_habitat_covs)])) > 0.25

# Optimize with NLM
inputTable <- cbind(transect_habitat_covs[,2:ncol(transect_habitat_covs)],detectionHist)
# re-scale our input variables
inputTable[,!grepl(names(inputTable),pattern="det|obs")] <- scale(inputTable[,!grepl(names(inputTable),pattern="det|obs")])
# make it something that singleSeasonOccupancy can digest
inputTable <- inputTable[,c('det',names(inputTable)[!grepl(names(inputTable),pattern="det")])]
  inputTable <- cbind(det=inputTable[,'det'],inputTable[,n[!grepl(n,pattern="perc_|det")]],inputTable[,n[grepl(n,pattern="perc_")]])
    inputTable <- inputTable[,!grepl(names(inputTable),pattern="obs")]

n <- names(inputTable)

# fit a test model
singleSeasonOccupancy(det_parameters=runif(n=1+length(n[!grepl(n,pattern="det|obs|perc_")])), occ_parameters=runif(n=1+length(n[grepl(n,pattern="perc_")])),table=inputTable)

# use R's built-in optimization to fit an optimal likelihood
m <- nlm(f=singleSeasonOccupancy,p=rep(0,13),
         det_parameters=3,
         occ_parameters=8,
         table=inputTable, hessian=TRUE)

# test combindations of input variables
combinations <- list()
vars <- paste(names(inputTable),c("a0","b0")) # tack our intercept terms on here.  They are handled internally by singleSeasonOccupancy()
for(i in 1:length(names(inputTable))){
  c <- utils::combn(vars,m=i)
  for(j in 1:ncol(c)){
    m <- nlm(f=singleSeasonOccupancy,p=rep(0,i),
             vars=c[,j],
             table=inputTable,hessian=TRUE)
    focal <- list()
      focal[[1]] <- m
        focal[[2]] <- c[,j]
    combinations[[length(combinations)+1]] <- focal
  }
  cat(".")
}; cat("\n");

permuted_aics <- unlist(lapply(combinations,FUN=function(x) AIC(x[[1]])))

which(permuted_aics ==min(permuted_aics))
