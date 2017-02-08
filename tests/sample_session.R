#
# MAIN
#

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
s <- imbcrTableToShapefile(filename=recursiveFindFile(name="RawData_PLJV_IMBCR_20161024.csv",root="/home/ktaylora/Incoming")[1])
# number of unique 1-km2 transects?
M <- length(as.character(unique(s$transectnum))) # number of transects

#
# parse our count observations and metadata into something we can use for est. p-detection
# and occupancy
#

transect_habitat_covs <- parseHabitatMetadataByTransect(s) # transect-level habitat covariates
        detectionHist <- parseStationCountsAsOccupancy(parseStationLevelMetadata(s,spp=spp[1])) # covariates for probability of detection

pairs(transect_habitat_covs,cex=0.7,pch=15,log=T)
abs(cor(transect_habitat_covs[,2:ncol(transect_habitat_covs)])) > 0.25

# Optimize with NLM
inputTable <- cbind(transect_habitat_covs[,2:ncol(transect_habitat_covs)],detectionHist)
# re-scale our input variables
inputTable[,!grepl(names(inputTable),pattern="det|obs")] <- scale(inputTable[,!grepl(names(inputTable),pattern="det|obs")])

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

m <- nlm(f=singleSeasonOccupancy,p=rep(0,10),
         vars=c("a0","tod","doy","intensity","b0","perc_ag","perc_grass", "perc_shrub", "perc_tree","perc_playa"),
         table=inputTable,hessian=TRUE)
