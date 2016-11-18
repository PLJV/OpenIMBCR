#
# Single-season occupancy modeling work for 2016 data
#
# Author: Kyle Taylor [2016]
#
# Much of this is based on earlier implementations described by David Pavlacky (Bird Conservancy of the Rockies)
# and by Andy Royle (USFWS).
#

require(raster)
require(rgdal)

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

AIC <- function(m) (m$minimum*2) + (2*length(m))
 SE <- function(m) sqrt(diag(solve(m$hessian)))

stripCommonName <- function(x) tolower(gsub(x,pattern=" |'|-",replacement=""))

recursiveFindFile <- function(name=NULL,root=Sys.getenv("HOME")){
  if(is.null(name)){
    return(NULL)
  } else {
    return(list.files(root,pattern=name,recursive=T,full.names=T))
  }
}

imbcrTableToShapefile <- function(filename=NULL,outfile=NULL,write=F){
  if(is.null(outfile) && write && !is.character(filename)){
    stop("cannot write shapefile to disk without an outfile= or filename.ext to parse")
  }
  # sanity-check to see if our output shapefile already exists
  s <- recursiveFindFile(outfile)[1]
  if(!is.null(s)){
    s <- landscapeAnalysis:::readOGRfromPath(s)
  } else {
    #
    # Parse ALL BCR raw data tables into a single table
    #
    if(class(filename) == "data.frame"){
      t <- filename
    } else {
      if(file.exists(filename)){
        t <- read.csv(filename)
      } else {
        t <- recursiveFindFile(name=filename)
          t <- lapply(t,read.csv)
            t <- do.call(what=rbind,t)
      }
    }
    names(t) <- tolower(names(t))
    #
    # iterate over each UTM zone in the table, creating SpatialPoints
    # projected to a focal UTM.  Then merge all of the zones together into
    # a single shapfile with an arbitrary CRS.
    #
    s <- list()
    for (zone in unique(na.omit(t$ptvisitzone))){
      s[[length(s) + 1]] <- na.omit(t[t$ptvisitzone == zone,])
      s[[length(s)]] <- SpatialPointsDataFrame(
                                  coords = data.frame(x = s[[length(s)]]$ptvisiteasting,
                                  y = s[[length(s)]]$ptvisitnorthing),
                                  data = s[[length(s)]],
                                  proj4string = CRS(projection(paste("+init=epsg:269", zone, sep = "")))
                        )
    }

    s <- lapply(s,FUN=spTransform,CRS(projection(s[[1]])))
      s <- suppressMessages(do.call(rbind,s))
        s$FID <- 1:nrow(s)
    # write to disk -- and allow some wiggle-room on filename conventions
    if(write){
      writeOGR(s,".",ifelse(is.null(outfile),gsub(filename,pattern=".csv",replacement=""),outfile),driver="ESRI Shapefile",overwrite=T)
    }
  }
  return(s)
}
#' fetch IMBCR Metadata to Landfire code conversions
#'
fetchImbcrToLandfireMetadata <-function(url="https://docs.google.com/spreadsheets/d/1wJ3Xwr67GTYYfim29cKQ9m3AJgoOQryYTHc9v5gJB2g/pub?gid=0&single=true&output=csv"){
  download.file(url,destfile="imbcr_meta_codes.csv",quiet=T);
    t <- read.csv("imbcr_meta_codes.csv")
      names(t) <- tolower(names(t))
  # clean-up and return
  file.remove("imbcr_meta_codes.csv")
  return(t[,c(1,2)])
}

#' accepts a SpatialPointsDataFrame of IMBCR data and parses observations into counts for a focal
#' species specified by the user.
#' @param spp character string specifying focal species common name
#' @export
parseDetectionsBySpecies <- function(s,spp=NULL){
  in_sample <- which(stripCommonName(s$common.name) %in% stripCommonName(spp[1]))
  out_sample <- which(!(1:length(s) %in% in_sample))

  t <- s@data
    t$common.name <- as.character(t$common.name)

  t[out_sample,]$common.name <- spp[1]
  t[out_sample,]$cl_count <- 0

  s@data <- t
  return(s)
}
#' parse relavent metadata for estimating p-detection and leave station-level count information in-place
#' return a list of M data.frames (one for each IMBCR transect) that can be post-processed later for occupancy or abundance
#' modeling.
#' @export
parseStationLevelMetadata <- function(s,spp=NULL){


        s_spp <- parseDetectionsBySpecies(s,spp=spp) # parse our SpatialPointsDataFrame for focal species (spp)
  n_occupancy <- length(unique(s_spp[s_spp$cl_count > 0,]$transectnum))/length(unique(s$transectnum))

  detectionHist <- list()
  for(t in as.character(unique(s_spp$transectnum))){

      counts   <- s_spp[s_spp$transectnum == t,]$cl_count
           det <- as.numeric(counts > 0)
      station  <- s_spp[s_spp$transectnum == t,]$point
      interval <- s_spp[s_spp$transectnum == t,]$timeperiod
           doy <- as.numeric(strftime(as.POSIXct(as.Date(as.character(s_spp[s_spp$transectnum == t,]$date), "%m/%d/%Y")),format="%j"))
           tod <- as.numeric(s_spp[s_spp$transectnum == t,]$starttime)
           obs <- as.character(s_spp[s_spp$transectnum == t,]$observer)

      distance_disqualifier <- s_spp[s_spp$transectnum == t,]$radialdistance
        distance_disqualifier <- (distance_disqualifier < 0 | distance_disqualifier > 400)

        detectionHist[[length(detectionHist)+1]] <- data.frame(
                                                      counts=counts,
                                                      detection=det,
                                                      station=station,
                                                      interval=interval,
                                                      doy=doy,
                                                      tod=tod,
                                                      obs=obs,
                                                      disqualified=distance_disqualifier
                                                    )
  }
  # return our list of data.frames to user
  return(detectionHist)
}
#' accepts a list of data.frames (as returned by parseStationLevelMetadata) and post-processes
#' the counts into binomial detections that can be used to fit occupancy.
#' @export
parseStationCountsAsOccupancy <- function(detectionHist,na.rm=F){
  # tod (HH:MM) -> Minutes since midnight
  hhmm_to_min <- function(x){
    (as.numeric(substr(as.character(x),1,nchar(as.character(x))-2))*60) + (as.numeric(substr(as.character(x),2,nchar(as.character(x)))))
  }
  for(i in 1:M){
    d <- as.numeric(aggregate(counts~station,detectionHist[[i]],function(x){sum(x>0)})$counts > 0)  # did we observe ANY birds across our 6 minute count for each station?
    if(!na.rm){
      det <- rep(".",16)
      det[which(as.numeric(unique(detectionHist[[i]]$station)) %in% 1:16)] <- d
        det <- paste(det,collapse="")
    } else {
      det <- d;
    }
    tod <- hhmm_to_min(round(median(detectionHist[[i]]$tod,na.rm=T)))
    doy <- round(median(detectionHist[[i]]$doy,na.rm=T))
    obs <- landscapeAnalysis:::Mode(detectionHist[[i]]$obs)
    intensity <- max(detectionHist[[i]]$station)
    detectionHist[[i]] <- data.frame(det=det,tod=tod,doy=doy,obs=obs,intensity=intensity)
  }
  do.call(rbind,detectionHist)
}
#'
#'
parseHabitatMetadataByTransect <- function(s){
  # fetch our source table of habitat associations
  t <- fetchImbcrToLandfireMetadata()
  # aggregate by IMBCR habitat codes
         AG  <- as.vector(t[t$to == "AG",1])
       TREES <- as.vector(t[t$to == "TR",1])
    WETLANDS <- as.vector(t[t$to == "WE",1])
  SHRUBLANDS <- as.vector(t[t$to == "SH",1])
       GRASS <- as.vector(t[t$to == "GR",1])
       OTHER <- as.vector(t[t$to == "XX",1])
  # build a data.frame and return to user
  transects <- unique(s$transectnum)
   habitat <- list()
  for(i in 1:length(transects)){
    # percent cover of our 1-km2 transect
    habSummary <- s[s$transectnum == transects[i],]@data$primaryhab
      habSummary <- as.data.frame(table(habSummary)*100/(length(habSummary)))
        names(habSummary) <- c("habitat","perc_cover")

    habSummary <-
    data.frame(perc_ag=sum(habSummary[habSummary$habitat %in% AG, 2]),
               perc_grass=sum(habSummary[habSummary$habitat %in% GRASS, 2]),
               perc_shrub=sum(habSummary[habSummary$habitat %in% SHRUBLANDS, 2]),
               perc_tree=sum(habSummary[habSummary$habitat %in% TREES,2]),
               perc_playa=habSummary[habSummary$habitat == "PL", 2],
               perc_wetland=sum(habSummary[habSummary$habitat %in% WETLANDS,2]),
               perc_urban=habSummary[habSummary$habitat == "UR", 2],
               perc_other=sum(habSummary[habSummary$habitat %in% OTHER,2])
               )
     # bind to a table
     habitat[[i]] <- cbind(transect=transects[[i]],habSummary)
     if(sum(habitat[[i]][1,2:ncol(habitat[[i]])]) < 100){
       warning(paste("transect:",transects[[i]]," % cover only adds up to ",round(sum(habitat[[i]][1,2:ncol(habitat[[i]])])),sep=""))
     }
  }
  do.call(rbind,habitat)
}
#' validate transect-level habitat metadata vs. LANDFIRE
#'
validateTransectMetadata <- function(s){
  focal <- s[s$transectnum == transect_habitat_covs$transect[1],][!duplicated(s[s$transectnum == transect_habitat_covs$transect[1],]$point),]
}
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

#
# Fit a basic model that estimates probability of detection using IMBCR metadata
# across 1-km transects and evaluate its accuracy using an 80/20 k-fold cross-validation
#

# fake some presence/absence data that is strongly correlated with "elevation" to inform
# building/testing our occupancy model
# intensity  <- matrix(abs(rnorm(n=M,mean=(M:1),sd=10)/M),ncol=1) # simulate getting kinda worse at sampling as we go, because we are tired
#          y <- rbinom(n=(M*6),size=1,prob=abs(((1:M)/M)*(1-intensity))) # there is a clear trend of increase associated with elevation, partially obscured by our sampling effort
#          y <- matrix(y,ncol=6) # format as six repeat visits per site
#
#  elevation <- matrix(rnorm(n=M,mean=1:M,sd=15),ncol=1)
# elevation2 <- matrix(rnorm(n=M,mean=(1:M)^2,sd=1),ncol=1)

#
# Fit a single-season occupancy model that allows for heterogeneity in detection probability
# across transects, using Andy Royle's (2012) model specification.
#

singleSeasonOccupancy <- function(parameters,vars=c("a0","intensity","b0","elevation","elevation2")){
  # name of all potential variables
  covarNames <- c("a0","intensity","b0","elevation","elevation2")
  intercept <- rep(1:M)
  # by default, set our coefficients = 0 for this run
  coeffs <- rep(0,length(covarNames))
    names(coeffs) <- covarNames
  # for those vars used in this run, set coefficients to the parameter value given by nlm()
  coeffs[vars] <- parameters
  # parameters on detection
  a0 <- coeffs[1];
  a1 <- coeffs[2];
  # parameters on occupancy
  b0 <- coeffs[3];
  b1 <- coeffs[4];
  b2 <- coeffs[5];
  # prediction for
  prob <- expit(a0*intercept + a1*intensity)
   psi <- expit(b0*intercept + b1*elevation + b2*elevation2)
  # solve for likelihood
  likelihood <- rep(NA,M)
  for(i in 1:M){
    detections <- y[i,] # individual detections for our repeat visits at site i
    na_det <- is.na(detections) # any NA values?
    nd <- sum(detections[!na_det]) # check for zero-detections
    p <- prob[i,] # what is the predicted probability of detections for this site, given our parameters?
    # calculate likelihood of occupancy, given our calculated probability of detection for the focal transect
    cp <- (p^detections)*((1-p)^(1-detections))
      cp[na_det] <- 1 # set any NA values to 1
    likelihood[i] <- log(prod(cp)*psi[i] + ifelse(nd==0,1,0)*(1-psi[i])) # joint probability across detections, e.g.: http://stats.stackexchange.com/questions/211848/likelihood-why-multiply
  }
  sum(-1*likelihood)
}

# Optimize with NLM
m <- nlm(singleSeasonOccupancy,c(0,0,0,0,0),vars=c("a0","intensity","b0","elevation","elevation2"),hessian=TRUE)
