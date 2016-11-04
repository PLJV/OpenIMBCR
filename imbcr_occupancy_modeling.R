#
# Single-season occupancy modeling work for 2016 data
#
# Author: Kyle Taylor
#
# Much of this is based on earlier implementations described by David Pavlacky (Bird Conservancy of the Rockies)
# and by Andy Royle (USFWS).
#

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

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

spp <-
c(start
  'Grasshopper sparrow',
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

s <- imbcrTableToShapefile(filename=recursiveFindFile(name="RawData_PLJV_IMBCR_20161024.csv",root="/home/ktaylora/Incoming")[1])

#
# Estimate naive occupancy for each species
#

      s_spp <- parseDetectionsBySpecies(s,spp=spp[1])
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


#
# Fit a basic model that estimates probability of detection using IMBCR metadata
# across 1-km transects and evaluate its accuracy using an 80/20 k-fold cross-validation
#

M <- 300 # number of sites

# fake some presence/absence data that is strongly correlated with "elevation" to inform
# building our occupancy model
intensity  <- matrix(abs(rnorm(n=M,mean=(M:1),sd=10)/M),ncol=1) # simulate getting kinda worse at sampling as we go, because we are tired
         y <- rbinom(n=(M*6),size=1,prob=abs(((1:M)/M)*(1-intensity))) # there is a clear trend of increase associated with elevation, partially obscured by our sampling effort
         y <- matrix(y,ncol=6) # format as six repeat visits per site

 elevation <- matrix(rnorm(n=M,mean=1:M,sd=15),ncol=1)
elevation2 <- matrix(rnorm(n=M,mean=(1:M)^2,sd=1),ncol=1)


singleScaleOccupancy <- function(parameters,vars=c("a0","intensity","b0","elevation","elevation2")){
  # name of all potential variables
  covarNames <- c("a0","intensity","b0","elevation","elevation2")
  intercept <- rep(1:M)
  # by default, set our coefficients = 0 for this run
  coeffs <- rep(0,length(covarNames))
    names(coeffs) <- covarNames
  # for those vars used in this run, set coefficients to the parameter value given by nlm()
  coeffs[vars] <- parameters
  a0 <- coeffs[1];
  a1 <- coeffs[2];

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
    likelihood[i] <- log(prod(cp)*psi[i] + ifelse(nd==0,1,0)*(1-psi[i])) # http://stats.stackexchange.com/questions/211848/likelihood-why-multiply
  }
  sum(-1*likelihood)
}


#
# Fit a single-season occupancy model that allows for heterogeneity in detection probability
# across transects, using Andy Royle's (2012) model specification.
#
