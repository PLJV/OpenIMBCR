#
# OpenIMBCR
#
# Author: Kyle Taylor [2016] (kyle.taylor@pljv.org)
#
# Much of this is based on earlier discussions and implementations described by David Pavlacky (Bird Conservancy of the Rockies)
# Rob Sparks (Bird Conservancy of the Rockies), and a lot of work by Andy Royle (USFWS).
#
# Please report bugs to Kyle Taylor
#

#' expit link function
#' @export
expit <- function(x) 1/(1+exp(-x))
#' logit link function
#' @export
logit <- function(x) log(x/(1-x))
#' calculate AIC from a likelihood function optimization procedure
#' @export
AIC <- function(m){
  if(as.character(class(m))=="unmarkedFitDS"){
    return(m@AIC)
  }
  # if it isn't an unmarked dataframe, assume we rolled our own optimization with nlm()
  return((m$minimum*2) + (2*length(m$estimate)))
}
#' calculate AICc from a likelihood function optimization procedure
#' @export
AICc <- function(m) AIC(m) + ((2*length(m$estimate)*(length(m$estimate)+1))/(m$size-length(m$estimate)-1))
#' calculated BIC from a likelihood function optimization procedure
#' @export
BIC <- function(m) (m$minimum*2) + (q*log(m$size))
#'
#' calculate SE from a likelihood function optimization procedure
#' @export
SE <- function(m) sqrt(diag(solve(m$hessian)))
#' recursively calculate all possible permutations of an input table n
permutations <- function(n){
   if(n==1){
       return(matrix(1))
   } else {
       sp <- permutations(n-1)
       p <- nrow(sp)
       A <- matrix(nrow=n*p,ncol=n)
       for(i in 1:n){
           A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
       }
       return(A)
   }
}
#' strip non-essential characters and spaces from a species common name in a field
stripCommonName <- function(x) tolower(gsub(x,pattern=" |'|-",replacement=""))
#' recursively find all files in folder (root) that match the pattern (name)
#' @export
recursiveFindFile <- function(name=NULL,root=Sys.getenv("HOME")){
  if(is.null(name)){
    return(NULL)
  } else {
    return(list.files(root,pattern=name,recursive=T,full.names=T))
  }
}
#' parse a source CSV (published by BCR) for IMBCR data and return as a SpatialPointsDataFrame
#' to the user
#' @export
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
      s[[length(s)]] <- sp::SpatialPointsDataFrame(
                                  coords = data.frame(x = s[[length(s)]]$ptvisiteasting,
                                  y = s[[length(s)]]$ptvisitnorthing),
                                  data = s[[length(s)]],
                                  proj4string = sp::CRS(raster::projection(paste("+init=epsg:269", zone, sep = "")))
                        )
      row.names(s[[length(s)]]) <- paste(letters[length(s)],row.names(s[[length(s)]]),sep=".") # based on : http://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r
    }

    s <- lapply(s,FUN=sp::spTransform,sp::CRS(raster::projection(s[[1]])))
      s <- do.call(rbind,s)
        s$FID <- 1:nrow(s)
    # write to disk -- and allow some wiggle-room on filename conventions
    if(write){
      rgdal::writeOGR(s,".",ifelse(is.null(outfile),gsub(filename,pattern=".csv",replacement=""),outfile),driver="ESRI Shapefile",overwrite=T)
    }
  }
  return(s)
}
#' fetch IMBCR Metadata to Landfire code conversions
#' @export
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

  detectionHist <- list()
  for(t in as.character(unique(s_spp$transectnum))){

        counts <- s_spp[s_spp$transectnum == t,]$cl_count
           det <- as.numeric(counts > 0)
      station  <- s_spp[s_spp$transectnum == t,]$point
      interval <- s_spp[s_spp$transectnum == t,]$timeperiod
          dist <- s_spp[s_spp$transectnum == t,]$radialdistance
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
                                                      dist=dist,
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
  for(i in 1:length(detectionHist)){
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
    intensity <- length(unique(detectionHist[[i]]$station))
    detectionHist[[i]] <- data.frame(det=det,tod=tod,doy=doy,obs=obs,intensity=intensity)
  }
  do.call(rbind,detectionHist)
}
#' generate plant community composition data from the associated IMBCR metadata
#'@export
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
#' extract (and optionally, summarize using the fun= argument) raster data across IMBCR transects
#' @export
extractByTransect <- function(s=NULL,r=NULL,fun=NULL){
  # sanity-check
  if(sum(is.null(list(s,r)))>0) stop("s= and r= arguments can't be null")
  # ensure consistent projections
  orig_crs <- sp::CRS(raster::projection(s))
  s <- sp::spTransform(s,sp::CRS(raster::projection(r)))
  # search for a meaningful transect field
  transect_field <- names(s)
    transect_field <- transect_field[grepl(transect_field, pattern="trns|transect")][1] # by convention, use the first field that has "transect" or "trns" in it
      if(length(transect_field) == 0) stop("couldn't find a meaningful IMBCR transect field in object s=")
  # iterate over our transect data, processing as we go
  transects <- unique(as.character(s@data[,transect_field]))
  for(i in 1:length(transects)){
    focal <- s@data[,transect_field] == transects[i]
      focal <- s[focal,]
        focal <- focal[!duplicated(focal$point),]
    # extract and assign our habitat values
    rVal <- extract(r,focal)
    # assign a summary statistic, specified by fun=function()
    if(is.function(fun)){
      s@data[s@data[,transect_field] == as.vector(unique(focal@data[,transect_field])),names(r)[1]] <- fun(rVal)
    }
    # assign our values point-by-point
    else {
      for(j in 1:length(focal$point)){
        s@data[s@data[,'point'] == focal$point[j],names(r)[1]] <-  rVal[j]
      }
    }
  }
  # restore our original CRS
  s <- sp::spTransform(s,orig_crs)
  return(s)
}
#' validate transect-level habitat metadata vs. LANDFIRE
#'
#' @export
validateTransectMetadata <- function(s){
  focal <- s[s$transectnum == transect_habitat_covs$transect[1],][!duplicated(s[s$transectnum == transect_habitat_covs$transect[1],]$point),]
}
#' Fit a single-season occupancy model  assumes a constant probability of species
#' detection across transects, which is probably inappropriate for IMBCR data and will
#' lead to inaccurate predictions of occupancy. This is "model m0" from the literature
#' and loosely follows Andy Royle's (2008) model specification. It is designed so that
#' parameter estimates can be derived through an optimization procedure like nlm() and
#' variables can be easily selected/de-selected, such as for calculating AIC.
#'
#' @param table a data.frame of detection histories ('det') and covariates on
#' occupancy ('b0') and probability of detection ('a0')
#' @param det_parameters proposed parameter values for detection, including 1 intercept (a0)
#' @param occ_parameters proposed paramter values for occupancy, including 1 intercept (b0)
#'
#' @export
singleSeasonOccupancy <- function(det_parameters=NULL, occ_parameters=NULL, p=NULL, table=NULL){
  if(is.null(table)){
    stop("table= argument requires an input table containing covariates and a 'det' field defining detection histories for your sites.")
  }
  # assign a vars vector containing a composite of our det_parameters and occ_parameters
  if(is.null(p)){
    vars <- c("a0", colnames(table[,2:(length(det_parameters)+1)]), "b0", colnames(table[,(length(det_parameters)+2):ncol(table)]))
  } else {
    vars <- c("a0",colnames(table[,2:(det_parameters+1)]), "b0", colnames(table[,(det_parameters+2):ncol(table)]))
  }
  # The p-parameter is typically assigned by nlm() and should include proposed values for our two intercept terms (in addition to our covariates)
  if(!is.null(p)){
    det_parameters = p[1:(det_parameters+1)]
    occ_parameters = p[(length(det_parameters)+1):length(p)]
  }
  if(length(det_parameters)+length(occ_parameters) > ncol(table)+2-1 ){ # accounting for intercept terms, but ignoring the first column (detections)
    stop("length of paramters is greater than the number of columns specified by table=")
  }
  # define the number of transects, number of stations per transect, and target variables we are considering for this iteration
            M <- nrow(table) # number of sites (IMBCR transects)
    nStations <- nchar(as.character(table[1,1])) # number of stations in each transect (should be 16)
   covarNames <- vars
  # re-build a consistent table (t) of detection histories, intercepts, and covariate data we can work with
  y <- table[,1] # assume our 'detection' field is always the first column
    y <- suppressWarnings(matrix(as.numeric(matrix(unlist(strsplit(as.character(y),split="")))),nrow=M,ncol=nStations))
  t <- matrix(rep(0,M*length(covarNames)),ncol=length(covarNames)) # build a table of zeros for our covariate data
    colnames(t) <- covarNames
     t[,"a0"] <- rep(1,M) # fill our intercept columns
     t[,"b0"] <- rep(1,M)
  # assign values for our focal table from user-specific source table
  t <- data.frame(t)
  table <- data.frame(table)
  t[,vars[!grepl(vars,pattern=0)]] <- table[,vars[!grepl(vars,pattern=0)]] # assign all covariate data (minus intercept data) from source table
  # by default, set our coefficients = 0 for this run
  coeffs <- rep(0,length(covarNames))
    names(coeffs) <- covarNames
  # assign proposed parameters for this optimization step, from nlm()
  tryCatch(coeffs[vars] <- c(det_parameters,occ_parameters),
           warning=function(w){
             stop("caught a warning assigning detection + occurrence parameters + intercept terms. Did you remember to generate proposed paramter values for the intercept terms?")
           })
  detection_coeffs <- 1:(which(grepl(vars,pattern="b0$"))-1)
  occupancy_coeffs <- which(grepl(vars,pattern="b0$")):length(vars)
  # parameters on detection
  # matrix operation : prob <- expit(a0*t[,'a0'] + a1*t[,"tod"] + a2*t[,"doy"] + a3*t[,"intensity"])
  prob <- sweep(t[,vars[detection_coeffs]], MARGIN=2, coeffs[detection_coeffs],`*`)
    prob <- expit(rowSums(prob))
  # matrix operation : psi <- expit(b0*t[,'b0'] + b1*t[,"perc_ag"] + b2*t[,"perc_grass"] + b3*t[,"perc_shrub"] + b4*t[,"perc_tree"] + b5*t[,"perc_playa"])
  psi <- sweep(t[,vars[occupancy_coeffs]], MARGIN=2, coeffs[occupancy_coeffs],`*`)
     psi <- expit(rowSums(psi))
  # solve for likelihood
  likelihood <- rep(NA,M)
  for(i in 1:M){
    detections <- y[i,] # individual detections for our repeat visits at site i
    na_det <- is.na(detections) # any NA values?
    nd <- sum(detections[!na_det]) # check for zero-detections
    p <- prob[i] # what is the predicted probability of detections for this site, given our parameters?
    # calculate likelihood of occupancy, given our calculated probability of detection for the focal transect
    cp <- (p^detections)*((1-p)^(1-detections))
      cp[na_det] <- 1 # set any NA values to 1
    likelihood[i] <- log(prod(cp)*psi[i] + ifelse(nd==0,1,0)*(1-psi[i])) # joint probability across detections, e.g.: http://stats.stackexchange.com/questions/211848/likelihood-why-multiply
  }
  sum(-1*likelihood)
}
#' fit the single-season, two-state occupancy model of Royle-Nichols (the "RN" model) to IMBCR data. This model can account for abundance-induced heterogeneity of detection and capitalize
#' on replicate observations of abundance within sampling units. Royle first implemented the algorithm for use with BBS data, but we can theoretically extend the model
#' within-season repeat observations made at IMBCR stations. Requires the user specify an upper_bound parameter for density within each 1-km transect. Like the single-season occupancy model (model m0), the RN model
#' does not allow for heterogeneity in detection within sites, which may result in biased estimates of detection in instances where an observer detects a species early in the sampling process.
#' Can be accomodated with a removal design.
#'
#' @export
singleSeasonRN <- function(det_parameters=NULL, occ_parameters=NULL, p=NULL, upper_bound=NULL, quasi_binom=FALSE, table=NULL){
  if(is.null(table)){
    stop("table= argument requires an input table containing covariates and a 'det' field defining detection histories for your sites.")
  }
  if(is.null(upper_bound)){
    stop("upper_bound= parameter (e.g., a rational prediction of an absolute maximum bird density per IMBCR transect) must be defined")
  }
  # assign a vars vector containing a composite of our det_parameters and occ_parameters
  if(is.null(p)){
    vars <- c("a0", colnames(table[,2:(length(det_parameters)+1)]), "b0", colnames(table[,(length(det_parameters)+2):ncol(table)]))
  } else {
    vars <- c("a0",colnames(table[,2:(det_parameters+1)]), "b0", colnames(table[,(det_parameters+2):ncol(table)]))
  }
  # The p-parameter is typically assigned by nlm() and should include proposed values for our two intercept terms (in addition to our covariates)
  if(!is.null(p)){
    det_parameters = p[1:(det_parameters+1)]
    occ_parameters = p[(length(det_parameters)+1):length(p)]
  }
  if(length(det_parameters)+length(occ_parameters) > ncol(table)+2-1 ){ # accounting for intercept terms, but ignoring the first column (detections)
    stop("length of paramters is greater than the number of columns specified by table=")
  }
  # define the number of transects, number of stations per transect, and target variables we are considering for this iteration
            M <- nrow(table) # number of sites (IMBCR transects)
    nStations <- nchar(as.character(table[1,1])) # number of stations in each transect (should be 16)
   covarNames <- vars
  # re-build a consistent table (t) of detection histories, intercepts, and covariate data we can work with
  y <- table[,1] # assume our 'detection' field is always the first column
    y <- suppressWarnings(matrix(as.numeric(matrix(unlist(strsplit(as.character(y),split="")))),nrow=M,ncol=nStations))
  t <- matrix(rep(0,M*length(covarNames)),ncol=length(covarNames)) # build a table of zeros for our covariate data
    colnames(t) <- covarNames
     t[,"a0"] <- rep(1,M) # fill our intercept columns
     t[,"b0"] <- rep(1,M)
  # assign values for our focal table from user-specific source table
  t <- data.frame(t)
  table <- data.frame(table)
  t[,vars[!grepl(vars,pattern=0)]] <- table[,vars[!grepl(vars,pattern=0)]] # assign all covariate data (minus intercept data) from source table
  # by default, set our coefficients = 0 for this run
  coeffs <- rep(0,length(covarNames))
    names(coeffs) <- covarNames
  # assign proposed parameters for this optimization step, from nlm()
  tryCatch(coeffs[vars] <- c(det_parameters,occ_parameters),
           warning=function(w){
             stop("caught a warning assigning detection + occurrence parameters + intercept terms. Did you remember to generate proposed paramter values for the intercept terms?")
           })
  detection_coeffs <- 1:(which(grepl(vars,pattern="b0$"))-1)
  occupancy_coeffs <- which(grepl(vars,pattern="b0$")):length(vars)
  # parameters on detection (the rate parameter of our poisson)
  # matrix operation : prob <- expit(a0*t[,'a0'] + a1*t[,"tod"] + a2*t[,"doy"] + a3*t[,"intensity"])
  r <- sweep(t[,vars[detection_coeffs]], MARGIN=2, coeffs[detection_coeffs],`*`)
    r <- expit(rowSums(r))
  # matrix operation : psi <- expit(b0*t[,'b0'] + b1*t[,"perc_ag"] + b2*t[,"perc_grass"] + b3*t[,"perc_shrub"] + b4*t[,"perc_tree"] + b5*t[,"perc_playa"])
  lambda <- sweep(t[,vars[occupancy_coeffs]], MARGIN=2, coeffs[occupancy_coeffs],`*`)
     psi <- exp(rowSums(lambda))
  # solve for likelihood
  likelihood <- rep(NA,M)
  for(i in 1:M){
    # determine a truncated (by upper_bound) probability for our Poisson count,
    gN <- dpois(0:upper_bound, lambda[i])
      gN <- gN/sum(gN)
    # individual detections for our repeat visits at site i (e.g., [1,1,1,3,2,1])
    detections <- y[i,]
    na_det <- is.na(detections) # any NA values?

    pmat <- 1 - outer((1-r[i,]),0:upper_bound,"^")

    focal <- t((pmat^detections)*(1-pmat)^(1-detections))
      focal[,na_det] <- 1
        focal <- apply(focal,1,prod)

    likelihood[i] <- sum(focal*gN)
 }
  sum(-1*likelihood)
}
#' fit a distance model to IMBCR station data (testing)
#'
singleSeasonDistance <- function(x=NULL, psi_params=NULL, sigma_params=NULL, as_radial_distance=T, link="half_normal"){
  # detection functions
      uniform <- function(u){ 1/u }
  half_normal <- function(u){ (2*pi*u) * (exp(-(u^2)/sigma2)) / (pi*x_max^2)  }
  hazard_rate <- function(u){ 1 - exp( -(u/a)^(-b) ) }
  nz<-500 # transects where species was not observed
  sigma2 <- exp(sigma_params[2])
  if(as_radial_distance){
    # x = r*sin(theta)
    x <- x * sin( (angle/360)*(2*pi) )
    x <- x/100 # format units as meters/100 for our detection function
  }
  x_max <- ceiling(max(x))
  nind<-length(x)
  y<-c(rep(1,nind),rep(0,nz))
  x<-c(x,rep(NA,nz)) # NA fill for our 0 capture

  lik <- function(parms){
    psi <- expit(parms[1])
    picap <- integrate(half_normal,0,x_max)$value # marginal probability of encounter
    #part1 <- sum(log(psi*exp(-(x[1:nind]^2)/sigma2) ) )
    part1 <- sum(log(half_normal(x[1:nind]))) # half-normal distance detection function (modified for point count survey)
    #part2 <-  nz*log( 1-psi*picap)
    part2 <- nz*log(1-picap)
    part3 <- 0
    -1*(part1+part2+part3) # log of a product is the sum of the log factors
  }

  out<-nlm(lik,c(logit(175/(nz+nind)),log(1.2) ),hessian=TRUE)
  print(out)

  psihat<- expit(out$estimate[1])
  N <-psihat*( nind+nz )
  D<- N/48
  cat("MLE Density: ",D,fill=TRUE)
}
