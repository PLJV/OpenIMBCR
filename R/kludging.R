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
  require(raster)
  require(rgdal)
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
      if(length(filename)==1){
        t <- read.csv(filename)
      # legacy support : rbind across multiple CSV files
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
      s <- do.call(sp::rbind.SpatialPointsDataFrame,s)
        s$FID <- 1:nrow(s)
    # write to disk -- and allow some wiggle-room on filename conventions
    if(write){
      rgdal::writeOGR(s,".",ifelse(is.null(outfile),gsub(filename,pattern=".csv",replacement=""),outfile),driver="ESRI Shapefile",overwrite=T)
    }
  }
  return(s)
}
#' fetch IMBCR Metadata to Landfire code conversions
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
#' @param attribute_all boolean flag indicating whether to take transect-wide metadata and attribute to points
#' @export
parseHabitatMetadataByTransect <- function(s, attribute_all=T){
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
     # merge our output and QA check our composition estimates
     habitat[[i]] <- cbind(transect=transects[[i]],habSummary)
     if(sum(habitat[[i]][1,2:ncol(habitat[[i]])]) < 100){
       warning(paste("transect:",transects[[i]]," % cover only adds up to ",
               round(sum(habitat[[i]][1,2:ncol(habitat[[i]])])),sep=""))
     }
  }
  # bind our summary statistics by transect identifier
  t_habitat_metadata <- do.call(rbind,habitat)
    n <- gsub(names(t_habitat_metadata),pattern="transect",replacement="transectnum")
      names(t_habitat_metadata) <- n
  if(attribute_all){
    # generalize our summary statistics across all point observations
    s@data <- merge(s@data,t_habitat_metadata,by="transectnum")
    return(s)
  } else {
    return(t_habitat_metadata)
  }
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
#' Re-build an IMBCR dataframe for a focal species so that all transect / stations
#' are represented with their distances or NA values so that the observations
#' are comprehensible by unmarked. This function is currently in testing.
#' @export
rebuildImbcrTable <- function(s=NULL, spp=NULL){
  # local function builds an NA-filled template of what a full
  # transect should look like
  make_transect_template <- function(x=NULL,spp=NULL,years=NULL){
    return(data.frame(
      common.name=spp,
      transectnum=x,
      radialdistance=rep(NA,6*16),
      point=unlist(lapply(1:16,FUN=function(x) rep(x,6))),
      timeperiod=rep(1:6,16),
      year=unlist(lapply(years,function(x){rep(x,6*16)}))
    ))
  }
  colnames(s@data) <- tolower(colnames(s@data))
  target_table <- do.call(
    rbind,
    lapply(
        as.character(unique(s$transectnum)),
        FUN=make_transect_template,spp=spp,years=unique(s$year)
      )
  )
  # merge produces x/y duplicates here, with x occuring first
  target_table <- merge(
    x=s@data[tolower(s@data$common.name) == tolower(spp),],
    y=target_table,
    by=c("transectnum","point","timeperiod","year","common.name","radialdistance"),
    all=T
  )
  target_table <- target_table[!duplicated(
    target_table[,c("transectnum","point","timeperiod","year")]),]
  return(target_table)
}
#' validate transect-level habitat metadata vs. LANDFIRE
#'
#' @export
validateTransectMetadata <- function(s){
  focal <- s[s$transectnum == transect_habitat_covs$transect[1],][!duplicated(s[s$transectnum == transect_habitat_covs$transect[1],]$point),]
}
#' accepts a named raster stack of covariates, an IMBCR SpatialPointsDataFrame,
#' and a species common name and returns a formatted unmarked distance data.frame
#' that can be used for model fitting with unmarked.
#' @export
buildUnmarkedDistanceDf <- function(r=NULL, s=NULL, spp=NULL,
                                    vars=c("doy","starttime"), #
                                    fun=mean,
                                    d=c(0,100,200,300,400,500,600,700,800)){
  # do our covariates in r=raster stack occur in our IMBCR data.frame object?
  if(sum(names(r) %in% names(s@data))<raster::nlayers(r)){
    s <- suppressWarnings(raster::extract(r,s,sp=T))
      s$doy <- as.numeric(strftime(as.POSIXct(as.Date(as.character(s$date), "%m/%d/%Y")),format="%j")) # convert date -> doy
        s@data <- s@data[,!grepl(names(s@data),pattern="FID")]
  }
  # kludging to select the covariates specified in s= that we will aggregate
  # and use at the transect level
  if(!is.null(vars)){
    vars <- append(names(r), vars)
  } else {
    vars <- names(r)
  }
  # parse our dataset for RNEP records
  t <- s[s$common.name == spp,]@data
  # build a distance table
  distances <- data.frame(distance=t$radialdistance,transect=t$transectnum)
  y <- unmarked::formatDistData(distances, distCol="distance",transectNameCol="transect",dist.breaks=d)
  # build a target matrix
  stateCovariates <- matrix(NA,ncol=length(vars),nrow=length(levels(t$transectnum))) # e.g., 300 total transects x n state covariates
    rownames(stateCovariates) <- levels(t$transectnum)
  # aggregate by field
  for(i in 1:length(vars)){
    stateCovariates[,i] <- aggregate(s@data[,vars[i]], by=list(Category=s@data$transectnum), FUN=fun, na.rm=T)[,2]
  }
  # specify covariate names
  colnames(stateCovariates) <- vars
    stateCovariates <- data.frame(stateCovariates) # unmarked expects this ahtos a data.frame
  # format our training data as umf
  return(unmarked::unmarkedFrameDS(y=as.matrix(y), siteCovs=stateCovariates, survey="point", dist.breaks=d, unitsIn="m"))
}
#' Validates a raw IMBCR data table for a user-specified species,
#' NA-filling temporal (minutes) and spatial (station) replicates across
#' transects where the species was not observed
validateTransectsForSpecies <- function(s=NULL, spp=NULL){
  return(NA)
}
#' use Monte Carlo randomization to test for heterogeneity in detection across
#' within-transect temporal replicates
monteCarloTemporalHeterogeneity <- function(s=NULL){
  return(NA)
}
#' use Monte Carlo randomization to test for heterogeneity in detection across
#' within-transect spatial replicates
monteCarloSpatialHeterogeneity <- function(s=NULL){
  return(NA)
}
