# Author: Kyle Taylor <kyle.taylor@pljv.org>
# Year : 2017
# Description : various tasks related to formatting and parsing raw IMBCR data
# into Spatial* and Unmarked* objects.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
# Define some useful local functions for manipulating IMBCR data
#
#' hidden function that greps for four-letter-codes
birdcode_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^bird")])
}
#' hidden function that greps for four-letter-codes
commonname_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^c.*.[.]n.*.")])
}
#' hidden function that greps for the distance field name
distance_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^rad")])
}
#' hidden function that greps for the transect field name
transect_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^tran")])
}
#' hidden function that greps for the timeperiod field name
timeperiod_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^tim")])
}
#' built-in (hidden) function that will accept the full path of a shapefile 
#' and parse the string into something that rgdal can understand (DSN + Layer).
parseLayerDsn <- function(x=NULL){
  path <- unlist(strsplit(x, split="/"))
    layer <- gsub(path[length(path)],pattern=".shp",replacement="")
      dsn <- paste(path[1:(length(path)-1)],collapse="/")
  return(c(layer,dsn))
}
#' built-in (hidden) function that will accept the full path of a shapefile and read using rgdal::ogr
#' @param path argument provides the full path to an ESRI Shapefile
readOGRfromPath <- function(path=NULL, verbose=F){
  landscapeAnalysis:::include('rgdal')
  path <- landscapeAnalysis:::parseLayerDsn(path)

  layer <- path[1]
    dsn <- path[2]

  return(rgdal::readOGR(dsn,layer,verbose=verbose))
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
#' parse a source CSV (published by BCR) for IMBCR data and return
#' as a SpatialPointsDataFrame to the user for post-processing or
#' conversion to unmarked data.frame objects for modeling
#' @export
imbcrTableToShapefile <- function(filename=NULL,outfile=NULL,
                                  write=F, calc_spatial_covs=T){
  if(is.null(outfile) && write && !is.character(filename)){
    stop("cannot write shapefile to disk without an outfile= or filename.ext to parse")
  }
  # sanity-check to see if our output shapefile already exists
  s <- recursiveFindFile(outfile)[1]
  if(!is.null(s)){
    s <- OpenIMBCR:::readOGRfromPath(s)
  # Parse ALL BCR raw data tables into a single table
  } else {
    if(inherits(filename,"data.frame")){
      t <- filename
    } else {
      if(length(filename)==1){
        t <- read.csv(filename)
      # legacy support : rbind across multiple CSV files
      } else {
        t <- lapply(recursiveFindFile(
              name=filename
            ),
            read.csv
          )
        t <- do.call(rbind, t)
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
        coords = data.frame(
            x = s[[length(s)]]$ptvisiteasting,
            y = s[[length(s)]]$ptvisitnorthing
          ),
        data = s[[length(s)]],
        proj4string = sp::CRS(raster::projection(
            paste("+init=epsg:269", zone, sep = "")
          ))
      )
      # based on : http://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r
      row.names(s[[length(s)]]) <-
        paste(
          letters[length(s)],
          row.names(s[[length(s)]]),
          sep="."
        )
    }
    # merge our segments and convert to a consistent
    # popular CRS
    s <- do.call(
        sp::rbind.SpatialPointsDataFrame,
        lapply(
          s,
          FUN=sp::spTransform,
          sp::CRS(raster::projection("+init=epsg:2163"))
        ))
    s$FID <- 1:nrow(s)
    # calculate spatial (lat/lon) covariates for each station?
    if(calc_spatial_covs){
     # calculate lat/lon covariates in WGS84
     coords <- sp::spTransform(s,"+init=epsg:4326")@coords
       coords <- cbind(coords, coords^2, log10(coords+361))
         colnames(coords) <- c("lon","lat","lon_2","lat_2","ln_lon","ln_lat")
     s@data <- cbind(s@data,coords)
       rm(coords)
    }
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
        distance_disqualifier <- (distance_disqualifier < 0 | distance_disqualifier > 800)

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
#' kludging to back-fill any transect stations in an imbcr data.frame
#' that were sampled, but where a focal species wasn't observed, with
#' NA values
#' @export
scrub_imbcr_df <- function(df,
                           allow_duplicate_timeperiods=F,
                           four_letter_code=NULL,
                           back_fill_all_na=F,
                           drop_na="none"){
  # throw-out any lurking 88 values, count before start values, and
  # -1 distance observations
  df <- df[!df@data[, timeperiod_fieldname(df)] == 88, ]
  df <- df[!df@data[, timeperiod_fieldname(df)] == -1, ]
  df <- df[!df@data[, distance_fieldname(df)]   == -1, ]
  # build a dataframe for our detections
  detected <- toupper(df@data[, birdcode_fieldname(df)]) ==
    toupper(four_letter_code)
  # define a pool of potential non-detections
  not_detected <- df[!detected, ]
  not_detected@data[, distance_fieldname(df)] <- NA
  not_detected@data[, birdcode_fieldname(df)] <- toupper(four_letter_code)
  not_detected@data[, commonname_fieldname(df)] <-
    as.character(df@data[which(detected == T)[1],commonname_fieldname(df)])
  not_detected@data[, 'cl_count']             <- 0 # not used, but stay honest
  # allow a single NA value for each station, but only keep the NA values
  # if we didn't observe the bird at that point -- by default, don't allow
  # duplicate NA's within a time-period
  not_detected <- not_detected[!duplicated(not_detected@data[,
                    c(transect_fieldname(not_detected), "year", "point",
                      if (allow_duplicate_timeperiods)
                        timeperiod_fieldname(not_detected)
                      else NULL
                    )
                  ]), ]
  # allow multiple detections at stations
  detected <- df[detected, ]
  # take the merge of detections and non-duplicated, non-detections as
  # our new data.frame
  transect_heuristic <- function(x=NULL){
    x <- x@data[, c(transect_fieldname(df), 'year', 'point')]
    return(round(sqrt(as.numeric(x[,1])) + sqrt(x[,2]) + sqrt(x[,3]),5))
  }
  not_detected <- not_detected[
      !transect_heuristic(not_detected) %in%
      transect_heuristic(detected),
    ]
  df <- rbind(y=detected, x=not_detected)
  # zero-inflation fix (1) : don't drop any transects
  if(grepl(tolower(drop_na), pattern="none")){
      df <- sp:::rbind.SpatialPointsDataFrame(detected, not_detected)
  }
  # zero-inflation fix (2) : drop transects without at least one detection
  else if(grepl(tolower(drop_na),pattern="some")){
    valid_transects <- unique(detected@data[,transect_fieldname(df)])
    df <- df[df@data[,transect_fieldname(df)] %in% valid_transects,]
  }
  # zero-inflation fix (3) : drop all NA values
  else if(grep(tolower(drop_na), pattern="all")){
    df <- detected
  }
  df[order(sqrt(as.numeric(df$transectnum))+sqrt(df$year)+sqrt(df$point)),]
}
#' accepts a formatted IMBCR SpatialPointsDataFrame and builds an
#' unmarkedFrameGDS data.frame that we can use for modeling with
#' the unmarked package. 
#' @export
build_unmarked_gds <- function(df=NULL,
                               numPrimary=1,
                               distance_breaks=NULL,
                               covs=NULL,
                               unitsIn="m",
                               summary_fun=median,
                               drop_na_values=T
                               ){
  if(inherits(df, "Spatial")){
    df <- df@data
  }
  # determine distance breaks / classes, if needed
  if(is.null(distance_breaks)){
    distance_breaks  = df$distance_breaks
    distance_classes = append(sort(as.numeric(unique(
                            df$dist_class))),
                            NA
                          )
  } else {
    distance_classes = append(1:length(distance_breaks)-1, NA)
  }
  # parse our imbcr data.frame into transect-level summaries
  # with unmarked::gdistsamp comprehension
  transects <- unique(df[,transect_fieldname(df)])
  # pool our transect-level observations
  transects <- do.call(rbind,
      lapply(
          transects,
          FUN=pool_by_transect_year,
          df=df, breaks=distance_breaks,
          covs=covs
        )
    )
  # bug fix : drop entries with NA values before attempting PCA or quantile pruning
  if(drop_na_values){
    transects <- transects[ !as.vector(rowSums(is.na(transects)) > 0) , ]
    transects <- transects[ ,!grepl(colnames(transects), pattern="_NA")]
  }
  # build our unmarked frame and return to user
  return(unmarked::unmarkedFrameGDS(
      # distance bins
      y=transects[,grepl(names(transects),pattern="distance_")],
      # covariates that vary at the site (transect) level
      siteCovs=transects[,!grepl(colnames(transects),pattern="distance_")],
      # not used (covariates at the site-year level)
      yearlySiteCovs=NULL,
      survey="point",
      unitsIn=unitsIn,
      dist.breaks=distance_breaks,
      numPrimary=numPrimary # should be kept at 1 (no within-season visits)
    ))
}
#' hidden function used to clean-up an unmarked data.frame (umdf) by dropping
#' any NA columns attributed by scrub_imbcr_df(), mean-center (scale) our site
#' covariates (but not sampling effort!), and do some optional quantile filtering
#' that drops covariates with low variance, which is a useful 'significance
#' pruning' precursor for principal components analysis. Prefer dropping the NA
#' bin here (rather than in scrub_imbcr_df), so that we still have an accurate
#' account of total sampling effort to attribute in scrub_unmarked_dataframe().
scrub_unmarked_dataframe <- function(x=NULL, normalize=T, prune_cutoff=NULL){
  row.names(x@y) <- NULL
  row.names(x@siteCovs) <- NULL
  x@y <- x@y[,!grepl(colnames(x@y), pattern="_NA")]
  x@obsToY <- matrix(x@obsToY[,1:ncol(x@y)],nrow=1)
  # do some quantile pruning of our input data, selectively dropping
  # an arbitrary number of variables based on a user-specified
  # low-variance threshold
  if(!is.null(prune_cutoff)){
    # e.g., what is the total variance for each cov across all sites?
    # drop those standardized variables with < prune_cutoff=0.05 variance
    effort_field <- ifelse(
        sum(grepl(colnames(x@siteCovs), pattern="effort")),
        "effort",
        NULL
      )
    vars_to_scale <- colnames(x@siteCovs)[
        !grepl(tolower(colnames(x@siteCovs)), pattern="effort")
      ]
    # bug-fix : only try to prune numeric variables
    is_numeric <- apply(
        x@siteCovs[1,vars_to_scale],
        MARGIN=2,
        FUN=function(x) !is.na(suppressWarnings(as.numeric(x)))
      )
    if(length(vars_to_scale)!=sum(is_numeric)){
      warning(paste(
          "the following input variables are not numeric and cannot be",
          "filtered by quantile and will not be pruned:",
          paste(
              vars_to_scale[!is_numeric],
              collapse=", "
            )
        ))
      vars_to_scale <- vars_to_scale[is_numeric]
    }
    # calculate relative variance across all sites for each variable (column)
    variance <- apply(
      x@siteCovs[,vars_to_scale],
      MARGIN=2,
      FUN=function(x) ( (x - min(x)) / (max(x)-min(x)) ) # quick min-max normalize
    )
    # min-max will return NA on no variance (e.g., divide by zero)
    variance[is.na(variance)] <- 0
    variance <- apply(
        variance,
        MARGIN=2,
        FUN=var
      )
    # drop variables that don't meet our a priori variance threshold
    dropped <- as.vector(variance < quantile(variance, p=prune_cutoff))
    if(sum(dropped)>0){
      warning(paste(
        "prune_cutoff dropped these variables due to very small variance: ",
        paste(colnames(x@siteCovs[,vars_to_scale])[dropped], collapse=", "),
        sep=""
      ))
      keep <- unique(c(
        names(is_numeric[!is_numeric]),
        effort_field,
        vars_to_scale[!dropped]
      ))
      x@siteCovs <- x@siteCovs[, keep]
    }
  }
  # normalize our site covariates?
  if(normalize){
    # don't try to normalize non-numeric values -- drop these as site covs
    x@siteCovs <-
      x@siteCovs[ , as.vector(unlist(lapply(x@siteCovs[1,], FUN=is.numeric)))]
    # don't normalize the "effort" field
    vars_to_scale <- colnames(x@siteCovs)[
        !grepl(tolower(colnames(x@siteCovs)), pattern="effort")
      ]
    # scaling call
    x@siteCovs[,vars_to_scale] <- as.data.frame(
        scale(x@siteCovs[,vars_to_scale])
      )
    # sanity check : do some variable pruning based on variance
    # from our normalization step -- drop variables with low variance
    # from consideration and report dropped variables to user
    dropped <- as.vector(unlist(lapply(x@siteCovs[1,], FUN=is.na)))
    if(sum(dropped)>0){
      warning(paste(
        "scale() dropped these variables due to very small variance: ",
        paste(colnames(x@siteCovs)[dropped], collapse=", "),
        sep=""
      ))
      x@siteCovs <- x@siteCovs[,!dropped]
    }
  }
  return(x)
}
