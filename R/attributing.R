# Author: Kyle Taylor <kyle.taylor@pljv.org>
# Year : 2017
# Description : various tasks related to attributing and joining spatial datasets
# for modeling with IMBCR datasets

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

#'
#' @export
calc_transect_effort <- function(df = NULL) {
  if(inherits(df,"Spatial")){
    df <- df@data
  }
  effort <- sapply(
    unique(as.character(df[,transect_fieldname(df)])),
    FUN = function(n){
      sum(length(unique(df[ df[,transect_fieldname(df)] == n, ]$point)), na.rm = T)
    }
  )
  return(effort)
}
#'
#' @export
calc_time_of_day <- function(df=NULL){
  df$tod <- as.numeric(df$starttime)
}
#'
#' @export
calc_day_of_year <- function(df=NULL){
  if(inherits(df,"Spatial")){
     s <- df
    df <- df@data
  }
  df$doy <- as.numeric(strftime(as.POSIXct(as.Date(as.character(
    df$date), "%m/%d/%Y")),format="%j"))
  if(exists("s")){
    s@data <- df
    return(s)
  } else {
    return(df)
  }
}
#' accepts a Spatial*DataFrame object with a distance observation
#' field name. Will attempt to arbitrarily bin the distances into
#' a number of distance bin intervals (specified by breaks). Will
#' reclass raw distances to bin identifier (e.g., distance class 3).
#' @export
calc_dist_bins <- function(df=NULL, p=0.95, breaks=10, force_shoulder=F){
  if(inherits(df,"Spatial")){
    df <- df@data
  }
  # Distance at which detections for all birds are considered suspect in IMBCR
  EMPIRICAL_DIST_CUTOFF <- 300 
  # what's our cut-off -- drop distance observations that are outside of the cut-off
  # if they are beyond our "difficult to detect" threshold
  cutoff <- quantile(df[,distance_fieldname(df)], p=p, na.rm=T)
  if (cutoff > EMPIRICAL_DIST_CUTOFF){
    df <- df[ is.na(df[,distance_fieldname(df)]) | (df[,distance_fieldname(df)] <= cutoff) , ]
  # if our cut-off isn't greater than our "difficult to detect" threshold then don't
  # drop any records
  } else {
    p <- 1
  }
  # define our bin intervals for a user-specified number of breaks
  if (length(breaks) == 1){
    if (force_shoulder){
      bin_intervals <- c(
        0,
        as.vector(quantile(
          df[,distance_fieldname(df)],
          probs = seq(0.5, 1, length.out = breaks-1), na.rm = T)
        )
      )
    } else {
      bin_intervals <- seq(
        from=0,
        to=quantile(df[,distance_fieldname(df)], p=p, na.rm=T),
        length.out=breaks
      )
    }
  # or did the user specify explicit breaks for us to use?
  } else {
    bin_intervals <- breaks
  }
  # calculate a detections matrix (by transect)
  df <- do.call(rbind, lapply(
    X = unique(as.character(df[,transect_fieldname(df)])),
    FUN = function(x){
      matrix(table(
          cut(
            df[ df[ , transect_fieldname(df)] == x, distance_fieldname(df) ],
            breaks = bin_intervals,
            labels = paste("dst_class_",1:(length(bin_intervals)-1),sep = "")
          )
        ),
        nrow=1
      )
    }
  ))
  # record names and return to user
  return(list(y=df, breaks=bin_intervals))
}
#' hidden function that calculates route centroids and attribute the source data
#' with number of detections for a bird of interest
#' @export
calc_route_centroids <- function(s=NULL, four_letter_code=NULL, use='detections'){
  transects <- as.vector(unique(
    s$transectnum
  ))
  pts <- lapply(
    X=transects, 
    FUN=function(transect){ 
      transect <- s[s$transectnum == transect, ] 
      detections <- transect$birdcode == four_letter_code
      # use sum of distances > 0 by default
      if(grepl(tolower(use), pattern="detections")){
        detections <- sum(transect@data[detections, distance_fieldname(s)] > 0, na.rm=T)
      # if anything else is specified, try to sum a data.frame by that variable
      } else {
        detections <- sum(transect@data[detections, use], na.rm=T)
      }
      pt <- rgeos::gCentroid(transect)
      pt$transect <- unique(transect$transectnum)
      pt@data[, use] <- detections
      return(pt)
    })
  return(do.call(rbind, pts))
}
#' calculate the latitude and longitude of feature centroids in an input dataset
#' and return an appended Spatial*DataFrame Object with lat/lon entries
calc_lat_lon <- function(s=NULL, proj_str="+init=epsg:4326"){
  coords <- as.data.frame(rgeos::gCentroid(sp::spTransform(s,proj_str), byid=T)@coords)
  colnames(coords) <- c("lon","lat")
  s@data <- cbind(s@data, coords)
  return(s)
}
#' hidden function that summarizes imbcr transect covariate data and metadata
#' by year (with list comprehension). This allows you to calculate covariates
#' at the IMBCR station level and then pool (summarize) the observations by
#' transect and year
pool_by_transect_year <- function(x=NULL, df=NULL, breaks=NULL, covs=NULL,
                                  summary_fun=median){
  breaks <- length(breaks)
  transect_year_summaries <- data.frame()
  # summarize focal_transect_year by breaking into counts within
  # distance classes and binding effort, year, and covs calculated
  # at the transect scale
  years <- sort(unique(df[df[,transect_fieldname(df)] == x, "year"]))
  for(year in years){
    focal_transect_year <- df[
      df[,transect_fieldname(df)] == x & df$year == year, ]
    # pre-allocate zeros for all bins
    distances <- rep(0,(breaks-1))
      names(distances) <- 1:(breaks-1)
    # build a pivot table of observed bins
    dist_classes <- sort(focal_transect_year$dist_class) # drop NA's
      dist_classes <- table(as.numeric(dist_classes))
    # merge pivot with pre-allocate table and add an NA bin
    distances[names(distances) %in% names(dist_classes)] <- dist_classes
      distances <- append(distances,
                          sum(is.na(focal_transect_year$dist_class)))
    distances <- as.data.frame(matrix(distances,nrow=1))
      names(distances) <- paste("distance_",c(1:(breaks-1),"NA"),sep="")
    # summarize each of the covs across the transect-year
    summary_covs <- matrix(rep(NA,length(covs)),nrow=1)
      colnames(summary_covs) <- covs
    for(cov in covs){
      # some covs are year-specific; filter accordingly
      cov_year <- names(focal_transect_year)[
          grepl(names(focal_transect_year),pattern=cov)
        ]
      if(length(cov_year)>1){
          cov_year <- cov_year[grepl(cov_year,pattern=as.character(year))]
        }
      summary_covs[,cov_year] <- summary_fun(
          focal_transect_year[,cov_year],
          na.rm=T
        )
    }
    # post-process pooled transect-year
    # keep most of our vars intact, but drop those that lack meaning at
    # the transect scale or that we have summarized above
    meta_vars <- colnames(df)[!colnames(df) %in%
              c(transect_fieldname(df), "year", "dist_class",
                distance_fieldname(df), "timeperiod", "point", "how",
                  "FID", "visual", "migrant", "cl_count", "cl_id",
                    "ptvisitzone", "ptvisiteasting", "ptvisitnorthing",
                      "rank", covs)]
    # build our summary transect-year data.frame
    focal_transect_year <- cbind(
        focal_transect_year[1, meta_vars],
        data.frame(transectnum=x, year=year),
        distances,
        summary_covs
      )
    # merge into our annual summary table
    transect_year_summaries <-
      rbind(transect_year_summaries,focal_transect_year)
  }
  return(transect_year_summaries)
}
#' shorthand vector extraction function that performs a spatial join attributes
#' vector features in x with overlapping features in y. Will automatically
#' reproject to a consistent CRS.
spatial_join <- function(x=NULL, y=NULL, drop=T){
  over <- sp::over(
      x = sp::spTransform(
          x,
          sp::CRS(raster::projection(y))
        ),
      y = y
    )
  x@data <- cbind(x@data, over)
  # drop non-overlapping features using the NA values
  # attributed to the first column of y@data
  if (drop) {
    x <- x[ !is.na(x@data[,colnames(y@data)[1]]) , ]
  }
  return(x)
}
#' generate a square buffer around a single input feature
buffer_grid_unit <- function(x=NULL, row=NULL, radius=1500){
  if (!is.null(row)){
    x <- x[row, ]
  } else {
    x <- unlist(x)[1, ]
  }
  # bogart a centroid for the enclosing grid square
  centroid <- rgeos::gCentroid(
      x,
      byid = F
    )@coords
  # build a single square polygon
  square <- sp::Polygon(
        coords = matrix(
          data=
            c(centroid[1]+radius, # x1
              centroid[1]-radius, # x2
              centroid[1]-radius, # x3
              centroid[1]+radius, # x4
              centroid[2]+radius, # y1
              centroid[2]+radius, # y2
              centroid[2]-radius, # y3
              centroid[2]-radius),# y4
          ncol=2
        ),
        hole=F
    )
    # return as a SpatialPolygons object that Raster can comprehend
    return(sp::SpatialPolygons(
        Srl=list(
          sp::Polygons(list(square),
          ID=ifelse(is.null(row), 1, row)
        )),
        proj4string=sp::CRS(raster::projection(x))
      ))
}
#' testing: generalization of the buffer_grid_units function that uses
#' lapply to buffer across all grid features in an input units
#' data.frame. Parallel implementations for this function are difficult,
#' because you have to push a copy of units onto each node
l_buffer_grid_units <- function(units=NULL, radius=1500){
  return(lapply(
      X=1:length(units@polygons),
      FUN=function(x) buffer_grid_unit(row=x, x=units, radius=radius)
    ))
}
#' testing: generalization of buffer_grid_unit with parallel comprehension
#' @export
par_buffer_grid_units <- function(units=NULL, radius=1500){
  if (!inherits(units, 'list')){
    units <- lapply(
        X=sp::split(units, 1:length(units@polygons)),
        FUN=methods::as,
        Class='SpatialPolygons'
      )
  }
  e_cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(
      cl=e_cl,
      varlist=c("buffer_grid_unit"),
      envir=environment()
    )
  parallel::clusterCall(
      e_cl,
      function(x) {
        library("rgeos");
        library("sp");
        library("raster");
      }
    )
  return(parallel::parLapply(
    cl=e_cl,
    X=units,
    fun=buffer_grid_unit,
    radius=radius
  ))
}
#' extract an input raster dataset using polygon feature(s). Bulk process list
#' objects using the parallel package. Avoid passing raster files that are
#' loaded on a network filesystem. Prefer local cached data, otherwise parallel
#' may fail unexpectedly.
#' @export
extract_by <- function(polygon=NULL, r=NULL){
  if (!inherits(polygon, 'list')){
    return(raster::crop(
        x=r,
        y=spTransform(polygon, sp::CRS(raster::projection(r)))
      ))
  }
  # simplify our input polygons list to save RAM
  if(inherits(polygon[[1]], 'SpatialPolygonsDataFrame')){
    polygon <- lapply(polygon, FUN = as, 'SpatialPolygons')
  }
  # default list comprehension: reproject to a consistent CRS
  # and then crop our input raster using a local cluster instance
  if(!raster::compareCRS(polygon[[1]],r)){
    warning(paste("input polygon= object(s) needed to be reprojected to",
    "the CRS of r= this may lead to RAM issues", sep = " "))
    polygon <- lapply(
        polygon,
        FUN=sp::spTransform,
        CRSobj=sp::CRS(raster::projection(r))
      )
  }
  e_cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl=e_cl, varlist=c("r"), envir=environment())
  ret <- parallel::parLapply(
      e_cl,
      X=polygon,
      fun=function(x){ raster::crop(x=r, y=x) }
    )
  parallel::stopCluster(cl=e_cl); rm(e_cl); gc();
  return(ret)
}
#' testing: chunk out the task of extracting buffers from a raster surface
#' using parallel in the most memory efficient way possible. Useful
#' for very large (>100,000) list of buffer polygons.
extract_by_large_df <- function(polygon=NULL, r=NULL){
  usng_extractions <- list();
  steps            <- round(seq(0, nrow(polygon), length.out=30));
  for(i in 1:(length(steps)-1)){
    usng_extractions <- append(
      usng_extractions,
      extract_by(polygon[(steps[i]+1):(steps[i+1]), ], r=r)
    )
  }
  rm(r, polygon, steps); gc();
  return(usng_extractions)
}
#' preform a binary reclassification on an input raster, with matching from
#' values substituted with 1's, and non-matching cells set to NA values. The
#' output rasters should be thematically consistent with the patch metric
#' calculations done with the SDMTools package.
binary_reclassify <- function(x=NULL, from=NULL, nomatch=NA){
  if (!inherits(x, 'list')){
    return(raster::match(x, table=from, nomatch=NA) >= 1)
  }
  return(lapply(lapply(
          x,
          FUN=raster::match,
          table=from,
          nomatch=nomatch
        ),
        FUN=raster::calc,
        fun=function(x,na.rm=F){x>=1}
      ))
}

#' shorthand total area calculation function
#' @export
calc_total_area <- function(x=NULL, area_of_cell = NULL){
   if(is.null(area_of_cell)){
     warning("no area units specified -- will use meters->square-kilometers")
     area_of_cell <- 10^-6
   }
   ret <- raster::cellStats(x, stat=sum, na.rm=T) *
          prod(raster::res(x)) * area_of_cell
   if (is.na(ret) | is.null(ret)){
     return(0)
   } else {
     return(ret)
   }
}
#' hidden function that will calculate interpatch distance using the raster::distance
#' function
#' @export
calc_interpatch_distance <- function(x=NULL, stat=mean){
  # if there are no habitat patches (issolation would theoretically be
  # very high), don't try to calc inter-patch distance because distance()
  # will throw an error
  if (sum(!is.na(raster::values(x))) == 0){
    return(9999)
  # if the whole raster is pure habitat (i.e., no inter-patches)
  } else if (sum(is.na(raster::values(x))) == raster::ncell(x)) {
    return(0)
  }
  else {
    ret <- try(stat(raster::values(raster::distance(x)), na.rm = T))
    if(class(ret) == "try-error"){
      return(NULL)
    } else {
      return(ret)
    }
  }
}
#' hidden function that will use SDMTools to calculate the mean patch area
#' of a given cell
#' @export
calc_mean_patch_area <- function(x=NULL, area_of_cell = NULL){
  if(is.null(area_of_cell)){
    warning("no area_of_cell argument supplied -- area will be calulated as square-kilometers")
    area_of_cell <- 10^-6
  }
  # if there are no habitat patches don't try to calc
  if (sum(!is.na(raster::values(x))) == 0){
    return(0)
  # if the whole raster is a giant patch
  } else if (sum(is.na(raster::values(x))) == raster::ncell(x)) {
    return(raster::ncell(x) * raster::prod(raster::res(x)) * area_of_cell)
  } else {
  # default units from SDMTools are m^2
    return(SDMTools::ClassStat(x)$mean.patch.area * area_of_cell)
  }
}
#' @export
calc_mean <- function(x=NULL, na.rm=T){
  # if there are no habitat patches don't try to calc
  if (sum(!is.na(raster::values(x))) == 0) {
    return(0)
  }
  # assume it's valid
  return( mean(raster::values(x), na.rm=na.rm) )
}
#' hidden function that will use SDMTools to calculate the number of unique patches
#' on a given raster object
#' @export
calc_patch_count <- function(x=NULL){
  # if there are no habitat patches don't try to calc
  if (sum(!is.na(raster::values(x))) == 0){
    return(0)
  # if the whole raster is a giant patch
  } else if (sum(is.na(raster::values(x))) == raster::ncell(x)) {
    return(1)
  } else {
  # default units from SDMTools are m^2
    return(SDMTools::ClassStat(x)$n.patches)
  }
}
#' shorthand function that applies an arbitrary, user-supplied summary
#' statistic function across an input list of raster objects, returning a
#' scalar value attributable to each grid unit. Will optionally
#' reclassify each input raster to a binary using binary_reclassify() if a
#' non-null value is passed to from=
#' @export
par_calc_stat <- function(X=NULL, fun=NULL, from=NULL, backfill_missing_w=0, ...){
  stopifnot(inherits(fun, 'function'))
  e_cl <- parallel::makeCluster(parallel::detectCores()-1)
  # assume our nodes will always need 'raster' and kick-in our
  # copy of the binary_reclassify shorthand while we are at it
  parallel::clusterExport(
      e_cl,
      varlist=c('binary_reclassify'),
      envir=environment()
    )
  parallel::clusterCall(
      e_cl,
      function(x) library("raster")
    )
  # if we have a non-null value for from, assume
  # that we want to do a re-classification
  ret <- unlist(parallel::parLapply(
      e_cl,
      # assume x is already a binary if from is NULL
      X = if(is.null(from)){
        X
      } else {
        # nested parallel call to reclassify our input list, if needed
        parallel::parLapply(
            cl=e_cl,
            X=X,
            fun=binary_reclassify,
            from=from
          )
      },
      fun = function(i, ...) { fun(i, ...) }
    ))
  # account for any null/na return values from our FUN statistic
  if (!is.null(backfill_missing_w) && length(ret < length(X))){
    null_ret_values <- seq(1, length(X))[
        !seq(1, length(X)) %in% as.numeric(names(ret))
      ]
    null_rets <- rep(backfill_missing_w, length(null_ret_values))
      names(null_rets) <- null_ret_values
    ret <- c(ret, null_rets)
      ret <- ret[order(as.numeric(names(ret)))]
  }
  # clean-up cluster and return FUN result to user as a vector
  parallel::stopCluster(cl=e_cl); rm(e_cl); gc();
  return(ret)
}
#' testing: std lapply implementation to benchmark against parLapply
l_calc_stat <- function(x, fun=NULL, from=NULL){
  return(lapply(
    X=if(!is.null(from)){
        lapply(x, FUN=binary_reclassify, from=from)
      } else {
        x
      },
    FUN=fun
  ))
}
