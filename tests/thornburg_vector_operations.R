#
# IMBCR Thornburg Modeling Precursor Work
#
# This project workspace is kludging for vector-based roving windows analyses
# and habitat fragmentation metric calculations. All pre-cursors to fitting
# an HDS model (Royle, 2004) from 1-km IMBCR transect data. The best parts
# of this will be unit-tested, cleaned-up, and rolled into the OpenIMBCR 'R'
# package.
#
# Authors: KT (kyle.taylor@pljv.org) [2017], DP, LG, AB, RS, AG
#

require(raster)
require(rgdal)
require(rgeos)
require(parallel)

#' generate a square buffer around a single input feature
buffer_grid_unit_by <- function(row=NULL, units=NULL, radius=1500){
  if (!is.null(row)){
    units <- units[row, ]
  }
  # bogart a centroid for the enclosing grid square
  centroid <- rgeos::gCentroid(
      units,
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
        proj4string=sp::CRS(raster::projection(units))
      ))
}
#' generalization of the buffer_grid_unit_by function that uses
#' lapply to buffer across all grid features in an input units
#' data.frame
l_buffer_grid_unit_by <- function(units=NULL, radius=1500){
  return(lapply(
      X=1:nrow(units),
      FUN=function(x) buffer_grid_unit_by(row=x, units=units, radius=radius)
    ))
}
#' extract an input raster dataset using polygon feature(s). Bulk process list
#' objects using the parallel package. Avoid passing raster files that are
#' loaded on a network filesystem. Prefer local cached data, otherwise parallel
#' may fail unexpectedly.
extract_by <- function(polygon=NULL, r=NULL){
  if (!inherits(polygon, 'list')){
    return(raster::crop(
        x=r,
        y=spTransform(polygon, sp::CRS(raster::projection(r)))
      ))
  }
  # default list comprehension: reproject to a consistent CRS
  # and then crop our input raster using a local cluster instance
  polygon <- lapply(
      polygon,
      sp::spTransform,
      sp::CRS(raster::projection(r))
    )
  e_cl <- parallel::makeCluster(parallel::detectCores()-1)
  clusterExport(cl=e_cl, varlist=c("r"))
  ret <- parLapply(
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
extract_by_large_df <- function(usng_buffers=NULL, units=NULL, r=NULL){
  usng_extractions <- list();
  steps            <- round(seq(0, nrow(units), length.out=30));
  for(i in 1:(length(steps)-1)){
    usng_extractions <- append(
      usng_extractions,
      extract_by(usng_buffers[(steps[i]+1):(steps[i+1])], r)
    )
  }
  rm(usng_buffers, r, units, steps); gc();
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
        FUN=calc,
        fun=function(x,na.rm=F){x>=1}
      ))
}
#' shorthand function that applies an arbitrary, user-supplied summary
#' statistic function across an input list of raster objects. Will optionally
#' reclassify each input raster to a binary using binary_reclassify() if a
#' non-null value is passed to from=
parLapply_calc_stat <- function(x, fun=NULL, from=NULL){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  ret <- unlist(parallel::parLapply(
      cl,
      # assume x is already a binary if from is NULL
      X=if(is.null(from)) x else binary_reclassify(x, from=from),
      fun=function(x) fun(x)
    ))
  # clean-up cluster and return FUN result to user as a vector
  parallel::stopCluster(cl=cl); rm(cl); gc();
  return(ret)
}
#' testing: std lapply implementation to benchmark against parLapply
#' this doesn't work -- currently returns the same value for every raster
l_calc_stat <- function(x, fun=NULL, from=NULL){
  return(lapply(
    X=if(is.null(from)) x else binary_reclassify(x, from=from),
    FUN=function(x) fun(x)
  ))
}
#
# MAIN
#

# read-in US national grid and our source landcover data
# subset our input units by a user-defined range, if possible
argv <- na.omit(as.numeric(commandArgs(trailingOnly = T)))

if(length(argv)>1){
  units <-
    readOGR(
      "/gis_data/Grids/",
      "1km_usng_pljv_region_v1.0",
      verbose=F
    )[(argv[1]+1):argv[2], ]
} else {
  units <-
    readOGR(
      "/gis_data/Grids/",
      "1km_usng_pljv_region_v1.0",
      verbose=F
    )
}

cat(" -- reading input raster/vector datasets\n")
r <- raster(paste("/gis_data/Landcover/PLJV_Landcover/LD_Landcover/",
    "PLJV_TX_MORAP_2016_CRP.img", sep=""
  ))

system.time(usng_buffers <- l_buffer_grid_unit_by(units))

# basic implementation for extracting that will use parallel by default,
# but fails if the grid units are large
system.time(usng_extractions <- extract_by(usng_buffers, r))

area_statistics <-
  data.frame(
      field_name=c(
        'sgp_area',
        'mgp_area',
        'ok_sg_area',
        'pl_area',
        'crp_area',
        'rd_area'
      ),
      src_raster_value=c(
        '75',
        '71',
        'c(85,87)',
        '12',
        '39',
        'c(44,41)'
      )
    )

# benchmarking parlapply implementation
# takes four minutes longer than lapply and won't use
# multiple cores on machines that are memory limited
for(i in 1:nrow(area_statistics)){
  units@data[, as.character(area_statistics[i, 1])] <-
    parLapply_calc_stat(
      # using our 3x3 buffered unit raster extractions
      usng_extractions,
      fun = function(x){
         # calculate units of total area in square-kilometers
         raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
      },
      # using these PLJV landcover cell values in the reclassification
      from = eval(parse(text=as.character(area_statistics[i, 2])))
    )
}

# benchmarking lapply implementation
# start <- Sys.time()
# for(i in 1:nrow(area_statistics)){
#   units@data[, as.character(area_statistics[i, 1])] <-
#     l_calc_stat(
#       # using our 3x3 buffered unit raster extractions
#       usng_extractions,
#       fun = function(x){
#          # calculate units of total area in square-kilometers
#          raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
#       },
#       # using these PLJV landcover cell values in the reclassification
#       from = eval(parse(text=as.character(area_statistics[i, 2])))
#     )
# }
# finish <- Sys.time()
# l_time <- finish-start

# within-unit patch metric calculations [Short, Mixed,
# shin oak/sand sage, CRP(?)]
habitat <- binary_reclassify(r, from=221:222)
# calculate a within-unit patch count
# calculate inter-patch distance
# calculate mean patch area
# do a PCA of our fragmentation metrics
# save to disk
