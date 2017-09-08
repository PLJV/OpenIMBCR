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
  # default list comprehension
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
  parallel::stopCluster(e_cl)
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

#
# MAIN
#

# read-in US national grid
cat(" -- reading input raster/vector datasets\n")
r <- raster(paste("/gis_data/Landcover/PLJV_Landcover/LD_Landcover/",
    "PLJV_TX_MORAP_2016_CRP.img", sep=""
  ))

units <- readOGR("/gis_data/Grids/","1km_usng_pljv_region_v1.0", verbose=F)

# subset our input units by a user-defined range, if possible
argv <- na.omit(as.numeric(commandArgs(trailingOnly=T)))
if(length(argv)>1){
  units <- units[(argv[1]+1):argv[2],]
}

# will take ~1.35 hours for 600,000 without threading
# system.time(usng_extractions <- lapply(
#     X=1:nrow(units),
#     FUN=buffer_grid_unit_by,
#     units=units
#   ))

# calculate total area composition metric over a 3x3 matrix (in vector space)
# step-wise implementation
# steps <- round(seq(0, nrow(units), length.out=30))
# steps <- lapply(
#     1:(length(steps)-1),
#     FUN=function(x){ units[(steps[x]+1):(steps[x+1]),] }
#   )

# this takes ~55.194 seconds for 10,000 units; ~ 0.9199 hours for 600,000
# steps <- round(seq(0, nrow(units), length.out=30))
# system.time(
#     usng_buffers <- unlist(lapply(steps, FUN=l_buffer_grid_unit_by))
#   )

system.time(usng_buffers <- l_buffer_grid_unit_by(units))

# buffer all grid units so that our area extractions are consistent with
# a 3x3 matrix -- try to do this in parallel to be efficient
# cl           <- parallel::makeCluster(6)
# usng_buffers <- list()
#
# for(step in steps){
#   clusterExport(
#       cl=cl,
#       varlist=c("buffer_grid_unit_by","l_buffer_grid_unit_by","step")
#     )
#   system.time(usng_buffers <- append(
#       usng_buffers,
#       parallel::parLapply(
#         cl,
#         X=steps,
#         fun=l_buffer_grid_unit_by,
#         simplify=T
#       )))
# }
#
# parallel::stopCluster(cl=cl)

# this parsing and for-looping is needed to deal with memory problems
# it's difficult to parallelize this task across a really big selection
# of grid units
  # usng_extractions <- list();
  # steps            <- round(seq(0, nrow(units), length.out=30));
  # for(i in 1:(length(steps)-1)){
  #   usng_extractions <- append(
  #     usng_extractions,
  #     extract_by(usng_buffers[(steps[i]+1):(steps[i+1])], r)
  #   )
  # }
  # rm(usng_buffers); gc();

# test a function wrapper
# usng_extractions <- extract_by_large_df(usng_buffers, units, r)

# basic implementation for extracting that will use parallel by default,
# but fails if the grid units are large
system.time(usng_extractions <- extract_by(usng_buffers, r))

# lapply equivalent -- not much faster
# steps <- lapply(
#     1:(length(steps)-1),
#     FUN=function(x){
#       usng_extractions <- append(usng_extractions, extract_by(
#         usng_buffers[(steps[x]+1):(steps[x+1])],
#         r
#       ))
#     }
#   )

# 284.368 seconds for 10,000 units; ~4.739 hours for 600,000
# system.time(usng_extractions <- extract_by(usng_buffers, r))

# 63.180 seconds for 10,000 units; ~1.053 hours for 600,000 units
cl <- parallel::makeCluster(parallel::detectCores()-1)
system.time(units$sgp_total_area <- unlist(parallel::parLapply(
    cl,
    X=binary_reclassify(usng_extractions, from=331),
    fun=function(x){
      # units of total area are in square-kilometers
      raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
    }
  )))
system.time(units$mgp_total_area <- unlist(parallel::parLapply(
    cl,
    X=binary_reclassify(usng_extractions, from=332),
    fun=function(x){
      # units of total area are in square-kilometers
      raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
    }
  )))
system.time(units$oak_sage_total_area <- unlist(parallel::parLapply(
    cl,
    X=binary_reclassify(usng_extractions, from=342:343),
    fun=function(x){
      # units of total area are in square-kilometers
      raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
    }
  )))
system.time(units$playas_total_area <- unlist(parallel::parLapply(
    cl,
    X=binary_reclassify(usng_extractions, from=121),
    fun=function(x){
      # units of total area are in square-kilometers
      raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
    }
  )))
system.time(units$crp_total_area <- unlist(parallel::parLapply(
    cl,
    X=binary_reclassify(usng_extractions, from=211),
    fun=function(x){
      # units of total area are in square-kilometers
      raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
    }
  )))
system.time(units$road_total_area <- unlist(parallel::parLapply(
    cl,
    X=binary_reclassify(usng_extractions, from=221:222),
    fun=function(x){
      # units of total area are in square-kilometers
      raster::cellStats(x, stat=sum) * prod(raster::res(x)) * 10^-6
    }
  )))
parallel::stopCluster(cl)
# within-unit patch metric calculations [Short, Mixed,
# shin oak/sand sage, CRP(?)]
habitat <- binary_reclassify(r, from=221:222)
# calculate a within-unit patch count
# calculate inter-patch distance
# calculate mean patch area
# do a PCA of our fragmentation metrics
# save to disk
