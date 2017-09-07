#
# IMBCR Thornberg Modeling Precursor Work
#
# This project workspace is kludging for vector-based roving windows analyses
# and habitat fragmentation metric calculations. All pre-cursors to fitting
# an HDS model (Royle, 2004) from 1-km IMBCR transect data.
#
# Authors: KT (kyle.taylor@pljv.org) [2017], DP, LG, AB, RS, AG, Ferris Bueller
#

require(raster)
require(rgdal)
require(rgeos)
require(parallel)

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

l_buffer_grid_unit_by <- function(units=NULL, radius=1500){
  return(lapply(
      X=1:nrow(units),
      FUN=function(x) buffer_grid_unit_by(row=x, units=units)
    ))
}

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
  e_cl <- parallel::makeCluster(11)
  clusterExport(cl=e_cl, varlist=c("r"))
  return(parLapply(
      e_cl,
      X=polygon,
      fun=function(x){ raster::crop(x=r, y=x) }
    ))
  # return(lapply(
  #     polygon,
  #     FUN=function(x){ raster::crop(x=r, y=x) }
  #   ))
}

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
    "PLJV_TX_MORAP_2016_CRP.img"
  )
units <- readOGR("/home/ktaylora/","1km_usng_pljv_region_v1.0", verbose=F)

# testing
units <- units[1:10000,]

# will take ~1.35 hours for 600,000 without threading
system.time(usng_extractions <- lapply(
    X=1:nrow(units),
    FUN=buffer_grid_unit_by,
    units=units
  ))

# calculate total area composition metric over a 3x3 matrix (in vector space)
steps <- round(seq(0, nrow(units), length.out=30))
steps <- lapply(
    1:(length(steps)-1),
    FUN=function(x){ units[(steps[x]+1):(steps[x+1]),] }
  )

# this takes ~55.194 seconds for 10,000 units; ~ 0.9199 hours for 600,000
system.time(
    usng_buffers <- unlist(lapply(steps, FUN=l_buffer_grid_unit_by))
  )
rm(steps)

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

# 284.368 seconds for 10,000 units; ~4.739 hours for 600,000
system.time(usng_extractions <- extract_by(usng_buffers, r))

# 63.180 seconds for 10,000 units; ~1.053 hours for 600,000 units
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

# within-unit patch metric calculations [Short, Mixed,
# shin oak/sand sage, CRP(?)]
habitat <- binary_reclassify(r, from=221:222)
# calculate a within-unit patch count
# calculate inter-patch distance
# calculate mean patch area
# do a PCA of our fragmentation metrics
# save to disk
