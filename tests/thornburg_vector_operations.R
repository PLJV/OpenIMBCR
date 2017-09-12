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

# We require clean POSIX threading and runtime argument
# comprehension to run -- I haven't tested any of this on
# MS Windows so let's stop by default on non-unix platforms

stopifnot(grepl(
    tolower(.Platform$OS.type), pattern = "unix"
  ))

require(raster)
require(rgdal)
require(rgeos)
require(parallel)

system("/usr/bin/clear")

#' generate a square buffer around a single input feature
buffer_grid_unit <- function(row=NULL, units=NULL, radius=1500){
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
#' testing: generalization of the buffer_grid_units function that uses
#' lapply to buffer across all grid features in an input units
#' data.frame
buffer_grid_units <- function(units=NULL, radius=1500){
  return(lapply(
      X=1:nrow(units),
      FUN=function(x) buffer_grid_unit(row=x, units=units, radius=radius)
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
  # simplify our input polygons list to save RAM
  if(inherits(polygon[[1]], 'SpatialPolygonsDataFrame')){
    polygon <- lapply(polygon, FUN = as, 'SpatialPolygons')
  }
  # default list comprehension: reproject to a consistent CRS
  # and then crop our input raster using a local cluster instance
  if(!raster::compareCRS(polygon[[1]],r)){
    warning(paste("input polygon= object(s) needed to be reprojected to",
    "the CRS of r= this may lead to RAM issues", sep = ""))
    polygon <- lapply(
        polygon,
        FUN=sp::spTransform,
        CRSobj=sp::CRS(raster::projection(r))
      )
  }
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
        FUN=calc,
        fun=function(x,na.rm=F){x>=1}
      ))
}
#' shorthand function that applies an arbitrary, user-supplied summary
#' statistic function across an input list of raster objects, returning a
#' scalar value attributable to each grid unit. Will optionally
#' reclassify each input raster to a binary using binary_reclassify() if a
#' non-null value is passed to from=
par_calc_stat <- function(x, fun=NULL, from=NULL, backfill_missing_w=0){
  stopifnot(!is.null(fun))
  e_cl <- parallel::makeCluster(parallel::detectCores()-1)
  # assume our nodes will always need 'raster' and kick-in our
  # copy of the binary_reclassify shorthand while we are at it
  parallel::clusterExport(e_cl, varlist=c('binary_reclassify'))
  parallel::clusterCall(
      e_cl,
      function(x) library("raster")
    )
  ret <- unlist(parallel::parLapply(
      e_cl,
      # assume x is already a binary if from is NULL
      X = if(is.null(from)){
        x
      } else {
        # nested parallel call to reclassify our input list, if needed
        parallel::parLapply(
            cl=e_cl,
            X=x,
            fun=binary_reclassify,
            from=from
          )
      },
      fun = fun
    ))
  # account for any null/na return values from our FUN statistic
  if (!is.null(backfill_missing_w) && length(ret < length(x))){
    null_ret_values <- seq(1, length(x))[
        !seq(1, length(x)) %in% as.numeric(names(ret))
      ]
    null_rets <- rep(backfill_missing_w, length(null_ret_values))
      names(null_rets) <- null_ret_values
    ret <- sort(c(ret, null_rets))
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
#
# MAIN
#

cat(" -- reading input raster/vector datasets\n")

# read-in US national grid and our source landcover data
# subset our input units by a user-defined range, if possible
argv <- na.omit(as.numeric(commandArgs(trailingOnly = T)))

r <- raster(paste("/gis_data/Landcover/PLJV_Landcover/LD_Landcover/",
    "PLJV_TX_MORAP_2016_CRP.img", sep=""
  ))

if(length(argv)>1){
  cat(
      " -- will process units chunkwise (",
      paste(argv, collapse = ":"),
      ")\n", sep = ""
    )
  units <-
    sp::spTransform(readOGR(
        "/gis_data/Grids/",
        "1km_usng_pljv_region_v1.0",
        verbose=F
      )[(argv[1]+1):argv[2], ],
      CRSobj=sp::CRS(raster::projection(r))
    )
} else {
  units <-
    sp::spTransform(readOGR(
        "/gis_data/Grids/",
        "1km_usng_pljv_region_v1.0",
        verbose=F
      ),
      CRSobj=sp::CRS(raster::projection(r))
    )
}



cat(" -- building 3x3 buffered grid units across project area\n")
system.time(usng_buffers <- buffer_grid_units(units))
gc();

# basic implementation for extracting that will use parallel by default,
# but fails if the grid units are large
cat(" -- extracting 3x3 buffered grid units across landcover raster\n")
system.time(usng_extractions <- extract_by(usng_buffers, r))
gc();

cat(" -- calculating habitat composition metrics\n")

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
    par_calc_stat(
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

configuration_statistics <- c(
    'pat_cnt',
    'mn_pat_ar',
    'in_pat_dist'
  )

cat(" -- extracting across our unbuffered grid units\n")

# ~ 1184 seconds per 24634 units
system.time(
  usng_extractions <- extract_by(
    polygon=unlist(split(units,1:nrow(units))),
    r=r
  )
)

rm(r); gc();

# within-unit patch metric calculations [Short, Mixed,
# shin oak/sand sage, CRP(?)]
cat(" -- building a habitat/not-habitat raster surfaces\n")

valid_habitat_values <- eval(parse(
    text=paste("c(",paste(area_statistics$src_raster_value[
      !grepl(area_statistics$field_name, pattern="rd_area")
    ], collapse = ","), ")", sep="")
  ))

cat(" -- calculating patch configuration metrics\n")

# ~200.377 seconds
# system.time(units@data[, as.character(configuration_statistics[1])] <-
#   l_calc_stat(
#     # using our using our un-buffered unit raster extractions
#     usng_extractions,
#     # parse the focal landscape configuration metric
#     fun = function(x){ SDMTools::ClassStat(x)$n.patches },
#     # using these PLJV landcover cell values in the supplemental
#     # reclassification
#     from = valid_habitat_values
#   ))
# ~ 49 seconds
system.time(units@data[, as.character(configuration_statistics[1])] <-
par_calc_stat(
    # using our using our un-buffered unit raster extractions
    usng_extractions,
    # parse the focal landscape configuration metric
    fun = function(x){ SDMTools::ClassStat(x)$n.patches },
    # using these PLJV landcover cell values in the supplemental
    # reclassification
    from = valid_habitat_values
  ))
system.time(units@data[, as.character(configuration_statistics[2])] <-
  par_calc_stat(
    # using our using our un-buffered unit raster extractions
    usng_extractions,
    # mean patch area function:
    fun = function(x){ SDMTools::ClassStat(x)$mean.patch.area },
    # using these PLJV landcover cell values in the supplemental
    # reclassification
    from = valid_habitat_values
  ))
system.time(units@data[, as.character(configuration_statistics[3])] <-
  par_calc_stat(
    # using our using our un-buffered unit raster extractions
    usng_extractions,
    # mean inter-patch distance function:
    fun = function(x){
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
          ret <- try(mean(raster::values(raster::distance(x)), na.rm = T))
          if(class(ret) == "try-error"){
            return(NULL)
          } else {
            return(ret)
          }
        }
    },
    # using these PLJV landcover cell values in the supplemental
    # reclassification
    from = valid_habitat_values,
    backfill_missing_w=9999
  ))

# do a PCA of our fragmentation metrics

# save to disk
cat(" -- finished: caching metrics to disk")
writeOGR(
    units,
    dsn=".",
    layer=paste(
        "units_attributed_",
        argv[1],"-",argv[2],
        sep=""
      ),
    overwrite=T,
    driver="ESRI Shapefile"
  )