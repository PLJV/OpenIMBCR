#
# IMBCR Thornburg Modeling Precursor Work
#
# This project workspace is kludging for vector-based roving windows analyses
# and habitat fragmentation metric calculations. All pre-cursors to fitting
# an HDS model (Royle, 2004) from 1-km IMBCR transect data.
#
# Accepts three arguments (last 2 are optional). 1.) path to gridded shapefile
# containing our IMBCR grid units; 2.) lower chunk; 3.) upper chunk. The chunks
# are literally lower and upper row numbers for the units shapefile that allow
# processing units in parallel on Unix.
#
# Authors: KT (kyle.taylor@pljv.org) [2017], DP, LG, AB, RS, AG
#

# We require clean POSIX threading and runtime argument
# comprehension to run -- I haven't tested any of this on
# MS Windows so let's stop by default on non-unix platforms

stopifnot(grepl(
    tolower(.Platform$OS.type), pattern = "unix"
  ))

#
# MAIN
#

cat(" -- reading input raster/vector datasets\n")

# read-in US national grid and our source landcover data
# subset our input units by a user-defined range, if possible
argv <- na.omit(suppressWarnings(as.numeric(commandArgs(trailingOnly = T))))

# r <- raster(paste("/gis_data/Landcover/PLJV_Landcover/LD_Landcover/",
#     "PLJV_TX_MORAP_2016_CRP.img", sep=""
#   ))

r <- raster::raster(paste("/gis_data/Landcover/NASS/Raster/",
    "2016_nass_crp_test_merge.tif", sep=""
  ))


if(length(argv)>1){
  cat(
      " -- will process units chunkwise (",
      paste(argv, collapse = ":"),
      ")\n", sep = ""
    )
  units <-
    sp::spTransform(rgdal::readOGR(
        "/gis_data/Grids/",
        "1km_usng_pljv_region_v2.0",
        verbose=F
      )[(argv[1]+1):argv[2], ],
      CRSobj=sp::CRS(raster::projection(r))
    )
} else {
  argv <- c(NULL,NULL)
  units <-
    sp::spTransform(rgdal::readOGR(
        "/gis_data/Grids/",
        "1km_usng_pljv_region_v2.0",
        verbose=F
      ),
      CRSobj=sp::CRS(raster::projection(r))
    )
}

cat(" -- building 3x3 buffered grid units across project area\n")
usng_buffers_9km <- OpenIMBCR:::par_buffer_grid_units(units)

# basic implementation for extracting that will use parallel by default,
# but fails if the grid units are large

cat(" -- extracting 3x3 buffered grid units across landcover raster\n")
usng_extractions_9km <- OpenIMBCR:::extract_by(usng_buffers_9km, r)

cat(" -- calculating habitat composition/configuration metrics\n")
area_statistics <-
  data.frame(
      field_name=c(
        'grass_ar',
        'shrub_ar',
        'wetland_ar'
      ),
      src_raster_value=c(
        '176',
        'c(64,152)',
        '195'
      )
    )

configuration_statistics <- c(
    'pat_ct',
    'mn_p_ar',
    'inp_dst'
  )

for(i in 1:nrow(area_statistics)){
  # units@data[, as.character(area_statistics[i, 1])] <-
  #   par_calc_stat(
  #     # using our 1 km unit raster extractions
  #     X=usng_extractions_1km,
  #     fun = calc_total_area,
  #     # using these PLJV landcover cell values in the reclassification
  #     from = eval(parse(text=as.character(area_statistics[i, 2])))
  #   )
  units@data[, as.character(area_statistics[i, 1])] <-
    OpenIMBCR::par_calc_stat(
      # using our 3x3 buffered unit raster extractions
      X=usng_extractions_9km,
      fun = OpenIMBCR::calc_total_area,
      # using these PLJV landcover cell values in the reclassification
      from = eval(parse(text=as.character(area_statistics[i, 2])))
    )
}

# within-unit patch metric calculations [NASS Grass]
cat(" -- building a habitat/not-habitat raster surfaces\n")
valid_habitat_values <- eval(parse(
    text=paste("c(",paste(area_statistics$src_raster_value[
      !grepl(area_statistics$field_name, pattern="rd_ar")
    ], collapse = ","), ")", sep="")
  ))
cat(" -- calculating patch configuration metrics\n")
units@data[, as.character(configuration_statistics[1])] <-
  OpenIMBCR::par_calc_stat(
      # using our using our un-buffered unit raster extractions
      usng_extractions_9km,
      # parse the focal landscape configuration metric
      fun = OpenIMBCR::calc_patch_count,
      # using these PLJV landcover cell values in the supplemental
      # reclassification
      from = valid_habitat_values
    )
units@data[, as.character(configuration_statistics[2])] <-
  OpenIMBCR::par_calc_stat(
    # using our using our un-buffered unit raster extractions
    usng_extractions_9km,
    # mean patch area function:
    fun = OpenIMBCR::calc_mean_patch_area,
    # using these PLJV landcover cell values in the supplemental
    # reclassification
    from = valid_habitat_values
  )
units@data[, as.character(configuration_statistics[3])] <-
  OpenIMBCR::par_calc_stat(
    # using our using our un-buffered unit raster extractions
    usng_extractions_9km,
    # mean inter-patch distance function:
    fun = OpenIMBCR::calc_interpatch_distance,
    # using these PLJV landcover cell values in the supplemental
    # reclassification
    from = valid_habitat_values,
    backfill_missing_w=9999
  )

# save to disk
cat(" -- finished: caching metrics to disk\n")
rgdal::writeOGR(
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
