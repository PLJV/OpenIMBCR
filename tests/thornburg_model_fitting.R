require(OpenIMBCR)
require(rgeos)

generate_imbcr_uniform_grid <- function(s=NULL){
  # use transects with a full 16 stations to align or uniform grid
  has_full_grid <- function(){
    transects <- unique(s$transectnum)
    return(data.frame(
        transect=transects,
        has_full_grid=unlist(lapply(transects, FUN=function(y){
          num_stations <- length(unique(
            s@data[
              s$transectnum == y, 'point']
            ))
          return(num_stations == 16)
        }))
    ))
  }
  # build our regional extent and find the centroid of each full transect
  e <- raster::extent(s)
  focal_transects <- has_full_grid()
    focal_transects <- focal_transects$transect[focal_transects$has_full_grid]
  s <- s[s$transectnum %in% focal_transects,]
  s <- do.call(
    rbind,
    lapply(focal_transects, FUN=function(x){
      centroid <- s[s$transectnum == x,]
        centroid <- rgeos::gCentroid(
            centroid[!duplicated(centroid@data[,'point']), ]
          )
      return(SpatialPointsDataFrame(
          centroid, data=data.frame(transectnum=x)
        ))
    })
  )
  # generate a uniform grid and return to user
  return(rasterize(
      s,
      raster(resolution=1000, crs=sp::CRS(raster::projection(s)), ext=e),
      # ymn=e@ymin, ymx=e@ymax+125, xmn=e@xmin+875, xmx=e@xmax),
      field=1,
      background=0
    ))
}
uniform_1km_grid <- generate_imbcr_uniform_grid(master_table)
writeRaster(uniform_1km_grid, "/home/ktaylora/uniform_imbcr_grid_test.tif", overwrite=T)
#
# main
#
master_table <- imbcrTableToShapefile(list.files(
  "/global_workspace/imbcr_number_crunching/",
  pattern="RawData_PLJV_IMBCR_20161201.csv",
  recursive=T,
  full.names=T
  ))

# 1.) calculate our detection parameters at the 1km transect scale

# 2.) build a data.frame from our station points and their respective
# USNG attributed grid cell

# 3.)
