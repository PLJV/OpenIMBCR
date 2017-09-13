#!/usr/bin/python3

import sys
import os
import glob
from osgeo import ogr


def printf(string=NULL):
    sys.stdout.write(str(string));
    sys.stdout.flush();

def run_r_thread(*kwargs):
  """ System wrapper for Rscript that calls our grid unit attribution script
  chunkwise """
  com = "Rscript"
  src = "/global_workspace/thornburg/thornburg_vector_operations.R"
  os.system(com + " " + src + " " + str(kwargs[0][0]) + " " + str(kwargs[0][1]))

def step_through_grid_units(step=24634, path="/gis_data/Grids/1km_usng_pljv_region_v1.0.shp"):
  """ Determine the number of segments in a fixed units shapefile and then
  process the units stepwise (using run_r_thread) """
  nFeatures = sp_count_features(path)
  segments = [(n, min(n+step, nFeatures)) for n in range(0, nFeatures, step)]
  for seg in segments:
    run_r_thread(seg)

def find_shapefile_segments(path='.', pattern="units_attributed_*.shp", **kwargs):
    """ shorthand for glob.glob file finding """
    return(glob.glob(pattern))

def sp_count_features(path=None):
    """ Returns the number of spatial features in a shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(path, 0)
    return(dataSource.GetLayer().GetFeatureCount())

def sp_merge_segments(**kwargs):
    """ List all shapefiles in the CWD and merge all features into singe file """

    outShapefile = "units_attributed.shp"
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)

    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer("units_arributed", geom_type=ogr.wkbPolygon)

    fileSegments = find_shapefile_segments()

    printf(" -- merging segments:");
    for file in fileSegments:
        printf(".")
        inDataSource = ogr.Open(file)
        inLayer = inDataSource.GetLayer()
        for feature in inLayer:
            outFeature = ogr.Feature(outLayer.GetLayerDefn())
            outFeature.SetGeometry(feature.GetGeometryRef().Clone())
            outLayer.CreateFeature(outFeature)
            outFeature = None
            outLayer.SyncToDisk()
    printf(" -- done\n");

if __name__ == "__main__" :
  step_through_grid_units()
  sp_merge_segments()
