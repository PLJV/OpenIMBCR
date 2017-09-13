#!/usr/bin/python

import sys, os
from osgeo import ogr

def sp_count_features(path):
    """ Returns the number of spatial features in a shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(path, 0)
    return(dataSource.GetLayer().GetFeatureCount())

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

def ls_shapefiles(path='.', **kwargs):
    """ this """
    pass
def sp_merge_segments(**kwargs):
    """ List all shapefiles in the CWD and merge all features into singe file """
    outShapefile = "units_arributed.shp"
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)

    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer("units_arributed", geom_type=ogr.wkbPolygon)
    # Add input Layer Fields to the output Layer
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # Get the output Layer's Feature Definition
    outLayerDefn = outLayer.GetLayerDefn()

    # Add features to the ouput Layer
    for i in range(0, inLayer.GetFeatureCount()):
        # Get the input Feature
        inFeature = inLayer.GetFeature(i)
        # Create output Feature
        outFeature = ogr.Feature(outLayerDefn)
        # Add field values from input Layer
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # Set geometry as centroid
        geom = inFeature.GetGeometryRef()
        centroid = geom.Centroid()
        outFeature.SetGeometry(centroid)
        # Add new feature to output Layer
        outLayer.CreateFeature(outFeature)

if __name__ == "__main__" :
  step_through_grid_units()
