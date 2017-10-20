#!/usr/bin/python3

import sys
import os
import getopt
import glob
from osgeo import ogr


def run_sys_thread(*args, **kwargs):
  """ System wrapper for Rscript that calls our R thread """
  if 'background' in kwargs:
    background = '&'
  else:
    background = ''
  os.system(' '.join(args) + background)


def find_files(path='.', pattern="*.rdata", **kwargs):
  """ shorthand for glob.glob file finding """
  return glob.glob(pattern, recursive=True)


def sp_count_features(path=None):
  """ Returns the number of spatial features in a shapefile"""
  driver = ogr.GetDriverByName("ESRI Shapefile")
  dataSource = driver.Open(path, 0)
  return (dataSource.GetLayer().GetFeatureCount())


def step_through_grid_units(step=24634, units_path="/gis_data/Grids/1km_usng_pljv_region_v3.0.shp", **kwargs):
  """ Determine the number of segments in a fixed units shapefile and then
  process the units stepwise (using run_sys_thread) """
  # script path handling
  if 'r_script_path' in kwargs:
    script_path = kwargs['r_script_path']
  else:
    script_path = "/global_workspace/thornburg/thornburg_vector_operations.R"
  nFeatures = sp_count_features(units_path)
  segments = [(n, min(n+step, nFeatures)) for n in range(0, nFeatures, step)]
  for seg in segments:
    run_sys_thread(
      "Rscript",
      script_path,
      str(seg[0]),
      str(seg[1]))


def step_through_training(*args, **kwargs):
  # script path handling
  if 'r_script_path' in kwargs:
    script_path = kwargs['r_script_path']
  else:
    script_path = "/home/ktaylora/Projects/OpenIMBCR/tests/thornburg_model_fitting.R"
  # script arguments handling
  if 'training_dataset_path' in kwargs:
    train_path = kwargs['training_dataset_path']
  else:
    train_path = "/global_workspace/thornburg/vector/units_attributed_training.shp"
  # process by four-letter code
  for birdcode in args[0]:
    run_sys_thread(
        "Rscript",
        script_path,
        train_path,
        birdcode)


def step_through_prediction(*args, **kwargs):
  if 'r_script_path' in kwargs:
    script_path = kwargs['r_script_path']
  else:
    script_path = "/home/ktaylora/Projects/OpenIMBCR/tests/thornburg_model_prediction.R"
  for file in args[0]:
    run_sys_thread(
      "Rscript",
      script_path,
      file)


def step_through_zip(*args):
  files = find_files(pattern="*.rdata")

def sp_merge_segments(**kwargs):
  """ List all shapefiles in the CWD and merge all features into singe file """

  outShapefile = "units_attributed.shp"
  outDriver = ogr.GetDriverByName("ESRI Shapefile")

  if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)

  outDriver = ogr.GetDriverByName('ESRI Shapefile')
  outDataSource = outDriver.CreateDataSource(outShapefile)
  outLayer = outDataSource.CreateLayer("units_arributed", geom_type=ogr.wkbPolygon)

  fileSegments = find_files(pattern="units_attributed_*.shp")

  print(" -- merging segments:")

  for file in fileSegments:
    print(".")
    inDataSource = ogr.Open(file)
    inLayer = inDataSource.GetLayer()
    # if this is the first file in the merge, create attribute table
    if outLayer.GetLayerDefn().GetFieldCount() == 0:
      for index in range(inLayer.GetLayerDefn().GetFieldCount()):
        srcField = inLayer.GetLayerDefn().GetFieldDefn(index)
        field = ogr.FieldDefn(srcField.GetName(), srcField.GetType())
        field.SetWidth(srcField.GetWidth())
        field.SetPrecision(srcField.GetPrecision())
        outLayer.CreateField(field)
    # create a feature geometry and populate all fields from source
    for feature in inLayer:
      outFeature = ogr.Feature(outLayer.GetLayerDefn())
      outFeature.SetGeometry(feature.GetGeometryRef().Clone())
      for index in range(feature.GetFieldCount()):
        outFeature.SetField(index, feature.GetField(index))
      outLayer.CreateFeature(outFeature)
      outFeature = None
    # cache the feature dump to disk and continue
    outLayer.SyncToDisk()

  print(" -- done\n")

def print_usage():
  print("usage:", sys.argv[0] "<flag(s)> <option(s)>")
  print(sys.argv[0], "-h --help")
  print(sys.argv[0], "-t --buildTrainingDataset",   ": build an IMBCR training dataset from R Vector Operations workflow")
  print(sys.argv[0], "-p --buildPredictionDataset", ": build an IMBCR prediction dataset from R Vector Operations workflow")
  print(sys.argv[0], "-b --buildModels", ": fit models and generate a shapefile of 1-km densities")
  print(sys.argv[0], "-c --codes", ": define four-letter bird codes used for model fitting")
  sys.exit(0)

if __name__ == "__main__":

  birdcodes = ['LBCU', 'MCLO', 'MOPL', 'BUOW', 'LOSH', 'CASP',
               'WITU', 'GRSP', 'AMGP', 'EAME', 'CCLO',
               'SCQU', 'SPPI', 'BEVI', 'NOBO', 'EUCD', 'MODO',
               'RINP', 'LASP', 'SWHA', 'STFL', 'CONI', 'DICK',
               'VESP', 'WEME', 'HOLA', 'LARB', 'COYE']

  try:
    options, remainder = getopt.getopt(
      sys.argv[1:],
      'c:t:p:h:b',
      ['codes',
       'buildTrainingDataset',
       'buildPredictionDataset',
       'help',
       'buildModels'
       ])

  except getopt.GetoptError as e:
    print(e)
    print("see:", sys.argv[0], "-h --help for help")
    sys.exit(1)

  if len(sys.argv) < 2 or any(d[0] in ('-h','--help') for d in options):
    print_usage()

  for key, value in options:
    if key in ('-c','--codes'):
      birdcodes = value.split(' ')
    elif key in ('-t','--buildTrainingDataset'):
      if value:
        dst_path = value
      else:
        dst_path = "/gis_data/Grids/units_attributed_training.shp"

      if not os.path.exists(path):
        step_through_grid_units()
        sp_merge_segments()
        os.system("mv units_attributed.* /gis_data/Grids/units_attributed_training.*")
    elif key in ('-p', '--buildPredictionDataset'):
      pass

    for code in birdcodes:
      step_through_training()
      step_through_prediction(find_files(pattern="*.rdata"))
      step_through_zip(find_files(pattern="*.rdata"))
