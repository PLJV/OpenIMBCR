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


def build_model_training_dataset(step=24634, units_path="/gis_data/Grids/1km_usng_pljv_region_v3.0.shp", **kwargs):
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


def step_through_build_models(*args, **kwargs):
  # script path handling
  if 'r_script_path' in kwargs:
    script_path = kwargs['r_script_path']
  else:
    script_path = "/home/ktaylora/Projects/OpenIMBCR/tests/thornburg_model_fitting.R"
  # script arguments handling
  if 'training_dataset_path' in kwargs:
    train_vector_path = kwargs['training_dataset_path']
  else:
    train_vector_path = "/global_workspace/thornburg/vector/units_attributed_training.shp"
  # process by four-letter code
  for birdcode in args[0]:
    run_sys_thread(
        "Rscript",
        script_path,
        train_vector_path,
        birdcode)


def step_through_make_model_predictions(*args, **kwargs):
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
  print("usage:", sys.argv[0], "<flag(s)> <option(s)>")
  print(sys.argv[0], "-h --help")
  print(sys.argv[0], "-t --buildTrainingDataset",   ": build an IMBCR training dataset from R Vector Operations workflow")
  print(sys.argv[0], "-p --buildPredictionDataset", ": build an IMBCR prediction dataset from R Vector Operations workflow")
  print(sys.argv[0], "-b --buildModels", ": fit models and generate a shapefile of 1-km densities")
  print(sys.argv[0], "-c --codes", ": define four-letter bird codes used for model fitting")
  print(sys.argv[0], "-z --compressSessionOutput", ": 7zip compress our session output")

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
      'c:s:d:t:p:h:b1:b2:b3:b4:b5:b6:z',
      ['codes',
       'src',
       'dst',
       'buildTrainingDataset',
       'buildPredictionDataset',
       'help',
       'build1kmDistanceModels',
       'build1kmMultinomialMixtureModels',
       'build1kmZIPModels',
       'build1kmPoissonGLMs',
       'build1kmNegBinomialGLMs',
       'build1KmDistanceInterceptKriging',
       'compressSessionOutput'
       ])

  except getopt.GetoptError as e:
    print(e)
    print("see:", sys.argv[0], "-h --help for help")
    sys.exit(-1)

  if len(sys.argv) < 2 or any(d[0] in ('-h','--help') for d in options):
    print_usage()

  src_path                 = '.'
  dst_path                 = '.'
  build_training_dataset   = False
  build_prediction_dataset = False
  fit_models               = False
  predict_models           = False
  compress_session_output  = False

  for key, value in options:
    if key in ('-s','--src'):
      src_path = str(value)
    elif key in ('-d','--dst'):
      dst_path = str(value)
    elif key in ('-c','--codes'):
      birdcodes = value.split(' ')
    elif key in ('-t','--buildTrainingDataset'):
      build_training_dataset = True
    elif key in ('-p', '--buildPredictionDataset'):
      build_prediction_dataset = True
    elif key in ('-b1', '--build1kmDistanceModels'):
      fit_models = True
  if build_training_dataset:
    try:
      build_model_training_dataset()
      sp_merge_segments()
      os.system("mv units_attributed.*" + dst_path + "/units_attributed_training.*")
    except Exception as e:
      print(e)
      sys.exit(-1)
  if build_prediction_dataset:
    pass
  if fit_models:
    try:
      for code in birdcodes:
        step_through_build_models(code)
    except Exception as e:
      print(e)
      sys.exit(-1)
  if predict_models:
    try:
      step_through_make_model_predictions(find_files(pattern="*.rdata"))
    except Exception as e:
      print(e)
      sys.exit(-1)
  if compress_session_output:
    step_through_zip(find_files(pattern="*.rdata"))