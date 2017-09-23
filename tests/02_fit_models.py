import sys
import os

def run_r_thread(*args):
  """ System wrapper for Rscript that calls our model-fitting script chunkwise """
  os.system("Rscript " + str(args[0]) + " " + str(args[1]) + " " + str(args[2]))


def step_through_birdcodes(*args):
  for birdcode in args[0]:
      run_r_thread(
        "/home/ktaylora/Projects/OpenIMBCR/tests/thornburg_model_fitting.R",
        "/global_workspace/thornburg/vector/units_attributed_training.shp",
        birdcode)


if __name__ == "__main__":
    birdcodes = ['LBCU', 'MCLO', 'MOPL', 'BUOW', 'LOSH', 'CASP',
                'WITU', 'GRSP', 'AMGP', 'EAME', 'CCLO',
                'SCQU', 'SPPI', 'BEVI', 'NOBO', 'EUCD', 'MODO',
                'RINP', 'LASP', 'SWHA', 'STFL', 'CONI', 'DICK',
                'VESP', 'WEME', 'HOLA', 'LARB', 'COYE']

    step_through_birdcodes(birdcodes)
