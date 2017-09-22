import sys
import os

def run_r_thread(com="Rscript",
                 src=None,
                 *kwargs):
  """ System wrapper for Rscript that calls our model-fitting script chunkwise """
  os.system(com + " " + src + " " + str(kwargs[0][0]) + " " + str(kwargs[0][1]))


def step_through_birdcodes(*kwargs):
  for birdcode in kwargs[0]:
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
