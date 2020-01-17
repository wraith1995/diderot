import sys
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
sys.path.append('/home/teocollin/gitcode/curvedMesh')
sys.path.append('isoNewton')

from firedrake import *
import ctypes as ct
import firedrake_build as fb
import os
import passing
import dill
import sys
import numpy as np
import argparse
import nrrd_utils as nu
import nrrd
import ctypes
import functools
from render import render
from load import buildFiredrakeTetMesh
from _ctypes import dlclose
intTy = ct.c_int32
floatTy = ct.c_double



def runTruth(outFile, pointsNrrd, femArgs,
              timeSteps=32, timeEps=0.0000001,
              stepSize=0.04, stepMax=32,
              errorMax=10000000, second=0, save=True):

    programNameArg = "evalProg"
    nameSpaceArg = "evalProg"
    outFileName = outFile
    library = ct.CDLL("./s3/" + programNameArg + ".so")
    program = passing.Library(library, nameSpace=nameSpaceArg)
    inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]],
            "timeSteps": [intTy(timeSteps)], "timeEps" : [floatTy(timeEps)],
              "stepSize": [floatTy(stepSize)], "stepMax": [intTy(stepMax)],
              "errorMax": [floatTy(errorMax)], "second": [intTy(second)]}
    if save:
        outputs = [("stream", 2, outFileName)]
    else:
        outputs = []
    namedInputs = {"startPoints": pointsNrrd}
    r = program.go(inputs, outputs, namedInputs=namedInputs, shutdown=False, verbose=False, time=True)
    dlclose(library._handle)
    return(r)
    
    # render(outFileName, outFileName, dim=dim)
    # render(outFileNamep, outFileNamep, dim=dim)
    # File(outFileName + ".pvd").write(f)
