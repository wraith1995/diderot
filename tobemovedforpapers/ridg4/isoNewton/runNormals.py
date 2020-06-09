import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
sys.path.append('/home/teocollin/gitcode/curvedMesh')
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
#from render import render
#from load import buildFiredrakeTetMesh


intTy = ct.c_int32
floatTy = ct.c_double


def getNormals(f, pointsNrrdFile):
    (_, femArgs) = fb.passAll(f, intTy, floatTy, geometric=True)
    programNameArg = "isoNewton/evalProg"
    nameSpaceArg = "evalProg"
    outFileName = "normals"
    outFileName1 = "stren"
    library = ct.CDLL("./" + programNameArg + ".so")
    program = passing.Library(library, nameSpace=nameSpaceArg)
    inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]]}
    namedInputs = {"ipos": pointsNrrdFile}
    outputs = [("ref", 1, outFileName), ("stren", 1, outFileName1)]
    program.go(inputs, outputs, namedInputs=namedInputs, verbose=True, workers=None, shutdown=False)

