import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
from firedrake import *
import ctypes as ct
import firedrake_build as fb
import os
import passing
import dill
import sys
import argparse

parser = argparse.ArgumentParser(description='Run a diderot program')
parser.add_argument("-s", action='store_true', default=False)
args = parser.parse_args()



# select the types
intTy = ct.c_int32
floatTy = ct.c_double

#Firedrake stuff would go here.
mesh = UnitCubeMesh(1,1,1)
space = FunctionSpace(mesh, "Lagrange", 2)
f = Function(space)
f = interpolate(Expression("x[0]*x[0] + x[1]*x[1] + x[2]*x[2]"), space)


# build json
jsonFile = "evalProg.json"
dataFile = "evalProg.dill"
if os.path.exists(jsonFile):
    pass
else:
    fb.spaceToJson(space, jsonFile, refCellDefault="simplex")
    exit(0)
# build data
(preFemArgs, femArgs) = fb.passAll(f, intTy, floatTy, geometric=not(args.s))
if args.s:
    with open(dataFile, "wb+") as f:
        dill.dump(preFemArgs, f)
    exit(0)
else:
    pass
# print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = []
outputs = []
program.go(inputs, outputs)

