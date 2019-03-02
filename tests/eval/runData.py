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
import numpy as np
import nrrd_utils as nu
parser = argparse.ArgumentParser(description='Run a diderot program')
parser.add_argument("-s", action='store_true', default=False)
parser.add_argument("--bot", default=0.0, type=float)
parser.add_argument("--top", default=2.0, type=float)
parser.add_argument("--n", default=5, type=int)
parser.add_argument("--d", default=5, type=int)
parser.add_argument("--num", default=100, type=int)
args = parser.parse_args()

# select the types
intTy = ct.c_int32
floatTy = ct.c_double



#Firedrake stuff would go here.
n = args.n
mesh = UnitCubeMesh(n, n, n)
space = FunctionSpace(mesh, "Lagrange", args.d)
f = Function(space)
numberOfDofs = len(f.dat.data)
# type to get the same
newDofs = np.random.uniform(low=args.bot, high=args.top, size=numberOfDofs)
#f.dat.data = newDofs

#generate points:
newXc = np.random.uniform(size=args.num)
newYc = np.random.uniform(size=args.num)
newZc = np.random.uniform(size=args.num)
points = map(list,zip(newXc, newYc, newZc))
pointsNrrd = "points.nrrd"
dataPoints = np.array(list(points))
kindString = "3-vector"
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)

#translates to make real results...
#get the points
#compare! yay!

# build json
jsonFile = "evalProg.json"
dataFile = "evalProg.dill"
#fb.spaceToJson(space, jsonFile, refCellDefault="simplex")
# build data
(preFemArgs, femArgs) = fb.passAll(f, intTy, floatTy, geometric=not(args.s))
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"cube": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]]}
namedInputs = {"points": pointsNrrd}
outputs = []
program.go(inputs, outputs, namedInputs=namedInputs)

