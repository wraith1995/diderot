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
import numpy as np
import argparse
import nrrd_utils as nu
import nrrd
parser = argparse.ArgumentParser(description='Run a diderot program')
parser.add_argument("-s", action='store_true', default=False)
args = parser.parse_args()



# select the types
intTy = ct.c_int32
floatTy = ct.c_double

#Firedrake stuff would go here.
mesh = UnitCubeMesh(4, 4, 4)
dim = mesh.geometric_dimension()
space = VectorFunctionSpace(mesh, "Lagrange", 2, dim=dim)
f = Function(space)
if dim == 2:
    x,y = SpatialCoordinate(mesh)
    f = interpolate(as_vector([-x,y]), space)
elif dim == 3:
    x,y,z = SpatialCoordinate(mesh)
    f = interpolate(as_vector([-x,y, z]), space)

# build json
jsonFile = "evalProg.json"
dataFile = "evalProg.dill"
if os.path.exists(jsonFile):
    pass
else:
    fb.spaceToJson(space, jsonFile, refCellDefault="simplex")
    exit(0)
# build data
(preFemArgs, femArgs) = fb.passAll(f, intTy, floatTy, geometric=True)

# xPoints = np.random.uniform(low = bot[0], high = top[0], size=sizes[0])
# yPoints = np.random.uniform(low = bot[1], high = top[1], size=sizes[1])
# zPoints = np.random.uniform(low = bot[2], high = top[2], size=sizes[2])
points = [[0.8, 0.1, 0.5]]#list(map(list, zip(*(xPoints, yPoints, zPoints)[0:dim]))) # list([[-2.42058914,  3.01164183, -7.03902682]])#
pointsNrrd = "points.nrrd"
dataPoints = np.array(points, dtype="float64")
kindString = "{0}-vector".format(dim)
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)
# print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]]}
outputs = [("stream", 2, "stream"), ("newStream", 2, "newStream")]
namedInputs = {"startPoints": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs)

