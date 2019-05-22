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
import ctypes
import functools
from render import render
parser = argparse.ArgumentParser(description='Run a diderot program')
parser.add_argument("-s", action='store_true', default=False)
args = parser.parse_args()



# select the types
intTy = ct.c_int32
floatTy = ct.c_double
myDim = 3
#Firedrake stuff would go here.
#
if myDim == 2:
    mesh = UnitSquareMesh(8,8)
elif myDim == 3:
    mesh = UnitCubeMesh(16, 16, 16)
dim = mesh.geometric_dimension()
space = VectorFunctionSpace(mesh, "Lagrange", 2, dim=dim)
f = Function(space)


if dim == 2:
    M = np.array([[-1, 1], [0, -1]], dtype="float64")
    center = np.array([0.5, 0.5], dtype="float64")
    Mtemp = np.asfarray(M, dtype="float64")
    Mptr = M.ctypes.data_as(ctypes.POINTER(floatTy))
    centerPtr = np.asfarray(center, dtype="float64").ctypes.data_as(ctypes.POINTER(floatTy))
    x,y = SpatialCoordinate(mesh)
    Mx = M.dot([x, y] - center)
    f = interpolate(as_vector(Mx), space)
elif dim == 3:
    M = np.array([[2, 1, 0], [0, 2, 1], [0, 0, 2]], dtype="float64")
    center = np.array([0.50, 0.50, 0.50], dtype="float64")
    Mtemp = np.asfarray(M, dtype="float64")
    Mptr = M.ctypes.data_as(ctypes.POINTER(floatTy))
    centerPtr = np.asfarray(center, dtype="float64").ctypes.data_as(ctypes.POINTER(floatTy))
    x,y,z = SpatialCoordinate(mesh)
    Mx = M.dot([x, y, z] - center)
    f = interpolate(as_vector(Mx), space)

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
#[0.8, 0.1, 0.5],
numPoints = 20*20*20
newXc = np.random.uniform(low=0, high=1, size=numPoints)
newYc = np.random.uniform(low=0, high=1, size=numPoints)
newZc = np.random.uniform(low=0, high=1, size=numPoints)
# if dim == 3:
#     points = [[0.8, 0.1, 0.5], [0.3, 0.3, 0.3]]#list(map(list, zip(*(xPoints, yPoints, zPoints)[0:dim]))) # list([[-2.42058914,  3.01164183, -7.03902682]])#
# elif dim == 2:
    #points = [[0.8, 0.2]]
#points = list(map(list, zip(*(newXc, newYc, newZc)[0:dim])))
points = [(0.05 * x + 0.025, 0.05 *y + 0.025, 0.05 * z + 0.025) for x in range(20) for y in range(20) for z in range(20)]
#points = [[0.225,0.525,0.575]]
pointsNrrd = "points.nrrd"
dataPoints = np.array(points, dtype="float64")
kindString = "{0}-vector".format(dim)
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)
# print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
outFileName = "stream"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]], "M" : [Mptr], "center": [centerPtr]}
outputs = [("stream", 2, outFileName), ("newStream", 2, "newStream")]
namedInputs = {"startPoints": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs)
render(outFileName, outFileName, dim=dim)
File(outFileName + ".pvd").write(f)
