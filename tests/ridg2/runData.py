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



eps = 0.01
rad = 0.5
iso = 0.0
fStren = 0.1
fBias = 0.01
zfuz = 0.15

#stren = 10 - 0.01 - zfuz = 0 - what we want almost. 

# select the types
intTy = ct.c_int32
floatTy = ct.c_double
myDim = 3
#Firedrake stuff would go here.
#

mesh = UnitCubeMesh(4, 4, 4)
dim = mesh.geometric_dimension()
space = FunctionSpace(mesh, "Lagrange", 2)
f = Function(space)

x,y,z = SpatialCoordinate(mesh)
f = interpolate((x-0.5) * (x-0.5) + (y-0.5)*(y-0.5) + (z-0.5) * (z-0.5), space)

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
numPoints = 1000
newXc = np.random.uniform(low=-4.0, high=4, size=numPoints)
newYc = np.random.uniform(low=-4.0, high=4, size=numPoints)
newZc = np.random.uniform(low=-4.0, high=4, size=numPoints)
# if dim == 3:
#     points = [[0.8, 0.1, 0.5], [0.3, 0.3, 0.3]]#list(map(list, zip(*(xPoints, yPoints, zPoints)[0:dim]))) # list([[-2.42058914,  3.01164183, -7.03902682]])#
# elif dim == 2:
    #points = [[0.8, 0.2]]
points = list(map(list, zip(*(newXc, newYc, newZc)[0:dim])))
#points = [(0.05 * x + 0.025, 0.05 *y + 0.025, 0.05 * z + 0.025) for x in range(20) for y in range(20) for z in range(20)]
#points = [[0.225,0.525,0.575]]
pointsNrrd = "points.nrrd"
dataPoints = np.array(points, dtype="float64")
kindString = "{0}-vector".format(dim)
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)
# print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
outFileName = "pos"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]], "eps" :[floatTy(eps)], "rad" :[floatTy(rad)], "v0" : [floatTy(iso)], "fStren": [floatTy(fStren)], "fBias": [floatTy(fBias)], "zfuz" : [floatTy(zfuz)]}
outputs = [("pos", 1, outFileName)]
namedInputs = {"ipos": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs, verbose=True, workers=None, shutdown=False)
render("pos", "normed", dim=3)
exit(0)
#render(outFileName, outFileName, dim=dim)
#File(outFileName + ".pvd").write(f)
