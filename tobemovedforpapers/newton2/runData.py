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


eps = 0.005
rad = 0.05
iso = 0.0 #0.25*0.25


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
space = FunctionSpace(mesh, "Lagrange", 4)
f = Function(space)

x,y,z = SpatialCoordinate(mesh)
xp = (x - 1.0/2.0)
yp = (y - 1.0/2.0)
zp = (z - 1.0/2.0)
steiner = xp*xp*yp*yp + yp*yp*zp*zp + xp*xp*zp*zp + xp*yp*zp
mitchel = 4* (xp*xp*xp*xp + (yp*yp + zp*zp) * (yp*yp + zp*zp)) + 17*(xp*xp * (yp*yp + zp*zp) - 20 * (xp*xp + yp*yp + zp*zp)) + 17 #stupid
torus = (xp*xp + yp*yp + zp*zp + 0.25*0.25 - 0.125*0.125)*(xp*xp + yp*yp + zp*zp + 0.25*0.25 - 0.125*0.125) - 4 * 0.25*0.25 * (xp * xp + yp * yp)
planeTest = xp+yp+zp #strange...
f = interpolate(planeTest, space)

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
newXc = np.random.uniform(low=0, high=1, size=numPoints)
newYc = np.random.uniform(low=0, high=1, size=numPoints)
newZc = np.random.uniform(low=0, high=1, size=numPoints)
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
outFileName = "ugg1"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]], "eps" :[floatTy(eps)], "rad" :[floatTy(rad)], "v0" : [floatTy(iso)]}
#("pos", 1, outFileName)
outputs = [ ("_pos", 1, "ugg1")]
namedInputs = {"ipos": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs, verbose=True, workers=None)
render(outFileName, outFileName, dim=dim)
#File(outFileName + ".pvd").write(f)
