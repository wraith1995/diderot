import sys
sys.path
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
from runNormals import getNormals

meshDatas = buildFiredrakeTetMesh("meshfiles/cs1.msh", 3, linearGmshBackup="lin_meshfiles/cs1.msh")
linMesh = meshDatas[-1]
curvedMesh = meshDatas[0]
mesh = curvedMesh

fStrTh = 24.0
fBias = 0.1
tipd = 0.1
fsEps = 0.1
geoEps = 0.1
mvmtEps = 0.1
rpcEps = 0.01
pcmvEps = 0.3

verb = 0
sfs = 0.5
hist = 0.5
pcp = 5



#stren = 10 - 0.01 - zfuz = 0 - what we want almost. 

# select the types
intTy = ct.c_int32
floatTy = ct.c_double
myDim = 3
#Firedrake stuff would go here.
#

dim = mesh.geometric_dimension()
space = FunctionSpace(mesh, "Lagrange", 6)
f = Function(space)

x,y,z = SpatialCoordinate(mesh)
#f = interpolate(y*y*x + z*z, space)
f = interpolate(z*z*(sin(x*x + y*y + z*z)), space)

# print(f.at([2.13336,1.92878,1.31913]))
# exit(0)

# import matplotlib.pyplot as plt
# plot(mesh, surface=True)
# plt.show()

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
numPoints =20000
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
#nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)
# print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
outFileName = "pos"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]],
          "fStrTh": [floatTy(fStrTh)], "fBias": [floatTy(fBias)], "fsEps": [floatTy(fsEps)],
          "geoEps": [floatTy(geoEps)], "mvmtEps": [floatTy(mvmtEps)], "rpcEps": [floatTy(rpcEps)],
          "pcmvEps": [floatTy(pcmvEps)], "tipd": [floatTy(tipd)]}
outputs = [("_pos", 1, outFileName)]
namedInputs = {"ipos": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs, verbose=True, workers=None, shutdown=False)
getNormals(f, "pos_0.nrrd")
render("pos", "normedp", dim=3, normalsFile="normals")
exit(0)
#render(outFileName, outFileName, dim=dim)
#File(outFileName + ".pvd").write(f)
