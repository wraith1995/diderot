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
m = Mesh("lin_meshfiles/cs1.msh")

meshDatas = buildFiredrakeTetMesh("meshfiles/cs1.msh", 3, linearGmshBackup="lin_meshfiles/cs1.msh")
linMesh = meshDatas[-1]
curvedMesh = meshDatas[0]
mesh = curvedMesh

fStrTh = 26.0
fBias = 0.1
tipd = 0.5
fsEps = 0.01
geoEps = 0.1
mvmtEps = 0.01
rpcEps = 0.01
pcmvEps = 0.3
# fStrTh = 26.0
# fBias = 0.1
# tipd = 0.5
# fsEps = 0.01
# geoEps = 0.1
# mvmtEps = 0.01
# rpcEps = 0.01
# pcmvEps = 0.3

verb = 0
sfs = 0.5
hist = 0.5
pcp = 5
intTy = ct.c_int32
floatTy = ct.c_double
myDim = 3


dim = mesh.geometric_dimension()
space = FunctionSpace(mesh, "Lagrange", 6)
f = Function(space)



x,y,z = SpatialCoordinate(mesh)
#f = interpolate(y*y*x + z*z, space)
f = interpolate(z*z*(sin(x*x + y*y + z*z)), space)
#f = interpolate(z*(sin(x*x + y*y + z*z)), space)
#getNormals(f, "pos_0.nrrd")
#render("pos", "normedp", dim=3, normalsFile="normals", scalarsFile="stren")
#exit(0)

#getNormals(f, "test1.nrrd")
#exit(0)
#render("pos", "normedpp", dim=3, normalsFile="normals", scalarsFile="stren")


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


pointsNrrd = "points.nrrd"
# numPoints = 20000
# newXc = np.random.uniform(low=-4.0, high=4, size=numPoints)
# newYc = np.random.uniform(low=-4.0, high=4, size=numPoints)
# newZc = np.random.uniform(low=-4.0, high=4, size=numPoints)
# points = list(map(list, zip(*(newXc, newYc, newZc)[0:dim])))
# dataPoints = np.array(points, dtype="float64")
# kindString = "{0}-vector".format(dim)
# nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)

# # print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
outFileName = "pos"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]],
          "fStrTh": [floatTy(fStrTh)], "fBias": [floatTy(fBias)], "fsEps": [floatTy(fsEps)],
          "geoEps": [floatTy(geoEps)], "mvmtEps": [floatTy(mvmtEps)], "rpcEps": [floatTy(rpcEps)],
          "pcmvEps": [floatTy(pcmvEps)], "tipd": [floatTy(tipd)], "hist" : [floatTy(hist)]}
outputs = [("_pos", 1, outFileName)]
namedInputs = {"ipos": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs, verbose=True, workers=None, shutdown=False)
getNormals(f, "pos_0.nrrd")
render("pos", "normedp", dim=3, normalsFile="normals", scalarsFile="stren")
exit(0)
#render(outFileName, outFileName, dim=dim)
#File(outFileName + ".pvd").write(f)
