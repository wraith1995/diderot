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


render("ugg1", "normed.vtk", dim=3, normalsFile="normals")
exit(0)
meshDatas = buildFiredrakeTetMesh("meshfiles/fuck3.msh", 3)
linMesh = meshDatas[-1]
curvedMesh = meshDatas[0]
mesh = curvedMesh
# T = curvedMesh.coordinates
# fs = T.function_space()
# m = T.function_space().mesh()
# V1 = VectorFunctionSpace(m, "Lagrange", 3)
# V0 = VectorFunctionSpace(m, "Lagrange", 1)

# zeroFuncSpace = VectorFunctionSpace(m, "DG", 0, dim=3)
# v = TestFunction(zeroFuncSpace)


# linT = interpolate(T, V0)
# deLinT = interpolate(linT, V1)
# err = interpolate(T - deLinT, V1)
# perCellErrr = inner(v, err)*dx(linMesh)
# errPerCell = assemble(perCellErrr)
# q =  [idx for (idx, v) in enumerate(errPerCell.dat.data) if np.linalg.norm(v) > 0.001]
# for idx in q:
#     for n in fs.cell_node_map().values[idx]:
#         print(np.linalg.norm(T.dat.data[n]))
    
# exit(0)



eps = 0.005
rad = 0.5
iso = 0.0
dim=3





# select the types
intTy = ct.c_int32
floatTy = ct.c_double

dim = mesh.geometric_dimension()
space = FunctionSpace(mesh, "Lagrange", 3)
f = Function(space)

x,y,z = SpatialCoordinate(mesh)
f = interpolate(x*x + y*y + z*z - 9.5*9.5, space)


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
numPoints = 10000
newXc = np.random.uniform(low=-10, high=10, size=numPoints)
newYc = np.random.uniform(low=-10, high=10, size=numPoints)
newZc = np.random.uniform(low=-10, high=10, size=numPoints)
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
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]], "eps" :[floatTy(eps)], "rad" :[floatTy(rad)], "v0" : [floatTy(iso)]}
#("pos", 1, outFileName)
outputs = [ ("_pos", 1, "ugg1")]
namedInputs = {"ipos": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs, verbose=True, workers=None, shutdown=False)
#render(outFileName, outFileName, dim=dim)
#File(outFileName + ".pvd").write(f)
getNormals(f, "ugg1_0.nrrd")
exit(0)
