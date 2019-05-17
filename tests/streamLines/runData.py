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
parser = argparse.ArgumentParser(description='Run a diderot program')
parser.add_argument("-s", action='store_true', default=False)
args = parser.parse_args()



# select the types
intTy = ct.c_int32
floatTy = ct.c_double

#Firedrake stuff would go here.
mesh = UnitCubeMesh(16, 16, 16)
dim = mesh.geometric_dimension()
space = VectorFunctionSpace(mesh, "Lagrange", 2, dim=dim)
f = Function(space)


if dim == 2:
    x,y = SpatialCoordinate(mesh)
    f = interpolate(as_vector([-x,y]), space)
elif dim == 3:
    M = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype="float64")
    Mtemp = np.asfarray(M, dtype="float64")
    Mptr = M.ctypes.data_as(ctypes.POINTER(floatTy))
    x,y,z = SpatialCoordinate(mesh)
    Mx = M.dot([x, y, z])
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
points = [[0.8, 0.1, 0.5], [0.3, 0.3, 0.3]]#list(map(list, zip(*(xPoints, yPoints, zPoints)[0:dim]))) # list([[-2.42058914,  3.01164183, -7.03902682]])#
pointsNrrd = "points.nrrd"
dataPoints = np.array(points, dtype="float64")
kindString = "{0}-vector".format(dim)
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)
# print(preFemArgs)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"meshData": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]], "M" : [Mptr]}
outputs = [("stream", 2, "stream"), ("newStream", 2, "newStream")]
namedInputs = {"startPoints": pointsNrrd}
program.go(inputs, outputs, namedInputs=namedInputs)

#rendering code:
def readScalarArary(fileName):
    readData, header = nrrd.read(fileName)
    return(readData)


def readLines(controlFile, seqFile, dim=3):
    seqData, seqHeader = nrrd.read(seqFile)
    size = functools.reduce(lambda x,y: x*y, seqData.shape, 1)
    seqData = seqData.T
    controlData, controlHeader = nrrd.read(controlFile)
    
    result = np.empty(controlData.shape[1:], dtype=object)
    #uhh---not general...?
    print("Shape:", controlData.shape)
    for j in range(controlData.shape[1]):
        for k in range(controlData.shape[2]):
            offset = controlData[0][j][k]
            limit = controlData[1][j][k]
            result[j][k] = []
            # if j != 0 or k != 0:
            #     continue
            for idx in range(0, limit):
                #  print(limit)
                result[j][k].append([seqData[offset + idx][i]
                                     for i in range(dim)])



    res = result.reshape(functools.reduce(lambda x, y: controlData.shape[1:]))
    return(res)
data = readLines("stream_0.nrrd", "stream_1.nrrd", dim=3)
writeToVtk(data, "test.vtk")

def writeToVtk(pointData, vtkFile):
    points = pointData
    lengths = map(len, points)
    totalLen = sum(lengths)
    vtkPoints = vtk.vtkPoints()
    vtkPoints.SetNumberOfPoints(totalLen)
    lines = vtk.vtkCellArray()
    vertices = vtk.vtkCellArray()
    j = 0
    for plList in points:
        numPoints = len(plList)
        lines.InsertNextCell(numPoints)
        for p in plList:
            vtkPoints.SetPoint(j, p[0], p[1], p[2])
            lines.InsertCellPoint(j)
            vertices.InsertNextCell(1)
            vertices.InsertCellPoint(j)
            j += 1

    polygon = vtk.vtkPolyData()
    polygon.SetPoints(vtkPoints)
    polygon.SetVerts(vertices)
    polygon.SetLines(lines)
    polygonMapper = vtk.vtkPolyDataMapper()
    polygonMapper.SetInputData(polygon)
    polygonMapper.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polygon)
    writer.SetFileName(vtkFile)
    writer.Update()
