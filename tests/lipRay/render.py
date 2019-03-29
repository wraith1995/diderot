import nrrd
import vtk
import numpy as np
import functools
#from paraview.simple import *

def readScalarArary(fileName):
    readData, header = nrrd.read(fileName)
    return(readData)


def readLines(controlFile, seqFile, dim=3):
    seqData, seqHeader = nrrd.read(seqFile)
    size = functools.reduce(lambda x,y: x*y, seqData.shape, 1)
    seqData = seqData.T
    controlData, controlHeader = nrrd.read(controlFile)
    
    result = np.empty(controlData.shape[1:], dtype=object)
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



    res = result.reshape(controlData.shape[1] * controlData.shape[2])
    return(res)

#  steps = readScalarArary("steps_0.nrrd")

points = readLines("spaceInter_0.nrrd", "spaceInter_1.nrrd")
lengths = map(len, points)
totalLen = sum(lengths)
vtkPoints = vtk.vtkPoints()
vtkPoints.SetNumberOfPoints(totalLen)
lines = vtk.vtkCellArray()
vertices = vtk.vtkCellArray()
j = 0
for plList in points:
    numPoints = len(plList)
    # print("umm", plList)

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
writer.SetFileName("test.vtk")
writer.Update()
