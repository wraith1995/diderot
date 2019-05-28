import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import nrrd
import numpy as np
import sympy as sp
import nrrd_utils as nu
import functools
import vtk
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
        offset = controlData[0][j]
        limit = controlData[1][j]
        result[j] = []
        for idx in range(0, limit):
            result[j].append(np.array([seqData[offset + idx][i]
                                 for i in range(dim)]))



    res = np.array([np.array(x) for x in result.reshape(controlData.shape[1])])
    return(res)

def make3d(x):
    if len(x) == 2:
        return([x[0], x[1], 0.0])
    elif len(x) == 3:
        return(x)
    else:
        raise Exception("Stupid")

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
            pp = make3d(p)
            vtkPoints.SetPoint(j, pp[0], pp[1], pp[2])
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
#data = readLines("stream_0.nrrd", "stream_1.nrrd", dim=2)
#writeToVtk(data, "test.vtk")

def render(fileNameIn, fileNameOut, dim=2):
    data = readLines(fileNameIn + "_0.nrrd", fileNameIn + "_1.nrrd", dim=dim)
    writeToVtk(data, fileNameOut + ".vtk")
    
