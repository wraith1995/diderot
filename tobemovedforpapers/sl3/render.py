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


def compareError(trueFile, errorFile, dim=3):

    trueData = readLines(trueFile + "_0.nrrd", trueFile + "_1.nrrd", dim=dim)
    errorData = readLines(errorFile + "_0.nrrd", errorFile + "_1.nrrd", dim=dim)
    lines = trueData.shape[0]
    #errors = (lenDiff, maxDiff, lastDiff, avDiff)
    errors = []
    for line in range(lines):
        trueLine = trueData[line]
        errorLine = errorData[line]
        maxLine = min(len(trueLine), len(errorLine))
        lenDiff = len(trueLine) - len(errorLine)
        error = abs(np.array(trueLine[0:maxLine]) - np.array(errorLine)[0:maxLine])
        maxError = np.max(error)
        lastDiff = np.max(error[-1])
        avDiff = np.mean(error)
        errors.append((lenDiff, maxError, lastDiff, avDiff))
    return(errors)



def compareErrorSingle(trueFile, errorFile):
    return(compareError(trueFile, errorFile)[0])

def compareToBackSolve(solFunc, timeStep, startPoint, line):
    steps = len(line)
    truth = np.array([solFunc(timeStep * n, startPoint) for n in range(steps)])
    errors = abs(lin - truth)
    return(np.max(errors), np.max(error[-1]), np.mean(error))

def compareToBackSolve(errorFile, pointsList, solFunc, timeStep):
    errorData = readLines(errorFile + "_0.nrrd", errorFile + "_1.nrrd", dim=dim)
    errors = [compareToBackSolve(solFunc, timeStep, point, line) for (point, line) in zip(pointsList, errorData)]
    return(errors)

def compareToBackSolveSingle(errorFile, pointsList, solFunc, timeStep):
    return(compareToBackSolve(errorFile, pointsList, solFunc, timeStep)[0])
    
    
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

def render(fileNameIn, fileNameOut, dim=3):
    data = readLines(fileNameIn + "_0.nrrd", fileNameIn + "_1.nrrd", dim=dim)
    writeToVtk(data, fileNameOut + ".vtk")
    
