import nrrd
import vtk
import numpy as np
import functools

# debug-able
# polyline-tube
# 2d option
# Particle/cube/whatever
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
            for idx in range(0, limit):
                result[j][k].append([seqData[offset + idx][i]
                                     for i in range(dim)])



    res = result.reshape(controlData.shape[1] * controlData.shape[2])
    return(res)

def read2dLines(controlFile, seqFile, z):
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
            for idx in range(0, limit):
                extra = [seqData[offset + idx][i]
                         for i in range(2)] + [z]
                result[j][k].append(extra)
    res = result.reshape(controlData.shape[1] * controlData.shape[2])
    return(res)


def readSegs(controlFile, seqFile, z=None):
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
            for idx in range(0, limit):
                extra1 = [seqData[offset + idx][i]
                          for i in range(3)]
                extra2 = [seqData[offset + idx + 3][i]
                          for i in range(3)]
                if z is not None:
                    extra1[2] = z
                    extra2[2] = z
                result[j][k].append([extra1, extra2])
    res = result.reshape(controlData.shape[1] * controlData.shape[2])
    resultp = []
    for r in res:
        resultp.extend(r)
        print(resultp)
    return(resultp)

def readPoints(controlFile, seqFile, z=None):
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
            for idx in range(0, limit):
                extra1 = [seqData[offset + idx][i]
                          for i in range(3)]
                if z is not None:
                    extra1[2] = z
                result[j][k].append(extra1)
    res = result.reshape(controlData.shape[1] * controlData.shape[2])
    resultp = []
    for r in res:
        resultp.append(r)
    return(resultp)

def renderLines(data, resultFile):
    lengths = map(len, data)
    totalLen = sum(lengths)
    vtkPoints = vtk.vtkPoints()
    vtkPoints.SetNumberOfPoints(totalLen)
    vLines = vtk.vtkCellArray()
    vertices = vtk.vtkCellArray()
    j = 0
    for plList in data:
        numPoints = len(plList)
        vLines.InsertNextCell(numPoints)
        for p in plList:
            vtkPoints.SetPoint(j, p[0], p[1], p[2])
            vLines.InsertCellPoint(j)
            vertices.InsertNextCell(1)
            vertices.InsertCellPoint(j)
            j += 1

    polygon = vtk.vtkPolyData()
    polygon.SetPoints(vtkPoints)
    polygon.SetVerts(vertices)
    polygon.SetLines(vLines)
    polygonMapper = vtk.vtkPolyDataMapper()
    polygonMapper.SetInputData(polygon)
    polygonMapper.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polygon)
    writer.SetFileName(resultFile)
    writer.Update()

def renderPoints(data, resultFile):
    vtkPoints = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    for (idxp, p) in enumerate(data):
        idx = vtkPoints.InsertNextPoint(p)
        vertices.InsertNextCell(idxp)
        vertices.InsertCellPoint(idx)
    result = vtk.vtkPolyData()
    result.SetPoints(vtkPoints)
    result.SetVerts(vertices)
    
    polygonMapper = vtk.vtkPolyDataMapper()
    polygonMapper.SetInputData(result)
    polygonMapper.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(result)
    writer.SetFileName(resultFile)
    writer.Update()

#Spheres: https://vtk.org/Wiki/VTK/Examples/Python/GeometricObjects/Display/Sphere
#https://gitlab.kitware.com/paraview/paraview/issues/17935
#The above all need renders
#Ask gordon about particles to do circles


def lines(name):
    controlFile = name + "_0.nrrd"
    seqFile = name + "_1.nrrd"
    result = readLines(controlFile, seqFile)
    renderLines(result, name + "_lines.vtk")
    return(result)


def lines2d(name):
    controlFile = name + "_0.nrrd"
    seqFile = name + "_1.nrrd"
    result = read2dLines(controlFile, seqFile)
    renderLines(result, name + "_2d_lines.vtk")
    return(result)


def segs(name):
    controlFile = name + "_0.nrrd"
    seqFile = name + "_1.nrrd"
    result = readSegs(controlFile, seqFile)
    renderLines(result, name + "_segs.vtk")
    return(result)


def segs2d(name, z):
    controlFile = name + "_0.nrrd"
    seqFile = name + "_1.nrrd"
    result = readSegs(controlFile, seqFile, z=z)
    renderLines(result, name + "_2d_segs.vtk")
    return(result)


def points(name, z=None):
    controlFile = name + "_0.nrrd"
    seqFile = name + "_1.nrrd"
    result = readPoints(controlFile, seqFile, z=z)
    renderPoints(result, name + "_points.vtk")
    return(result)
