import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import nrrd
import numpy as np
import sympy as sp
import nrrd_utils as nu
import functools
import vtk

def make3d(x):
    if len(x) == 2:
        return([x[0], x[1], 0.0])
    elif len(x) == 3:
        return(x)
    else:
        raise Exception("Stupid")

def writePointsToVtk(data, vtkFile, normals=None, scalars=None):
    vtkPoints = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    numPoints = len(data)
    pointNormalsArray = vtk.vtkDoubleArray()
    scalarArray = vtk.vtkDoubleArray()
    scalarArray.SetNumberOfComponents(1)
    pointNormalsArray.SetNumberOfComponents(3)
    countAdded = 0
    for (idx, point) in enumerate(data):
        pointId = vtkPoints.InsertNextPoint(point)
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(pointId)
        countAdded += 1
    polygon = vtk.vtkPolyData()
    polygon.SetPoints(vtkPoints)
    polygon.SetVerts(vertices)
    if normals is not None:
        pointNormalsArray.SetNumberOfTuples(countAdded)
        for (idx, normal) in enumerate(normals):
            #z is 2nd, x is 3rd
            #z is 1st
            pointNormalsArray.SetTuple(idx, tuple(normal))
            
        polygon.GetPointData().SetNormals(pointNormalsArray)
    if scalars is not None:
        scalarArray.SetNumberOfTuples(countAdded)
        for (idx, datam) in enumerate(scalars):
            scalarArray.SetTuple(idx, tuple([datam]))
        polygon.GetPointData().SetScalars(scalarArray)
    polygonMapper = vtk.vtkPolyDataMapper()
    polygonMapper.SetInputData(polygon)
    polygonMapper.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polygon)
    writer.SetFileName(vtkFile)
    writer.Update()
    
    
#data = readLines("stream_0.nrrd", "stream_1.nrrd", dim=2)
#writeToVtk(data, "test.vtk")

def render(fileNameIn, fileNameOut, dim=2, normalsFile=None, scalarsFile=None):
    (data, _) = nrrd.read(fileNameIn + "_0.nrrd")
    reformulatedData = nu.expectedOrder(data)
    normals = None
    scalars= None
    if normalsFile is not None:
        (ndata, _) = nrrd.read(normalsFile + "_0.nrrd")
        (ndata1, _) = nrrd.read(scalarsFile + "_0.nrrd")
        normReform = nu.expectedOrder(ndata)
        normals = normReform
        scalars = ndata1
    writePointsToVtk(reformulatedData, fileNameOut + ".vtk", normals=normals, scalars=scalars)
    

#render("ugg1", "ugg")
