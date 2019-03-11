import ctypes as ct
from ctypes import POINTER, c_int, c_double, c_void_p, c_float, c_int32, c_uint
import numpy as np


def makeFiniteArray(array, ty):
    return ((ty * len(array))(*array))

def makeArrayType(array, ty):
    # array.flatten()
    arrayp = array  # array.flatten()
    if ty == c_float or ty == c_double:
        arrayp = np.asfarray(arrayp, dtype=ty)
    dataType = arrayp.ctypes.data_as(POINTER(ty))
    return(dataType)

def makeVoidPointerTuple(val):
    return(val,ct.cast(ct.pointer(val), ct.c_void_p))

def makeMeshType(ctylesInt, ctylesFloat, extraNamesAndTypes=[]):
    #  make extra fields
    sortedExtraFields = sorted(extraNamesAndTypes, key=lambda x: x[0])
    names = [x[0] for x in sortedExtraFields]
    print(sortedExtraFields)
    class _CFunction(ct.Structure):
        """C struct that represents a base mesh"""
        _fields_ = [
            ("indexMap", POINTER(ctylesInt)),
            ("coordMap", POINTER(ctylesFloat)),
            ("dim", ctylesInt),
            ("mapDim", ctylesInt),
            ("numCells", ctylesInt),
            ("index", c_void_p),
            ("con", POINTER(c_int))
        ] + sortedExtraFields
    def build(indexMap, coordMap, dim, mapDim, numCells, index, con,
              extraValues=[]):
        if len(extraValues) != len(names):
            raise Exception("mesh is mixing extra arguments")
        
        ty = _CFunction()
        ty.indexMap = makeArrayType(indexMap, ctylesInt)
        ty.coordMap = makeArrayType(coordMap, ctylesFloat)
        ty.dim = dim
        ty.mapDim = mapDim
        ty.numCells = numCells
        ty.index = ct.cast(index, ct.c_void_p)
        ty.con = makeArrayType(con, ctylesInt)
        for (name, val) in zip(names, extraValues):
            setattr(ty, name, val)
        return(makeVoidPointerTuple(ty))

    return(_CFunction, build)

def makeSpaceType(meshTy, ctylesInt):
    class _CFunction(ct.Structure):
        """C struct that represents a base mesh"""
        _fields_ = [
            ("indexMap", POINTER(ctylesInt)),
            ("mesh", meshTy)
        ]
    def build(indexMap, spaceDim, mesh):
        ty = _CFunction()
        ty.indexMap = makeArrayType(indexMap, ctylesInt)
        ty.mesh = mesh
        return(makeVoidPointerTuple(ty))
        
    return(_CFunction, build)


def makeFuncType(spaceTy, ctylesFloat):
    class _CFunction(ct.Structure):
        """C struct that represents a base mesh"""
        _fields_ = [
            ("coordMap", POINTER(ctylesFloat)),
            ("space", spaceTy)
        ]
    def build(coordMap, space):
        ty = _CFunction()
        ty.coordMap = makeArrayType(coordMap, ctylesFloat)
        ty.space = space
        return(makeVoidPointerTuple(ty))
    return(_CFunction, build)


def makeAllTypes(ctylesInt, ctylesFloat,
                 extraTys=[]):
    (meshTy, buildMesh) = makeMeshType(ctylesInt, ctylesFloat,
                                       extraNamesAndTypes=extraTys)
    (spaceTy, buildSpace) = makeSpaceType(meshTy, ctylesInt)
    (funcTy, buildFunc) = makeFuncType(spaceTy, ctylesFloat)
    
    def buildAll(meshIndexMap, meshCoordMap,
                 dim, meshMapDim, numCell,
                 sIndex, con,
                 spaceIndexMap,
                 spaceDim,
                 funcCoordMap, extraVals=[]):
        (meshVal, meshValPtr) = buildMesh(meshIndexMap, meshCoordMap, dim,
                                          meshMapDim,
                                          numCell, sIndex, con,
                                          extraValues=extraVals)
        (spaceVal, spaceValPtr) = buildSpace(spaceIndexMap, spaceDim, meshVal)
        (_, funcValPtr) = buildFunc(funcCoordMap, spaceVal)
        return((meshValPtr, spaceValPtr, funcValPtr))
    return(buildMesh, buildSpace, buildFunc, buildAll)




