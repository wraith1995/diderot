import ctypes as ct
from ctypes import POINTER, c_int, c_double, c_void_p, c_float, c_int32, c_uint
import numpy as np


def makeArrayType(array, ty):
    array.flatten()
    if ty == c_float or ty == c_double:
        array = np.asfarray(array, dtype=ty)
    dataType = array.ctypes.data_as(POINTER(ty))
    return(dataType)

def makeVoidPointerTuple(val):
    return(val,ct.cast(ct.pointer(val), ct.c_void_p))

def makeMeshType(ctylesInt, ctylesFloat):
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
        ]
    def build(indexMap, coordMap, dim, mapDim, numCells, index, con):
        ty = _CFunction()
        ty.indexMap = makeArrayType(indexMap, ctylesInt)
        ty.coordMap = makeArrayType(coordMap, ctylesFloat)
        ty.dim = dim
        ty.mapDim = mapDim
        ty.numCells = numCells
        ty.index = ct.cast(index, ct.c_void_p)
        ty.con = makeArrayType(con, ctylesInt)
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


def makeAllTypes(ctylesInt, ctylesFloat):
    (meshTy, buildMesh) = makeMeshType(ctylesInt, ctylesFloat)
    (spaceTy, buildSpace) = makeSpaceType(meshTy, ctylesInt)
    (funcTy, buildFunc) = makeFuncType(spaceTy, ctylesFloat)
    
    def buildAll(meshIndexMap, meshCoordMap,
                 dim, meshMapDim, numCell,
                 sIndex, con,
                 spaceIndexMap,
                 spaceDim,
                 funcCoordMap):
        (meshVal, meshValPtr) = buildMesh(meshIndexMap, meshCoordMap, dim, meshMapDim,
                            numCell, sIndex, con)
        (spaceVal, spaceValPtr) = buildSpace(spaceIndexMap, spaceDim, meshVal)
        (funcVal, funcValPtr) = buildFunc(funcCoordMap, spaceVal)
        return((meshValPtr, spaceValPtr, funcValPtr))
    return(buildMesh, buildSpace, buildFunc, buildAll)




