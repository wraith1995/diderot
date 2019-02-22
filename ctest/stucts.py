import ctypes as ct
from ctypes import POINTER, c_int, c_double, c_void_p, c_float, c_int32, c_uint


floatString = "float32"
ctylesFloat = c_float

intString = "int32"
ctylesInt = c_int32

def makeArrayType(array, ty):
    array.flatten()
    dataType = array.ctypes.data_as(POINTER(ty))
    return(dataType)

def makeMeshType():
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
    return(_CFunction)

def makeSpaceType(meshTy):
    class _CFunction(ct.Structure):
        """C struct that represents a base mesh"""
        _fields_ = [
            ("indexMap", POINTER(ctylesInt)),
            ("mesh", meshTy)
        ]
    return(_CFunction)


def makeFuncType(spaceTy):
    class _CFunction(ct.Structure):
        """C struct that represents a base mesh"""
        _fields_ = [
            ("coordMap", POINTER(ctylesFloat)),
            ("space", spaceTy)
        ]
    return(_CFunction)



