import numpy as np
import numpy.ctypeslib as npct
import ctypes as c
from ctypes import POINTER, c_int, c_double, c_void_p, c_float, c_int32,

floatString = "float64"
ctylesFloat = c_double

intString = "int32"
ctylesInt = c_int32



def makeMeshStruct():
    class MeshStruct(ctypes.Structure):
        """C struct that represents a base mesh"""
        _fields_ = [
            ("indexMap", POINTER(ctylesInt)),
            ("coordMap", POINTER(ctylesFloat)),
            ("dim", ctylesInt),
            ("mapDim", ctylesInt),
            ("numCells", ctylesInt)
        ]
    return(makeMeshStruct)


meshType = makeMeshStruct()
