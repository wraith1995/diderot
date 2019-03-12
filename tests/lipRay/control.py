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
from itertools import repeat
import nrrd_utils as nu

intTy = ct.c_int32
floatTy = ct.c_double

#Firedrake stuff would go here.
mesh = Mesh("test.msh")
space = FunctionSpace(mesh, "Lagrange", 2)
f = Function(space)
ue = Expression("1.0-(x[0]*x[0]+x[1]*x[1])",degree=2)
f = interpolate(ue, space)
fw = File("test.pvd")
fw.write(f)
u = TrialFunction(space)
v = TestFunction(space)
fp = interpolate(Expression("4.0",degree=2), space)
a = inner(grad(u),grad(v))*dx 
L = fp*v*dx
u0 = Constant(0.0)
bc = DirichletBC(space, u0, "on_boundary")
u = Function(space)
solve(a == L, u, bc )
fw = File("test1.pvd")
fw.write(u)
g = Function(space)
g = interpolate(u-f, space)
fw = File("test2.pvd")
fw.write(g)


mesh = UnitSquareMesh(1,1)
space = FunctionSpace(mesh,"Lagrange",2)
f = Function(space)
g = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5))", degree=2)
f = interpolate(g,space)
g = f
fw = File("test3.pvd")
fw.write(g)


#
camEyePy = [-1, -1, 2.0]
camAtPy = [1.0, 1.0, 0.5]
camUpPy = [0.5, 0.5, 1.0]
camFovPy = 30
iresUPy = 500
iresVPy = 500
initDict = {"camEye": [(floatTy*3)(*camEyePy)],
            "camAt": [(floatTy*3)(*camAtPy)],
            "camUp": [(floatTy*3)(*camUpPy)],
            "camFOV": [floatTy(camFovPy)],
            "iresU": [intTy(iresUPy)],
            "iresV": [intTy(iresVPy)]}

(_, femArgs) = fb.passAll(g, intTy, floatTy)

# seg program:
programNameArgSeg = "seg/evalProg"
nameSpaceArgSeg = "evalProg"
inputs = {"a" : [femArgs[0]], **initDict}
outputs = [("rayCellInter", 2, "rayCellInter"), ("spaceInter", 2, "spaceInter"), ("cellInter", 2, "cellInter")]
segLib = ct.CDLL("./" + programNameArgSeg + ".so")
segProgram = passing.Library(segLib, nameSpace=nameSpaceArgSeg)
segProgram.go(inputs, outputs)


numCell = mesh.coordinates.dat.data.shape[0]
valsLArray = np.array(list(repeat(8.0, numCell)))
valsL = np.asfarray(valsLArray).ctypes.data_as(ct.POINTER(floatTy))
valsG = valsL

(_, femArgsp) = fb.passAll(g, intTy, floatTy,
                        extraData=[("L", valsL, ct.POINTER(floatTy)),
                                   ("G", valsG, ct.POINTER(floatTy))])
# load indecies
# load times
loadedIndeices = nu.get("rayCellInter_0.nrrd")
loadedIndeices = loadedIndeices.reshape(2 * loadedIndeices.shape[1] * loadedIndeices.shape[2])
indexStr = "index.nrrd"
nu.writeSequence(loadedIndeices, indexStr)
s = "cellInter_1.nrrd"


inputs = {"a": [femArgsp[0]], "b": [femArgsp[1]], "c": [femArgsp[2]], **initDict}
namedInputs = {"times": "rayCellInter_1.nrrd", "indecies": indexStr, "cells": s}
outputs = [("newCells", 2, "newCells"), ("intervals", 2, "intervals")]

programNameArgSub = "sub/evalProg"
nameSpaceArgSub = "evalProg"
subLib = ct.CDLL("./" + programNameArgSub + ".so")
subProgram = passing.Library(subLib, nameSpace=nameSpaceArgSub)
subProgram.go(inputs, outputs, namedInputs=namedInputs)




programNameArgfind = "find/evalProg"
nameSpaceArgfind = "evalProg"
findLib = ct.CDLL("./" + programNameArgfind + ".so")
findProgram = passing.Library(findLib, nameSpace=nameSpaceArgfind)
loadedIndeicesP = nu.get("intervals_0.nrrd")
loadedIndeicesP = loadedIndeicesP.reshape(2 * loadedIndeicesP.shape[1] * loadedIndeicesP.shape[2])
indexStrP = "index_0.nrrd"
nu.writeSequence(loadedIndeicesP, indexStrP)
inputs = {"a": [femArgsp[0]], "b": [femArgsp[1]], "c": [femArgsp[2]], **initDict}
namedInputs = {"times": "intervals_1.nrrd", "indecies": indexStrP,
               "cells": "newCells_1.nrrd"}
outputs = [("rgba", 1, "rgba"), ("out", 1, "out")]
findProgram.go(inputs, outputs, namedInputs=namedInputs)


# run that...
# do the same thing


