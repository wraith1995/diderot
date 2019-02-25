from firedrake import *
import ctypes as ct
import firedrake_build as fb
import os
import passing
import dill

#Firedrake stuff:
mesh = UnitSquareMesh(1,1)
space = FunctionSpace(mesh, "Lagrange", 4)
f = Function(space)
f = interpolate(Expression("x[0]*x[0] + x[1]*x[1]"), space)

#Diderot stuff:
#choose types...
intTy = ct.c_int32
floatTy = ct.c_double
# build json
jsonFile = "test.json"
dataFile = "ugg.dill"
if os.path.exists(jsonFile):
    pass
else:
    fb.spaceToJson(space, "test.json", refCellDefault="simplex")
# build data
(preFemArgs, femArgs) = fb.passAll(f, intTy, floatTy)
# print(preFemArgs)
programNameArg = "justTypes"
nameSpaceArg = "justTypes"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)

femArgsNamed = {"a": [femArgs[0]], "b": [femArgs[1]], "c": [femArgs[2]]}
outputs = ["pos"]
program.go(femArgsNamed, outputs)

