from firedrake import *
import ctypes as ct
import firedrake_build as fb
import os
import passing
import dill
import sys
sys.path
sys.path.append('/path/to/the/example_file.py')
#Firedrake stuff:
mesh = UnitCubeMesh(1,1,1)
space = FunctionSpace(mesh, "Lagrange", 2)
f = Function(space)
f = interpolate(Expression("x[0]*x[0] + x[1]*x[1] + x[2]*x[2]"), space)


#Diderot stuff:
#choose types...
intTy = ct.c_int32
floatTy = ct.c_double

# build json
saveToDill = False
jsonFile = "test1.json"
dataFile = "ugg1.dill"
if os.path.exists(jsonFile):
    pass
else:
    fb.spaceToJson(space, jsonFile, refCellDefault="simplex")
# build data
(preFemArgs, femArgs) = fb.passAll(f, intTy, floatTy, geometric=not(saveToDill))
if saveToDill:
    with open(dataFile, "wb+") as f:
        dill.dump(preFemArgs, f)
    exit(0)
else:
    pass
# print(preFemArgs)
programNameArg = "justTypes"
nameSpaceArg = "justTypes"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)

femArgsNamed = {"a": [femArgs[0]], "b": [femArgs[1]], "c": [femArgs[2]]}
outputs = ["pos"]
program.go(femArgsNamed, outputs)

