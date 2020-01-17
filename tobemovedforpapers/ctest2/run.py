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
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')

#Firedrake stuff:
mesh = UnitCubeMesh(1,1,1)
space = FunctionSpace(mesh, "Lagrange", 2)
f = Function(space)
f = interpolate(Expression("(x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5)"), space)
#outfile = File("output.pvd")
#outfile.write(f)

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
femArgsNamed = {"a": [femArgs[0]], "b" : [femArgs[1]], "c" : [femArgs[2]]}
#femArgsNamed = {"a": [femArgs[0]], }
outputs = [("rgba", 1, "rgba")] #[("rayCellInter", 2, "rayCellInter"), ("spaceInter", 2, "spaceInter")]
program.go(femArgsNamed, outputs)

