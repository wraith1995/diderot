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
import argparse
import numpy as np
import nrrd_utils as nu
parser = argparse.ArgumentParser(description='Run a diderot program')
parser.add_argument("-s", action='store_true', default=False)
parser.add_argument("--bot", default=0.0, type=float)
parser.add_argument("--top", default=2.0, type=float)
parser.add_argument("--n", default=5, type=int)
parser.add_argument("--d", default=5, type=int)
parser.add_argument("--num", default=100, type=int)
args = parser.parse_args()

# select the types
intTy = ct.c_int32
floatTy = ct.c_double



#Firedrake stuff would go here.
n = args.n
mesh = UnitCubeMesh(n, n, n)
space = FunctionSpace(mesh, "Lagrange", args.d)
spacep = VectorFunctionSpace(mesh, "Lagrange", args.d - 1, dim=3)
f = Function(space)
numberOfDofs = len(f.dat.data)
# type to get the same
newDofs = np.array(np.random.uniform(low=args.bot, high=args.top, size=numberOfDofs),  dtype=floatTy)
f = Function(space, val=newDofs)
gradF = project(grad(f), spacep)

f = interpolate(Expression("(x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5)"), space)
gradF = interpolate(grad(f), spacep)
gradF = interpolate(Expression(("2*x[0] - 1.0", "2*x[1] - 1.0", "2*x[2] - 1.0")), spacep)

# floating point b/c sparsity or something else
# basis funcs b/c sparsity?
#generate points:
newXc = np.random.uniform(size=args.num)
newYc = np.random.uniform(size=args.num)
newZc = np.random.uniform(size=args.num)
points = map(list, zip(newXc, newYc, newZc))
pointsNrrd = "points.nrrd"
dataPoints = np.array(list(points))
kindString = "3-vector"
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)

def makeResult(point):
    cell = mesh.locate_cell(point)
    if cell is None:
        cell = -1
        world = float("nan")
        vec = [world, world, world]
    else:
        world = f.at(point)
        vec = gradF.at(point)
    return((cell, world, vec))

fResults = list(map(makeResult, dataPoints))

   

#translates to make real results...
#get the points
#compare! yay!

# build json
jsonFile = "evalProg.json"
dataFile = "evalProg.dill"
#fb.spaceToJson(space, jsonFile, refCellDefault="simplex")
# build data
(preFemArgs, femArgs) = fb.passAll(f, intTy, floatTy, geometric=not(args.s))
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = {"cube": [femArgs[0]], "space": [femArgs[1]], "data": [femArgs[2]]}
namedInputs = {"points": pointsNrrd}
outputs = [("c", 1, "c"), ("translated", 1, "translated"),
           ("result", 1, "result"), ("gradResult", 1, "gradResult")]
program.go(inputs, outputs, namedInputs=namedInputs)

cells = nu.get("c_0.nrrd")
refVecs = nu.get("translated_0.nrrd")
resultValues = nu.get("result_0.nrrd")
gradResult = nu.expectedOrder(nu.get("gradResult_0.nrrd"))

ourResults = zip(cells, resultValues, gradResult)

eps = 0.00000001
tests = zip(fResults, ourResults)
#print(list(fResults), list(ourResults))
cellsT = 0
valsT = 0
gradsT = 0 
for (idx, test) in enumerate(tests):
    print(idx)
    if test[0][0] != test[1][0]:
        print("Cells are different at {0}".format(idx))
        print("Got fcell {0} and cell {1}".format(test[0][0], test[1][0]))
        cellsT+=1
    if abs(test[0][1] - test[1][1]) > eps:
        print("Values are too far aparent at {0}".format(idx))
        print("Got fvalue {0} and value {1}".format(test[0][1], test[1][1]))
        valsT+=1
    fgrad = test[0][2]
    ograd = test[1][2]
    errs = [abs(a - b) for (a, b) in zip(fgrad, ograd)]

    print(errs, fgrad, ograd)
    tests = [e > eps for e in errs]
    if any(tests):
        gradsT+=1
        print("Grads are to far appart at {0}".format(idx))
        print("Got fgrad {0} and grad {1}".format(fgrad, ograd))
        for (idxp, (err, errTest)) in enumerate(zip(errs, tests)):
            if errTest:
                print("At idxp, the {0}th entry has error {1}".format(idxp, err))


print("cell errors: {0}".format(cellsT))
print("value errors: {0}".format(valsT))
print("grad errors: {0}".format(gradsT))
