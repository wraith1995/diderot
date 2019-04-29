import sys
sys.path.append("../symb")
import loadSc
import sympy as sp
import ctypes as ct
import os
import passing
import pickle
import sys
import argparse
import stucts as st
import numpy as np
parser = argparse.ArgumentParser(description='Run a diderot program')
import math
xVar = loadSc.x
yVar = loadSc.y
zVar = loadSc.z

def findFunctions(p):
    grad = [sp.diff(p, xVar), sp.diff(p, yVar), sp.diff(p, zVar)]
    print("grad:", grad)
    H = np.array([[float(sp.diff(g, xVar)), float(sp.diff(g, yVar)), float(sp.diff(g, zVar))] for g in grad])
    print("H:",H)
    eigenValues, eigenVectors = np.linalg.eig(H)
    print(eigenValues)
    idx = eigenValues.argsort()
    eigenValues = eigenValues[idx]
    print(eigenValues)
    print(eigenVectors)
    eigenVectors = eigenVectors[:,idx]
    print(eigenVectors[0])
    print(eigenVectors[1])
    print(np.dot(eigenVectors[0], grad))
    print(np.dot(eigenVectors[1], grad))

a = 1.0
b = 2.0
c = 3.0
theta = 30.0
A0 = np.identity(3)
Az = np.array([[math.cos(theta), -math.sin(theta), 0], [math.sin(theta), math.cos(theta), 0], [0, 0, 1]])
Ay = np.array([[math.cos(theta), 0, math.sin(theta)], [0, 1, 0], [-math.sin(theta), 0, math.cos(theta)]])
Ax = np.array([[1.0, 0, 0.0], [0, math.cos(theta), -math.sin(theta)], [0, math.sin(theta), math.cos(theta)]])
A = Ax
y = A.dot([loadSc.x, loadSc.y, loadSc.z])
#y = [loadSc.x, loadSc.y, loadSc.z]
p = -0.5 * (a*loadSc.x**2 + b*loadSc.y**2 + c*loadSc.z**2)
print("P:",p)
print("Y:",y)
print("A",A)
pPrime = p.xreplace({loadSc.x : y[0], loadSc.y: y[1], loadSc.z: y[2]})
print(pPrime)
findFunctions(pPrime)
oneCoeff = loadSc.mono_coeffs(sp.Poly(pPrime), loadSc.monos)

# select the types
intTy = ct.c_int32
floatTy = ct.c_double

pklfile = "scarray.pkl"
data = pickle.load(open(pklfile, "rb"))



builder = st.makeAllTypes(intTy, floatTy)
z = data[0][9].shape
reset = list(data[0])
print(z)
reset[9] = np.zeros(z)
for i in range(len(oneCoeff)):
    reset[9][i] = oneCoeff[i]
args = builder[-1](*reset)

#main selection


# python camera control:
camEyePy = [0.01, -3.1, 0.01]
camAtPy = [0.0, 0.0, 0.0]
camUpPy = [0.0, 0.0, 1.0]
lightVsp = camEyePy  
camFovPy = 15
iresUPy = 300
iresVPy = 300
camNearPy = -2
camFarPy = 4.0

initDict = {"camEye": [(floatTy*3)(*camEyePy)],
            "camAt": [(floatTy*3)(*camAtPy)],
            "camUp": [(floatTy*3)(*camUpPy)],
            "camFOV": [floatTy(camFovPy)],
            "camNear": [floatTy(camNearPy)],
            "camFar": [floatTy(camFarPy)],
            "iresU": [intTy(iresUPy)],
            "iresV": [intTy(iresVPy)]}


#python feature control
fStrThPy = 0.0
thickPy = 0.01
featureDict = {"fStrTh" : [floatTy(fStrThPy)], "thick" : [floatTy(thickPy)]}

cellChoice = 0 #int(sys.argv[1])
print("Cell {0}".format(cellChoice))
inputs = {"block": [args[0]],
          "space": [args[1]],
          "U": [args[2]],
          "rayStep": [floatTy(0.005)],
          "cellChoice": [intTy(cellChoice)],
          **initDict, **featureDict}
outputs = [("rgba", 1, "rgba{0}".format(cellChoice))]
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
program.go(inputs, outputs, verbose=True)
