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
import nrrd
parser = argparse.ArgumentParser(description='Run a diderot program')
import math

xVar = loadSc.x
yVar = loadSc.y
zVar = loadSc.z
a = 1.0
b = 2.0
c = 3.0
globalP = -(2 * xVar**3 + 5*yVar**3 + 2* zVar**3 + 2 * xVar * zVar)
#(loadSc.z*loadSc.x*loadSc.y)**3  - 0.5 * (a*loadSc.x**2 + b*loadSc.y**2 + c*loadSc.z**2) - 10 * loadSc.y**2*loadSc.x - 5 * loadSc.z**3 * loadSc.y

# select the types
intTy = ct.c_int32
floatTy = ct.c_double

def angTranslate(x):
    return(math.pi/180.0 * x)

def samplePolyToNrrd(poly, dataNrrd, samples=300, orig=[-0.6, -0.6, -0.6], end=[0.6, 0.6, 0.6]):
    f = sp.lambdify([xVar, yVar, zVar], poly)
    basis1 = ((end[0] - orig[0]) / (samples - 1), 0.0, 0.0)
    basis2 = (0.0, (end[1] - orig[1]) / (samples - 1), 0.0)
    basis3 = (0.0, 0.0, (end[2] - orig[2]) / (samples - 1))
    data = np.zeros((samples, samples, samples), dtype="float64")
    for i in range(samples):
        x = orig[0] + basis1[0] * i
        for j in range(samples):
            y = orig[1] + basis2[1] * j
            for k in range(samples):
                z = orig[2] + basis3[2] * k
                data[i][j][k] = float(f(x, y, z))
    options = {"kinds": ["space", "space", "space"],
               "encoding": "raw",
               "space dimension": 3,
               "space directions": [basis1, basis2, basis3],
               "space origin": orig,
               "centerings": ["node", "node", "node"]}
    nrrd.write(dataNrrd, data, options)
    text = """NRRD0006\n# Complete NRRD file format specification at:\n# http://teem.sourceforge.net/nrrd/format.html
type: double
dimension: 3
sizes: {0} {0} {0}
kinds: space space space
centers: cell cell cell
endian: little
encoding: raw
space dimension: 3
space directions: {1} {2} {3}
space origin: {4}
data file: {5}
byte skip: -1
    """.format(samples, tuple(basis1), tuple(basis2), tuple(basis3), tuple(orig), dataNrrd)
    # with open(outNrrd, "w+") as ffile:
    #     ffile.write(text)
                


def findFunctions(p):
    print("Looking at function:",p)
    grad = [sp.diff(p, xVar), sp.diff(p, yVar), sp.diff(p, zVar)]
    H = np.array([[(sp.diff(g, xVar)), (sp.diff(g, yVar)), (sp.diff(g, zVar))] for g in grad])
    print("Grad:", grad)
    print("Hess:",H)
    # eigenValues, eigenVectors = np.linalg.eig(H)
    # idx = eigenValues.argsort()
    # eigenValues = eigenValues[idx]
    # eigenVectors = eigenVectors[:,idx]
    # print("Looking for zeros given by:")
    # print(np.dot(eigenVectors[0], grad))
    # print(np.dot(eigenVectors[1], grad))
    # print("eigenvalues are:", eigenValues)


def makeData(xRot=0.0, yRot=0.0, zRot=0.0, datafile="data.nrrd"):
    a = 1.0
    b = 2.0
    c = 3.0
    xRot = angTranslate(xRot)
    yRot = angTranslate(yRot)
    zRot = angTranslate(zRot)
    Ax = np.array([[1.0, 0, 0.0], [0, math.cos(xRot), -math.sin(xRot)], [0, math.sin(xRot), math.cos(xRot)]])
    Ay = np.array([[math.cos(yRot), 0, math.sin(yRot)], [0, 1, 0], [-math.sin(yRot), 0, math.cos(yRot)]])
    Az = np.array([[math.cos(zRot), -math.sin(zRot), 0], [math.sin(zRot), math.cos(zRot), 0], [0, 0, 1]])
    A = Ax.dot(Ay).dot(Az)
    y = A.dot([loadSc.x, loadSc.y, loadSc.z])
    p = globalP
    pPrime = p.xreplace({loadSc.x : y[0], loadSc.y: y[1], loadSc.z: y[2]})
    print("Loading poly:", pPrime)
    samplePolyToNrrd(pPrime, datafile)
    findFunctions(pPrime)
    oneCoeff = loadSc.mono_coeffs(sp.Poly(pPrime), loadSc.monos)

    pklfile = "scarray.pkl"
    data = pickle.load(open(pklfile, "rb"))

    
    z = data[0][9].shape
    reset = list(data[0])
    reset[9] = np.zeros(z)
    print("Coeffs:",oneCoeff)
    for i in range(len(oneCoeff)):
        reset[9][i] = oneCoeff[i]
    #args = builder[-1](*reset)
    return(reset, A)

#main selection
builder = st.makeAllTypes(intTy, floatTy)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
#write nrrds
commandFile = open("commands.txt", "w+")
def compare(xRot=0.0, yRot=0.0, zRot=0.0):
    (reset, rot) = makeData(xRot=xRot, yRot=yRot, zRot=zRot)
    args = builder[-1](*reset)
    # python camera control:
    camEyePy = [0.0, -3.1, 0.0]
    camAtPy = [0.0, 0.0, 0.0]
    camUpPy = [0.0, 0.0, 1.0]
    lightVsp = camEyePy  
    camFovPy = 15
    iresUPy = 300
    iresVPy = 300
    camNearPy = -2
    camFarPy = 4.0
    rotArg = (floatTy * 9)(*(list(rot.flatten())))
    
    initDict = {"camEye": [(floatTy*3)(*camEyePy)],
                "camAt": [(floatTy*3)(*camAtPy)],
                "camUp": [(floatTy*3)(*camUpPy)],
                "camFOV": [floatTy(camFovPy)],
                "camNear": [floatTy(camNearPy)],
                "camFar": [floatTy(camFarPy)],
                "iresU": [intTy(iresUPy)],
                "iresV": [intTy(iresVPy)],
                "rot": [rotArg]}
    #python feature control
    fStrThPy = 0.0
    thickPy = 0.005
    anal = False
    featureDict = {"fStrTh" : [floatTy(fStrThPy)], "thick" : [floatTy(thickPy)]}

    cellChoice = 0 #int(sys.argv[1])
    # print("Cell {0}".format(cellChoice))
    inputs1 = {"block": [args[0]],
               "space": [args[1]],
               "U": [args[2]],
               "rayStep": [floatTy(0.005)],
               "cellChoice": [intTy(cellChoice)],
               "analytical": [ct.c_bool(False)],
               **initDict, **featureDict}
    inputs2 = {"block": [args[0]],
               "space": [args[1]],
               "U": [args[2]],
               "rayStep": [floatTy(0.005)],
               "cellChoice": [intTy(cellChoice)],
               "analytical": [ct.c_bool(True)],
               **initDict, **featureDict}
    ourNrrd = "ourResult_{0}_{1}_{2}".format(xRot, yRot, zRot)
    theirNrrd = "theirResult_{0}_{1}_{2}".format(xRot, yRot, zRot)
    outputs1 = [("rgba", 1, ourNrrd)]
    outputs2 = [("rgba", 1, theirNrrd)]
    program = passing.Library(library, nameSpace=nameSpaceArg)
    program.go(inputs1, outputs1, verbose=False, shutdown=False)
    print("ROt is:", rot)
    program.go(inputs2, outputs2, verbose=False, shutdown=False)
    n1 = "/home/teocollin/gitcode/diderot/tests/pagot/vr/" + ourNrrd + "_0"
    ourData, _ = nrrd.read(n1 + ".nrrd")
    n2 = "/home/teocollin/gitcode/diderot/tests/pagot/vr/" + theirNrrd + "_0"
    theirData, _ = nrrd.read(n2 + ".nrrd")
    renderString = "unu quantize -b 16 -i {0}.nrrd -o {0}.png"
    diffString = "unu 2op - {0}.nrrd {1}.nrrd -o {0}_prime.nrrd"
    appendString = "convert {0}.png {1}.png {2}.png -append {0}_diff.png"
    doString = "{0}; {1}; {2}; {3}; {4}".format(
        renderString.format(n1),
        renderString.format(n2),
        diffString.format(n1, n2),
        renderString.format(n1 + "_prime"),
        appendString.format(n1, n2, n1 + "_prime")
    )
    print("Command: {0}".format(doString))
    commandFile.write(doString)
    diff = ourData - theirData
    s = diff.shape
    total = s[0] * s[1] * s[2]
    diffp = diff.reshape(total)
    err = np.linalg.norm(diffp, ord=float("inf"))
    avErr = np.linalg.norm(diffp) / total
    print("With files {0} and {1}, there errors are ({2}, {3})".format(ourNrrd,
                                                                       theirNrrd, err, avErr))
    if (err > 0.001 and avErr > 0.000001 ):
        print("Error occured with file {0}".format(ourNrrd))
        
#compare(zRot=math.pi * (10.0/180.0))

angs = [x * 15 for x in range(0,5)]
for xRot in angs:
    for yRot in angs:
        for zRot in angs:
            print("Doing ({0}, {1}, {2})".format(xRot, yRot, zRot), file=sys.stderr)
            compare(xRot=xRot, yRot=yRot, zRot=zRot)
            #exit(0)
#compare()
commandFile.close()
