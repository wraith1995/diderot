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
<<<<<<< HEAD
from mpmath import mpf, mp
=======
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
xVar = loadSc.x
yVar = loadSc.y
zVar = loadSc.z

# select the types
intTy = ct.c_int32
floatTy = ct.c_double


def samplePolyToNrrd(poly, dataNrrd, samples=300, orig=[-0.6, -0.6, -0.6], end=[0.6, 0.6, 0.6]):
    f = sp.lambdify([xVar, yVar, zVar], poly)
    basis1 = ((end[0] - orig[0]) / (samples - 1), 0.0, 0.0)
    basis2 = (0.0, (end[1] - orig[1]) / (samples - 1), 0.0)
    basis3 = (0.0, 0.0, (end[2] - orig[2]) / (samples - 1))
<<<<<<< HEAD
    data = np.zeros((samples, samples, samples), dtype="float64")
=======
    data = np.zeros((samples, samples, samples), dtype="float32")
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
    for i in range(samples):
        x = orig[0] + basis1[0] * i
        for j in range(samples):
            y = orig[1] + basis2[1] * j
            for k in range(samples):
                z = orig[2] + basis3[2] * k
<<<<<<< HEAD
                data[i][j][k] = float(f(x, y, z))
    options = {"kinds": ["space", "space", "space"],
               "encoding": "raw",
               "space dimension": 3,
               "space directions": [basis1, basis2, basis3],
               "space origin": orig,
               "centerings": ["node", "node", "node"]}
=======
                data[i][j][k] = f(x, y, z)
    options = {"kinds": ["space", "space", "space"],
               "space dimension": 3,
               "space directions": [basis1, basis2, basis3],
               "space origin": orig,
               "centerings": ["cell", "cell", "cell"]}
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
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
    H = np.array([[float(sp.diff(g, xVar)), float(sp.diff(g, yVar)), float(sp.diff(g, zVar))] for g in grad])
<<<<<<< HEAD
    print("G:", grad)
    print("H:", H)
=======
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
    eigenValues, eigenVectors = np.linalg.eig(H)
    idx = eigenValues.argsort()
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
<<<<<<< HEAD
    print(eigenVectors)
=======
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
    print("Looking for zeros given by:")
    print(np.dot(eigenVectors[0], grad))
    print(np.dot(eigenVectors[1], grad))
    print("eigenvalues are:", eigenValues)


def makeData(xRot=0.0, yRot=0.0, zRot=0.0, a=1.0, b=2.0, c=3.0, datafile="data.nrrd"):
    a = 1.0
    b = 2.0
    c = 3.0
    Ax = np.array([[1.0, 0, 0.0], [0, math.cos(xRot), -math.sin(xRot)], [0, math.sin(xRot), math.cos(xRot)]])
    Ay = np.array([[math.cos(yRot), 0, math.sin(yRot)], [0, 1, 0], [-math.sin(yRot), 0, math.cos(yRot)]])
    Az = np.array([[math.cos(zRot), -math.sin(zRot), 0], [math.sin(zRot), math.cos(zRot), 0], [0, 0, 1]])
    A = Ax.dot(Ay).dot(Az)
    y = A.dot([loadSc.x, loadSc.y, loadSc.z])
    p = -0.5 * (a*loadSc.x**2 + b*loadSc.y**2 + c*loadSc.z**2)
    pPrime = p.xreplace({loadSc.x : y[0], loadSc.y: y[1], loadSc.z: y[2]})
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
<<<<<<< HEAD
    return(reset, A)
=======
    return(reset)
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315

#main selection
builder = st.makeAllTypes(intTy, floatTy)
programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
#write nrrds
commandFile = open("commands.txt", "w+")
def compare(xRot=0.0, yRot=0.0, zRot=0.0):
<<<<<<< HEAD
    (reset, rot) = makeData(xRot=xRot, yRot=yRot, zRot=zRot)
    args = builder[-1](*reset)
=======
    reset = makeData(xRot=xRot, yRot=yRot, zRot=zRot)
    args = builder[-1](*reset)
    
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
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
<<<<<<< HEAD

    #make rot array
    rotArg = (floatTy * 9)(*(list(rot.flatten())))
=======
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
    
    initDict = {"camEye": [(floatTy*3)(*camEyePy)],
                "camAt": [(floatTy*3)(*camAtPy)],
                "camUp": [(floatTy*3)(*camUpPy)],
                "camFOV": [floatTy(camFovPy)],
                "camNear": [floatTy(camNearPy)],
                "camFar": [floatTy(camFarPy)],
                "iresU": [intTy(iresUPy)],
<<<<<<< HEAD
                "iresV": [intTy(iresVPy)],
                "rot": [rotArg]}
=======
                "iresV": [intTy(iresVPy)]}
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
    
    
    #python feature control
    fStrThPy = 0.0
    thickPy = 0.01
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
<<<<<<< HEAD
    print("ROt is:", rot)
=======
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
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
        
<<<<<<< HEAD
=======
    
>>>>>>> 3b62777ac36b781b70d4cdf090b32d975878b315
angs = [x * 10 for x in range(0,18)]
for xRot in angs:
    for yRot in angs:
        for zRot in angs:
            compare(xRot=xRot, yRot=yRot, zRot=zRot)
    
#compare()
commandFile.close()
