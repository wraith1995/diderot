import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import ctypes as ct
import os
import passing
import pickle
import sys
import argparse
import stucts as st
import numpy as np
parser = argparse.ArgumentParser(description='Run a diderot program')


# select the types
intTy = ct.c_int32
floatTy = ct.c_double

pklfile = "scarray.pkl"
data = pickle.load(open(pklfile, "rb"))

builder = st.makeAllTypes(intTy, floatTy)
args = builder[-1](*data[0])

#main selection


# python camera control:
camEyePy = [0.0, -3.1, 0.0]
camAtPy = [0.0, 0.0, 0.0]
camUpPy = [0.0, 0.0, 1.0]
lightVsp = camEyePy  
camFovPy = 20
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

featureDict = {"fStrTh" : [floatTy(fStrThPy)]}

for cellChoice in range(119):
    print("Cell {0}".format(cellChoice))
    inputs = {"block": [args[0]],
              "space": [args[1]],
              "U": [args[2]],
              "rayStep": [floatTy(0.001)],
              "cellChoice": [intTy(cellChoice)],
              **initDict, **featureDict}
    outputs = [("rgba", 1, "rgba{0}".format(cellChoice))]
    programNameArg = "evalProg"
    nameSpaceArg = "evalProg"
    library = ct.CDLL("./" + programNameArg + ".so")
    program = passing.Library(library, nameSpace=nameSpaceArg)
    program.go(inputs, outputs, verbose=True)
