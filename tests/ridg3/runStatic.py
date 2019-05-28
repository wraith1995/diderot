import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import ctypes as ct
import os
import passing
import dill
import stucts as st

import argparse
parser = argparse.ArgumentParser(description='Run a diderot given static data')
args = parser.parse_args

intTy = ct.c_int32
floatTy = ct.c_double

#load data in.
dataFile = "evalProg.dill"
with open(dataFile, "rb") as f:
    data = dill.loads(f.read())
  
datap = list(data)
datap[5] = ct.c_void_p(0)
datapp = tuple(datap)

builders = st.makeAllTypes(intTy, floatTy)
buildAll = builders[-1]
femArgs = buildAll(*datapp)


\programNameArg = "evalProg"
nameSpaceArg = "evalProg"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)
inputs = []
outputs = []
program.go(femArgsNamed, outputs)
