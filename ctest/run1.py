import ctypes as ct
import os
import passing
import dill
import stucts as st


intTy = ct.c_int32
floatTy = ct.c_double
# build json
jsonFile = "test.json"
dataFile = "ugg.dill"

with open(dataFile, "rb") as f:
    data = dill.loads(f.read())
    
datap = list(data)
datap[5] = ct.c_void_p(0)
datapp = tuple(datap)

builders = st.makeAllTypes(intTy, floatTy)
buildAll = builders[-1]
femArgs = buildAll(*datapp)

programNameArg = "justTypes"
nameSpaceArg = "justTypes"
library = ct.CDLL("./" + programNameArg + ".so")
program = passing.Library(library, nameSpace=nameSpaceArg)

femArgsNamed = {"a": [femArgs[0]], "b": [femArgs[1]], "c": [femArgs[2]]}
outputs = ["pos"]
program.go(femArgsNamed, outputs)

