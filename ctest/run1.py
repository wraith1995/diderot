import ctypes as ct
import os
import passing
import dill
import stucts as st


intTy = ct.c_int32
floatTy = ct.c_float
# build json
jsonFile = "test.json"
dataFile = "ugg.dill"

with open(dataFile, "rb") as f:
    data = dill.loads(f.read())
    

print(data)
datap = list(data)
print(datap[5])
datap[5] = ct.c_void_p(0) #ct.cast(ct.POINTER(datap[5]), ct.c_void_p)
datapp = tuple(datap)
print(datapp)
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

