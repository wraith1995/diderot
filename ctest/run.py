import ctypes as ct
import firedrake_build as fb
import structs as st

programNameArg = "justTypes"
nameSpaceArg = "justTypes"

library = ct.CDLL("./" + programNameArg + ".so")
