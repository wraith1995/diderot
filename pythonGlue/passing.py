import numpy as np
import numpy.ctypeslib as npct
import ctypes as c
from ctypes import POINTER, c_int, c_double, c_void_p, c_float, c_int32, c_uint
import ctypes as ct

# These sorts of things should be arguments to the library, which can also supply help converting arugments to the proper form
# automation of mesh and other types of inputs will require more trickery - maybe modifying generated header file to have those types.



class Library:
    def __init__(self, lib, nameSpace="diderot"):
        self.lib = lib
        self.nameSpace = nameSpace
        self.world = 0

    def errorCheck(self, name="", returnCheck=0):
        any_errors_func = self.lib.__getattr__(self.nameSpace + "_any_errors")
        get_errors_func = self.lib.__getattr__(self.nameSpace + "_get_errors")
        if (self.world == 0):
            raise Exception("World Pointer is void")
        if returnCheck != 0:
            print("Warning, it is possible an error occured:{0}".format(name))
            if (any_errors_func(self.world)):
                errors = get_errors_func(self.world)
                errorBtypes = ct.c_char_p(errors).value
                raise Exception("Error occured: {0} {1}".format(name, errorBtypes))

    def get_func(self, name):
        funcName = self.nameSpace + "_" + name
        return(self.lib.__getattr__(funcName))

    def create_world(self):
        world_func = self.get_func("new_world")
        world = world_func()
        self.world = ct.c_void_p(world)
        self.errorCheck(name="create_world")
        return(world)

    def init_world(self):
        initFunc = self.get_func("init_world")
        world = self.world
        ret = initFunc(self.world)
        self.errorCheck(name="create_world", returnCheck=ret)

    def set_input(self, inputName, args):
        argsPrime = [self.world] + args
        argFunc = self.get_func("input_set_" + inputName)
        ret = argFunc(*argsPrime)
        self.errorCheck(name="input_set_" + inputName, returnCheck=ret)

    def create_strands(self):
        strandFunc = self.get_func("create_strands")
        strandRet = strandFunc(self.world)
        self.errorCheck(name="create_strands", returnCheck=strandRet)

    def runStrands(self):
        runFunc = self.get_func("run")
        print("run")
        steps = runFunc(self.world, ct.c_uint(0))
        print("end")
        return(steps)

    def shutDown(self):
        runFunc = self.get_func("shutdown")
        runFunc(self.world)
        self.world = 0

    def saveOutput(self, fileName, outName, num=1):
        nrrdNew = self.lib.__getattr__("nrrdNew")
        nrrdSave = self.lib.__getattr__("nrrdSave")
        nrrdNuke = self.lib.__getattr__("nrrdNuke") # use this.
        nrrdNewLambda = lambda x: nrrdNew()
        nrrds = [nrrdNewLambda(x) for x in range(num)]
        nrrdTests = [x==0 for x in nrrds]
        if (any(nrrdTests)):
            raise Exception("Unable to allocate Nrrd(s) for ouput {0}".format(outName))
        outFunc = self.get_func("output_get_" + outName)
        outArgs = [self.world] + nrrds
        outRet = outFunc(*outArgs)
        self.errorCheck(name="outputing " + outName, returnCheck=outRet)
        for (idx, nrrd) in enumerate(nrrds):
            fileNameIdx = fileName + "_" + str(idx) + ".nrrd"
            saveArgs = [ct.c_char_p(str.encode(fileNameIdx)), nrrd, 0]
            saveRet = nrrdSave(*saveArgs)
            print("Saving:"+fileNameIdx) #c_char_p
            self.errorCheck(name="Saving to nrrd: {0}".format(fileName),
                            returnCheck=saveRet)

    def go(self, inputs, outputs):
        self.create_world()
        self.init_world()
        print("Running")
        for name in inputs:
            self.set_input(name, inputs[name])
        print("Set inputs")
        self.create_strands()
        print("Create Strands")
        steps = self.runStrands()
        print("Ran steps: {0}".format(steps))
        for output in outputs:
            fileName = output[2]
            self.saveOutput(fileName, output[0], num=output[1])
        self.shutDown()



# meshType = makeMeshType()
# mesh =  meshType()
# spaceType = makeSpaceType(meshType)
# space = spaceType()
# funcType = makeFuncType(spaceType)
# func = funcType()


# mesh.dim = 2
# mesh.mapDim = 4
# mesh.numCells = 2
# c_int_p = POINTER(ctylesInt)
# indexMap = np.array([[0, 1, 2, 3], [4, 5, 6, 7]], dtype=intString)
# indexMap.flatten()
# coordMap = np.array([[1.0, 2.0], [3.0, 5.0], [7.0, 11.0], [13.0, 17.0], [19.0, 23.0], [27.0, 29.0],[31.0, 37.0], [41.0, 43.0]], dtype=floatString)
# coordMapPrime = np.array([2.0,3.0,5.0,7.0,11.0,13.0,17.0,19.0,23.0,29.0,31.0,37.0,41.0, 43.0, 47.0, 53.0, 59.0, 61.0, 68.0, 71.0, 73.0, 79.0, 83.0, 89.0, 97.0, 101.0, 103.0, 107.0, 109.0, 113.0, 127.0, 131.0], dtype=floatString)
# print(coordMapPrime.reshape(8,2,2))
# coordMap.flatten()

# mesh.indexMap = indexMap.ctypes.data_as(c_int_p)
# mesh.coordMap = coordMap.ctypes.data_as(POINTER(ctylesFloat))

# space.indexMap = indexMap.ctypes.data_as(c_int_p)
# space.mesh = mesh
# func.coordMap = coordMapPrime.ctypes.data_as(POINTER(ctylesFloat))
# func.space = space

# meshPtr = ct.cast(ct.pointer(mesh), ct.c_void_p)
# spacePtr = ct.cast(ct.pointer(space), ct.c_void_p)
# funcPtr = ct.cast(ct.pointer(func), ct.c_void_p)

# inputs = {"a" : [meshPtr], "b" : [spacePtr], "c" : [funcPtr]}
# outputs = ["z"]

# doingSometing = Library(library, nameSpace=nameSpaceArg)
# doingSometing.go(inputs, outputs)
