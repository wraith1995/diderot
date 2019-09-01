import signal
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
    def runTime(self):
        f = self.get_func("get_runTime")
        return(f())
    def create_world(self):
        world_func = self.get_func("new_world")
        world = world_func()
        self.world = ct.c_void_p(world)
        self.errorCheck(name="create_world")
        return(world)
    def setVerobse(self, bl):
        verboseFunc = self.get_func(name="set_verbose")
        verboseFunc(self.world, ct.c_bool(bl))
    def setWorkers(self, workers):
        workerFunc = self.get_func(name="set_num_workers")
        workerFunc(self.world, ct.c_uint(workers))
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
    def set_input_by_name(self, inputName, string):
        argsPrime = [self.world, ct.c_char_p(str.encode(string))]
        argFunc = self.get_func("input_set_by_name_" + inputName)
        ret = argFunc(*argsPrime)
        self.errorCheck(name="input_set_by_name_" + inputName, returnCheck=ret)
    def create_strands(self):
        strandFunc = self.get_func("create_strands")
        strandRet = strandFunc(self.world)
        self.errorCheck(name="create_strands", returnCheck=strandRet)

    def runStrands(self, steps=0, time=False):
        if not(time):
            runFunc = self.get_func("run")
            print("run")
            steps = runFunc(self.world, ct.c_uint(steps))
            print("end")
            return(steps)
        else:
            time =-10.0
            timeVal = ct.c_double(time)
            ugg = ct.byref(timeVal)
            hmm = (ct.c_double * 1)()
            runFunc = self.get_func("run_with_time")
            print("run")
            steps = runFunc(self.world, ct.c_uint(steps),
                            hmm)
            print("end")
            return(steps, hmm[0])

    def shutDown(self):
        runFunc = self.get_func("shutdown")
        runFunc(self.world)
        self.world = 0

    def saveOutput(self, fileName, outName, num=1, snapshot=False):
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
        if snapshot:
            outFunc = self.get_func("snapshot_" + outName)
        outRet = outFunc(*outArgs)
        self.errorCheck(name="outputing " + outName, returnCheck=outRet)
        for (idx, nrrd) in enumerate(nrrds):
            fileNameIdx = fileName + "_" + str(idx) + ".nrrd"
            saveArgs = [ct.c_char_p(str.encode(fileNameIdx)), nrrd, 0]
            saveRet = nrrdSave(*saveArgs)
            print("Saving:"+fileNameIdx) #c_char_p
            self.errorCheck(name="Saving to nrrd: {0}".format(fileName),
                            returnCheck=saveRet)

    def setup(self, inputs, namedInputs=[], verbose=False, workers=None):
        self.create_world()
        print("Created World")
        self.init_world()
        print("Running")
        self.setVerobse(verbose)
        if workers is not None:
            self.setWorkers(workers)
        for name in inputs:
            print("Set input:{0}".format(name))
            self.set_input(name, inputs[name])
        print("Set inputs")
        for name in namedInputs:
            print("Setting input: {0}".format(name))
            self.set_input_by_name(name, namedInputs[name])
        print("Set named inputs")
        self.create_strands()
        print("Create Strands")

    def snapshotGo(self, inputs, outputs, namedInputs=[], verbose=False, workers=None, stepsize=1):
        print("setting up:")
        self.setup(inputs, namedInputs=namedInputs, verbose=verbose, workers=workers)
        print("snapshot run")
        itter = 0
        while True:
            test = self.runStrands(steps=stepsize, time=None)
            print("Snapshot ran ({0}):{1}".format(itter, test))
            for output in outputs:
                fileName = str(itter) + output[2]
                self.saveOutput(fileName, output[0], num=output[1], snapshot=True)
            itter+=1
            if test == 0:
                return(0)

    def go(self, inputs, outputs, namedInputs=[], verbose=False, shutdown=True,
           workers=None, steps=0, time=False):
        self.create_world()
        print("Created World")
        self.init_world()
        print("Running")
        self.setVerobse(verbose)
        if workers is not None:
            self.setWorkers(workers)
        for name in inputs:
            print("Set input:{0}".format(name))
            self.set_input(name, inputs[name])
        print("Set inputs")
        for name in namedInputs:
            print("Setting input: {0}".format(name))
            self.set_input_by_name(name, namedInputs[name])
        print("Set named inputs")
        self.create_strands()
        print("Create Strands")
        result = self.runStrands(steps=steps, time=None)
        print("Ran steps: {0}".format(result))
        for output in outputs:
            fileName = output[2]
            self.saveOutput(fileName, output[0], num=output[1])
        if shutdown:
            try:
                self.shutDown()
            except Exception:
                print("Time out in shutdown")
        return(result)
