import sys
import csv
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
sys.path.append('/home/teocollin/gitcode/curvedMesh')
sys.path.append('s1')
sys.path.append('s2')
sys.path.append('s3')
from firedrake import *
import ctypes as ct
import firedrake_build as fb
import os
import passing
import dill
import sys
import numpy as np
import argparse
import nrrd_utils as nu
import nrrd
import ctypes
import functools
from render import render, compareError, compareErrorSingle, compareToBackSolveSingle
from load import buildFiredrakeTetMesh
from runGuided import runGuided
from runGuidedErr import runGuidedErr
from runTruth import runTruth
#get the damn mesh
meshDatas = buildFiredrakeTetMesh("meshfiles/gsc5.msh", 3, linearGmshBackup="lin_meshfiles/gsc5.msh")
linMesh = meshDatas[-1]
curvedMesh = meshDatas[0]
mesh = curvedMesh

# select the types
intTy = ct.c_int32
floatTy = ct.c_double
dim = 3

#build the function
beta = 0.1
M = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 0]], dtype="float64")
center = np.array([0.0, 0.0, beta], dtype="float64")
Mtemp = np.asfarray(M, dtype="float64")
Mptr = M.ctypes.data_as(ctypes.POINTER(floatTy))
centerPtr = np.asfarray(center, dtype="float64").ctypes.data_as(ctypes.POINTER(floatTy))
x,y,z = SpatialCoordinate(mesh)
Mx = M.dot([x, y, z]) + center
space = VectorFunctionSpace(mesh, "Lagrange", 2, dim=3)
f = interpolate(as_vector(Mx), space)



(_, femArgs) = fb.passAll(f, intTy, floatTy, geometric=True)
#points = [[-0.8, 0.1, 0.2]]
#points = [[0.85, 0.0, 0.2]]
points = [[-0.875, 0.1, 0.2]]
pointsNrrd = "points.nrrd"
dataPoints = np.array(points, dtype="float64")
kindString = "{0}-vector".format(dim)
nu.writeSequence(dataPoints, pointsNrrd, dataKind=kindString)

# dumb outright/fix validity issues...
#loop over, calc errors and vtks
#calc error via time step/time slowdown
#invalid position creation...
#runtime calc
#PVO IS MORE INTERESTING BECAUSE ISO
#other mesh.
#perfoamnce.
# runTruth("umm", pointsNrrd, femArgs, stepSize=0.5, stepMax=500)
# render("umm", "umm", dim=dim)
# runTruth("umm1", pointsNrrd, femArgs, stepSize=0.05, stepMax=500)
# render("umm1", "umm1", dim=dim)
# # exit(0)

#add another mesh
#add timing calculate
#add error calculate
#yay!
#write this section
#writhe methods + conclusion

stepSizes = [0.5, 0.05, 0.005, 0.0005]
stepsMaxes = [5000]
#common to newtons
timeStepses = [32]
secondes = [0, 1]
timeEpses = [0.0000001]
#specific to error:
errorMaxes = [0.1, 0.01, 0.001, 0.0001, 0.0000001]
prefixNrrd = "nrrdResults/"
prefixVtk = "vtkResults/"
errFileDict = dict()
errorDict = dict()
timingDict = dict()
allOtherErrors = []
for stepSize in stepSizes:
    for stepsMax in stepsMaxes:
        errorComps = []
        first = "run_{0}_{1}".format(str(stepSize), stepsMax)
        firstNrrd = prefixNrrd + first + "_truth"
        timeTrue = runTruth(firstNrrd, pointsNrrd, femArgs, stepSize=stepSize,
                            stepMax=stepsMax)
        render(firstNrrd, prefixVtk + first + "_truth")
        timingDict[(stepSize, stepsMax)] = (timeTrue[1], dict())
        errorDict[(stepSize, stepsMax)] = dict()
        for timeStep in timeStepses:
            for timeEps in timeEpses:
                for second in secondes:
                    guide = first + "_{0}_{1}_{2}_guide".format(str(timeStep), str(timeEps), str(second))
                    guideNrrd = prefixNrrd + guide
                    timeGuide = runGuided(guideNrrd, pointsNrrd, femArgs,
                                          stepSize=stepSize, stepMax=stepsMax,
                                          timeSteps=timeStep, timeEps=timeEps,
                                          second=second)
                    timingDict[(stepSize, stepsMax)][1][(timeStep, timeEps, second)] = (timeGuide[1], dict())
                    errorDict[(stepSize, stepsMax)][(timeStep, timeEps, second)] = (compareErrorSingle(firstNrrd, guideNrrd), dict())
                    render(guideNrrd, prefixVtk + guide)
                    errorComps.append(guideNrrd)
                    for maxError in errorMaxes:
                        guideErr = guide + "_{0}_checked".format(str(maxError))
                        guideErrNrrd = prefixNrrd + guideErr
                        timeErrGuide = runGuidedErr(guideErrNrrd, pointsNrrd, femArgs,
                                                    stepSize=stepSize, stepMax=stepsMax,
                                                    timeSteps=timeStep, timeEps=timeEps,
                                                    second=second, errorMax=maxError)
                        timingDict[(stepSize, stepsMax)][1][(timeStep, timeEps, second)][1][maxError] = timeErrGuide[1]
                        errorDict[(stepSize, stepsMax)][(timeStep, timeEps, second)][1][maxError] = compareErrorSingle(firstNrrd, guideErrNrrd)
                        render(guideErrNrrd, prefixVtk + guideErr)
                        errorComps.append(guideErrNrrd)
        print(errorComps)
        errFileDict[(stepSize, stepsMax)] = errorComps

def makeTiminingList(timings):
    titles = ["stepSize", "stepsMax", "guideTimeStep", "guideTimeEps", "dervApproxType", "maxErrorParam", "time"]
    data = []
    for entry in timings.keys():
        firstTwo = entry
        print(firstTwo)
        data.append([firstTwo[0], firstTwo[1], float("nan"), float("nan"), float("nan"), float("nan"), timings[entry][0]])
        tempTimings = timings[entry][1]
        for entryNext in tempTimings.keys():
            data.append([firstTwo[0], firstTwo[1], entryNext[0], entryNext[1], entryNext[2], float("inf"), tempTimings[entryNext][0]])
            tempTimingsPrime = tempTimings[entryNext][1]
            for entryFin in tempTimingsPrime.keys():
                data.append([firstTwo[0], firstTwo[1], entryNext[0], entryNext[1], entryNext[2], entryFin, tempTimingsPrime[entryFin]])
    return(titles, data)

def makeErrorsList(timings):
    titles = ["stepSize", "stepsMax", "guideTimeStep", "guideTimeEps", "dervApproxType", "maxErrorParam", "stepDiffs", "maxError", "averageError"]
    data = []
    for entry in timings.keys():
        firstTwo = entry
        print(firstTwo)
        tempTimings = timings[entry]
        for entryNext in tempTimings.keys():
            errors = [tempTimings[entryNext][0][0], tempTimings[entryNext][0][1], tempTimings[entryNext][0][3]]
            data.append([firstTwo[0], firstTwo[1], entryNext[0], entryNext[1], entryNext[2], float("inf"), *errors])
            tempTimingsPrime = tempTimings[entryNext][1]
            for entryFin in tempTimingsPrime.keys():
                otherErrors = [tempTimingsPrime[entryFin][0], tempTimingsPrime[entryFin][1], tempTimingsPrime[entryFin][3]]
                data.append([firstTwo[0], firstTwo[1], entryNext[0], entryNext[1], entryNext[2], entryFin, *otherErrors])
    return(titles, data)

def makeTimingsCsv(fileName, timings):
    (titles, data) = makeTiminingList(timings)
    with open(fileName, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(titles)
        for datam in data:
            writer.writerow(datam)

def makeErrorsCsv(fileName, timings):
    (titles, data) = makeErrorsList(timings)
    with open(fileName, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(titles)
        for datam in data:
            writer.writerow(datam)

def makeCombined(fileName, timings, errors):
    (titles, dataTimings) = makeTiminingList(timings)
    (titlesp, dataErrors) = makeErrorsList(errors)
    combinedTitles = titles + ["stepDiffs", "maxError", "averageError"]
    with open(fileName, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(combinedTitles)
        for timing in dataTimings:
            if not(timing[2] == timing[2]):
                writer.writerow(timing + [0.0, 0.0, 0.0])
                continue
            if len(dataErrors) > 0:
                topDataErr = dataErrors[0]
                dataErrors.remove(topDataErr)
                writer.writerow(timing + topDataErr[-3:])
                continue
                
    

            
makeTimingsCsv("timings.cvs", timingDict)
makeErrorsCsv("errors.cvs", errorDict)
makeCombined("combined.cvs", timingDict, errorDict)

#paraview default settings -one done - but others ugg
#backsolve error... for every nrrd, for each time step, compute correct value -> yay..***
#make another mesh
