import sys
import pickle
import numpy as np
import sympy as sp
import meshio
from enum import Enum
from math import log, ceil
from mpmath import *
#mp.dps = 100
#print(mp)

import ctypes as ct
import time
class Axes(Enum):
    X = 1
    Y = 2
    Z = 3


affineLib = ct.CDLL("/home/teocollin/gitcode/intervals/affinePoly/ap.so")
testPolyBool = affineLib.__getattr__("testPolyBool")
testPolyBoolExp = affineLib.__getattr__("testPvoExplt")
testPolyEval = affineLib.__getattr__("testPolyEval")

def testPolyEvalFunction(mins, coords):
    xmin = ct.c_double(mins[0])
    ymin = ct.c_double(mins[1])
    zmin = ct.c_double(mins[2])
    array = np.asfarray(coords, dtype="float64")
    arrayptr = array.ctypes.data_as(ct.POINTER(ct.c_double))
    return((testPolyEval(xmin, ymin, zmin, arrayptr)))

def polyTestFunction(mins, maxes, coords):
    xmin = ct.c_double(mins[0])
    ymin = ct.c_double(mins[1])
    zmin = ct.c_double(mins[2])
    
    xmax = ct.c_double(maxes[0])
    ymax = ct.c_double(maxes[1])
    zmax = ct.c_double(maxes[2])

    array = np.asfarray(coords, dtype="float64")
    arrayptr = array.ctypes.data_as(ct.POINTER(ct.c_double))
    return(bool(testPolyBool(xmin, ymin, zmin, xmax, ymax, zmax, arrayptr)))


def polyTestFunctionExp(mins, maxes, grad, hess):
    xmin = ct.c_double(mins[0])
    ymin = ct.c_double(mins[1])
    zmin = ct.c_double(mins[2])
    
    xmax = ct.c_double(maxes[0])
    ymax = ct.c_double(maxes[1])
    zmax = ct.c_double(maxes[2])

    grads = [np.asfarray(g, dtype="float64") for g in grad]
    gradPtrs = [array.ctypes.data_as(ct.POINTER(ct.c_double)) for array in grads]
    hesses = [np.asfarray(g, dtype="float64") for g in hess]

    hessPtrs = [array.ctypes.data_as(ct.POINTER(ct.c_double)) for array in hesses]

    return(bool(testPolyBoolExp(xmin, xmax, ymin, ymax, zmin, zmax,
                             gradPtrs[0], gradPtrs[1], gradPtrs[2],
                             hessPtrs[0], hessPtrs[1], hessPtrs[2],
                             hessPtrs[4], hessPtrs[5], hessPtrs[8])))

# print(polyTestFunction((-0.1, -0.1, -0.1), (0.1, 0.1, 0.1), 2*np.ones(680)))
# exit(0)

x,y,z = sp.Symbol("x"),sp.Symbol("y"), sp.Symbol("z")

"""
Here is how we understand rectangles:

0 - 1 - 2 - 3

4 - 5 - 6 - 7


"""


def topology(rect):
    #bottom segs:
    seg1 = [rect[0], rect[1]]
    seg2 = [rect[1], rect[2]]
    seg3 = [rect[2], rect[3]]
    seg4 = [rect[3], rect[0]]
    
    # legs:
    seg5 = [rect[0], rect[4]]
    seg6 = [rect[1], rect[5]]
    seg7 = [rect[2], rect[6]]
    seg8 = [rect[3], rect[7]]


    #roof:

    seg9 = [rect[4], rect[5]]
    seg10 = [rect[5], rect[6]]
    seg11 = [rect[6], rect[7]]
    seg12 = [rect[7], rect[4]]

    #bottom square
    bsquare = [rect[0], rect[1], rect[2], rect[3]]
    #top square
    tsquare = [rect[4], rect[5], rect[6], rect[7]]

    forwardSquare = [rect[0], rect[1], rect[5], rect[4]]
    rightSquare = [rect[1], rect[2], rect[6], rect[5]]
    backSquare = [rect[2], rect[3], rect[7], rect[6]]
    leftSquare = [rect[3], rect[0], rect[4], rect[7]]

    top = {"bot": bsquare, "top": tsquare,
           "for": forwardSquare, "back": backSquare,
           "right": rightSquare, "left": leftSquare}

    return(top)

def quadDivide(square):
    """
    We assume 4  3 
              1  2
    and produce as such.
    """
    sq1 = square[0]
    sq2 = square[1]
    sq3 = square[2]
    sq4 = square[3]

    mid1 = (sq1 + sq2)/2.0 #np arrays
    mid2 = (sq2 + sq3)/2.0
    mid3 = (sq3 + sq4)/2.0
    mid4 = (sq4 + sq1)/2.0
    """
       mid3
    mid4   mid2
       mid1
    """

    com = (sq1 + sq2 + sq3 + sq4)/4.0

    sub1 = [sq1, mid1, com, mid4]
    sub2 = [mid1, sq2, mid2, com]
    sub3 = [com, mid2, sq3, mid3]
    sub4 = [mid4, com, mid3, sq4]
    """
    sub4  sub3
    sub1  sub2
    """
    return([sub1, sub2, sub3, sub4])



##This is a dumb way to do it - the real way to do it is to represent cubes by points with radius in max norm.
##Supply con

def pointToRect(p, rad):
    e1 = np.array([rad, 0, 0])
    e2 = np.array([0, rad, 0])
    e3 = np.array([0, 0, rad])
    p1 = p - e1 - e2 - e3
    p2 = p - e1 + e2 - e3
    p3 = p - e1 + e2 + e3
    p4 = p - e1 - e2 + e3
    p5 = p + e1 - e2 - e3    
    p6 = p + e1 + e2 - e3
    p7 = p + e1 + e2 + e3
    p8 = p + e1 - e2 + e3
    return([p1, p2, p3, p4, p5, p6, p7, p8])

def pointToQuadCenters(p, rad):
    e1 = np.array([rad, 0, 0])
    e2 = np.array([0, rad, 0])
    e3 = np.array([0, 0, rad])

    p1 = p - e1
    p2 = p + e1
    p3 = p - e2
    p4 = p + e2
    p5 = p - e3
    p6 = p + e3
    return([(p1, Axes.X), (p2, Axes.X), (p3, Axes.Y), (p4, Axes.Y), (p5, Axes.Z), (p6, Axes.Z)])

def pointTo3Interval(p, rad):
    e1 = np.array([rad, 0, 0])
    e2 = np.array([0, rad, 0])
    e3 = np.array([0, 0, rad])

    pmin = p - e1 - e2 - e3
    pmax = p + e1 + e2 + e3
    #intervals = list(zip(pmin, pmax))
    return((pmin, pmax))

def quadBasis(axes, rad):
    e1 = np.array([rad, 0, 0])
    e2 = np.array([0, rad, 0])
    e3 = np.array([0, 0, rad])
    w1 = None
    w2 = None
    if axes == Axes.X:
        w1 = e2
        w2 = e3
    elif axes == Axes.Y:
        w1 = e1
        w2 = e3
    elif axes == Axes.Z:
        w1 = e1
        w2 = e2
    else:
        raise Exception("Invalid axes failure")
    return((w1, w2))
def pointTo2Interval(p, axes, rad):
    (e1, e2) = quadBasis(axes, rad)
    pmin = p - e1 - e2
    pmax = p + e1 + e2
    intervals = list(zip(pmin, pmax))
    return(intervals)


def subdivide4(p, axes, rad):
    (w1, w2) = quadBasis(axes, rad)
    x1 = w1 / 2.0
    x2 = w2 / 2.0

    e1 = x1 + x2
    e2 = x1 - x2

    p1 = p + e1
    p2 = p - e1
    p3 = p + e2
    p4 = p - e2
    return([(point, axes) for point in [p1, p2, p3, p4]])

def subDivide8(p, rad):
    e1 = np.array([rad, 0, 0])/2.0
    e2 = np.array([0, rad, 0])/2.0
    e3 = np.array([0, 0, rad])/2.0

    w1 = e1 + e2
    w2 = e1 - e2

    q1 = p + w1
    q2 = p + w2
    q3 = p - w1
    q4 = p - w2

    p1 = q1 + e3
    p2 = q1 - e3
    p3 = q2 + e3
    p4 = q2 - e3
    p5 = q3 + e3
    p6 = q3 - e3
    p7 = q4 + e3
    p8 = q4 - e3
    return([p1, p2, p3, p4, p5, p6, p7, p8])



def writeToVtk(fileName, pointQuads):
    points = [tuple(point) for quad in pointQuads for point in quad]
    uniqPoints = list(set(points))
    quadInts = []
    for quad in pointQuads:
        tempIdxes = []
        for p in quad:
            idx = uniqPoints.index(tuple(p))
            tempIdxes.append(idx)
        quadInts.append(tempIdxes)
    cells = {
        "hexahedron" : np.array(quadInts, dtype="int32")
    }
    writePoints = np.array(uniqPoints)
    meshio.write_points_cells(fileName, writePoints, cells)

# test
# z = np.array([0.0, 0.0, 0.0])
# r = pointToRect(z, 0.5)
# newQuads = subDivide8(z, 0.5)
# newPoints = [pointToRect(hexa, 0.25) for hexa in newQuads]
# newQuadsAgain = [zzz for zz in newQuads for zzz in subDivide8(zz, 0.25)]
# newPointsAgain = [pointToRect(hexa, 0.25/2.0) for hexa in newQuadsAgain]

# writeToVtk("ugg.vtk", newPointsAgain) #, file_format="vtk-ascii")
#subdivide, interval all, test all, return allright

def buildWgH(i):
    pklfile = "gH.pkl"
    data = pickle.load(open(pklfile, "rb"))
    datam = data[i]
    return(datam["grad"], datam["hess"])
def buildW(i):
    pklfile = "sc.pkl"
    data = pickle.load(open(pklfile, "rb"))
    datam = data[i]
    polyText = datam["poly"]
    prePoly = eval(polyText)
    poly = prePoly
    #print(sp.Poly(poly).expand())
    gradPoly = [sp.diff(poly, x), sp.diff(poly, y), sp.diff(poly, z)]
    gradPolyArray = np.array(gradPoly)
    hessianPoly = [[sp.diff(g, x), sp.diff(g, y), sp.diff(g, z)]
                   for g in gradPoly]
    doted = np.array(hessianPoly).dot(gradPolyArray)
    wfield = np.cross(gradPolyArray, doted)
    lf = sp.lambdify([x, y, z], wfield, "sympy")
    #print(wfield[0].expand())
    #print(wfield[1].expand())
    #print(wfield[2].expand())
    pklfilePvo = "pvo.pkl"
    dataPvo = pickle.load(open(pklfilePvo, "rb"))
    datamPvo = dataPvo[i]
    #print(datamPvo)
    #print(len(datamPvo))
    
    
    return(wfield, lf, datamPvo)





#(p, f) = buildW(10)
#load in actual data
#call the damn C function
#yay!
def subdivideOct(i, epsilon):
    (p, f, pvoCoeffs) = buildW(i)
    (grad, hess) = buildWgH(i)
    Oi = ceil(log(1.0/epsilon, 2))
    center = np.array([0.0, 0.0, 0.0])
    currentRadius = 0.5
    hexes = [center]
    itter = 0
    # print(pvoCoeffs)
    # print(p[0].expand())
    # print(p[1].expand())
    # print(p[2].expand())
    internalStart = time.time()
    while itter < Oi:
        #print(itter, "<", Oi)
        newHexes = []
        for hexa in hexes:
            (mins, maxes) = pointTo3Interval(hexa, currentRadius)
            intervals = [mpi((a, b)) for (a, b) in zip(mins, maxes)]
            #print(mins, maxes)
            #test0 = polyTestFunction(mins, maxes, pvoCoeffs[0])
            #test1 = polyTestFunction(mins, maxes, pvoCoeffs[1])
            #test2 = polyTestFunction(mins, maxes, pvoCoeffs[2])o
            if polyTestFunctionExp(mins, maxes, grad, hess):
                #print("Eliniated something")
                continue
            else:
                newHexes.extend(subDivide8(hexa, currentRadius * 0.5))
        hexes = newHexes
        currentRadius = currentRadius * 0.5
        itter += 1
    internalEnd = time.time()
    print("Internal took:", (internalEnd - internalStart), " seconds")
    return(hexes)
# when match depth reached, quadize
# repeat.

#newton
for i in range(119):
    start = time.time()
    l = len(subdivideOct(i, 0.1))
    print("Cell:", i)
    print("Remaining:", l)
    end = time.time()
    print("Took:", (end - start), " seconds")
