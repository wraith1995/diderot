import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import pickle
import numpy as np
import sympy as sp
from itertools import repeat, product
from jsonbuilder import buildPythonic
from json import dumps
x,y,z = sp.Symbol("x"),sp.Symbol("y"), sp.Symbol("z")
def prod(s):
    z = [q[0] for q in s]
    ret =  np.product(z)
    return(ret)


def mono_index_(var_dim, degree):
    length = var_dim
    total = degree
    if length == 1:
        yield (total,)
        return

    for i in range(total + 1):
        for t in mono_index_(length - 1, total - i):
            yield (i,) + t

def mono_index(var_dim, degree):
    stuff = []
    for deg in range(1,degree+1):
        stuff += list(mono_index_(var_dim, deg))
    non_sorted = stuff
    non_sorted.append(tuple(repeat(0 ,var_dim)))
    result = sorted(non_sorted, key = sum)
    return(result)

def gen_monos(syms, idxes):
    "Creates a monomial basis for polynomials in [[vars]] of degree degree"
    syms_prime = syms[0]
    len_test = len(idxes[0])
    var_dim = len(syms_prime)
    if len_test != var_dim:
        raise Exception("Dimensions of idxes and syms incompatible")
    idx = []
    for opt in idxes:
        temp = product([a**b  for (a,b) in  zip(syms_prime, opt)])
        mono = prod(temp)
        idx.append(mono)
    return(idx)

def mono_coeffs(p, monos):
    "Transform a basis p to a monomial basis"
    coeffs = []
    for mon in monos:
        coeff = p.coeff_monomial(mon)
        coeffs.append(coeff)
    return(coeffs)


mono_idxes = mono_index(3, 6)
mono_len = len(mono_idxes)
monos = gen_monos([[x, y, z]], mono_idxes)

verts = [[-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
         [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5]]
topology = dict()
topology[0] = dict() #ignore
topology[1] = {0: [0, 1], 1: [1, 2], 2: [2, 3], 3: [3, 0],
               4: [4, 5], 5: [5, 6], 6: [6, 7], 7: [7, 0],
               8: [0, 4], 9: [1, 5], 10: [2, 6], 11: [3, 7]}
topology[2] = {0: [0, 1, 2, 3], 1: [4, 5, 6, 7],
               2: [0, 4, 1, 5],
               3: [1, 5, 2, 6],
               4: [2, 6, 3, 7],
               5: [3, 7, 0, 1]}
topology[3] = [] #ignore
# json = buildPythonic(3, [x, y, z], [1, x, y, z],
#                      monos, verts, topology, [0,0,0], insideFile="checkInside.in")
# with open("pagotFile.json", "w+") as f:
#     txt = dumps(json, indent=4)
#     f.write(txt)
#Time to build shit...
#Do a test of evaluation
#Alright.


hexPoints = np.random.uniform(low=-0.5, high=0.5, size=(100,3))

evalResults = []
evalVecResults = []

pklfile = "sc.pkl"
data = pickle.load(open(pklfile, "rb"))

def buildData():
    meshIndexMap = []
    meshCoordMap = []
    dim = 3
    meshMapDim = 4
    numCells = 0
    sIndex = None
    con = np.array([])
    spaceIndexMap = []
    spaceDim = mono_len
    funcCoordMap = []

    for datam in data:
        meshIndexMap.append(list(range(numCells * 4, numCells * 4 + 4)))
        com = datam["com"]
        print(com)
        meshCoordMap.append(list(com)) #1
        meshCoordMap.append(list(np.array([1, 0, 0]))) #x
        meshCoordMap.append(list(np.array([0, 1, 0]))) #y
        meshCoordMap.append(list(np.array([0, 0, 1]))) #z

        spaceIndexMap.append(list(range(numCells * mono_len,
                                        numCells * mono_len + mono_len)))
        polyText = datam["poly"]
        prePoly = eval(polyText)
        poly = sp.poly(prePoly)
        coeffs = mono_coeffs(poly, monos)
        lambd = sp.lambdify([x,y,z], prePoly)
        evalResults.append([lambd(*x) for x in hexPoints])
        lambdx = sp.lambdify([x,y,z], sp.diff(prePoly, x))
        lambdy = sp.lambdify([x,y,z], sp.diff(prePoly, y))
        lambdz = sp.lambdify([x,y,z], sp.diff(prePoly, z))
        evalVecResults.append([[lambdx(*x), lambdy(*x), lambdz(*x)] for x in hexPoints])
        funcCoordMap += coeffs
        numCells += 1

    meshIndexMap = np.array(meshIndexMap, dtype="int32")
    print("index:",meshIndexMap)
    meshCoordMap = np.array(meshCoordMap, dtype="float64")
    print("Coords:", meshCoordMap)
    spaceIndexMap = np.array(spaceIndexMap, dtype="int32")
    funcCoordMap = np.array(funcCoordMap, dtype="float64")

    eval1 = np.array(evalResults, dtype="float64")
    eval2 = np.array(evalVecResults, dtype="float64")
    print(eval1.shape, eval2.shape)
    ret = (meshIndexMap, meshCoordMap,
           dim, meshMapDim, numCells,
           sIndex, con,
           spaceIndexMap,
           spaceDim,
           funcCoordMap)
    print(eval2)
    return((ret, hexPoints, eval1, eval2))

result = buildData()
pklfilePrime = "scarray.pkl"
txt = pickle.dumps(result)
with open(pklfilePrime, "wb+") as f:
    f.write(txt)
