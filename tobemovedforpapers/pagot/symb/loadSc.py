import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import pickle
import numpy as np
import sympy as sp
from itertools import repeat, product
from jsonbuilder import buildPythonic
from json import dumps
from sympy.printing.ccode import ccode
from sympy.codegen.rewriting import create_expand_pow_optimization
x,y,z = sp.Symbol("x"),sp.Symbol("y"), sp.Symbol("z")
def prod(s):
    z = [q[0] for q in s]
    ret =  np.product(z)
    return(ret)

def prodPrime(s):
    z = [q[0] for q in s]
    ret =  sp.Mul(z, evaluate=False)
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

def npow(a, b):
    print("a,b:",a,b)
    if b == 0:
        return(1)
    fuck = [a for j in range(b)]
    temp = sp.Mul(*fuck, evaluate=False)
    return(temp)
def gen_monos_power(syms, idxes):
    "Creates a monomial basis for polynomials in [[vars]] of degree degree"
    syms_prime = syms[0]
    len_test = len(idxes[0])
    var_dim = len(syms_prime)
    if len_test != var_dim:
        raise Exception("Dimensions of idxes and syms incompatible")
    idx = []
    for opt in idxes:
        termps = [npow(a, b) for (a, b) in zip(syms_prime, opt)]
        temp = sp.Mul(*termps, evaluate=False)
        print(opt, " to ", temp)
        idx.append(temp)
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
monos =  gen_monos_power([[x, y, z]], mono_idxes)

mono_idxes14 = mono_index(3, 14)
mono_len14 = len(mono_idxes14)
monos14 = gen_monos_power([[x, y, z]], mono_idxes14)

def pow_to_mul(expr):
    """
    Convert integer powers in an expression to Muls, like a**2 => a*a.
    """
    pows = list(expr.atoms(sp.Pow))
    syms = list(expr.atoms(sp.Symbol))
    if any(not e.is_Integer for b, e in (i.as_base_exp() for i in pows)):

        raise ValueError("A power contains a non-integer exponent")
    repl = list(zip(pows, (sp.Mul(*[b]*e,evaluate=False) for b,e in (i.as_base_exp() for i in pows))))
    return(sp.Mul(*[t[1] for t in repl], evaluate=False))

# def multiTermToPoly(expr):
#     args = list(expr.args)
#     multed = []
#     for exp in args:
#         if len(exp) == 1:
#             multed.append()
#         if type(exp[1]) == int:
#             multed.append(sp.Pow(exp[0], exp[1]))
#             if len(exp) != 2:
#                 raise Exception("Big ops")
#         else:
#             multed.append(exp)
#     multedPrime = [pow_to_mul(exp) for exp in multed]
#     fin = sp.Mul(*multed, evaluate=False)
#     print("From {0} to {1} with {2}".format(expr, fin, args))
#     return(fin)

def multiTermToPoly(expr):
    exprp = expr.replace(lambda x: x.is_Pow and x.exp > 0, lambda x: sp.Mul(*[x.base]*x.exp, evaluate=False))
    print(expr, " to ", exprp)
    return(exprp)
def fMul(x): return(sp.Mul(*x, evaluate=False))

def findLowestPower(var, npow, cseDict):
    for i in reversed(range(1, npow+1)):
        
    

def writeCseEvaluation(monos, type="AAF"):
    extraVarCounter = 0
    cseDict = dict()
    
    for (idx, mono) in enumerate(monos):
        subMuls = mono.expr(sp.Mul)
        if len(subMuls) == 0:
            subMuls = [mono]
        
    #ugg

coeffs = np.array([sp.Symbol("c[{0}]".format(x)) for x in range(mono_len)])
r = coeffs.dot(monos)
#print(r.simplify())
#print(sp.polys.polyfuncs.horner(sp.Poly(r), wrt=x))
cstrs = " + ".join([ "c[{0}] * ".format(idx) + str((m)) for (idx, m) in enumerate(monos)])
with open("ugg.txt", "w+") as f:
    f.write(cstrs)
exit(0)

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


def buildGradAndHess(fileName):
    arrays = []
    for (idx, datam) in enumerate(data):
        print("Idx:", idx)
        polyText = datam["poly"]
        poly = eval(polyText)
        gradPoly = [sp.diff(poly, x), sp.diff(poly, y), sp.diff(poly, z)]

        gradPolyArray = np.array(gradPoly)
        hessianPoly = [[sp.diff(g, x), sp.diff(g, y), sp.diff(g, z)]
                       for g in gradPoly]
        grads = []
        hess = []
        for g in gradPolyArray:
            grads.append(mono_coeffs(sp.Poly(g), monos))
        for grad in hessianPoly:
            for g in grad:
                hess.append(mono_coeffs(sp.Poly(g), monos))
        arrays.append({"grad": grads, "hess": hess})
    text = pickle.dumps(arrays)
    with open(fileName, "wb+") as f:
        f.write(text)


def buildPVOData(fileName):
    arrays = []
    for (idx, datam) in enumerate(data):
        print("Idx:", idx)
        polyText = datam["poly"]
        poly = eval(polyText)
        gradPoly = [sp.diff(poly, x), sp.diff(poly, y), sp.diff(poly, z)]

        gradPolyArray = np.array(gradPoly)
        hessianPoly = [[sp.diff(g, x), sp.diff(g, y), sp.diff(g, z)]
                   for g in gradPoly]
        doted = np.array(hessianPoly).dot(gradPolyArray)
        wfield = np.cross(gradPolyArray, doted)
        p1 = sp.Poly(wfield[0])
        p2 = sp.Poly(wfield[1])
        p3 = sp.Poly(wfield[2])
        coeffs1 = mono_coeffs(p1, monos11)
        coeffs2 = mono_coeffs(p2, monos11)
        coeffs3 = mono_coeffs(p3, monos11)
        print(p1.total_degree())
        print(p2.total_degree())
        print(p3.total_degree())
        arrays.append([coeffs1, coeffs2, coeffs3])
    text = pickle.dumps(arrays)
    with open(fileName, "wb+") as f:
        f.write(text)
        

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

#buildPVOData("pvo.pkl")
# result = buildData()
# pklfilePrime = "scarray.pkl"
# txt = pickle.dumps(result)
# with open(pklfilePrime, "wb+") as f:
#     f.write(txt)
