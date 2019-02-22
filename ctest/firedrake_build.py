from firedrake import *
import tsfc
import ufl
import sympy as sp
import json
import numpy as np

EPSILON = 0.000000001

def computeInverseIndex(shape):
    subArraySizes = [product(shape[0:x+1]) for x in range(len(shape))]
    subArraySizes.reverse()
    subArraySizesLeft = subArraySizes[1:]
    arraySize = subArraySizes[0]
    
    def tab(x):
        collect = []
        for y in subArraySizesLeft:
            mod = x % y
            div = x // y
            collect.append(div)
            x = mod
        collect.append(x)
        return(collect)
    idxs = [tab(i) for i in range(arraySize)]
    return(idxs)

def cleanEpsilon(x):
    xp = float(x)
    if abs(xp) < EPSILON:
        return(0.0)
    else:
        return(xp)

def monos(degree, var):
    # We order our nonomials based on the following indexing algorithim
    dim = len(var)
    shape = [degree + 1 for i in range(dim)]
    idxes = computeInverseIndex(shape)
    def makeVar(idxes):
        temp = [a**b for (a, b) in zip(var, idxes)]
        result = product(temp)
        return(result)
    return([makeVar(idx) for idx in idxes])
def processSympyPoly(poly, var):
    maxDegree = poly.degree()
    polyVars = var  # [sp.Symbol("x{0}".format(i)) for i in range(dim)]
    monoPolys = monos(maxDegree, polyVars)
    def getCoeff(x,p):
        try:
            return(x.coeff_monomial(p))
        except:
            return(0.0)
                
    expansion = [cleanEpsilon(getCoeff(poly,mono)) for mono in monoPolys]
    result = {"coeffs": expansion, "degree": maxDegree}
    return((result, maxDegree))

def elementToBasis(elem):
    fiatElem = tsfc.fiatinterface.create_element(elem)
    dim  = elem.cell().topological_dimension()
    # check an assert with top and goem dim
    if isinstance(elem, ufl.MixedElement): # we don't have access to fiat classes...
        firstFiatElem = fiatElem.elements()[0]
        if all([firstFiatElem == x for x in fiatElem.elements()]):
            fiatElem = firstFiatElem
        else:
            raise Exception("We don't handle non-uniform mixed elements - use .split()")

    if isinstance(elem, ufl.TensorProductElement):
        raise Exception("We don't handle tensor product elements yet")
    if isinstance(elem, ufl.EnrichedElement):
        raise Exception("We don't handle direct sum elements yet")
    if isinstance(elem, ufl.RestrictedElement):
        raise Exception("We don't handle restricted elements yet")
    var = [sp.Symbol("x{0}".format(i)) for i in range(dim)]
    basisFunctions = fiatElem.tabulate(0, var)
    basisFunctions = basisFunctions[tuple([0 for i in range(dim)])]
    print(basisFunctions)
    # should be a straightup array
    if len(basisFunctions.shape) != 1:
        raise Exception("Don't know how to handle non-vector basis shape")
    basis = [processSympyPoly(sp.Poly(x), var) for x in basisFunctions]
    return((basis, dim))
    
#accelerate function:
#newton function
#refCell

def makeAccelerate():
    a = {"insert": "spat.in",
         "conservative": True,
         "includes": [
             "<spatialindex/capi/sidx_api.h>",
             "<spatialindex/capi/sidx_impl.h>"
         ],
         "linkDirs": [],
         "includeDirs": [],
         "libs": [
             "lspatialindex_c",
             "lspatialindex"
         ]}
    return(a)

def makeRefCell(geometry, epsilon=0.000000001, itters=16,
                newtonTol=0.00000001):
    a = {"type": "other",
         "epsilon": epsilon,
         "newtonParams": {"contraction": True, "itters": itters,
                          "killAfterTol": False, "newtonTol": newtonTol},
         "geometry": geometry}
    return(a)

def buildGeometry(elem):
    fiatElement = tsfc.fiatinterface.create_element(elem)
    refElem = fiatElement.get_reference_element()
    vertices = list(refElem.get_vertices())
    toNode = list(range(len(vertices)))
    topology = refElem.get_topology()
    dims = list(topology.keys())[1:-1]  # higher than points, but not itself.
    objects = {}
    for x in dims:
        key = "object" + str(x)
        entities = topology[x]
        entity = []
        print(entities)
        for e in entities:
            ep = entities[e]
            obj = {"entity": ep, "plane": ep[0:x + 1]}
            entity.append(obj)
        objects[key] = entity
    objects["vertices"] = vertices
    objects["toNode"] = toNode
    return(objects)
#glue it all togeather now.

def makeConstant(name, ty, value):
    a = {"name": name,
         "type": ty,
         "value": value}
    return(a)

def processSpace(V):
    mesh = V.mesh()
    meshElem = mesh.ufl_coordinate_element()
    spaceElem = V.ufl_element()
    (meshBasis, meshDim) = elementToBasis(meshElem)
    (spaceBasis, meshDimP) = elementToBasis(spaceElem)

    geometry = buildGeometry(meshElem)
    refCell = makeRefCell(geometry)
    accelerate = makeAccelerate()
    meshMapDim = len(meshBasis)
    meshBasisP = [x for (x,y) in meshBasis]
    spaceBasisP = [x for (x,y) in spaceBasis]
    meshMapDegree = max([y for (x,y) in meshBasis])
    print(meshMapDegree)
    mesh = {
        "basis": {"polys": meshBasisP},
        "constants": [
            makeConstant("dim", "int", meshDim),
            makeConstant("meshMapDim", "int", meshMapDim),
            makeConstant("degree", "int", meshMapDegree)
        ],
        "refCell": refCell,
        "accelerate": accelerate
    }
    spaceMapDim = len(spaceBasis)
    spaceMapDegree = max([y for (x,y) in spaceBasis])
    shape = list(V.shape)
    shapeLen = len(shape)
    space = {
        "basis" : {"polys" : spaceBasisP},
        "constants" : [
            makeConstant("dim", "int", meshDimP),
            makeConstant("spaceMapDim", "int", spaceMapDim),
            makeConstant("shape", "int[" + (str(shapeLen)) + "]", shape),
            makeConstant("degree", "int", spaceMapDegree)
        ]
    }
    function = {"constants": []}
    result = {
        "constants": [],
        "mesh": mesh,
        "space": space,
        "function": function
    }
    return(result)

def spaceToJson(V, jsonFile):
    dictJson = processSpace(V)
    with open(jsonFile, "w+") as f:
        dumped = json.dumps(dictJson, indent=4)
        f.write(dumped)


mesh = UnitSquareMesh(5,5)
space = FunctionSpace(mesh, "Lagrange", 4)
spaceToJson(space, "test.json")
