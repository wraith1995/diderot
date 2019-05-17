from firedrake import *
import tsfc
import ufl
import sympy as sp
import json
import numpy as np
import ctypes as ct
import stucts as st
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
    maxDegree = poly.total_degree()
    polyVars = var  # [sp.Symbol("x{0}".format(i)) for i in range(dim)]
    monoPolys = monos(maxDegree, polyVars)
    def getCoeff(x, p):
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

    # should be a straightup array
    if len(basisFunctions.shape) != 1:
        raise Exception("Don't know how to handle non-vector basis shape")
    basis = [processSympyPoly(sp.Poly(x), var) for x in basisFunctions]
    return((basis, dim))

# would be nice if this were just pure c
def buildCellConnections(mesh, ty):
    mesh.init() # we need mesh topology now.
    numCell = mesh.num_cells()
    facetsPercell = mesh.ufl_cell().num_facets()
    result = (-1) * np.ones((numCell, facetsPercell, 2), dtype=ty)
    interiorFacets = mesh.interior_facets
    facetCellMap = interiorFacets.facet_cell_map.values
    localFacetData = interiorFacets.local_facet_dat.data
    num_facets = localFacetData.shape[0]
    for x in range(num_facets):
        [c1, c2] = facetCellMap[x]
        [l1, l2] = localFacetData[x]
        result[c1][l1][0] = c2
        result[c1][l1][1] = l2
        result[c2][l2][0] = c1
        result[c2][l2][1] = l1
    print(result.shape)
    return(result) # all exterior facets are (-1, -1)

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

def makeRefCell(geometry, start, refCellDefault="other", epsilon=0.000000001, itters=16,
                newtonTol=0.00000001):
    a = {"type": refCellDefault,
         "epsilon": epsilon,
         "newtonParams": {"contraction": True, "itters": itters,
                          "killAfterTol": False, "newtonTol": newtonTol, "start": start},
         "geometry": geometry}
    return(a)

def buildGeometry(elem):
    fiatElement = tsfc.fiatinterface.create_element(elem)
    refElem = fiatElement.get_reference_element()
    vertices = list(refElem.get_vertices())
    start = list(np.average(vertices, axis=0))
    toNode = list(range(len(vertices)))
    topology = refElem.get_topology()
    dims = list(topology.keys())[1:-1]  # higher than points, but not itself.
    objects = {}
    for x in dims:
        key = "object" + str(x)
        entities = topology[x]
        entity = []

        for e in entities:
            ep = entities[e]
            obj = {"entity": ep, "plane": ep[0:x + 1]}
            entity.append(obj)
        objects[key] = entity
    objects["vertices"] = vertices
    objects["toNode"] = toNode
    return(objects,start)
#glue it all togeather now.

def makeConstant(name, ty, value):
    a = {"name": name,
         "type": ty,
         "value": value}
    return(a)

def processSpace(V, refCellDefault="other"):
    mesh = V.mesh()
    meshElem = mesh.ufl_coordinate_element()
    spaceElem = V.ufl_element()
    (meshBasis, meshDim) = elementToBasis(meshElem)
    (spaceBasis, meshDimP) = elementToBasis(spaceElem)

    (geometry, start) = buildGeometry(meshElem)
    refCell = makeRefCell(geometry, start, refCellDefault=refCellDefault)
    accelerate = makeAccelerate()
    meshMapDim = len(meshBasis)
    meshBasisP = [x for (x,y) in meshBasis]
    spaceBasisP = [x for (x,y) in spaceBasis]
    meshMapDegree = max([y for (x,y) in meshBasis])
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

def spaceToJson(V, jsonFile, refCellDefault="other"):
    dictJson = processSpace(V,refCellDefault=refCellDefault)
    with open(jsonFile, "w+") as f:
        dumped = json.dumps(dictJson, indent=4)
        f.write(dumped)

def passMeshHelper(mesh, intTy, floatTy, geometric=True):
    intPtr = ct.POINTER(intTy)
    floatPtr = ct.POINTER(floatTy)
    meshCoords = mesh.coordinates
    meshIndexMap = meshCoords.cell_node_map().values
    meshCoordsMap = meshCoords.dat.data

    dim = mesh.topological_dimension()
    meshMapDim = meshCoords.cell_node_map().values.shape[1]
    numCell = meshCoords.cell_node_map().values.shape[0]
    if geometric:
        sIndex = mesh.spatial_index and mesh.spatial_index.ctypes
        # copied from firedrake code;why though?
    else:
        sIndex = 0
    conBuild = buildCellConnections(mesh, intTy)
    con = conBuild
    return((meshIndexMap, meshCoordsMap, dim, meshMapDim, numCell, sIndex, con))

def passSpaceHelper(space, intTy):
    intPtr = ct.POINTER(intTy)
    spaceIndexMap = space.cell_node_map().values

    spaceDim = space.cell_node_map().values.shape[1]
    return(spaceIndexMap, spaceDim)

def passFuncHelper(func, floatTy):
    floatPtr = ct.POINTER(floatTy)
    funcCoordMap = func.dat.data
    return((funcCoordMap,))

def passAll(func, intTy, floatTy, geometric=True, extraData=[]):
    space = func.function_space()
    mesh = space.mesh()
    meshTuple = passMeshHelper(mesh, intTy, floatTy, geometric=geometric)
    spaceTuple = passSpaceHelper(space, intTy)
    funcTuple = passFuncHelper(func, floatTy)
    combined = meshTuple + spaceTuple + funcTuple
    (_, _, _, buildAll) = st.makeAllTypes(intTy, floatTy, extraTys=list(map(lambda x: (x[0], x[2]), extraData)))
    return((combined, buildAll(*combined,
                               extraVals=list(map(lambda x: x[1], extraData)))))
            
