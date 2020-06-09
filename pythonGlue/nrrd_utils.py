# requires pynrrd: pip install pynrrd
import nrrd
import numpy as np
import functools

def writeSequence(data, fileName, dataKind=""):
    count = len(data.shape) - 1
    endShape = list(range(1, count + 1))
    targetShape = tuple(endShape + [0])
    result = np.transpose(data, targetShape)
    if dataKind != "":
        header = {"kinds": [dataKind, "list"]}
    else:
        header = {"kinds": ["list"]}
    nrrd.write(fileName, result, header=header)
def get(fileName):
    readData, _ = nrrd.read(fileName)
    return(readData)
def expectedOrder(data):
    shape = data.shape
    count = len(shape) - 1
    startShape = list(range(1, count + 1))
    endShape = tuple(startShape + [0])
    result = np.transpose(data, endShape)
    return(result)


def readLines(controlFile, seqFile, dim=3):
    seqData, _ = nrrd.read(seqFile)
    size = functools.reduce(lambda x,y: x*y, seqData.shape, 1)
    seqData = seqData.T
    controlData, _ = nrrd.read(controlFile)
    result = np.empty(controlData.shape[1:], dtype=object)
    for j in range(controlData.shape[1]):
        for k in range(controlData.shape[2]):
            offset = controlData[0][j][k]
            limit = controlData[1][j][k]
            result[j][k] = []
            for idx in range(0, limit):
                if dim > 1:
                    result[j][k].append([seqData[offset + idx][i]
                                         for i in range(dim)])
                else:
                    result[j][k].append(seqData[offset + idx])
                                        



    res = result
    return(res)
