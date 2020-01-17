import sys
sys.path
sys.path.append('/home/teocollin/gitcode/diderot/pythonGlue')
import nrrd_utils as nu
import nrrd

data = nrrd.read("ugg_0.nrrd")[0]
dataOrdered = nu.expectedOrder(data)
for datam in dataOrdered:
    print(sum(datam))
