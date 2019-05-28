import nrrd
import numpy as np
(a,_) = nrrd.read("ugg1_0.nrrd")
for b in a.T:
    print(np.linalg.norm(b))
