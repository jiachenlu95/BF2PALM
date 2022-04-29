import numpy as np
import BF2PALM as bm

def selectRegion(strList):
    cc = []
    for i in range(np.shape(strList)[0]):
        tmp = [bm.convert(strList[i][1]),bm.convert(strList[i][0])]
        cc.append(tmp)
    return(cc)