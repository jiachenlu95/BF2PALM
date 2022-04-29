import numpy as np
def cal_lf(cccc,dxdy):
    tt = 0
    for x in range(cccc.shape[0]):
        tmp = cccc[x,:]
        if tmp.any()!=0: # In case encountered a zero array
            cc = np.nonzero(tmp)
            cc0 = cc[0][0]
            slice1 = []
            for i in range(np.size(cc)):
                if i+1 == np.size(cc):
                    slice1.append(tmp[cc0:cc[0][-1]+1]) 
                elif(cc[0][i]+1)!=(cc[0][i+1]): # pick up a contiune substance
                    cc1 = cc[0][i]
                    slice1.append(tmp[cc0:cc1+1])
                    cc0 = cc[0][i+1]
            for k in range(np.shape(slice1)[0]):
                tt += np.max(slice1[k])
    return(tt*dxdy/(cccc.size*dxdy*dxdy))