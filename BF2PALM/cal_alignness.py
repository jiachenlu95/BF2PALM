import numpy as np
names = globals()

def cal_alignness(cccc):
    maxRow = np.zeros([cccc.shape[0]])
    
    # the idea is to pick up all continuous (x>=2) empty space, document their index, and evaluate the "alignness" of this building array  
    for x in range(cccc.shape[0]):
        tmp = cccc[x,:]
        names['slice'+str(x)] = []
        if tmp.any()!=0: # In case encountered a zero array
            cc = np.where(tmp==0)[0]
            cc0 = cc[0]
            first = True
            for i in range(cc.size):
                try:
                    if (cc[i]+1)!=(cc[i+1]):
                        cc1 = cc[i]
                        names['slice'+str(x)].append(np.linspace(cc0,cc1,cc1-cc0+1,dtype='int'))
                        if (cc1-cc0+1)>maxRow[x]:
                            maxRow[x] = cc1-cc0+1
                            if first:
                                tmp1 = cc1-cc0+1
                                first = False
                        cc0 = cc[i+1]
                        
                except:
                    tmp = cc[-1] - cc0 + 1 + tmp1
                    if tmp > maxRow[x]:
                        maxRow[x] = tmp
        else:
            maxRow[x] = cccc.shape[1]
    alignness = np.round(np.mean(maxRow)/cccc.shape[1],4)
    return(alignness)