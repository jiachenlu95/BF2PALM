import numpy as np
import math



def count_edge(xy,lineO,lineL,res):
    lineO = lineO
    lineL = lineL
    for line in range(np.shape(xy)[0]-1):
        lineO.append(np.rad2deg(np.arctan((xy[line][1]-xy[line+1][1])/(xy[line][0]-xy[line+1][0]+1e-10))))
        lineL.append(math.dist(xy[line], xy[line+1])*res)
    return(lineO,lineL)
        
        
def cal_entropy(lineO,lineL):
    angle = np.zeros(np.shape(lineO)[0])
    
    for i in range(np.shape(lineO)[0]):
        angle[i] = np.round(lineO[i]/10)*10

    x,y=zip(*sorted(zip(angle.flatten(),np.array(lineL).flatten())))

    weighted = np.zeros(np.unique(angle).size)
    i = 0
    idx = 0

    while i < angle.size-1:
        if x[i] == x[i+1]:
            weighted[idx] = weighted[idx] + y[i]
            i += 1
        else:

            weighted[idx] = weighted[idx] + y[i]
            i += 1
            idx += 1
    weighted[idx] += y[-1]
    angle = np.unique(angle)
    
    if (angle[0] == -90) and (angle[-1] == 90):
        print('Merged two angle data')
    #weighted.size == 19: # in case 19 directions are recorded
        weighted[-1] += weighted[0]
        weighted = np.delete(weighted, 0)
        angle = np.delete(angle, 0)


    
    ################
    x = np.deg2rad(angle)
    y = weighted
    #################
    
    
    return(x,y)


