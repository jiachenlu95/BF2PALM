import numpy as np
def noHoles(domain,ths):
    # 1 1 1
    # 1 x 1
    # 1 1 1
    tmp = domain
    ori = np.zeros(domain.shape)
    for i in range(domain.shape[0]-1):
        for j in range(domain.shape[1]-1):
            ori[i,j] = domain[i,j]
            t = 0
            total = 0
            if (domain[i,j] != domain[i,j+1]):
                t+=1
                total += domain[i,j+1]

            if (domain[i,j] != domain[i,j-1]):
                t+=1
                total += domain[i,j-1]
            if (domain[i,j] != domain[i+1,j]):
                t+=1
                total += domain[i+1,j]

            if (domain[i,j] != domain[i-1,j]):
                t+=1
                total += domain[i-1,j]

            if (domain[i,j] != domain[i+1,j+1]): 
                t+=1
                total += domain[i+1,j+1]

            if (domain[i,j] != domain[i-1,j-1]):
                t+=1
                total += domain[i-1,j-1]

            if (domain[i,j] != domain[i+1,j-1]): 
                t+=1
                total += domain[i+1,j-1]
            if (domain[i,j] != domain[i-1,j+1]):
                t+=1
                total += domain[i-1,j+1]
            if t > ths: 
                #print(domain[i,j])
                tmp[i,j] = 16 #int(total/10000)*16

                #print('point fixed')

    return(tmp,ori)