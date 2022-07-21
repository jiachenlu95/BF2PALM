import numpy as np
names = globals()


def cal_lp(topo):
    t = 0
    for i in range(topo.shape[0]):
        for j in range(topo.shape[1]):
            if topo[i,j] != 0:
                t+=1
    return(t/topo.size)

def cal_lf(domain,dxdy):
    tt = 0
    for x in range(domain.shape[0]):
        tmp = domain[x,:]
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
    return(tt*dxdy/(domain.size*dxdy*dxdy))
    

def cal_lw(domain,dxdy):
    # Construct a new topo with empty space equals nan
    cc = np.zeros(domain.shape)
    for i in range(domain.shape[0]):
        for j in range(domain.shape[1]):
            if domain[i,j] != 0:
                cc[i,j] = domain[i,j]
            else:
                cc[i,j] = np.nan
    return((Lw(cc,dxdy)+Lw(np.transpose(cc),dxdy))/(cc.size*dxdy*dxdy))


    

def cal_alignness(domain): 
    gamma_x = np.zeros([domain.shape[0]])
    
    # the idea is to pick up all continuous (x>=2) empty space, document their index, and evaluate the "alignness" of this building array
    
    for x in range(domain.shape[0]):
        row_x = domain[x,:]
        names['slice'+str(x)] = [] # list of continued streets
        if row_x.any()!=0: # In case encountered a zero array
            empty_i = np.where(row_x==0)[0] # collect the index of empty grids
            BC = False
            first = False
            if (empty_i[0] == 0) and (empty_i[-1] == domain.shape[1]-1): # if both ends is open.
                BC = True
                first = True # Boolean for the periodic BC
            empty_0 = empty_i[0] # The first empty grid for each row in the alongwind direction.
            
            for i in range(empty_i.size): # Loop over every empty grid
                try:
                    if (empty_i[i]+1)!=(empty_i[i+1]): # if encounters incontinued empty grids - streets discontinued
                        empty_1 = empty_i[i] # the end point of the continued street
                        names['slice'+str(x)].append(np.linspace(empty_0,empty_1,empty_1-empty_0+1,dtype='int')) # record the index for this continued street
                        if (empty_1-empty_0+1)>gamma_x[x]: # if this street is longer than the longest record.
                            gamma_x[x] = (empty_1-empty_0+1) # normalized by building height ahead of this street section.
                            if first: # if it is the first the continued street
                                tmp1 = empty_1-empty_0+1 # record this value for the following periodic BC consideration
                                first = False 
                        empty_0 = empty_i[i+1] # update the starting point for the continued street.  
                        
                except: # for the last continued street, the second phase of considering periodic BC - adding 
                    if BC:
                        tmp = empty_i[-1] - empty_0 + 1 + tmp1
                    else:
                        tmp = empty_i[-1] - empty_0 + 1 
                        
                    if tmp > gamma_x[x]:
                            gamma_x[x] = tmp
                        
        else:
            gamma_x[x] = domain.shape[1]#/bldH # empty street.
    
    alignness = np.round(np.mean(np.array(gamma_x))/domain.shape[1],4)
    
    return(alignness, np.array(gamma_x)/domain.shape[1])

def Lw(topo,dxdy):  # input the topography and output the pressure on windward and leeward face. 
 
    def pressureDefPre(topo): # Count how much pb and pf was sampled from topo
        pfNN = np.zeros(topo.shape[1])
        pbNN = np.zeros(topo.shape[1])

        for j in range(topo.shape[1]):
            for i in range(topo.shape[0]):
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface
                        pfNN[j] +=1
                        #print('cao')
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face
                        pbNN[j] +=1
                except:
                    0  
        return(pfNN,pbNN)

    pfNN,pbNN = pressureDefPre(topo)
    
    pf = []
    pb = []
    pfO = []
    pbO = []
    distO = []
    pfN = 0
    pbN = 0
    dist = []
    o1 = 0; o2 = 0; o3 = 0
    for j in range(topo.shape[1]):
        
        if ~np.isnan(topo[0,j]) and ~np.isnan(topo[-1,j]): # no B grid in the first and the end - normal row
            for i in range(topo.shape[0]):
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface
                        pf.append(topo[i-1,j])
                        pfN += 1
                        itmp = i
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face
                        pb.append(topo[i+1,j])
                        pbN += 1
                        dist.append(i-itmp+1) # index of the paired pf and pb
                except:
                    0
            
        if ~np.isnan(topo[0,j]) and np.isnan(topo[-1,j]): # no B grid in the first but the end
            count = 0
            #print('Outlier 1 found at' + str(j) + ' th row')
            o1+=1
            for i in range(topo.shape[0]):

                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface normal except for the last
                        if  count == pfNN[j]-1: # if that's the last frontal face 
                            
                            pf.append(topo[i-1,j])
                            pb.append(topo[0,j])
                            pfN += 1
                            pbN += 1
                            dist.append(topo.shape[0]-i) # should be the length of the last continued building grids
                            #print(topo.shape[0]-i)
                            
                        else:
                            pf.append(topo[i-1,j])
                            pfN += 1
                            itmp = i
                            count += 1              
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face normal
                        pb.append(topo[i+1,j])
                        pbN += 1
                        dist.append(i-itmp+1) # index of the paired pf and pb
                        
                except:
                    1
                          
        if np.isnan(topo[0,j]) and ~np.isnan(topo[-1,j]): # no B grid in the end but the first
            
            first = True
            #print('Outlier 2 found at' + str(j) + ' th row')
            o2+=1
            
            for i in range(topo.shape[0]):
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface normal
                        #print('sampling pf'+str(i)+str(j))
                        pf.append(topo[i-1,j])
                        pfN += 1
                        itmp = i
                        
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face normal except for the first
                        if first: 
                            #print('sampling pb'+str(i)+str(j))
                            pb.append(topo[i+1,j])
                            pbN += 1
                            dist.append(i) # grid count to the end - distance fixed
                            first = False
                        else:
                            pb.append(topo[i+1,j])
                            #print('sampling pb'+str(i)+str(j))
                            pbN += 1
                            dist.append(i-itmp+1) # index of the paired pf and pb     
                except:
                    2

        if np.isnan(topo[0,j]) and np.isnan(topo[-1,j]): # B grid in the end and the first, record them in a seperate array
            count = 0
            first = True
            o3+=1
            #print('Outlier 3 found at' + str(j) + ' th row')
            
            for i in range(topo.shape[0]):
                
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface normal except for the last
                        
                        if  count == pfNN[j]-1: # if that's the last frontal face
                            pfO.append(topo[i-1,j])
                            upSize = topo.shape[0]-i
                            
                        else:
                            pf.append(topo[i-1,j])
                            pfN += 1
                            itmp = i
                            count += 1
                except:
                    3
                            
            for i in range(topo.shape[0]):
                
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face normal except for the first
                        
                        if first:
                            #print('qppen')
                            pbO.append(topo[i+1,j])
                            lowSize = i
                            distO.append(lowSize+upSize+1) # index of the paired pf and pb

                            first = False
                            
                        else:
                            pb.append(topo[i+1,j])
                            pbN += 1
                            dist.append(i-itmp+1) # index of the paired pf and pb
                        
                except:
                    3
                    
            
    pf.extend(pfO)
    pb.extend(pbO)
    dist.extend(distO)
                    
    return((np.array(pf)+np.array(pb)).sum()*dxdy)

# Add the alignedness profile later.



