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
    
    alignness = np.round(np.mean(gamma_x)/domain.shape[1],4)
    
    return(alignness, gamma_x/domain.shape[1])



# Add the alignedness profile later.



