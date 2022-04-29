def pressureDefT(topo):   
    
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
    topo = np.transpose(topo)
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
    return(np.array(pf),np.array(pb),np.array(dist))