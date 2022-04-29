import numpy as np

def cal_WGSdist(Gx1,Gx2,Gy1,Gy2): # Calculate real distance of two points with longitude and latitude
    R = 6371 
    x1 = Gx1
    x2 = Gx2
    y1 = Gy1
    y2 = Gy2
    dLon = np.deg2rad(x2-x1)
    dLat = np.deg2rad(y2-y1)

    a = np.sin(dLat/2)**2 + np.cos(np.deg2rad(y1))*np.cos(np.deg2rad(y2))*np.sin(dLon/2)**2
    
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    dist = R*c*1000 # distance in m
    return(dist)