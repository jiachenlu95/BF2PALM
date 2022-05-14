import numpy as np
import re
import math


names = globals()

def Parse_GE(bbox):
    """
    
    Parse the longitude/latitude coordinates extracted from GE to osmnx input
    
    Parameters
    -----------
    bbox : list
    List of LL coordinates extracted from Google Earth
    
    Returns
    -------
    bbox : list
    List of LL coordinates
        
    """
    
    if bbox[-1] in ['N', 'E']:
        multiplier = 1 
    else:
        multiplier = -1
 
    return multiplier * sum(float(x) / 60 ** n for n, x in enumerate(re.split('°|\'|\"', bbox[:-2])))

def cal_WGSdist(Gx1,Gx2,Gy1,Gy2): 
    '''
    Calculate real distance of two points with longitude and latitude
    '''
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

def shear(angle,x,y):
    
    """
    Rasterize 
    
    |1  -tan(𝜃/2) |  |1        0|  |1  -tan(𝜃/2) | 
    |0      1     |  |sin(𝜃)   1|  |0      1     |
    Parameters
    -----------
    bbox :
    
    res : 
    
    bldH :  
    
    Returns
    -------
    
    
    is_polygon : bool
        True if the tags are for a polygon type geometry
        
    """
        
    # shear 1
    tangent=math.tan(angle/2)
    new_x=round(x-y*tangent)
    new_y=y
    
    #shear 2
    new_y=round(new_x*math.sin(angle)+new_y)      #since there is no change in new_x according to the shear matrix

    #shear 3
    new_x=round(new_x-new_y*tangent)              #since there is no change in new_y according to the shear matrix
    
    return new_y,new_x



def rotate(image,angle):
    """
    Rasterize 
    
    Parameters
    -----------
    bbox :
    
    res : 
    
    bldH :  
    
    Returns
    -------
    
    
    is_polygon : bool
        True if the tags are for a polygon type geometry
        
    """
    angle=angle               # Ask the user to enter the angle of rotation

    # Define the most occuring variables
    angle=math.radians(angle)                               #converting degrees to radians
    cosine=math.cos(angle)
    sine=math.sin(angle)

    height=image.shape[0]                                   #define the height of the image
    width=image.shape[1]                                    #define the width of the image

    # Define the height and width of the new image that is to be formed
    new_height  = round(abs(image.shape[0]*cosine)+abs(image.shape[1]*sine))+1
    new_width  = round(abs(image.shape[1]*cosine)+abs(image.shape[0]*sine))+1

    # define another image variable of dimensions of new_height and new _column filled with zeros
    output=np.zeros([new_height,new_width])
    image_copy=output.copy()


    # Find the centre of the image about which we have to rotate the image
    original_centre_height   = round(((image.shape[0]+1)/2)-1)    #with respect to the original image
    original_centre_width    = round(((image.shape[1]+1)/2)-1)    #with respect to the original image

    # Find the centre of the new image that will be obtained
    new_centre_height= round(((new_height+1)/2)-1)        #with respect to the new image
    new_centre_width= round(((new_width+1)/2)-1)          #with respect to the new image


    for i in range(height):
        for j in range(width):
            #co-ordinates of pixel with respect to the centre of original image
            y=image.shape[0]-1-i-original_centre_height                   
            x=image.shape[1]-1-j-original_centre_width 

            #Applying shear Transformation                     
            new_y,new_x=shear(angle,x,y)
            new_y=new_centre_height-new_y
            new_x=new_centre_width-new_x

            output[new_y,new_x]=image[i,j]                          #writing the pixels to the new destination in the output image


    plt.imshow(output)
    plt.axis('equal')
    return(output)


def select_region(bbox):
    """
    
    Project the longitude/latitude coordinates to input of osmnx
    
    Parameters
    -----------
    bbox : list
    List of LL coordinates extracted from Google Earth
    
    Returns
    -------
    bbox_osm : array
    
    bbox_osm_ : list
        
    """
    bbox_osm = []
    for i in range(np.shape(bbox)[0]):
        tmp = [Parse_GE(bbox[i][1]),Parse_GE(bbox[i][0])]
        bbox_osm.append(tmp)
        
    ymin = np.min([sum(bbox_osm, [])[0],sum(bbox_osm, [])[2],sum(bbox_osm, [])[4],sum(bbox_osm, [])[6]])
    ymax = np.max([sum(bbox_osm, [])[0],sum(bbox_osm, [])[2],sum(bbox_osm, [])[4],sum(bbox_osm, [])[6]])
    xmin = np.min([sum(bbox_osm, [])[1],sum(bbox_osm, [])[3],sum(bbox_osm, [])[5],sum(bbox_osm, [])[7]])
    xmax = np.max([sum(bbox_osm, [])[1],sum(bbox_osm, [])[3],sum(bbox_osm, [])[5],sum(bbox_osm, [])[7]])  
    
    bbox_osm_ = np.array([xmax,xmin,ymin,ymax])
    
    return(bbox_osm,bbox_osm_)


def projection(bbox_cor,resx,resy):
    
    """
    Calculate the projection mutiplier  
    
    Parameters
    -----------
    bbox :
    
    
    Returns
    -------
    x_mul : float
    y_mul : float

    Multiplier in x and y direction
    
    
    is_polygon : bool
        True if the tags are for a polygon type geometry
        
    """
    
    # Collect all values together
    [y_max,y_min,x_min,x_max] = bbox_cor
    
    
    dist_x = cal_WGSdist(x_min,x_max,y_min,y_min)/resx
    
    dist_y = cal_WGSdist(x_min,x_min,y_min,y_max)/resy
    # multiplier on x and y coordinates
    
    x_mul = dist_x / (x_max - x_min)
    y_mul = dist_y / (y_max - y_min)
    
    
    return(x_mul,y_mul) 

def centroid(vertexes):
    x_list = [vertex [0] for vertex in vertexes]
    y_list = [vertex [1] for vertex in vertexes]
    len = len(vertexes)
    x = sum(x_list)/_len
    y = sum(y_list)/_len
    
    return(x, y)  


def noHoles(domain,ths):
    """
    Eliminate small and unconnected air space -- ths ~ 3

    """
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


def pressureDefT(topo):   
    """
    Identifiy windward and leeward grid for realistic irregular geometry.
    """
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
    
    # Initialization
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
