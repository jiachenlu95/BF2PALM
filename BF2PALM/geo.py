import numpy as np
import re


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



def rotate(tmp,angle):
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
    image = tmp             # Load the image
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

