import numpy as np
import osmnx as ox
from shapely.geometry import Point,Polygon
import matplotlib.pyplot as plt
import math
import re
import pandas as pd
import BF2PALM as bm
names = globals()


def extractDomainOSM(bbox,res,bldH):

    resx = res
    resy = res
    bbox = bbox
    bb = bm.selectRegion(bbox)
    yy = np.min([sum(bb, [])[0],sum(bb, [])[2],sum(bb, [])[4],sum(bb, [])[6]]),np.max([sum(bb, [])[0],sum(bb, [])[2],sum(bb, [])[4],sum(bb, [])[6]]),
    xx = np.min([sum(bb, [])[1],sum(bb, [])[3],sum(bb, [])[5],sum(bb, [])[7]]),np.max([sum(bb, [])[1],sum(bb, [])[3],sum(bb, [])[5],sum(bb, [])[7]]),
    B = ox.geometries_from_bbox(xx[1], xx[0], yy[0], yy[1], {"building": True})# new
    B.reset_index(level=0, inplace=True)
    bounds = B.total_bounds
    y_min = bounds[0]
    x_min = bounds[1]
    y_max = bounds[2]
    x_max = bounds[3]
    geo = B['geometry'].values
    cent = []
    for idxP in range(geo.shape[0]):
        cent.append(geo[idxP].centroid)
    data = {'geometry':geo,'cent':cent}
    df = pd.DataFrame(data)
    bb = bm.selectRegion(bbox)
    xq = []
    yq = []
    for i in range(np.shape(bbox)[0]):
        xq.append(bb[i][0])
        yq.append(bb[i][1])

    x_min = np.min(xq)
    x_max = np.max(xq)
    y_min = np.min(yq)
    y_max = np.max(yq)
    dist_y = bm.cal_WGSdist(x_min,x_min,y_min,y_max)/resy
    dist_x = bm.cal_WGSdist(x_min,x_max,y_min,y_min)/resx
    x_mul = dist_x / (x_max - x_min)
    y_mul = dist_y / (y_max - y_min)
    aa = bbox
    d = {}
    i = 0
    for idxP in range(df.shape[0]):
        try:
            d[i] = {'geometry': df['geometry'][idxP].exterior.coords.xy, 'cent': df['cent'][idxP].coords.xy}
            i+=1
        except:
            i = i
    df = pd.DataFrame.from_dict(d, "index")
    qq = df
    bb = np.zeros(np.shape(aa))
    xxx = []
    yyy = []
    for i in range(np.shape(aa)[0]):
        bb[i][0] = int(np.round((bm.selectRegion(aa)[i][0] - x_min)*x_mul))
        bb[i][1] = int(np.round((bm.selectRegion(aa)[i][1] - y_min)*y_mul))
        xxx.append(int(np.round((bm.selectRegion(aa)[i][0] - x_min)*x_mul)))
        yyy.append(int(np.round((bm.selectRegion(aa)[i][1] - y_min)*y_mul)))

    domain = np.zeros([np.max(xxx)-np.min(xxx),np.max(yyy)-np.min(yyy)]) 
    Gx = math.ceil(dist_x)
    Gy = math.ceil(dist_y)
    lineO = []
    lineL = []
    area = []
    for idxP in range(qq.shape[0]):
        if Point(qq['cent'][idxP][0][0],qq['cent'][idxP][1][0]).within(Polygon(bm.selectRegion(aa))): 
            # poly.bounds
            nx = []
            ny = []
            xy = []
            Pol = []
            for i in range(np.shape(qq['geometry'][idxP])[1]):
                xx = int(np.round((qq['geometry'][idxP][0][i] - x_min)*x_mul))
                yy = int(np.round((qq['geometry'][idxP][1][i] - y_min)*y_mul))
                nx.append(xx)
                ny.append(yy)
                xy.append([xx,yy])
                Pol.append((xx-np.min(xxx),yy-np.min(yyy)))

            xmin = np.min(nx) - np.min(xxx)
            xmax = np.max(nx) - np.min(xxx)
            ymin = np.min(ny) - np.min(yyy)
            ymax = np.max(ny) - np.min(yyy)

        # Construct edge
        for ii in range(np.size(nx)):
            try:
                domain[nx[ii]-np.min(xxx),ny[ii]-np.min(yyy)] = bldH # Fill the outline with building height
            except:
                continue

        pol = Polygon(Pol)
        area.append(pol.area)
        # Construct interior by filling grids with building height
        for x in range (xmin,xmax):
            for y in range(ymin,ymax):
                p = Point(x,y)
                if p.within(pol) or p.intersects(pol):
                    try:
                        domain[x,y] = bldH  # Fill inside with building height
                    except:
                        continue

            for q in range(np.shape(xy)[0]-1):
                lineO.append(np.rad2deg(np.arctan((xy[q][1]-xy[q+1][1])/(xy[q][0]-xy[q+1][0]+1e-10))))
                lineL.append(math.dist(xy[q], xy[q+1])*res) 

    plt.imshow(domain)
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('domain.pdf')

         

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


    xtick_font = {
        "family": "DejaVu Sans",
        "size": 10,
        "weight": "bold",
        "alpha": 1.0,
        "zorder": 3,
    }

    color="#003366"
    edgecolor="k"
    linewidth=0.5
    alpha=0.7

    
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

    fig, ax = plt.subplots(figsize=(6,6), subplot_kw={"projection": "polar"})


    ax.set_theta_zero_location("N") # 
    ax.set_theta_direction("clockwise")
    
    ax.set_ylim(top=y.max())

    # configure the y-ticks and remove their labels
    ax.set_yticks(np.linspace(0, y.max(), 5))
    ax.set_yticklabels(labels="")

    # configure the x-ticks and their labels
    xticklabels = ["N", "", "E", "", "S", "", "W", ""]
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=xticklabels, fontdict=xtick_font)
    ax.tick_params(axis="x", which="major", pad=-2)

    ax.bar(
        x,
        height=y,
        width=0.1,
        align="center",
        bottom=0,
        zorder=2,
        color=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
    )
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    plt.tight_layout()
    plt.savefig('polarO.png',dpi=300)

    dist = 2 # Edge evenly distributed to just 2 orientations
    hMin = 1/dist*np.log(1/dist)*dist

    dist = 18 # Edge evenly distributed to all 18 orientations
    hMax = 1/dist*np.log(1/dist)*dist

    H_w = np.nansum(weighted/sum(weighted)*1*np.log(weighted/sum(weighted)*1))
    phi = 1-(hMax-H_w)/(hMax-hMin)


    print('$\\phi=$'+str(phi))

    print(domain.shape)

    return(domain,angle,weighted,phi,np.array(area))