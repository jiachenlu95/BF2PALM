import osmnx as ox
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point,Polygon
import math
import re
import pandas as pd



from . import topo
from . import entropy
from . import rasterization
from . import results
from . import evalgeo


def extract_domain(res,bldH,LLbbox,name):
	"""
	Main
	"""
	# Resolution on the x and y direction is set to equal.
	resx = res
	resy = res
	bldH = bldH


	# Convert bbox to format that OSM can read and download the building footprint as geopandas dataframe
	bbox,bbox_cor = topo.select_region(LLbbox)
	# Download data from OSM with bbox as input
	gdf = ox.geometries_from_bbox(bbox_cor[0],bbox_cor[1],bbox_cor[2],bbox_cor[3], {"building": True})
	gdf.reset_index(level=0, inplace=True)
	# Extract the geometry and building level
	geom = gdf['geometry'].values
	bHeight = gdf['building:levels'].values

	# Calculate the centriod of each urban structure
	# If centriod of one structure falls in the selected polygon, include it in the domain
	centroid = []
	for idxP in range(geom.shape[0]):
	    centroid.append(geom[idxP].centroid)
	data = {'geometry':geom,'cent':centroid}
	# reproduce the dataframe with centeriod
	gdf = pd.DataFrame(data)

	x_mul,y_mul = topo.projection(bbox_cor,resx,resy)

	domain,gdf,x_start,y_start = rasterization.initialize(bbox,gdf,x_mul,y_mul,bbox_cor)

	lineO = [] # edge orientation
	lineL = [] # edge length
	area = [] # total area
	areaH = [] # area with realistic building height


	for idxP in range(gdf.shape[0]):# if centriod of the structure fall within the bbox
	    if Point(gdf['cent'][idxP][0][0],gdf['cent'][idxP][1][0]).within(Polygon(bbox)):
	        [domain,areaH,area] = rasterization.process_building(lineO,lineL,domain,idxP,gdf,bHeight,x_mul,y_mul,x_start,y_start,bbox_cor,bldH,area,areaH,res)

	angle,weighted = entropy.cal_entropy(lineO,lineL)

	results.print_domain(domain) 
	phi = results.print_entropy(angle,weighted) 
	results.print_results(phi,domain,area,areaH,bldH)

	angleRotate = 0
	results.showDiagram(domain,angle,weighted,phi,angleRotate,np.array(area),name)

	return(domain)


def extract_domainKML(res,bldH,LLbbox,name):
	"""
	Main
	"""
	# Resolution on the x and y direction is set to equal.
	resx = res
	resy = res
	bldH = bldH


	# Convert bbox to format that OSM can read and download the building footprint as geopandas dataframe
	#bbox,bbox_cor = topo.select_region(LLbbox)

	bbox = LLbbox


	ymin = np.min([sum(bbox, [])[0],sum(bbox, [])[2],sum(bbox, [])[4],sum(bbox, [])[6]])
	ymax = np.max([sum(bbox, [])[0],sum(bbox, [])[2],sum(bbox, [])[4],sum(bbox, [])[6]])
	xmin = np.min([sum(bbox, [])[1],sum(bbox, [])[3],sum(bbox, [])[5],sum(bbox, [])[7]])
	xmax = np.max([sum(bbox, [])[1],sum(bbox, [])[3],sum(bbox, [])[5],sum(bbox, [])[7]])
	bbox_cor = np.array([xmax,xmin,ymin,ymax])


	# Download data from OSM with bbox as input
	gdf = ox.geometries_from_bbox(bbox_cor[0],bbox_cor[1],bbox_cor[2],bbox_cor[3], {"building": True})

	gdf.reset_index(level=0, inplace=True)
	# Extract the geometry and building level
	geom = gdf['geometry'].values
	bHeight = gdf['building:levels'].values

	# Calculate the centriod of each urban structure
	# If centriod of one structure falls in the selected polygon, include it in the domain
	centroid = []
	for idxP in range(geom.shape[0]):
	    centroid.append(geom[idxP].centroid)
	data = {'geometry':geom,'cent':centroid}
	# reproduce the dataframe with centeriod
	gdf = pd.DataFrame(data)

	x_mul,y_mul = topo.projection(bbox_cor,resx,resy)

	domain,gdf,x_start,y_start = rasterization.initialize(bbox,gdf,x_mul,y_mul,bbox_cor)

	lineO = [] # edge orientation
	lineL = [] # edge length
	area = [] # total area
	areaH = [] # area with realistic building height


	for idxP in range(gdf.shape[0]):# if centriod of the structure fall within the bbox
	    if Point(gdf['cent'][idxP][0][0],gdf['cent'][idxP][1][0]).within(Polygon(bbox)):
	        [domain,areaH,area] = rasterization.process_building(lineO,lineL,domain,idxP,gdf,bHeight,x_mul,y_mul,x_start,y_start,bbox_cor,bldH,area,areaH,res)

	angle,weighted = entropy.cal_entropy(lineO,lineL)

	results.print_domain(domain) 
	phi = results.print_entropy(angle,weighted) 
	results.print_results(phi,domain,area,areaH,bldH)

	angleRotate = 0
	results.showDiagram(domain,angle,weighted,phi,angleRotate,np.array(area),name)

	return(domain)