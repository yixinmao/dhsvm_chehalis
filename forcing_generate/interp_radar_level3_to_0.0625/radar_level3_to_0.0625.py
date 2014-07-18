#!/usr/local/bin/python

from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pdb
import datetime as dt
import os.path
from shapely.geometry import Polygon
from shapely.geometry import Point


####################################################################
########################## User defined ############################
####################################################################
liv_match_list_path = '/usr1/ymao/dhsvm_chehalis/forcing_generate/interp_radar_level3_to_0.0625/Livneh_points_latlon'  # station_ID; lat; lon
radar_level3_path = '/raid/ymao/CEE599_HW3_radar/NEXRAD_level3_3hourly_nc'
output_dir_base = '/usr1/ymao/dhsvm_chehalis/forcing_generate/interp_radar_level3_to_0.0625/output'
output_dir_prec_result = '/usr1/ymao/dhsvm_chehalis/forcing_generate/interp_radar_level3_to_0.0625/output/precip_interp_from_radar_level3_m'
Liv_vic_output_path = '/raid/ymao/dhsvm_chehalis/vic_generate_forcing/vic_results_from_Livneh_updated_1990-Apr2014'  # Livneh vic output path (3 hourly)
Liv_vic_prec_column = 5  # number of column of precipitation in Livneh vic output files; start from 1
Liv_vic_start_time = start_time = dt.datetime(year=1990, month=1, day=1, hour=0, minute=0)

cell_size = 0.0625  # size of each Livneh grid cell, deg
radar_lat = 47.1158
radar_lon = -124.1069

nazi = 360
ngate = 115
dazi = 1
dgate = 2000  

start_time = dt.datetime(year=2011, month=9, day=8, hour=1, minute=30)
end_time = dt.datetime(year=2014, month=4, day=13, hour=22, minute=30)

# calculate
liv_match_list = np.loadtxt(liv_match_list_path)
nstn = np.shape(liv_match_list)[0]
duration = end_time - start_time
ntime = (duration.days*24 + duration.seconds//3600)/3+1

##########################################################################################
####### Decide which radar points are in which 1/16 cell and calculate distance ##########
########################################################################################## 
print 'Calculting which radar points are in which Livneh grid cell...' 
m = Basemap(llcrnrlon=-124, llcrnrlat=46, urcrnrlon=-122, urcrnrlat=48, resolution='i', projection='tmerc', lon_0=radar_lon, lat_0=radar_lat)
radar_x, radar_y = m(radar_lon, radar_lat)

stn_radar_list = []

for stn in range(nstn):
	stn_ID = liv_match_list[stn, 0]
	stn_lat = liv_match_list[stn, 1]
	stn_lon = liv_match_list[stn, 2]
	stn_radar_list.append([])

	# calculate xy for the stn and four corners of the grid cell
	stn_x, stn_y = m(stn_lon, stn_lat)  # station xy
	corner_x = np.empty(4)
	corner_y = np.empty(4)
	corner_x[0], corner_y[0] = m(stn_lon-cell_size/2, stn_lat-cell_size/2)  # lower left
	corner_x[1], corner_y[1] = m(stn_lon+cell_size/2, stn_lat-cell_size/2)  # lower right
	corner_x[2], corner_y[2] = m(stn_lon-cell_size/2, stn_lat+cell_size/2)  # upper left
	corner_x[3], corner_y[3] = m(stn_lon+cell_size/2, stn_lat+cell_size/2)  # upper right

	# calculate the range of distance and angle between radar center and the four corners
	corner_dis = np.empty(4)
	corner_azi = np.empty(4)
	for i in range(4):
		corner_dis[i] = np.sqrt(np.power((corner_x[i]-radar_x),2) + np.power((corner_y[i]-radar_y),2))
		corner_azi[i] = 90 - np.arctan((corner_y[i]-radar_y)/(corner_x[i]-radar_x))/np.pi*180.0

	# buffer a little; calculate which radar points are in which 1/16 grid cell
	min_dis_ind = int(round(np.min(corner_dis)/dgate) - 1)  # index start from 0
	max_dis_ind = int(round(np.max(corner_dis)/dgate) + 1)
	min_azi_ind = int(round(np.min(corner_azi/dazi)) - 1)  # index start from 0
	max_azi_ind = int(round(np.max(corner_azi/dazi)) + 1) 
	poly = Polygon(((corner_x[0],corner_y[0]),(corner_x[1],corner_y[1]),(corner_x[3],corner_y[3]),(corner_x[2],corner_y[2])))
	for azi_ind in range(min_azi_ind, max_azi_ind+1):
		for dis_ind in range(min_dis_ind, max_dis_ind+1):
			azi = azi_ind * dazi
			dis = dis_ind * dgate
			x = radar_x + dis * np.sin(azi/180.0*np.pi)
			y = radar_y + dis * np.cos(azi/180.0*np.pi)
			point = Point(x,y)
			if point.within(poly)==True:  # if this radar point is inside this 1/16 grid cell
				distance = np.sqrt(np.power((x-stn_x),2) + np.power((y-stn_y),2))
				stn_radar_list[stn].append([azi_ind, dis_ind, distance]) 


####################################################################
##################### Interpolate data #############################
####################################################################
# load level 3 data
print 'Loading radar level 3 data...'
data_l3 = np.empty([ntime, nazi, ngate])
data_l3[:] = np.nan
f = open('%s/missing_times_list' %output_dir_base, 'w') # year; month; day; hour; minute
times_exist = np.empty(ntime)  # 1 for exist; 0 for missing
for i in range(ntime):
	time = start_time + dt.timedelta(hours=i*3)
	print time
	filename = '%s/level3_3hourly_%s' %(radar_level3_path, time.strftime("%Y%m%d%H%M"))
	if os.path.isfile(filename)==1:   # if the file for this time exists
		times_exist[i] = 1
		nc = Dataset(filename)
		data_l3[i] = nc.variables['Precip3hr'][:] * 25.4  # unit: inch -> mm
		nc.close()
		for azi in range(nazi):
			for gate in range(ngate):
				if np.isnan(data_l3[i,azi,gate])==True:  # if is nan, no precipitation
					data_l3[i,azi,gate] = 0
	else:   # if the file for this time does NOT exist
		times_exist[i] = 0
		f.write('%d %d %d %d %d\n' %(time.year, time.month, time.day, time.hour, time.minute))
f.close()


# interpolation
print 'Interpolating...'
data_interp_0625 = np.empty([ntime, nstn])
for stn in range(nstn):
	print 'Station %d' %(stn+1)
	flag_overlap = -1 

	# calculate weights
	n_radar_points = len(stn_radar_list[stn])
	w = np.empty(n_radar_points)  # weight
	for i in range(n_radar_points):
		if np.abs(stn_radar_list[stn][i][2])<0.01:  # if a radar point overlaps the station
			flag_overlap = i
			break
		else:
			w[i] = 1.0 / (stn_radar_list[stn][i][2] * stn_radar_list[stn][i][2])
	w_sum = np.sum(w)

	# interpolation
	if flag_overlap>=0:
	# if a radar point overlaps the station, directly give its data to the station
		liv_data = np.loadtxt('%s/fluxes_%.5f_%.5f' %(Liv_vic_output_path, liv_match_list[stn,1], liv_match_list[stn,2]))
		for t in range(ntime):
			if times_exist[t]==1:  # if radar data is available at this time step
				data_interp_0625[t, stn] = data_l3[t, stn_radar_list[stn][flat_overlap][0], stn_radar_list[stn][flag_overlap][1]]
			else:  # if radar data is missing at this time step, use Livneh's precip data
				tdelta = time - (Liv_vic_start_time + dt.timedelta(hours=1.5))
				row = (tdelta.days*24 + tdelta.seconds//3600)/3+1  # number of row of this time step in liv_data; start from 1
				data_interp_0625[t, stn] = liv_data[row-1, Liv_vic_prec_column-1]
	else: 
	# if no radar points overlap the station, interpolating
		liv_data = np.loadtxt('%s/fluxes_%.5f_%.5f' %(Liv_vic_output_path, liv_match_list[stn,1], liv_match_list[stn,2]))
		for t in range(ntime):
			time = start_time + dt.timedelta(hours=t*3)
			if times_exist[t]==1:   # if radar data is available at this time step
				prec = 0
				for i in range(n_radar_points):
					prec = prec + w[i] * data_l3[t, stn_radar_list[stn][i][0], stn_radar_list[stn][i][1]]
				data_interp_0625[t, stn] = prec / w_sum
			else:   # if radar data is missing at this time step, use Livneh's precip data
				tdelta = time - (Liv_vic_start_time + dt.timedelta(hours=1.5))
				row = (tdelta.days*24 + tdelta.seconds//3600)/3+1  # number of row of this time step in liv_data; start from 1
				data_interp_0625[t, stn] = liv_data[row-1, Liv_vic_prec_column-1]


####################################################################
###################### Writing results #############################
####################################################################
print 'Writng results...'
for stn in range(nstn):
	f = open('%s/stn.%d_%.5f_%.5f' %(output_dir_prec_result,liv_match_list[stn,0],liv_match_list[stn,1],liv_match_list[stn,2]), 'w')
	for i in range(ntime):
		time = start_time + dt.timedelta(hours=i*3) - dt.timedelta(hours=1.5)  # convert 1:30 to 0:00
		f.write('%d %d %d %d %.6f\n' %(time.year,time.month,time.day,time.hour,data_interp_0625[i,stn]/1000))
f.close()





