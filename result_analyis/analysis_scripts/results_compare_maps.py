#!/usr/local/bin/python

################# compare two 2D maps from model output from 2 runs #################

from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import pdb
import datetime as dt

#################### user defined ####################
result_path_1 = '/usr1/ymao/dhsvm_chehalis/result_analyis/model_results/results_v3.1.2_20140529'
result_path_2 = '/usr1/ymao/dhsvm_chehalis/result_analyis/model_results/results_20140625_v3.1.2_modify_interp'
output_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/output/20140625_compare_before_after_modify_interp_v3.1.2/maps'

map_type = 1; # 1 for SWE; 2 for soil moist; 3 for ET

if map_type==1:  # SWE
	map_name = 'Map.Snow.Swq.asc'
	cbar_title = 'SWE [m]'
	nrow = 918
	ncol = 736
	nmap = 6
	map_title = ['SWE, 04/01/2006', 'SWE, 4/01/2007', 'SWE, 04/01/2008', 'SWE, 04/01/2009', 'SWE, 04/01/2010', 'SWE, 04/01/2011']
	fig_name = 'Map.SWE.v3.1'
	vmin = 0
	vmax = 10
elif map_type==2:  # soil moisture
	map_name = 'Map.1.Soil.Moist.asc'
	cbar_title = 'Soil Moisture [m]'
	nrow = 918
	ncol = 736
	nmap = 12
	map_title = ['Soil moisture, 03/01/2006', 'Soil moisture, 08/01/2006', 'Soil moisture, 03/01/2007', 'Soil moisture, 08/01/2007', 'Soil moisture, 03/01/2008', 'Soil moisture, 08/01/2008', 'Soil moisture, 03/01/2009', 'Soil moisture, 08/01/2009', 'Soil moisture, 03/01/2010', 'Soil moisture, 08/01/2010', 'Soil moisture, 03/01/2011', 'Soil moisture, 08/01/2011']
	fig_name = 'Map.Soil.Moist.v3.1'
	vmin = 0
	vmax = 0.5
elif map_type==3:  # ET
	map_name = 'Map.Evap.ETot.asc'
	cbar_title = 'Total ET [mm/3hours]'
	nrow = 918
	ncol = 736
	nmap = 12
	map_title = ['ET, 03/01/2006', 'ET, 08/01/2006', 'ET, 03/01/2007', 'ET, 08/01/2007', 'ET, 03/01/2008', 'ET, 08/01/2008', 'ET, 03/01/2009', 'ET, 08/01/2009', 'ET, 03/01/2010', 'ET, 08/01/2010', 'ET, 03/01/2011', 'ET, 08/01/2011']
	fig_name = 'Map.tot.ET.v3.1'
	vmin = 0
	vmax = 0.5

################### load data ####################
data_1 = np.loadtxt('%s/%s' %(result_path_1, map_name))
data_2 = np.loadtxt('%s/%s' %(result_path_2, map_name))
map_1 = np.empty([nmap, nrow, ncol])
for i in range(nmap):
	map_1[i] = data_1[i*nrow:(i+1)*nrow, :]
map_2 = np.empty([nmap, nrow, ncol])
for i in range(nmap):
	map_2[i] = data_2[i*nrow:(i+1)*nrow, :]


################### plot map difference #####################
for i in range(nmap):
	fig = plt.figure()
	if map_type==1:  # SWE
		map_mask_1 = np.ma.masked_less(map_1[i], 0.01)  # for SWE
		map_mask_2 = np.ma.masked_less(map_2[i], 0.01)  # for SWE
	elif map_type==2:  # soil moisture
		map_mask_1 = np.ma.masked_less(map_1[i], 0.0)  # for soil moisture
		map_mask_2 = np.ma.masked_less(map_2[i], 0.0)  # for soil moisture
	elif map_type==3:  # ET
		map_mask_1 = np.ma.masked_less(map_1[i]*1000, 0.0000001)  # for ET
		map_mask_2 = np.ma.masked_less(map_2[i]*1000, 0.0000001)  # for ET
	cs = plt.pcolormesh(map_mask_2-map_mask_1, cmap='RdBu', vmax=0.05, vmin=-0.05)
	plt.gca().invert_yaxis()
	plt.axis('equal')
	plt.xlim(1,ncol)
	plt.ylim(nrow, 1)
	cbar = plt.colorbar(cs, orientation='horizontal', fraction=0.05, pad=0.06)
	cbar.set_label('Difference [mm]', fontsize=14)
	plt.title(map_title[i])
	plt.savefig('%s/%s.diff.abs.%d.png' %(output_path, fig_name, i), format='png')


#plt.show()





