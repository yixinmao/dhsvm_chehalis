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
result_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/model_results/results_v2.1_20140529'
output_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/output/compare_v2.1_v2.1_20140529/maps_test'

map_type = 3; # 1 for SWE; 2 for soil moist; 3 for ET


if map_type==1:  # SWE
	map_name = 'Map.Snow.Swq.asc'
	cbar_title = 'SWE [m]'
	nrow = 918
	ncol = 736
	nmap = 6
	map_title = ['SWE, 04/01/2006', 'SWE, 4/01/2007', 'SWE, 04/01/2008', 'SWE, 04/01/2009', 'SWE, 04/01/2010', 'SWE, 04/01/2011']
	fig_name = 'Map.SWE.v2.1'
	vmin = 0
	vmax = 10
elif map_type==2:  # soil moisture
	map_name = 'Map.1.Soil.Moist.asc'
	cbar_title = 'Soil Moisture [m]'
	nrow = 918
	ncol = 736
	nmap = 12
	map_title = ['Soil moisture, 03/01/2006', 'Soil moisture, 08/01/2006', 'Soil moisture, 03/01/2007', 'Soil moisture, 08/01/2007', 'Soil moisture, 03/01/2008', 'Soil moisture, 08/01/2008', 'Soil moisture, 03/01/2009', 'Soil moisture, 08/01/2009', 'Soil moisture, 03/01/2010', 'Soil moisture, 08/01/2010', 'Soil moisture, 03/01/2011', 'Soil moisture, 08/01/2011']
	fig_name = 'Map.Soil.Moist.v2.1'
	vmin = 0
	vmax = 0.5
elif map_type==3:  # ET
	map_name = 'Map.Evap.ETot.asc'
	cbar_title = 'Total ET [mm/3hours]'
	nrow = 918
	ncol = 736
	nmap = 12
	map_title = ['ET, 03/01/2006', 'ET, 08/01/2006', 'ET, 03/01/2007', 'ET, 08/01/2007', 'ET, 03/01/2008', 'ET, 08/01/2008', 'ET, 03/01/2009', 'ET, 08/01/2009', 'ET, 03/01/2010', 'ET, 08/01/2010', 'ET, 03/01/2011', 'ET, 08/01/2011']
	fig_name = 'Map.tot.ET.v2.1'
	vmin = 0
	vmax = 0.5

################### load data ####################
data = np.loadtxt('%s/%s' %(result_path, map_name))
map = np.empty([nmap, nrow, ncol])
for i in range(nmap):
	map[i] = data[i*nrow:(i+1)*nrow, :]


################### plot map #####################
for i in range(nmap):
	fig = plt.figure()
	if map_type==1:  # SWE
		map_mask = np.ma.masked_less(map[i], 0.01)  # for SWE
	elif map_type==2:  # soil moisture
		map_mask = np.ma.masked_less(map[i], 0.0)  # for soil moisture
	elif map_type==3:  # ET
		map_mask = np.ma.masked_less(map[i]*1000, 0.0000001)  # for ET
	cs = plt.pcolormesh(map_mask, cmap='jet_r', vmin=vmin, vmax=vmax)
	plt.gca().invert_yaxis()
	plt.axis('equal')
	plt.xlim(1,ncol)
	plt.ylim(nrow, 1)
	cbar = plt.colorbar(cs, orientation='horizontal', fraction=0.045, pad=0.06)
	cbar.set_label(cbar_title, fontsize=14)
	plt.title(map_title[i])
	plt.savefig('%s/%s.%d.png' %(output_path, fig_name, i), format='png')


#plt.show()





