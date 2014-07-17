#!/usr/local/bin/python

######### compare results from two runs (same running time period) ###########


from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import pdb
import datetime as dt
from matplotlib.dates import YearLocator

################# user defined ################
result_path_1 = '/usr1/ymao/dhsvm_chehalis/result_analyis/model_results/results_20140709_1990-2011'
result_path_2 = '/usr1/ymao/dhsvm_chehalis/result_analyis/model_results/results_20140712_1990-2013'
label1 = 'Old Livneh data'
label2 = 'Updated Livneh data'
output_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/output/compare_20140709_20140712'
stream_gauge_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/stream_gauge_data/recent_gauge_daily_data_2006_2011'
sim_start_time_1 = dt.datetime(year=1990, month=1, day=1, hour=0)
sim_end_time_1 = dt.datetime(year=2011, month=12, day=31, hour=21)
sim_start_time_2 = dt.datetime(year=1990, month=1, day=1, hour=0)
sim_end_time_2 = dt.datetime(year=2013, month=12, day=31, hour=21)
start_time = dt.datetime(year=1990, month=1, day=1, hour=0)  # plot start time
end_time = dt.datetime(year=2011, month=12, day=31, hour=21)  # plot end time
dtime = 3  # hours

duration = end_time - start_time
ntime = (duration.days*24 + duration.seconds//3600)/dtime+1
dates = []
for i in range(ntime):
	dates.append(start_time + dt.timedelta(hours=i*dtime))

plot_year_interval = 5

################## Aggregated.Values ################
#f = open('%s/Aggregated.Values' %result_path_1, 'r')
#titles1 = f.readline().rstrip("\n").split()
#f.close()
#f = open('%s/Aggregated.Values' %result_path_2, 'r')
#titles2 = f.readline().rstrip("\n").split()
#f.close()
#
#agg1 = np.loadtxt('%s/Aggregated.Values' %result_path_1, skiprows=1, usecols=range(1,len(titles1)))
#agg2 = np.loadtxt('%s/Aggregated.Values' %result_path_2, skiprows=1, usecols=range(1,len(titles2)))
#
#print "Titles in file 1 not existing in file 2:"
#for i in range(len(titles1)-1):  # for every column in file1, find the column with the same title in file 2 -> plot; i in terms of agg
#	if titles1[i+1]=="Precip":   # in v2.1, Precip; in v3.1.2, Precip(m)
#		j = titles2.index("Precip(m)") - 1
#	elif titles1[i+1]=="EvapTot":  
#		j = titles2.index("TotEvap") - 1
#	elif titles1[i+1]=="SoilMoist0":  
#		j = titles2.index("SoilMoist1") - 1
#	elif titles1[i+1]=="SoilMoist1":  
#		j = titles2.index("SoilMoist2") - 1
#	elif titles1[i+1]=="SoilMoist2":  
#		j = titles2.index("SoilMoist3") - 1
#	elif titles1[i+1]=="Perc0":  
#		j = titles2.index("Perc1") - 1
#	elif titles1[i+1]=="Perc1":  
#		j = titles2.index("Perc2") - 1
#	elif titles1[i+1]=="Perc2":  
#		j = titles2.index("Perc3") - 1
#	elif titles1[i+1]=="OverSnow":  
#		j = titles2.index("SnowCover") - 1
#	elif (titles1[i+1] in titles2) == False:   # if this title in file1 does NOT exist in file2
#		print titles1[i+1]
#		continue
#	else:  # if this title in file 1 exists in file 1
#		j = titles2.index(titles1[i+1]) - 1
#	fig = plt.figure()
#	ax = plt.axes()
#	ax.plot_date(dates, agg1[0:np.shape(agg1)[0]-1,i], 'b-', label=label1)
#	ax.plot_date(dates, agg2[0:np.shape(agg2)[0]-1,j], 'r-', label=label2)
#	plt.legend()
#	plt.title(titles2[j+1])
#	fig.savefig('%s/Aggregated.Values.%s.png' %(output_path, titles2[j+1]), format='png')
#
#print "\nTitles in file 2 not existing in file 1:"
#for i in range(len(titles2)-1): 
#	if titles2[i+1]=="Precip(m)":   
#		continue
#		j = titles2.index("TotEvap") - 1
#	elif (titles2[i+1] in titles1) == False:  
#		print titles2[i+1]
#		continue


################## Mass.Balance ################
#titles = ['Total runoff', 'Total amount of water in the canopy', 'Total amount of water in the soil', 'Total amount of SWE', 'Total amount of sat. subsurface flow', 'Total amount of water intercepted by channels', 'Total amount of water intercepted by roads', 'Total water returned by culverts to the land surface', 'Total amount of ET', 'Total amount of precipitation', 'Total sublimation from snow on the ground', 'Total sublimation from snow in the canopy', 'Total water during the previous time step', 'Total flow from culverts to the channel', 'Total surface flow to the channel', 'Total mass balance error for the current time step'] 
#agg1 = np.loadtxt('%s/Mass.Balance' %result_path_1, skiprows=0, usecols=range(1,17))
#agg2 = np.loadtxt('%s/Mass.Balance' %result_path_2, skiprows=0, usecols=range(1,17))
#
#for i in range(16):
#	fig = plt.figure()
#	ax = plt.axes()
#	ax.plot_date(dates, agg1[:,i], 'b-', label=label1)
#	ax.plot_date(dates, agg2[:,i], 'r-', label=label2)
#	plt.legend()
#	plt.title(titles[i])
#	fig.savefig('%s/Mass.Balance.%d.png' %(output_path, i+1), format='png')

################### Streamflow - compare both data with stream gauge ################
#f = open('%s/Streamflow.Only' %result_path_1, 'r')
#titles = f.readline().rstrip("\n").split()  # date; each gague
#agg1 = np.loadtxt('%s/Streamflow.Only' %result_path_1, skiprows=2, usecols=range(1,len(titles)))
#agg2 = np.loadtxt('%s/Streamflow.Only' %result_path_2, skiprows=2, usecols=range(1,len(titles)))   # m3/s
#
#
#for i in range(len(titles)-1):
#	# load stream gague data
#	date_gauge = []   # YYYY-MM-DD
#	flow = []  # cfs
#	if titles[i+1]!='12019310' and titles[i+1]!='12021800':
#		f = open('%s/%s.txt' %(stream_gauge_path, titles[i+1]), 'r')
#		for line in f.readlines():
#			date_gauge.append(dt.datetime.strptime(line.split()[2], '%Y-%m-%d'))
#			if len(line.split())<=3:
#				flow.append(np.nan)
#			else:
#				flow.append(line.split()[3])
#
#	# plot
#	fig = plt.figure()
#	ax = plt.axes()
#	ax.plot_date(dates[1:ntime], agg1[:,i]*35.31463782/3600/3, 'b-', label=label1) # converted to cfs
#	ax.plot_date(dates[1:ntime], agg2[:,i]*35.31463782/3600/3, 'r-', label=label2) # converted to cfs
#	plt.xlim(dt.datetime(year=2009,month=1,day=1,hour=0), end_time)
#	years = YearLocator()
#	ax.xaxis.set_major_locator(years)
#	if titles[i+1]!='12019310' and titles[i+1]!='12021800':
#		ax.plot_date(date_gauge, flow, 'k-', label='Gague observation')
#	plt.legend()
#	plt.title(titles[i+1])
#	plt.ylabel('Streamflow (cfs)')
#	fig.savefig('%s/Streamflow.Only.%s.png' %(output_path, titles[i+1]), format='png')
##	f = open('%s/%s.txt' %(output_path, titles[i+1]), 'w')
##	for t in range(ntime-1):
##		f.write('%d %d %d %d %f %f\n' %(dates[t+1].year, dates[t+1].month, dates[t+1].day, dates[t+1].hour, agg1[t,i], agg2[t,i]))
##	f.close()


################## Streamflow - compare 2 runs ################
f = open('%s/Streamflow.Only' %result_path_1, 'r')
titles = f.readline().rstrip("\n").split()  # date; each gague
agg1 = np.loadtxt('%s/Streamflow.Only' %result_path_1, skiprows=2, usecols=range(1,len(titles)))
agg2 = np.loadtxt('%s/Streamflow.Only' %result_path_2, skiprows=2, usecols=range(1,len(titles)))   # m3/s


for i in range(len(titles)-1):
#	# load stream gague data
#	date_gauge = []   # YYYY-MM-DD
#	flow = []  # cfs
#	if titles[i+1]!='12019310' and titles[i+1]!='12021800':
#		f = open('%s/%s.txt' %(stream_gauge_path, titles[i+1]), 'r')
#		for line in f.readlines():
#			date_gauge.append(dt.datetime.strptime(line.split()[2], '%Y-%m-%d'))
#			if len(line.split())<=3:
#				flow.append(np.nan)
#			else:
#				flow.append(line.split()[3])

	# plot
	fig = plt.figure()
	ax = plt.axes()
	start_lag_1 = start_time - sim_start_time_1
	start_ind_1 = (start_lag_1.days*24 + start_lag_1.seconds//3600)/dtime
	end_lag_1 = end_time - sim_start_time_1
	end_ind_1 = (end_lag_1.days*24 + end_lag_1.seconds//3600)/dtime
	start_lag_2 = start_time - sim_start_time_2
	start_ind_2 = (start_lag_2.days*24 + start_lag_2.seconds//3600)/dtime
	end_lag_2 = end_time - sim_start_time_2
	end_ind_2 = (end_lag_2.days*24 + end_lag_2.seconds//3600)/dtime

	ax.plot_date(dates[1:ntime], agg1[start_ind_1:end_ind_1,i]*35.31463782/3600/3, 'b-', label=label1) # converted to cfs
	ax.plot_date(dates[1:ntime], agg2[start_ind_2:end_ind_2,i]*35.31463782/3600/3, 'r-', label=label2) # converted to cfs
	ax.plot_date(dates[1:ntime], (agg2[start_ind_1:end_ind_1,i]-agg1[start_ind_2:end_ind_2,i])*35.31463782/3600/3, 'k-', label='Difference') # converted to cfs
	plt.xlim(start_time, end_time)
	years = YearLocator(plot_year_interval)
	ax.xaxis.set_major_locator(years)
#	if titles[i+1]!='12019310' and titles[i+1]!='12021800':
#		ax.plot_date(date_gauge, flow, 'k-', label='Gague observation')
	print '%s' %titles[i+1]
	plt.legend()
	plt.title(titles[i+1])
	plt.ylabel('Streamflow (cfs)')
	fig.savefig('%s/Streamflow.Only.diff.%s.png' %(output_path, titles[i+1]), format='png')


