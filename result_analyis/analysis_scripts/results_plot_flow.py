#!/usr/local/bin/python

######### plot streamflow results and compare with gauges ###########

from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import pdb
import datetime as dt
from matplotlib.dates import YearLocator

################# user defined ################
result_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/model_results/results_20140712_1990-2013'
label = 'Simulated'
output_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/output/results_20140712_1990-2013'
stream_gauge_path = '/usr1/ymao/dhsvm_chehalis/result_analyis/stream_gauge_data/recent_gauge_daily_data_all_record'
start_time = dt.datetime(year=2011, month=1, day=1, hour=0)  # plot start time
end_time = dt.datetime(year=2013, month=12, day=31, hour=21) # plot end time
sim_start_time = dt.datetime(year=1990, month=1, day=1, hour=0) # simulation start time
sim_end_time = dt.datetime(year=2013, month=12, day=31, hour=21) # simulation end time
dtime = 3  # hours

duration = end_time - start_time
ntime = (duration.days*24 + duration.seconds//3600)/dtime+1
dates = []
for i in range(ntime):
	dates.append(start_time + dt.timedelta(hours=i*dtime))

plot_year_interval = 5

################### Plot streamflow - compare simulated data with stream gauge ################
f = open('%s/Streamflow.Only' %result_path, 'r')
titles = f.readline().rstrip("\n").split()  # date; each gague
agg1 = np.loadtxt('%s/Streamflow.Only' %result_path, skiprows=2, usecols=range(1,len(titles)))

print np.shape(agg1)
print '%d' %ntime

for i in range(len(titles)-1):
	# load stream gague data
	date_gauge = []   # YYYY-MM-DD
	flow = []  # cfs
	if titles[i+1]!='12021800':
		f = open('%s/%s.txt' %(stream_gauge_path, titles[i+1]), 'r')
		for line in f.readlines():
			date_gauge.append(dt.datetime.strptime(line.split()[2], '%Y-%m-%d'))
			if len(line.split())<=3:
				flow.append(np.nan)
			else:
				flow.append(line.split()[3])

	# plot
	fig = plt.figure(figsize=(12,8))
	ax = plt.axes()
	start_lag = start_time - sim_start_time
	start_ind = (start_lag.days*24 + start_lag.seconds//3600)/dtime
	end_lag = end_time - sim_start_time
	end_ind = (end_lag.days*24 + end_lag.seconds//3600)/dtime
	ax.plot_date(dates[1:ntime], agg1[start_ind:end_ind,i]*35.31463782/3600/3, 'r-', label=label) # converted to cfs
	plt.xlim(start_time, end_time)
	years = YearLocator(plot_year_interval)
	ax.xaxis.set_major_locator(years)
	if titles[i+1]!='12021800':
		print '%s' %titles[i+1]
		ax.plot_date(date_gauge, flow, 'k--', label='Gague observation')
	plt.legend()
	plt.title(titles[i+1])
	plt.ylabel('Streamflow (cfs)')
	fig.savefig('%s/Streamflow.Only.2011-2013.%s.png' %(output_path, titles[i+1]), format='png')
#	f = open('%s/%s.txt' %(output_path, titles[i+1]), 'w')
#	for t in range(ntime-1):
#		f.write('%d %d %d %d %f %f\n' %(dates[t+1].year, dates[t+1].month, dates[t+1].day, dates[t+1].hour, agg1[t,i], agg2[t,i]))
#	f.close()


