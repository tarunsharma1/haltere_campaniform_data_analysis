import copy
import matplotlib.pyplot as plt
import random
import numpy as np
import pickle
import glob
import scipy.interpolate as spi
import scipy.fftpack
import scipy.signal as sps
from scipy.optimize import leastsq
import seaborn as sns
import matplotlib as mpl
from scipy import stats
import will_amplitude_and_phase_calc
from scipy.integrate import cumulative_trapezoid
from sklearn.metrics import r2_score, mean_squared_error
from will_amplitude_and_phase_calc import estimate_phase, estimate_amplitude
import pandas as pd


########################
# this script reads the saved pickle files for flies on a genotype (one per fly averaged over 3 trials), 
# resamples them to a common time stamp, averages and then plots everything on one plot i.e
# plots the average and also each flies trace in lighter colors. It does this for all the speeds
# i.e new plot per rotation speed. The trace for each fly represents its left - right baseline subtracted WBA.

########################


def sign_with_positive_zero(arr):
    return np.sign(arr) + (arr == 0)


###################################################################################################
# Adjust Spines (Dickinson style, thanks to Andrew Straw)
###################################################################################################

# NOTE: smart_bounds is disabled (commented out) in this function. It only works in matplotlib v >1.
# to fix this issue, try manually setting your tick marks (see example below) 
def adjust_spines(ax,spines, spine_locations={}, smart_bounds=True, xticks=None, yticks=None, linewidth=1, spineColor='black'): # ivo: spineColor
    if type(spines) is not list:
        spines = [spines]
        
    # get ticks
    if xticks is None:
        xticks = ax.get_xticks()
    if yticks is None:
        yticks = ax.get_yticks()
        
	#spine_locations_dict = {'top': 10, 'right': 10, 'left': 10, 'bottom': 10}
    spine_locations_dict = {'top': 6, 'right': 6, 'left': 6, 'bottom': 6}
    for key in spine_locations.keys():
        spine_locations_dict[key] = spine_locations[key]
        
    if 'none' in spines:
        for loc, spine in ax.spines.items():
        #for loc, spine in ax.spines.iteritems(): #ivo:  was used in python 2.7
            spine.set_color('none') # don't draw spine
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])
        return
    
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward',spine_locations_dict[loc])) # outward by x points
            spine.set_linewidth(linewidth)
            spine.set_color(spineColor) #ivo
            ax.tick_params(colors=spineColor) #ivo
            ax.tick_params(length=linewidth*4) #ivo
            ax.tick_params(pad=linewidth*4) #ivo
            ax.tick_params(direction="in") #ivo
        else:
            spine.set_color('none') # don't draw spine
            
    # smart bounds, if possible
    if int(mpl.__version__[0]) > 0 and smart_bounds: 
        for loc, spine in ax.spines.items():
            if loc in ['left', 'right']:
                ticks = yticks
            if loc in ['top', 'bottom']:
                ticks = xticks
            spine.set_bounds(ticks[0], ticks[-1])

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    if 'top' in spines:
        ax.xaxis.set_ticks_position('top')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])    
    
    if 'left' in spines or 'right' in spines:
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, rotation=0)
    if 'top' in spines or 'bottom' in spines:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, rotation=0)
    	
    for line in ax.get_xticklines() + ax.get_yticklines():
        #line.set_markersize(6)
        line.set_markeredgewidth(linewidth)


def resample(x1, y1, x2, kind, extrapolate=True):
	# helper function
	if kind == 'spline':
		spline = spi.CubicSpline(x1, y1)
		y2 = spline.__call__(x2, extrapolate=extrapolate)
	else:
		fill_value = "extrapolate" if extrapolate else []
		#interp = spi.interp1d(x1, y1, kind=kind, bounds_error=False, fill_value=fill_value)
		interp = spi.interp1d(x1, y1, kind=kind, bounds_error=True)
		
		y2 = interp(x2)
	return y2


def fourier_transform_to_get_cutoff(signal, t):

	### one alternative suggested by email to reduce noise was to take the second half of the signal,invert it as the first half and use this as the cleaned signal
	# N = signal.shape[0]
	# signal_second_half = signal[int(N/2)::]
	# ## move the signal down to 0 to fix asymmetry
	# signal_second_half = signal_second_half - signal_second_half[0]
	# signal_new = np.zeros_like(signal)
	# signal_new[0:(int(N%2) + int(N/2))] = -1.0*signal_second_half
	# signal_new[int(N/2)::] = signal_second_half

	# signal = signal_new



	# helper function for me to visualize and see 
	# fourier transform to freq domain and look at signal to see what cutoff value to use for low pass filter

	N = signal.shape[0]
	dt = t[2] - t[1]

	#print ('sampling freq:', 1.0/dt)

	fft = 1.0/N * np.fft.fft(signal)
	fft = fft[:N//2]
	fftfreq = np.fft.fftfreq(N, dt)[:N//2] 
	
	#plt.plot(fftfreq, np.abs(fft))
	#plt.show()
	
	# the first highest fft value is going to be a line ..look at fun with fft if you set the signal as 5 + sinx biggest spike is a line then is the one for the sine curve
	# so we can safely ignore the first fft value.. and we know its going to be at 0 freq i.e at fft[0]

	idx = 1 + np.argmax(np.abs(fft[1:]))
	#print ('index:', idx)

	amp = np.abs(fft[idx]) * 2
	#print ('amp:', amp)

	phase =  np.angle(fft[idx]) * 180/np.pi
	#print ('phase in degrees:',phase)

	freq = np.abs(fftfreq[idx])

	## add offset np.mean(signal) just for visualizing the fits better
	#phase = 90
	signal_reconstucted =  amp*np.cos(2*np.pi*freq*t + phase*np.pi/180) + np.mean(signal)


	### lets try doing an optimization here using scipy to increase fit quality
	optimize_func = lambda x: x[0]*np.cos(2*np.pi*freq*t+x[1]*np.pi/180) + np.mean(signal) - signal
	est_amp, est_phase, est_mean = leastsq(optimize_func, [amp, phase, np.mean(signal)])[0]
	#print ('##### scipy optimized #####')
	#print ('amp:', est_amp)
	#print ('phase', est_phase)
	#print ('frequency:', est_freq)
	optimized_signal_reconstructed =  est_amp*np.cos(2*np.pi*freq*t + est_phase*np.pi/180) + np.mean(signal)

	plt.clf()
	#plt.ylim(-8,8)
	return optimized_signal_reconstructed, est_amp, est_phase, signal
	#return signal_reconstucted, amp, phase, signal

	######################## Will's method of amp and phase estimation #################################
	# num_cycle = 5
	# if t[-1] - t[0] < 100: ### speed 2
	# 	frequency = 1/18.01
	# else:
	# 	frequency = 1/25.84
	# phase_rad_est = estimate_phase(signal, dt, frequency, num_cycle, t)
	# phase_deg_est = np.rad2deg(phase_rad_est)
	# phase_deg_est = (phase_deg_est - 90)%360
	# amp_estimate = estimate_amplitude(signal, dt, frequency, num_cycle, t)
	# signal_reconstucted =  amp_estimate*np.cos(2*np.pi*frequency*t + phase_deg_est*np.pi/180) + np.mean(signal)
	# return signal_reconstucted, amp_estimate, phase_deg_est, signal


	
def calculate_average_over_cycles_per_fly(average_left_minus_right_og, t, resampled_trajectories, speed, head, store_baselines=False, left=False, right=False):
	### calculate average over 5 sinusoid cycles per fly.. they are already resampled to have the same t
	list_of_amps_per_fly = []
	list_of_phases_per_fly = []
	list_of_signals_per_fly = []
	list_of_signals_fitted_per_fly = []
	list_of_optimized_signals_fitted_per_fly = []
	list_of_vertical_offsets_per_fly = []


	num_cycles = 5
	# t1 = int(round(t.shape[0]/float(num_cycles)))
	# resampling_dt = np.mean(np.diff(t[0:t1]))
	# t2 = np.arange(len(t[0:t1])) * resampling_dt
	
	
	plt.clf()

	average_left_minus_right = copy.deepcopy(average_left_minus_right_og)
	

	for fly_number,WBA in enumerate(average_left_minus_right):
		### this loop is per fly.. I am fitting a sinusoid and calculating the amp and phase per cycle and then taking the mean in the end as the value per fly

		### going to use zero crossing points from the mean trajectory to split up into 5 cycles
		zero_crossing_idxs = np.where(np.abs(np.diff(np.sign(mean_trajectory)))==2)[0]
		## add 0 to the start 
		zero_crossing_idxs = np.insert(zero_crossing_idxs, 0,0)

		### sometimes due to noise there maybe be consequtive indices like [0, 900, 1800,1801, 1802, 1803, 2700...]. We want to group them and take maximum of each group
		split_indices = np.where(np.diff(zero_crossing_idxs) > 100)[0] + 1
		groups = np.split(zero_crossing_idxs, split_indices)
		zero_crossing_idxs = [group.max() for group in groups]


		trajectory_cycles = []
		
		cycles = []
		for x in range(0, (num_cycles-1)*2, 2):
			#temp.append(resample(t[x*t1:(x+1)*t1] - t[x*t1] , WBA[x*t1:(x+1)*t1], t2, kind='linear'))
			cycles.append(WBA[zero_crossing_idxs[x]:zero_crossing_idxs[x+2]])
			trajectory_cycles.append(resampled_trajectories[fly_number][zero_crossing_idxs[x]:zero_crossing_idxs[x+2]])
		## add last cycle
		cycles.append(WBA[zero_crossing_idxs[x+2]::])
		trajectory_cycles.append(resampled_trajectories[fly_number][zero_crossing_idxs[x+2]::])


		### now I just take the minimum cycle and cut off everything to that
		min_cycle_len = min(trajectory_cycles[x].shape[0] for x in range(num_cycles))
		

		for c in range(num_cycles):
			cycles[c] = cycles[c][0:min_cycle_len]
			trajectory_cycles[c] = trajectory_cycles[c][0:min_cycle_len]
			
		mean_trajectory_cycle = np.mean(trajectory_cycles, axis=0)
		mean_WBA_cycle = np.mean(cycles, axis=0)

		if left or right:
			## baseline subtract the first value so they start at 0
			mean_WBA_cycle = mean_WBA_cycle - mean_WBA_cycle[0]

		fitted_curve, amp, phase,_ = fourier_transform_to_get_cutoff(mean_WBA_cycle, np.arange(0, min_cycle_len)* np.diff(t)[0])
		#print (f'fourier phase {phase}')
		_,a, phase_trajectory,_ = fourier_transform_to_get_cutoff(mean_trajectory_cycle, np.arange(0,min_cycle_len) * np.diff(t)[0])
		
		
			
		print ('phases per cycle for fly :' + str(fly_number) + ': ')
		print (phase_trajectory - phase )
		
		
		phase = phase - phase_trajectory
		average_over_sinusoid_cycles = mean_WBA_cycle

		optimized_fitted_sin_curve = fitted_curve
		t2 = np.arange(0, min_cycle_len)* np.diff(t)[0]
		
		
		list_of_amps_per_fly.append(amp)
		list_of_phases_per_fly.append(phase)
		
		list_of_signals_per_fly.append(average_over_sinusoid_cycles)
		list_of_optimized_signals_fitted_per_fly.append(optimized_fitted_sin_curve)
		list_of_vertical_offsets_per_fly.append(np.mean(average_over_sinusoid_cycles))

	
	### now plot everything
	# plt.figure(figsize=(5,5))
	
	# for i in range(0,len(list_of_amps_per_fly)):
	# 	plt.subplot(5,5,i+1)
	# 	if head==0:
	# 		plt.ylim(-40,40)
	# 	if head==2:
	# 		plt.ylim(-0.1,0.1)
	# 	else:
	# 		plt.ylim(-20,20)
	# 	plt.title('amp: ' + str(round(list_of_amps_per_fly[i],1)) + ', phase:' + str(round(list_of_phases_per_fly[i],1)))
	# 	#plt.title('phase:' + str(round(list_of_phases_per_fly[i],1)) + ' rmse: ' + str(round(np.sqrt(mean_squared_error(list_of_signals_per_fly[i], list_of_optimized_signals_fitted_per_fly[i])), 2)))
	# 	if head==0:
	# 		plt.ylabel('L-R WBA (degrees)')
	# 	elif head==1:
	# 		plt.ylabel('Head roll angle (degrees)')
	# 	elif head==2:
	# 		plt.ylabel('Head pitch (pixels)')
	# 	elif head==3:
	# 		plt.ylabel('Head yaw angle (degrees)')

	# 	plt.plot(t2, list_of_signals_per_fly[i] , 'g', alpha=0.5)
		
	# 	plt.plot(t2, list_of_optimized_signals_fitted_per_fly[i], 'r', alpha=0.5)
	# 	plt.plot(t2, 0.02* resampled_trajectories[i][0:t2.shape[0]], 'b', alpha=0.5)
	# 	plt.xticks([])
	
	
	# ##plt.title('fit for fly :' + str(1+fly_number) + ' amplitude :' + str(amp))
	# ##plt.title('all fits')
	# plt.show()



	if head:
		## save to pickle file for plotting all genotypes on one plot
		if head==1:
			pickle.dump(list_of_amps_per_fly, open('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_scatter_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'-head-roll.p', 'wb'))
		elif head==2:
			pickle.dump(list_of_amps_per_fly, open('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_scatter_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'-head-pitch.p', 'wb'))
		elif head==3:
			print ('####### saving head ################')
			pickle.dump(list_of_amps_per_fly, open('./data_for_scatter_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'-head-yaw.p', 'wb'))
			pickle.dump(list_of_phases_per_fly, open('./data_for_scatter_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'-head-yaw-phases.p', 'wb'))

	else:
		## save to pickle file for plotting all genotypes on one plot
		if left:
			pickle.dump(list_of_amps_per_fly, open('./data_for_scatter_plot_new/'+subdir[:-1]+'-left-speed-'+str(speed+1)+'.p', 'wb'))

		elif right:
			pickle.dump(list_of_amps_per_fly, open('./data_for_scatter_plot_new/'+subdir[:-1]+'-right-speed-'+str(speed+1)+'.p', 'wb'))

		else:
			print ('####### saving wings ################')
			pickle.dump(list_of_amps_per_fly, open('./data_for_scatter_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'.p', 'wb'))
			pickle.dump(list_of_phases_per_fly, open('./data_for_scatter_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'-phases.p', 'wb'))
			
			if store_baselines:
				pickle.dump([list_of_amps_per_fly, avg_baseline_left[str(speed)], avg_baseline_right[str(speed)]], open('./data_for_baseline_correlation_plot_new/'+subdir[:-1]+'-speed-'+str(speed+1)+'.p', 'wb'))
				#pickle.dump([list_of_vertical_offsets_per_fly, mean_baseline_left[str(speed)], mean_baseline_right[str(speed)]], open('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_baseline_correlation_plot/'+subdir.split('-')[0]+'-'+subdir.split('-')[1]+'-speed-'+str(speed+1)+'-offsets.p', 'wb'))
				
			if show_plots:
				plt.scatter([0]*len(list_of_amps_per_fly),list_of_amps_per_fly)
				plt.ylim(0,30)
				plt.show()

	return list_of_signals_per_fly


## flag to indicate whether head tracking has been done already on this genotype or not
head_data = False
show_plots = False

parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8, 'svg.fonttype':'none'}
plt.rcParams.update(parameters)

subdir = 'UXJ96yawpos/'
all_flies = glob.glob('../hdf5/'+subdir+'/*')
pickle_files = []

for i in all_flies:
	pickle_files.append(glob.glob('../hdf5/'+subdir+i.split('/')[-1]+'/avg_left_minus_right*.p'))
	

#pickle_files = [pickle_files[6]]

print (pickle_files)

#pickle_files = [['/home/tarun/catkin_ws/src/trajectories_autostep_ros/bagfiles/hdf5/KirJ00yawpos/Nov-1-fly9/avg_left_minus_right_Nov-1-fly9.p']]

num_speeds = 2


avg_baseline_left = {}
avg_baseline_right = {}
avg_baseline_head = {}
avg_baseline_left_trace = {}
avg_baseline_right_trace = {}
avg_baseline_head_trace = {}
avg_baseline_time = {}
avg_baseline_trajectory = {}

for i in range(1, num_speeds):
	trajectories = []
	trajectories_time = []

	avg_baseline_left[str(i)] = []
	avg_baseline_right[str(i)] = []
	avg_baseline_left_trace[str(i)] = []
	avg_baseline_right_trace[str(i)] = []
	avg_baseline_head_trace[str(i)] = []
	avg_baseline_time[str(i)] = []
	avg_baseline_trajectory[str(i)] = []
	avg_baseline_head[str(i)] = []

	# new plot for each speed
	if show_plots:
		figures, axes = plt.subplots()

	# will be used for resampling
	time_and_WBA_for_all_flies = []
	
	fig, axes = plt.subplots()
	fig.set_size_inches(2.2,0.7)
	
	for f in pickle_files:
		# skip empty entries in list
		print (f)

		if f==[]:
			continue
		dict_for_pickle = pickle.load(open(f[0], 'rb'))

		## each key for the dict is the speed..the value is a list of [t, WBA]
		## ok so for a speed, the time for each fly starts at 0 and ends at the exact same (to first decimal) value
		## the number of points in time (and hence WBA) are different though fly to fly..so will have to resample

		time = dict_for_pickle[i][0]
		WBA = dict_for_pickle[i][1]
		head_yaw = dict_for_pickle[i][2]

		left_WBA = dict_for_pickle['mean_left_'+str(i)]
		right_WBA = dict_for_pickle['mean_right_'+str(i)]

		#WBA = WBA - np.mean(WBA)
		#left_WBA = left_WBA - np.mean(left_WBA)

		if head_data:
			head_roll = dict_for_pickle[i][2]
			head_pitch = dict_for_pickle[i][3]
			head_yaw = dict_for_pickle[i][4]
		
		if 'trajectory_1' and 'trajectory_1_time' in dict_for_pickle.keys():	
			trajectories.append(dict_for_pickle['trajectory_' + str(i)][0])
			trajectories_time.append(dict_for_pickle['trajectory_' + str(i) + '_time'][0])

		if 'mean_baseline_left_1' in dict_for_pickle.keys():

			avg_baseline_left_trace[str(i)].append(dict_for_pickle['mean_baseline_left_trace_'+str(i)])### in rkf_analysis, we resampled all baselines across all flies to common timescale manually.
			avg_baseline_right_trace[str(i)].append(dict_for_pickle['mean_baseline_right_trace_'+str(i)])
			avg_baseline_head_trace[str(i)].append(dict_for_pickle['mean_baseline_head_trace_'+str(i)])
			avg_baseline_time[str(i)].append(dict_for_pickle['baseline_time_'+str(i)])
			avg_baseline_trajectory[str(i)].append(dict_for_pickle['mean_baseline_trajectory_' + str(i)])

		if head_data:
			time_and_WBA_for_all_flies.append([time, WBA, left_WBA, right_WBA, head_roll, head_pitch, head_yaw])
		else:
			#if abs(dict_for_pickle['mean_baseline_left_'+str(i)] - dict_for_pickle['mean_baseline_right_'+str(i)]) < 10:
			time_and_WBA_for_all_flies.append([time, WBA, left_WBA, right_WBA, head_yaw])

		# first lets just plot all of them on one plot
		
		#axes.set_ylim(-40,40)
		if show_plots:
			r = random.random()
			b = random.random()
			g = random.random()
			## setting the color to the same number for all r,g,b makes it a shade of gray
			color = (0.5, 0.5, 0.5)

			#axes.set_xlabel('Time (seconds)')
			#axes.set_ylabel('WBA yaw (degrees)')
			baseline_time =  np.linspace(-15,0,len(avg_baseline_right_trace[str(i)][-1]))
			
			## wings
			#axes.plot(np.concatenate([baseline_time,time]), np.concatenate([avg_baseline_left_trace[str(i)][-1] - avg_baseline_right_trace[str(i)][-1], WBA]) , c=color, linewidth=0.25)
			## head
			#axes.plot(np.concatenate([baseline_time,time]), np.concatenate([avg_baseline_head_trace[str(i)][-1], head_yaw]) , c=color, linewidth=0.25 )


	
	### resample to a common timestamp across flies for baseline and also for stimulus motion - baseline first ##

	# t_baseline = np.arange(0, min([avg_baseline_time[str(i)][x][-1] for x in range(len(time_and_WBA_for_all_flies))]), 1/72.0)
	#for k,fly in enumerate(time_and_WBA_for_all_flies):
	# 	avg_baseline_left_trace[str(i)][k] = resample(avg_baseline_time[str(i)][k], avg_baseline_left_trace[str(i)][k], t_baseline, kind='linear')
	# 	avg_baseline_right_trace[str(i)][k] = resample(avg_baseline_time[str(i)][k], avg_baseline_right_trace[str(i)][k], t_baseline, kind='linear')
	# 	avg_baseline_head_trace[str(i)][k] = resample(avg_baseline_time[str(i)][k], avg_baseline_head_trace[str(i)][k], t_baseline, kind='linear')
		
	# 	avg_baseline_trajectory[str(i)][k] = resample(avg_baseline_time[str(i)][k], avg_baseline_trajectory[str(i)][k], t_baseline, kind='linear')
	# 	avg_baseline_time[str(i)][k] = t_baseline



	# mean_baseline_left_trace = np.mean(avg_baseline_left_trace[str(i)], axis=0)
	# mean_baseline_right_trace = np.mean(avg_baseline_right_trace[str(i)], axis=0)
	# mean_baseline_head_trace = np.mean(avg_baseline_head_trace[str(i)], axis=0)
	# mean_baseline_trajectory = np.mean(avg_baseline_trajectory[str(i)], axis=0)

	#### we don't actually care about the baseline trace though just the mean values of the baseline responses
	for k in range (0, len(time_and_WBA_for_all_flies)):
		avg_baseline_left[str(i)].append(np.mean(avg_baseline_left_trace[str(i)][k][900::]))
		avg_baseline_right[str(i)].append(np.mean(avg_baseline_right_trace[str(i)][k][900::]))
		avg_baseline_head[str(i)].append(np.mean(avg_baseline_head_trace[str(i)][k][900::]))




	########## resampling across flies for stimulus motion ###
	
	t = np.arange(0, min([time_and_WBA_for_all_flies[x][0][-1] for x in range(len(time_and_WBA_for_all_flies))]), 1/72.0)


	average_left_minus_right = []
	average_head_roll = []
	average_head_pitch = []
	average_head_yaw = []
	average_left = []
	average_right = []
	resampled_trajectories = []



	for k,fly in enumerate(time_and_WBA_for_all_flies):
		WBA = resample(fly[0], fly[1], t, kind='linear')
		left_WBA = resample(fly[0], fly[2], t, kind='linear')
		right_WBA = resample(fly[0], fly[3], t, kind='linear')
		yaw = resample(fly[0], fly[4], t, kind='linear')

		resampled_trajectories.append(resample(trajectories_time[k], trajectories[k], t, kind='linear'))

		### mean subtract every trace so that they all start roughly at 0 offset on the Y axis
		#average_left_minus_right.append(WBA - np.mean(WBA))
		average_left_minus_right.append(WBA)
		average_left.append(left_WBA)
		average_right.append(right_WBA)
		average_head_yaw.append(yaw)

		if head_data:
			### mean subtract every trace so that they all start roughly at 0 offset on the Y axis
			roll = resample(fly[0], fly[4], t, kind='linear')
			pitch = resample(fly[0], fly[5], t, kind='linear')
			yaw = resample(fly[0], fly[6], t, kind='linear')
			
			average_head_roll[k] = roll
			average_head_pitch[k] = pitch
			average_head_yaw[k] = yaw


		
	std = np.std(average_left_minus_right, axis = 0)
	mean = np.mean(average_left_minus_right, axis = 0)
	mean_left = np.mean(average_left, axis=0)
	mean_right = np.mean(average_right, axis=0)
	mean_yaw = np.mean(average_head_yaw, axis = 0)

	mean_trajectory = np.mean(resampled_trajectories, axis=0)

	# angular_position = cumulative_trapezoid(mean_trajectory, t, initial=0)
	
	# mean_norm = (mean_yaw - np.mean(mean_yaw)) / np.std(mean_yaw)
	# trajectory_norm = (angular_position - np.mean(angular_position)) / np.std(angular_position)

	# plt.plot(t, mean_norm, 'g')
	# plt.plot(t, trajectory_norm, 'b')
	# plt.show()
	

	### now that we have the average trajectory and behavioral data across all flies on the same time vector, 
	### we remove the divets (delays) at the start and at the end by determining the indices at the start and end when angular velocity is greater than 1 deg/sec.
	thresholded_start_index = np.where(mean_trajectory <= -1)[0][0] - 1
	thresholded_end_index = np.where(mean_trajectory >= 1)[0][-1] + 1

	# optimized_fitted_sin_curve, amp, phase, average_over_sinusoid_cycles  = fourier_transform_to_get_cutoff(mean,t)
	# print (f'phase of mean response {phase}')
	# optimized_fitted_sin_curve, amp, phase, average_over_sinusoid_cycles  = fourier_transform_to_get_cutoff(mean_trajectory,t)
	# print (f'phase of mean trajectories {phase}')

	
	# fig, ax1 = plt.subplots()
	# ax1.plot(t, mean_left, 'g')
	# ax1.set_ylim(-40,40)
	# adjust_spines(ax1,['bottom', 'left'],xticks=[0,20],  yticks=[-20,20], linewidth=0.6,spineColor='k')

	# ax2 = ax1.twinx()
	# ax2.plot(t, mean_right, 'r')
	# ax2.set_ylim(-40,40)
	# adjust_spines(ax2,['bottom', 'left'],xticks=[0,20],  yticks=[-20,20], linewidth=0.6,spineColor='k')
	
	# # ax3 = ax1.twinx()
	# # ax3.plot(t, mean, 'k')
	# # ax3.set_ylim(-50,50)

	# ax4 = ax1.twinx()
	# ax4.plot(t, mean_trajectory, 'b')
	# ax4.set_ylim(-530,530)
	# adjust_spines(ax4,['bottom', 'right'],xticks=[0,20],  yticks=[-200,200], linewidth=0.6,spineColor='k')

	# plt.show()

	mean_trajectory = mean_trajectory[thresholded_start_index:thresholded_end_index]
	t = t[thresholded_start_index:thresholded_end_index] - t[thresholded_start_index]
	mean_left = mean_left[thresholded_start_index:thresholded_end_index]
	mean_right = mean_right[thresholded_start_index:thresholded_end_index]
	mean = mean[thresholded_start_index:thresholded_end_index]
	mean_yaw = mean_yaw[thresholded_start_index:thresholded_end_index]

	for k in range(0,len(average_left_minus_right)):
		average_left_minus_right[k] = average_left_minus_right[k][thresholded_start_index:thresholded_end_index]
		average_left[k] = average_left[k][thresholded_start_index:thresholded_end_index]
		average_right[k] = average_right[k][thresholded_start_index:thresholded_end_index]
		average_head_yaw[k] = average_head_yaw[k][thresholded_start_index:thresholded_end_index]
		resampled_trajectories[k] = resampled_trajectories[k][thresholded_start_index:thresholded_end_index]

		
	# optimized_fitted_sin_curve, amp, phase, average_over_sinusoid_cycles  = fourier_transform_to_get_cutoff(mean,t)
	# print (f'phase of mean response {phase}')
	# optimized_fitted_sin_curve, amp, phase, average_over_sinusoid_cycles  = fourier_transform_to_get_cutoff(mean_trajectory,t)
	# print (f'phase of mean trajectories {phase}')


	######################## plots of mean traces with confidence interval using seaborn ################################
	
	# fig, ax1 = plt.subplots(figsize=(6,2.0))
	
	# ax1.set_ylim(-60,60)
	# adjust_spines(ax1,['bottom', 'left'],xticks=[0,20],  yticks=[-20,20], linewidth=0.6,spineColor='k')

	# data = {'t': list(t)*len(average_left_minus_right), 
    #     'wba': np.array(average_left_minus_right).flatten()}
	# df_multi = pd.DataFrame(data)

	# for k in range(0,len(average_left_minus_right)):
	# 	ax1.plot(t, average_left_minus_right[k] - np.mean(average_left_minus_right[k]),c='g', alpha=0.1)

	# ax1.plot(t, mean, c='g')
	# #sns.lineplot(x='t', y='wba', data=df_multi, axes=ax1, color='r')
	
	# ax4 = ax1.twinx()
	# ax4.plot(t, mean_trajectory, 'b')
	# ax4.set_ylim(-530,530)
	# adjust_spines(ax4,['bottom', 'right'],xticks=[0,20],  yticks=[-200,200], linewidth=0.6,spineColor='k')


	# plt.show()


	##################################################################################################################### 


	# mean_norm = (mean - np.mean(mean)) / np.std(mean)
	# trajectory_norm = (mean_trajectory - np.mean(mean_trajectory)) / np.std(mean_trajectory)

	# plt.plot(t, mean, 'g')
	# plt.plot(t, trajectory_norm, 'b')
	# plt.show()

	#mean_baseline_left_trace = np.mean(np.array(avg_baseline_left_trace[str(i)]), axis=0)
	#mean_baseline_left_trace = np.array(avg_baseline_left_trace[str(i)][0])
	#mean_baseline_right_trace = np.mean(np.array(avg_baseline_right_trace[str(i)]), axis=0)
	#mean_baseline_right_trace = np.array(avg_baseline_right_trace[str(i)][0])
	#mean_baseline_head_trace = np.mean(np.array(avg_baseline_head_trace[str(i)]), axis=0)
	
	# plt.clf()
	# plt.ylim(-40,40)
	# plt.plot(np.linspace(0,15,mean_baseline_right_trace.shape[0]), mean_baseline_head_trace, 'k')
	# plt.plot(np.linspace(0,15,mean_baseline_right_trace.shape[0]),[0]*mean_baseline_right_trace.shape[0],'b' )
	# plt.show()

	if head_data:
		mean_roll = np.mean(average_head_roll, axis = 0)
		mean_pitch = np.mean(average_head_pitch, axis = 0)
		mean_yaw = np.mean(average_head_yaw, axis = 0)

	
	#a1 = np.linspace(-15,0,mean_baseline_right_trace.shape[0])
	# if show_plots:
		
		
	# 	### plot average across all ###
	# 	## wings
	# 	axes.plot(np.concatenate([a1, t+a1[-1]]), np.concatenate([mean_baseline_left_trace, mean_left]), c='k', label='mean wba yaw', linewidth=0.6)
	# 	axes2 = axes.twinx()
	# 	axes2.plot(np.concatenate([a1, t+a1[-1]]), np.concatenate([mean_baseline_right_trace, mean_right]), c='r', label='mean wba yaw', linewidth=0.6)

	# 	### head
	# 	#axes.plot(np.concatenate([a1, t+a1[-1]]), np.concatenate([mean_baseline_head_trace, mean_yaw]), c='k', label='mean wba yaw', linewidth=0.6)

	# 	#axes.set_ylim(-50,50)
	# 	axes2.set_ylim(0,80)
	# 	axes.set_ylim(0,80)
	# 	axes.set_xlim(-15, a1[-1]+t[-1])
	# 	#axes.set_ylabel('L- R wba yaw mean across flies (deg) - Red')

	
	# have to resample even the rockafly trajectory
	#resampled_trajectory = resample(trajectories_time[i], trajectories[i], t, kind='linear')
	
	## store the traces of individual flies for head-wing phase analysis
	pickle.dump([average_left_minus_right, average_head_yaw, resampled_trajectories, t], open('./data_for_wing_head_phase_new/'+subdir[:-1]+'-speed-'+str(i+1)+'.p', 'wb'))
	

	if show_plots:
		axes3 = axes.twinx()
		#axes3.plot(t, resampled_trajectory,c='b')
		axes3.plot(np.concatenate([a1, t+a1[-1]]), np.concatenate([[0]*mean_baseline_right_trace.shape[0],resampled_trajectory]), c='b', linewidth=0.6)
		#axes3.set_ylabel('')

		adjust_spines(axes,['bottom', 'left'],xticks=[-15,5],  yticks=[30,70], linewidth=0.6,spineColor='k')
				
		adjust_spines(axes3,['bottom', 'right'],xticks=[-15,5], yticks=[-200, 200], linewidth=0.6,spineColor='k')
		adjust_spines(axes2,['bottom'],xticks=[-15,5], linewidth=0.6,spineColor='k')

		
		#adjust_spines(axes3,['right'],yticks=[-200,0,200],linewidth=0.8,spineColor='k')

		
		#plt.title(subdir + ' n = '+str(average_left_minus_right.shape[0]) + ' speed ' + str(i+1))
		axes.set_xlabel('20 s')
		axes.set_ylabel(' Wing \n response \n (degs)')
		axes3.set_ylabel('Angular \n velocity \n' + r'(deg $s^{-1}$)')
		

		plt.show()
		# plt.clf()
		# plt.figure(figsize=(3,0.7))		
		# plt.ylim(-40,40)
		# plt.xlim(-15, t[-1]+a1[-1])
		# plt.xlabel('Time (s)', fontsize=8, font='Arial')
		# plt.ylabel('Wing response \n (degs)', fontsize=8, font='Arial')
		# plt.plot(np.concatenate([a1, t+a1[-1]]), np.concatenate([mean_baseline_head_trace, mean_yaw]), 'k')
		
		# #plt.plot(t, mean_yaw, c='k', label='yaw')
		# #plt.plot(t, resampled_trajectory, c='b')

		# plt.show()

		###############		
		# if head_data:
		# 	plt.clf()		
		# 	plt.ylim(-30,30)
		# 	plt.xlabel('Time (seconds)')
		# 	plt.ylabel('Head roll angle (degrees) ')
		# 	plt.plot(t, mean_roll, c='r', label='left - right (avg)')
		# 	plt.show()





	## we want, for each fly, the average over 5 sinusoid cycles and the value of the fit to that in order to plot on a scatter plot
	## the average_left_minus_right containts WBA for all flies for a particular speed, resampled already to a common 't'
	calculate_average_over_cycles_per_fly(average_left_minus_right, t, resampled_trajectories,i, head=0, store_baselines=True)
	

	list_of_signals_per_fly_left = calculate_average_over_cycles_per_fly(average_left, t, resampled_trajectories, i, head=0, store_baselines=False, left=True)
	list_of_signals_per_fly_right = calculate_average_over_cycles_per_fly(average_right, t, resampled_trajectories, i, head=0, store_baselines=False, left=False, right=True)

	calculate_average_over_cycles_per_fly(average_head_yaw, t, resampled_trajectories, i,  head=3)

	## do the same for the head
	# if head_data:
	# 	calculate_average_over_cycles_per_fly(average_head_roll, t, resampled_trajectories, i, head=1)
	# 	calculate_average_over_cycles_per_fly(average_head_pitch, t, resampled_trajectories, i, head=2)
	# 	calculate_average_over_cycles_per_fly(average_head_yaw, t, resampled_trajectories, i, head=3)

	############ ipsi-contra analysis ###################################
	'''
	Basically we want to test if across the genotypes with campaniforms silenced, if the increase in amplitude of one wing is reduced while decreasing the
	amplitude on the other side?  
	We want to see if the slope of the line fitted through amplitude of the inside wing (decreasing in amplitude) plotted against number of campaniforms silenced
	is significant and same for outside (increasing in amplitude) wing. We want to do this for each fly.

	For yaw-pos, outside and inside need to be swapped thats all. 

	'''
	inside_wing_values = []
	outside_wing_values = []

	for z in range(0, len(list_of_signals_per_fly_left)):
		left_trace = list_of_signals_per_fly_left[z]
		right_trace = list_of_signals_per_fly_right[z]

		## taking the mean of a window at quarter period and three quarter period. Window length being 1second (72 indices)
		left_inside_quater = np.mean(left_trace[int(len(left_trace)*0.25) - 72:int(len(left_trace)*0.25) + 72])
		left_outside_three_quater = np.mean(left_trace[int(len(left_trace)*0.75) - 72:int(len(left_trace)*0.75) + 72])

		right_inside_three_quater = np.mean(right_trace[int(len(right_trace)*0.75) - 72:int(len(right_trace)*0.75) + 72])
		right_outside_quater = np.mean(right_trace[int(len(right_trace)*0.25) - 72:int(len(right_trace)*0.25) + 72])

		## take the mean of those two per fly
		inside_wing_values.append((left_inside_quater + right_inside_three_quater)/2)
		outside_wing_values.append((left_outside_three_quater + right_outside_quater)/2)


	if subdir[-4:]=='pos/':
		print ('##### pos : Swapping inside and outside ##')
		swap = inside_wing_values 
		inside_wing_values = outside_wing_values
		outside_wing_values = swap

	pickle.dump([inside_wing_values, outside_wing_values], open('./data_for_ipsi_contra_new/'+subdir[:-1]+'-inside-outside-speed-'+str(i+1)+'.p', 'wb'))



	###################################################################









	continue


	#### also plot the average over the 5 sinusoid cycles of the mean trace (over all flies)
	num_cycles = 5
	

	temp = []
	temp_head = []
	t1 = int(round(t.shape[0]/float(num_cycles)))
	
	# we want to resample the 5 sinusoids to length t1 because the total number of indexes might not be divisible by 5.

	# for z in range(0,average_left_minus_right.shape[0]):
	# 	plt.plot(range(0,mean.shape[0]), average_left_minus_right[z,:])
	# plt.show()


	resampling_dt = np.mean(np.diff(t[0:t1]))
	t2 = np.arange(len(t[0:t1])) * resampling_dt

	for x in range(0, num_cycles):
		temp.append(resample(t[x*t1:(x+1)*t1] - t[x*t1] , mean[x*t1:(x+1)*t1], t2, kind='linear'))

		temp_head.append(resample(t[x*t1:(x+1)*t1] - t[x*t1] , mean_yaw[x*t1:(x+1)*t1], t2, kind='linear'))
		
		if head_data:
			temp_head.append(resample(t[x*t1:(x+1)*t1] - t[x*t1] , mean_roll[x*t1:(x+1)*t1], t2, kind='linear'))

	average_over_sinusoid_cycles = np.mean(temp, axis=0)
	#average_over_sinusoid_cycles = temp[0]
	
	std_dev = np.std(temp, axis=0)

	average_over_sinusoid_cycles_yaw = np.mean(temp_head, axis=0)
	std_dev_yaw = np.std(temp_head, axis=0)

	if head_data:
		average_over_sinusoid_cycles_roll = np.mean(temp_head, axis=0)	
		std_dev_roll = np.std(temp_head, axis=0)

	#### now we also want to fit a sine curve to this average_over_sinusoid_cycles 
	optimized_fitted_sin_curve, amp, phase, _ = fourier_transform_to_get_cutoff(average_over_sinusoid_cycles,t2)
	#import ipdb;ipdb.set_trace()
	#if show_plots:
	plt.clf()
	figure, axes = plt.subplots()
	## plot everything
	
	axes.plot(t2, average_over_sinusoid_cycles, 'g')
	axes.plot(t2, optimized_fitted_sin_curve, 'r')
	
	axes.set_xlabel('Time (seconds)')
	axes.set_ylabel(' L - R WBA (degrees)')
	axes.set_ylim(-40,40)
	#plt.ylim(-25,25)
	##axes.plot(t[0:t1], fitted_sin_curve, 'r')
	
	#plt.fill_between(t[0:t1], average_over_sinusoid_cycles-std_dev, average_over_sinusoid_cycles+std_dev,color='green',alpha=0.2)
	axes2 = axes.twinx()
	axes2.plot(t2, mean_trajectory[0:t1], 'b')
	axes2.set_ylabel('Stimulus angular velocity (degrees/second)')
	#plt.xlabel('Time (seconds)')
	#plt.ylabel('L - R WBA (degrees)')
	plt.title(subdir + ' n = '+str(len(average_left_minus_right)) + ' collapsed ' + ' speed ' + str(i+1) + 'amp:' + str(amp) + 'phase:' + str(phase))
	plt.show()
	plt.clf()


	optimized_fitted_sin_curve, amp, phase, _ = fourier_transform_to_get_cutoff(average_over_sinusoid_cycles_yaw,t[0:t1])
	#import ipdb;ipdb.set_trace()
	## plot everything
	figure, axes = plt.subplots()
	# plot everything
	axes.plot(t[0:t1], average_over_sinusoid_cycles_yaw, 'g')
	axes.set_xlabel('Time (seconds)')
	axes.set_ylabel(' Head yaw (degrees)')
	axes.set_ylim(-20,20)
	#plt.ylim(-25,25)
	##axes.plot(t[0:t1], fitted_sin_curve, 'r')
	
	#plt.fill_between(t[0:t1], average_over_sinusoid_cycles-std_dev, average_over_sinusoid_cycles+std_dev,color='green',alpha=0.2)
	axes2 = axes.twinx()
	axes2.plot(t[0:t1], mean_trajectory[0:t1], 'b')
	axes2.set_ylabel('Stimulus angular velocity (degrees/second)')
	#plt.xlabel('Time (seconds)')
	#plt.ylabel('L - R WBA (degrees)')
	plt.title(subdir + ' head yaw n = '+str(len(average_left_minus_right)) + ' collapsed ' + ' speed ' + str(i+1) )
	plt.show()
	plt.clf()


	# if head_data:
	# 	#optimized_fitted_sin_curve, amp, phase, _ = fourier_transform_to_get_cutoff(average_over_sinusoid_cycles_roll,t[0:t1])
	# 	## plot everything
	# 	plt.plot(t[0:t1], average_over_sinusoid_cycles_roll, 'g')
	# 	plt.fill_between(t[0:t1], average_over_sinusoid_cycles_roll-std_dev_roll, average_over_sinusoid_cycles_roll+std_dev_roll,color='green',alpha=0.2)
	# 	plt.ylim(-25,25)
	# 	plt.plot(t[0:t1], optimized_fitted_sin_curve, 'r')			
	# 	plt.xlabel('Time (seconds)')
	# 	plt.ylabel('Head roll angle (degrees)')
	# 	plt.title('head roll ' + subdir + ' n = '+str(average_left_minus_right.shape[0]) + ' collapsed ' + ' speed ' + str(i+1) )
	# 	plt.show()
	# 	plt.clf()


	########### make nice seaborn plots which automatically gives confidence intervals ##############################
	
	num_cycles = 5
	temp = []
	temp_head = []
	t1 = int(round(t.shape[0]/float(num_cycles)))
	
	# we want to resample the 5 sinusoids to length t1 because the total number of indexes might not be divisible by 5.
	resampling_dt = np.mean(np.diff(t[0:t1]))
	t2 = np.arange(len(t[0:t1])) * resampling_dt

	

	for f in range(0, len(average_left_minus_right)):
		f_WBA = average_left_minus_right[f]
		f_head_yaw = average_head_yaw[f]

		for x in range(0, num_cycles):
			temp.append(resample(t[x*t1:(x+1)*t1] - t[x*t1] , f_WBA[x*t1:(x+1)*t1], t2, kind='linear'))

			temp_head.append(resample(t[x*t1:(x+1)*t1] - t[x*t1] , f_head_yaw[x*t1:(x+1)*t1], t2, kind='linear'))
			
		
	
	# # ########################### I want to this plot where I have cycle wise collapsed traces in order to measure drift ########
	# figure, axes = plt.subplots()
	# cmap = plt.get_cmap('jet', 5)
	
	# for outer in range(0,5):
	# 	cycle_list = []
	# 	for f in range(outer, len(temp_head), 5):
	# 		cycle_list.append(temp_head[f])
			
	# 	sns.lineplot(list(t2)*len(cycle_list), list(np.array(cycle_list).flatten()),ax = axes, color=cmap(outer), alpha=1-(0.1*outer))
	# 	axes.set_xlim([0,t2[-1]])
	# 	axes.set_ylim([-20,20])
	# 	axes.set_xlabel('Time (seconds)')
	# 	axes.set_ylabel('Head yaw  (degrees)')

	# norm = mpl.colors.Normalize(vmin=0, vmax=5)
	  
	# # creating ScalarMappable
	# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	# sm.set_array([])
	  
	# plt.colorbar(sm, ticks=np.linspace(1, 5, 5))

	# plt.show()
	# plt.clf()



	##############################################################################################################################

	########### adaptation : I want to make a plot with cycle_number on the x axis and amp (per fly per cycle) on the y axis. Fit a line and test if the line is significant. ##########3
	# fits_per_fly_per_cycle = {1:[], 2:[], 3:[], 4:[], 5:[]}
	
	# for outer in range(0, num_cycles):
	# 	for f in range(outer, len(temp), 5):
	# 		_, fit,_, _ = fourier_transform_to_get_cutoff(temp[f], t2)
	# 		fits_per_fly_per_cycle[outer+1].append(fit)

	# adaptation_x_axis = []
	# adaptation_y_axis =[]
	# num_flies = len(fits_per_fly_per_cycle[1])
	# for k in fits_per_fly_per_cycle.keys():
	# 	adaptation_x_axis.extend([k]*num_flies)
	# 	adaptation_y_axis.extend(fits_per_fly_per_cycle[k])

	# adaptation_x_axis = np.array(adaptation_x_axis)
	# adaptation_y_axis = np.array(adaptation_y_axis)
	# slope, intercept, r, p, se = stats.linregress(adaptation_x_axis, adaptation_y_axis)
	# plt.plot(adaptation_x_axis, adaptation_x_axis*slope + intercept, 'r')
	# print ('slope is', slope)
	# print ('p value is ', p)
	
	# plt.title('wing adaptation slope: ' + str(round(slope,2)) + ' p value: '+ str(round(p,3)))
	# plt.xlabel('cycle number')
	# plt.ylabel('amplitude of fit to response')
	# plt.scatter(adaptation_x_axis, adaptation_y_axis)
	# plt.show()	


	# plt.clf()
	# fits_per_fly_per_cycle = {1:[], 2:[], 3:[], 4:[], 5:[]}
	
	# for outer in range(0, num_cycles):
	# 	for f in range(outer, len(temp_head), 5):
	# 		_, fit,_, _ = fourier_transform_to_get_cutoff(temp_head[f], t2)
	# 		fits_per_fly_per_cycle[outer+1].append(fit)

	# adaptation_x_axis = []
	# adaptation_y_axis =[]
	# num_flies = len(fits_per_fly_per_cycle[1])
	# for k in fits_per_fly_per_cycle.keys():
	# 	adaptation_x_axis.extend([k]*num_flies)
	# 	adaptation_y_axis.extend(fits_per_fly_per_cycle[k])

	# adaptation_x_axis = np.array(adaptation_x_axis)
	# adaptation_y_axis = np.array(adaptation_y_axis)
	# slope, intercept, r, p, se = stats.linregress(adaptation_x_axis, adaptation_y_axis)
	# print ('slope is', slope)
	# print ('p value is ', p)
	# plt.plot(adaptation_x_axis, adaptation_x_axis*slope + intercept, 'r')
	
	# plt.title('head adaptation slope: ' + str(round(slope,2)) + ' p value: '+ str(round(p,3)))
	# plt.xlabel('cycle number')
	# plt.ylabel('amplitude of fit to response')
	# plt.scatter(adaptation_x_axis, adaptation_y_axis)
	# plt.show()


	########################################### collapsed with confidence interval in red ########################################################################

	# plt.clf()
	# figure, axes = plt.subplots()
	# figure.set_size_inches(2,1)

	# ## write data to pickle file for other plotting code
	# if i==1:
	# 	print ('writing wing and head data to pickle file for 2nd speed')
	# 	list_for_pickle = [t2, temp, temp_head, resampled_trajectory[0:t1]]
	# 	pickle.dump( list_for_pickle, open( subdir[:-1] + '_data_for_lineplot.p', 'wb' ) )
	
	# sns.lineplot(list(t2)*len(temp), list(np.array(temp).flatten()),ax = axes, sort=False,color='k')
	# axes.set_xlim([0,t2[-1]])
	# axes.set_ylim([-40,40])
	# adjust_spines(axes,['bottom', 'left'],xticks=[0,20],yticks=[-20,0,20], linewidth=0.8,spineColor='k')
	# axes2 = axes.twinx()
	# axes2.plot(t2, resampled_trajectory[0:t1], 'b' )
	# adjust_spines(axes2,['bottom'],xticks=[0,20], linewidth=0.8,spineColor='k')
	
	# axes.set_xlabel('Time (s)')
	# axes.set_ylabel('Wing response \n (deg)')
	

	# #axes2.set_ylabel('Stimulus angular velocity (degrees/second)')
	# #plt.title(subdir + ' WBA yaw n = '+str(average_left_minus_right.shape[0]) + ' collapsed ' + ' speed ' + str(i+1))
	#plt.show()


	# figure, axes = plt.subplots()
	
	# sns.lineplot(list(t2)*len(temp_head), list(np.array(temp_head).flatten()),ax = axes, color='r')
	# axes.set_xlim([0,t2[-1]])
	# axes.set_ylim([-20,20])
	# axes.set_xlabel('Time (seconds)')
	# axes.set_ylabel('head yaw grand mean (deg) - Red')

	# axes2 = axes.twinx()
	# axes2.plot(t2, resampled_trajectory[0:t1], 'b' )
	# #axes2.set_ylabel('Stimulus angular velocity (degrees/second)')
	# plt.title(subdir + ' head yaw n = '+str(average_left_minus_right.shape[0]) + ' collapsed ' + ' speed ' + str(i+1))


	# plt.show()