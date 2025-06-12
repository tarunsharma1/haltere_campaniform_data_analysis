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
from sklearn.preprocessing import StandardScaler
import math
import itertools
import scatter_plot

def bootstrap(data, n):
	## randomly resample n times and take the mean each time
	bootstrapped_data = np.zeros(n)
	for i in range(0,n):
		sample = np.random.choice(data, size=len(data))
		bootstrapped_data[i] = np.mean(np.array(sample))
	return bootstrapped_data


def confidence_interval(data):
	## get the 95% confidence interval by getting the 2.5th and 97.5th percentile of the data
	conf_interval = np.percentile(data,[2.5,97.5])
	print (conf_interval)
	return conf_interval[0], conf_interval[1]



def resample(x1, y1, x2, kind, extrapolate=True):
	# helper function
	if kind == 'spline':
		spline = spi.CubicSpline(x1, y1)
		y2 = spline.__call__(x2, extrapolate=extrapolate)
	else:
		fill_value = "extrapolate" if extrapolate else []
		interp = spi.interp1d(x1, y1, kind=kind, bounds_error=False, fill_value=fill_value)
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
	#print ('frequency:',freq)

	## add offset np.mean(signal) just for visualizing the fits better
	#phase = 90
	signal_reconstucted =  amp*np.cos(2*np.pi*freq*t + phase*np.pi/180) + np.mean(signal)


	### lets try doing an optimization here using scipy to increase fit quality
	optimize_func = lambda x: x[0]*np.cos(2*np.pi*x[1]*t+x[2]*np.pi/180) + np.mean(signal) - signal
	est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [amp, freq, phase, np.mean(signal)])[0]
	#print ('##### scipy optimized #####')
	#print ('amp:', est_amp)
	#print ('phase', est_phase)
	#print ('frequency:', est_freq)
	optimized_signal_reconstructed =  est_amp*np.cos(2*np.pi*est_freq*t + est_phase*np.pi/180) + np.mean(signal)

	#plt.clf()
	#plt.ylim(-8,8)
	#return optimized_signal_reconstructed, est_amp, est_phase, signal
	return signal_reconstucted, amp, phase, signal


num_cycles = 5
genotype_names = ['Control','R14B04', 'R17G01', 'R28C05', 'R60B12','R58F02', 'R31A09','R22E04','R86H12']
list_of_genotypes = ['KirJ00yawbothdirections', 'KirJ70yawbothdirections', 'KirJ73yawbothdirections', 'Kir28C05yawbothdirections','KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections',  'KirJ75yawbothdirections', 'KirJ96yawbothdirections']
 

data_for_scatter_wings = []

data_for_scatter_head_ = []
data_for_scatter_x = []


x_list = []
y_list = []
x_labels = []
y_conf_min = []
y_conf_max = []
bootstrapped_means = []


speed = 2




for k,genotype in enumerate(list_of_genotypes):
	print (genotype)
	
	avg_wba_, avg_head_yaw_, trajectory_, t_ = pickle.load(open('./data_for_wing_head_phase/'+genotype+'-speed-'+str(speed)+'.p', 'rb'))
	wing_phases = []
	head_phases = []
	wing_trajectory_difference = []
	head_trajectory_difference = []
	wing_head_difference = []
	grand_mean_wings = []
	grand_mean_ = []


	

	## for each fly in the collection avg_wba and avg_head_yaw, take average over 5 sinusoid cycles and calculate phase

	for fly in range(0,len(avg_wba_)):
		
		wba_fly_ = avg_wba_[fly]
		head_yaw_fly_ = avg_head_yaw_[fly]
	
		temp = []
		temp_head = []
		resampling_dt = np.mean(np.diff(t_))
		
		tx = np.arange(len(t_)) * resampling_dt
		t1 = int(round(tx.shape[0]/float(num_cycles)))

		# we want to resample the 5 sinusoids to length t1 because the total number of indexes might not be divisible by 5.

		t2 = tx[0:t1]
			
		for x in range(0, num_cycles):
			temp.append(resample(t_[x*t1:(x+1)*t1] - t_[x*t1] , wba_fly_[x*t1:(x+1)*t1], t2, kind='linear'))

			temp_head.append(resample(t_[x*t1:(x+1)*t1] - t_[x*t1] , head_yaw_fly_[x*t1:(x+1)*t1], t2, kind='linear'))
			

		average_over_sinusoid_cycles_ = np.mean(temp, axis=0)
		std_dev = np.std(temp, axis=0)

		average_over_sinusoid_cycles_yaw_ = np.mean(temp_head, axis=0)
		std_dev_yaw = np.std(temp_head, axis=0)
		
		grand_mean_.append(average_over_sinusoid_cycles_yaw_)
		grand_mean_wings.append(average_over_sinusoid_cycles_)


		#####################################################################
		## get phase for each fly
		optimized_fitted_sin_curve, amp, w_phase, _ = fourier_transform_to_get_cutoff(average_over_sinusoid_cycles_,t2)
		wing_phases.append(w_phase)
		

		optimized_fitted_sin_curve, amp, h_phase, _ = fourier_transform_to_get_cutoff(average_over_sinusoid_cycles_yaw_,t2)
		head_phases.append(h_phase)
	

	y_list_temp = np.array(head_phases)

	bootstrapped_data = bootstrap(y_list_temp, 10000)
	conf_min, conf_max = confidence_interval(bootstrapped_data)
	y_conf_min.append(conf_min)
	y_conf_max.append(conf_max)

	x_list.append([k]*len(y_list_temp))
	y_list.append(y_list_temp)
	bootstrapped_means.append(np.mean(bootstrapped_data))

	x_labels.append(genotype_names[k])


fig, axes = plt.subplots()
#fig.set_size_inches(2.2,1.5)
fig.set_size_inches(2.0,1.3)

#fig.set_size_inches(0.2,0.5)

x_axis = list(itertools.chain.from_iterable(x_list))
## add some jitter for visaulization purposes

x_axis = x_axis + np.random.uniform(low=-0.11, high=0.11, size=len(x_axis))

axes.scatter(x_axis, list(itertools.chain.from_iterable(y_list)), marker='.', c=[[0.5,0.5,0.5]],s=4)

## plot the means and conf intervals
#plt.scatter(range(0,len(list_of_genotypes)), [np.mean(y) for y in y_list],c='r')

## convert the confidence interval values into mean - value so as to use matplotlib errbar function
y_err_lower = []
y_err_upper = []
for i,y in enumerate(y_list):
	#mean = np.mean(y)
	mean = bootstrapped_means[i]
	y_err_lower.append(mean - y_conf_min[i])
	y_err_upper.append(y_conf_max[i] - mean)


axes.errorbar(range(0,len(list_of_genotypes)), bootstrapped_means, yerr=[y_err_lower, y_err_upper], fmt=".", c='k', ecolor='k', elinewidth=1)

		
#adjust_spines(axes,['bottom', 'left'],xticks=range(0,len(list_of_genotypes)), yticks=[8,15], linewidth=0.6,spineColor='k')
#adjust_spines(axes,['left'], yticks=[20,30], linewidth=0.6,spineColor='k')
scatter_plot.adjust_spines(axes,['bottom', 'left'],xticks=range(0,len(list_of_genotypes)), yticks=[100, 130], linewidth=0.6,spineColor='k')

#axes.set_xticks(range(0,len(list_of_genotypes)), x_labels)
axes.set_xticklabels(x_labels, rotation=90)
#axes.set_xlim(-0.12,0.12)
#axes.set_ylim(0,70)
#plt.ylim(0,4)
axes.set_ylabel('Head \n phase \n (degs)')
#plt.ylabel('Ratio of sinusoid fits to wing and head yaw')
#plt.xlabel('genotypes',fontsize=8, font='Arial')
#plt.title('  Baseline left + right of wba UX  ')
plt.show()
plt.clf()






	##### calculate phase difference 95% CIs for wings and head #############
# 	wing_phases = []
# 	head_phases = []
# 	wing_head_differences = []
# 	for x in range (0,len(grand_mean_wings)):
# 		wing_set = np.array(grand_mean_wings[x])
# 		head_set = np.array(grand_mean_[x])
# 		optimized_fitted_sin_curve_wings, amp, wings_phase, _ = fourier_transform_to_get_cutoff(wing_set ,np.arange(0,wing_set.shape[0]))
# 		wing_phases.append(wings_phase)
# 		optimized_fitted_sin_curve_wings, amp, head_phase, _ = fourier_transform_to_get_cutoff(head_set ,np.arange(0,head_set.shape[0]))
# 		head_phases.append(head_phase)
# 		wing_head_differences.append(wings_phase - head_phase)


# 	print ('len of wing phases is:', len(wing_phases))
# 	bootstrapped_wings = bootstrap(wing_phases, 10000)
# 	print ('conf interval is', confidence_interval(bootstrapped_wings))
# 	print ('bootstrapped mean is ', np.mean(np.array(bootstrapped_wings)))
# 	bootstrapped_diff = bootstrap(wing_head_differences, 10000)
# 	print ('conf_interval of diff is ', confidence_interval(bootstrapped_diff))




# 	########################################################################

# 	#fig, axes = plt.subplots()

# 	mean2 = np.mean(np.array(grand_mean_),axis=0)
# 	mean_wings = np.mean(np.array(grand_mean_wings), axis=0)

	
# 	#axes2.plot(np.arange(0, mean2.shape[0]), mean2, 'r')
# 	optimized_fitted_sin_curve_, amp, w_phase, _ = fourier_transform_to_get_cutoff(mean2 ,np.arange(0,mean2.shape[0]))
# 	max_index_ = np.argmax(optimized_fitted_sin_curve_)
# 	print ('phase DLC :', w_phase)
# 	data_for_scatter_head_.append(w_phase)


# 	optimized_fitted_sin_curve_wings, amp, wings_phase, _ = fourier_transform_to_get_cutoff(mean_wings ,np.arange(0,mean_wings.shape[0]))
# 	max_index_wings = np.argmax(optimized_fitted_sin_curve_wings)
# 	print ('phase wings :', wings_phase)
# 	data_for_scatter_wings.append(wings_phase)

# 	#axes2.plot(t2, mean_wings, 'r')
# 	#axes2.plot(t2,optimized_fitted_sin_curve_wings)


# 	print ('time difference wings and DLC:', t2[max_index_wings] - t2[max_index_])

# 	data_for_scatter_x.append(k)
	

	

# plt.ylim(85,130)

# plt.scatter(data_for_scatter_x, data_for_scatter_head_)
# plt.scatter(data_for_scatter_x, data_for_scatter_wings)
# plt.xlabel(subdirs, rotation=45)
# plt.show()