import glob
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import scipy.fftpack
import scipy.signal as sps
import pickle
from sklearn.metrics import r2_score

def sign_with_positive_zero(arr):
    return np.sign(arr) + (arr == 0)


class RkfAnalysis:
	def __init__(self, number_of_speeds=2):

		## flag to indicate if top view camera data is available
		self.head_data = False

		# 4 for roll , theres more for yaw
		self.number_of_speeds = number_of_speeds
		# these two are only for the title for plots
		self.list_of_periods_used = [25.84, 18.01]
		self.list_of_peak_velocities_used = [2*np.pi*1/t*1440 for t in self.list_of_periods_used]

		#self.number_of_cycles_per_speed = 5
		
		self.kinefly_data = {}
		self.trajectory_data = {}
		self.command_data = {}
		self.top_camera_data = {}

		# this will contain the kinefly data and trajectory data for each speed (period - all cycles)
		self.kinefly_data_separated_on_speed = {}
		self.kinefly_data_separated_on_speed['left'] = []
		self.kinefly_data_separated_on_speed['right'] = []
		self.kinefly_data_separated_on_speed['head'] = []
		self.kinefly_data_separated_on_speed['aux'] = []
		self.kinefly_data_separated_on_speed['timestamp'] = []
		
		self.trajectory_data_separated_on_speed = {}
		self.trajectory_data_separated_on_speed['position'] = []
		self.trajectory_data_separated_on_speed['timestamp'] = []

		self.trajectory_data_separated_on_speed_baseline = {}
		self.trajectory_data_separated_on_speed_baseline['position'] = []
		self.trajectory_data_separated_on_speed_baseline['timestamp'] = []

		
		# this will contain the kinefly data for each speed baseline (15 second sleep) and its timestamp (these two are not so important)
		self.kinefly_data_separated_on_speed_baseline = {}
		self.kinefly_data_separated_on_speed_baseline['left'] = []
		self.kinefly_data_separated_on_speed_baseline['right'] = []
		self.kinefly_data_separated_on_speed_baseline['head'] = []
		self.kinefly_data_separated_on_speed_baseline['aux'] = []
		self.kinefly_data_separated_on_speed_baseline['timestamp'] = []
		
		self.IMU_data = {}
		self.IMU_data['timestamp'] = []
		self.IMU_data['gyro_x'] = []
		self.IMU_data['gyro_y'] = []
		self.IMU_data['gyro_z'] =[]
		self.IMU_data['acc_x'] = []
		self.IMU_data['acc_y'] = []
		self.IMU_data['acc_z'] = []

		self.top_camera_data_separated_on_speed = {}
		self.top_camera_data_separated_on_speed['roll'] = []
		self.top_camera_data_separated_on_speed['pitch'] = []
		self.top_camera_data_separated_on_speed['timestamp'] = []

		self.top_camera_data_separated_on_speed_baseline = {}
		self.top_camera_data_separated_on_speed_baseline['roll'] = []
		self.top_camera_data_separated_on_speed_baseline['pitch'] = []
		self.top_camera_data_separated_on_speed_baseline['timestamp'] = []


		self.trajectory_speed_segregation_indices = []

		self.counter = 0

		# this is made as a class variable because it is used in the 3D headtracking code
		self.a3 = []
		self.h3 = []

	def read_hdf5(self, filename):
		hf = h5py.File(filename, 'r')
		print (hf.keys())

		if self.head_data==False:
			#self.kinefly_data['head'] = np.array(hf.get('kinefly_data_head'))
			head_yaw_file = filename.split('_compressed')[0] + '_head_yaw.p'
			self.kinefly_data['head'] = np.array(pickle.load(open(head_yaw_file,'rb')))

		self.kinefly_data['left'] = np.array(hf.get('kinefly_data_left'))
		self.kinefly_data['right'] = np.array(hf.get('kinefly_data_right'))
		self.kinefly_data['aux'] = np.array(hf.get('kinefly_data_aux'))
		self.kinefly_data['timestamp'] = np.array(hf.get('kinefly_data_timestamp')) 
		

		#self.trajectory_data['position'] = np.array(hf.get('trajectory_data'))
		#self.trajectory_data['timestamp'] = np.array(hf.get('trajectory_data_timestamp'))
		
		self.command_data['command'] = np.array(hf.get('command_data'))
		self.command_data['timestamp'] = np.array(hf.get('command_data_timestamp'))
		
		self.IMU_data['timestamp'] = np.array(hf.get('IMU_time'))
		self.IMU_data['gyro_x'] = np.array(hf.get('IMU_gyro_x_data'))
		self.IMU_data['gyro_y'] = np.array(hf.get('IMU_gyro_y_data'))
		self.IMU_data['gyro_z'] = np.array(hf.get('IMU_gyro_z_data'))
		self.IMU_data['acc_x'] = np.array(hf.get('IMU_acc_x_data'))
		self.IMU_data['acc_y'] = np.array(hf.get('IMU_acc_y_data'))
		self.IMU_data['acc_z'] = np.array(hf.get('IMU_acc_z_data'))

		#import ipdb;ipdb.set_trace()

		if self.head_data:
			head_roll_file = filename.split('_compressed')[0] + '_head_roll.p'
			self.top_camera_data['roll'] = np.array(pickle.load(open(head_roll_file,'rb')))
			self.top_camera_data['timestamp'] = np.array(hf.get('top_camera_frames_timestamp'))

			head_pitch_file = filename.split('_compressed')[0] + '_head_pitch.p' 
			self.top_camera_data['pitch'] = np.array(pickle.load(open(head_pitch_file,'rb')))

			head_yaw_file = filename.split('_compressed')[0] + '_head_yaw.p'
			self.kinefly_data['head'] = np.array(pickle.load(open(head_yaw_file,'rb')))

		### for some reason (output by DLC), the head data is usually one index less than the timestamp data..to avoid issues, I am going to add the last value again to the head data
		if self.kinefly_data['head'].shape[0] != self.kinefly_data['timestamp'].shape[0]:
			self.kinefly_data['head'] = np.append(self.kinefly_data['head'], self.kinefly_data['head'][-1]) 
		#import ipdb;ipdb.set_trace()


		
	def find_trajectory_speed_change_points(self):
		########################################################################################################################################################
		## the point of this function is to get the indices of the times at which the sinusoid motions start and end ..so time indices (from the trajectory, in this case IMU data) of [start1, end1, start2, end2]
		## this will be used to separate kinefly data into different segments
		## grbl: as I dont have a trajectory, I am going to use IMU_gyro_z data and IMU_timestamp for this.
		## command data commands are published when move is starting..so they are [-1,25.9,-2,-1,18,-2] for sleeping,running,returning,sleeping,running,returning
		idxs_move_start = list(np.intersect1d(list(np.where(self.command_data['command']!=-1)[0]), list(np.where(self.command_data['command']!=-2)[0])))
		## add the end of the moves so the indices + 1
		idxs = [[i,i+1] for i in idxs_move_start]
		motion_idxs = [item for sublist in idxs for item in sublist]

		motion_times = [self.command_data['timestamp'][i] for i in motion_idxs]

		## the motion_time_idxs are got from the IMU timestamp
		self.trajectory_speed_segregation_indices = []
		
		for i in range(0,len(motion_times), 2):
			start = np.where(self.IMU_data['timestamp'] >= motion_times[i])[0][0]
			end = np.where(self.IMU_data['timestamp'] <= motion_times[i+1])[0][-1]
			
			
			# window_len = 30
			# for start_idx in range(max(0,start), start + window_len):
			# 	window = self.IMU_data['gyro_z'][start_idx: start_idx + window_len]
				
			# 	## within this window check for number of zero crossings.. if there is noise like in the small dashed section, number of zero crossings would be high.
			# 	### once the sinusoid starts, num_zero_crossings should drop to 0
			# 	num_zero_crossings = np.where(np.abs(np.diff(sign_with_positive_zero(window))) == 2)[0].shape[0]

			# 	if num_zero_crossings == 0:
			# 		break
			# print (f'old start idx was {start} and new start idx is {start_idx}')
			# assert start_idx - start < 30

			# ### repeat the same thing for the extra bit at the end as well
			# for end_idx in range(end-60, end-window_len):
			# 	window = self.IMU_data['gyro_z'][end_idx: end_idx + window_len]

			# 	## when the sinusoid is over and the noise starts, number of zero crossings will be high. During sinusoid they will be zero. So we want the index of the first time its nonzero.
			# 	num_zero_crossings = np.where(np.abs(np.diff(sign_with_positive_zero(window))) == 2)[0].shape[0]

			# 	if num_zero_crossings!=0:
			# 		break

			# end_idx += window_len
			# print (f'old end idx was {end} and new end idx is {end_idx}')

			# assert end - end_idx < 30

			# resampling_dt = 72
			# t = np.arange(self.IMU_data['timestamp'][start_idx], self.IMU_data['timestamp'][end_idx-1], 1/resampling_dt)
			# motion_trace_resampled = self.resample(self.IMU_data['timestamp'][start_idx:end_idx], self.IMU_data['gyro_z'][start_idx:end_idx], t, kind='linear')
			
			# signal_reconstucted, amp, phase, signal = self.fourier_transform_to_get_cutoff(motion_trace_resampled, t - t[0])

			# ### lets go in increments shifting the start index and calculating r2_score to find the start index with the lowest r2 score. Then do the same for the end after
			# for new_start_idx in range(0, 30):
			# 	trace = motion_trace_resampled[new_start_idx::]
			# 	r2 = r2_score(trace[0:20], signal_reconstucted[0:20])
			# 	if r2 <=0:
			# 		break
			# 	#signal_reconstucted, amp, phase, signal = self.fourier_transform_to_get_cutoff(trace, t[new_start_idx::] - t[new_start_idx])
			# 	print (f'shifted index {new_start_idx} r2 score {r2}')

			# time_of_new_start_idx = t[new_start_idx]
			# ## index of original signal (before resampling)

			# og_idx = np.where(self.IMU_data['timestamp'] >= time_of_new_start_idx)[0][0]
			# print (f'index is {og_idx}')
			

			self.trajectory_speed_segregation_indices.append(start)
			self.trajectory_speed_segregation_indices.append(end)

		print (self.trajectory_speed_segregation_indices)
		

	


	def adjust_spines(self, ax_handle, spines = ["bottom", "left"]):
	    # helper function for plotting nicely

	    """ ax_handle is just the name of the axes object you want to manipulate
	    and spines is a list of strings referring to the spine positions to include """
	    for loc, spine in ax_handle.spines.items():
	        if loc in spines:
	            spine.set_position(('outward', 10))  # outward by 10 points
	        else:
	            spine.set_color('none')  # don't draw spine
	    # turn off ticks where there is no spine
	    if 'left' in spines:
	        ax_handle.yaxis.set_ticks_position('left')
	    else:
	        # no yaxis ticks
	        ax_handle.yaxis.set_ticks([])
	    if 'bottom' in spines:
	        ax_handle.xaxis.set_ticks_position('bottom')
	    else:
	        # no xaxis ticks
	        ax_handle.xaxis.set_ticks([])

	def visualize_autostep_motion_and_kinefly_data(self):
		# function to plot the autostep position against time and velocity, acceleration against time.. also to overlay L - R traces over autostep data
		
		## plotting autostep data
		for i in range(0,self.number_of_speeds):
			#plt.plot(self.trajectory_data_separated_on_speed['timestamp'][i], self.trajectory_data_separated_on_speed['position'][i])
			#plt.ylim(-40,40)
			#plt.plot(self.kinefly_data_separated_on_speed['timestamp'][i] - self.kinefly_data_separated_on_speed['timestamp'][i][0], self.kinefly_data_separated_on_speed['left'][i] - self.kinefly_data_separated_on_speed['right'][i])
			plt.clf()
			font = {'family' : 'arial',
			        'weight' : 'normal',
			        'size'   : 12}
			matplotlib.rc('font', **font)
			plt.rcParams['svg.fonttype'] = 'none'

			plt.rcParams['xtick.labelsize']=12
			plt.rcParams['ytick.labelsize']=12
			plt.rcParams['xtick.major.width'] = 0.75
			plt.rcParams['xtick.minor.width'] = 0.75
			plt.rcParams['ytick.major.width'] = 0.75
			plt.rcParams['ytick.minor.width'] = 0.75
			plt.rcParams['axes.linewidth']    = 0.75

			fig = plt.figure(figsize= (5,5))
			ax = plt.subplot(111)
			ax.set_ylim(-25, 25)
			ax.set_xlim(0, 160)
			plt.xticks(np.linspace(0,160, 9, endpoint = True))
			ax.spines['bottom'].set_bounds(0, 160)
			ax.spines['left'].set_bounds(-25, 25)
			plt.yticks(np.linspace(-25, 25, 5, endpoint = True))
			ax.set_xlabel('time (seconds)')
			ax.set_ylabel('WBA (degrees)')

			self.adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])

			### Note in the run_functions function i have a method called ax.twinx() which lets you plot the trajectory and the wing responses (two different Y axis on same X axis)

			#plt.title('peak velocity: '+str(self.list_of_peak_velocities_used[i]))
			#plt.xlabel('time (seconds)')
			#plt.ylabel('WBA (radians)')
			plt.plot(self.kinefly_data_separated_on_speed['timestamp'][i] - self.kinefly_data_separated_on_speed['timestamp'][i][0], self.kinefly_data_separated_on_speed['left'][i], c='g', label='left wing')
			
			plt.plot(self.kinefly_data_separated_on_speed['timestamp'][i] - self.kinefly_data_separated_on_speed['timestamp'][i][0], self.kinefly_data_separated_on_speed['right'][i], c='b', label='right wing')
			plt.legend()

			import ipdb;ipdb.set_trace()
			plt.show()
		



	def resample_kinefly_and_autostep_data(self):
		# do the same thing that kellan did, that is resample the kinefly data in certain timebins and resample autostep data using interpolation and the new kinefly times
		## the fast fourier transform (for sine fit) needs constant time intervals ..so can do this later for now without resampling ill just plot L-R vs responses
		## i can resample now using the segregated_by_speed function
		
		# for each speed, the first time for kinefly and trajectory are almost the same (different only in the 100s of milliseconds)
		for i in range(0, self.number_of_speeds):
			desired_start_time = max(self.trajectory_data_separated_on_speed['timestamp'][i][0], self.kinefly_data_separated_on_speed['timestamp'][i][0])
			desired_end_time = min(self.trajectory_data_separated_on_speed['timestamp'][i][-1], self.kinefly_data_separated_on_speed['timestamp'][i][-1])

			#resampling_dt = np.mean(np.diff(self.trajectory_data_separated_on_speed['timestamp'][i]))
			#t = np.arange(len(self.trajectory_data_separated_on_speed['timestamp'][i])) * resampling_dt
			
			resampling_dt = 72
			t = np.arange(desired_start_time, desired_end_time, 1/resampling_dt)
			
			# resample kinefly data to this new evenly spaced time vector
			for key in ['left','right','head']:
				#self.kinefly_data_separated_on_speed[key][i] = self.resample(self.kinefly_data_separated_on_speed['timestamp'][i]- self.kinefly_data_separated_on_speed['timestamp'][i][0], self.kinefly_data_separated_on_speed[key][i], t, kind='linear')
				self.kinefly_data_separated_on_speed[key][i] = self.resample(self.kinefly_data_separated_on_speed['timestamp'][i], self.kinefly_data_separated_on_speed[key][i], t, kind='linear')
			
			### still force time to start from 0 though for averaging across trials/flies
			self.kinefly_data_separated_on_speed['timestamp'][i] = t - t[0]
			
			## resample trajectory data to this new evenly spaced time vector
			#self.trajectory_data_separated_on_speed['position'][i] = self.resample(self.trajectory_data_separated_on_speed['timestamp'][i]-self.trajectory_data_separated_on_speed['timestamp'][i][0], self.trajectory_data_separated_on_speed['position'][i], t, kind='linear')
			self.trajectory_data_separated_on_speed['position'][i] = self.resample(self.trajectory_data_separated_on_speed['timestamp'][i], self.trajectory_data_separated_on_speed['position'][i], t, kind='linear')
			self.trajectory_data_separated_on_speed['timestamp'][i] = t - t[0]

			## resample head roll to trajectory timestamp as well
			if self.head_data:
				self.top_camera_data_separated_on_speed['roll'][i] = self.resample(self.top_camera_data_separated_on_speed['timestamp'][i]- self.top_camera_data_separated_on_speed['timestamp'][i][0], self.top_camera_data_separated_on_speed['roll'][i], t, kind='linear')
				self.top_camera_data_separated_on_speed['pitch'][i] = self.resample(self.top_camera_data_separated_on_speed['timestamp'][i]- self.top_camera_data_separated_on_speed['timestamp'][i][0], self.top_camera_data_separated_on_speed['pitch'][i], t, kind='linear')
				
				self.top_camera_data_separated_on_speed['timestamp'][i] = t
			
			#plt.plot(self.kinefly_data_separated_on_speed['timestamp'][i], self.kinefly_data_separated_on_speed['left'][i] - self.kinefly_data_separated_on_speed['right'][i],'b')
			#plt.plot(self.kinefly_data_separated_on_speed['timestamp'][i], self.kinefly_data_separated_on_speed['head'][i],'r')
			#plt.show()

	def resample(self, x1, y1, x2, kind, extrapolate=True):
		# helper function
		if kind == 'spline':
			spline = spi.CubicSpline(x1, y1)
			y2 = spline.__call__(x2, extrapolate=extrapolate)
		else:
			fill_value = "extrapolate" if extrapolate else []
			#interp = spi.interp1d(x1, y1, kind=kind, bounds_error=False, fill_value=fill_value)
			interp = spi.interp1d(x1, y1, kind=kind, bounds_error=True)
			y2 = interp(x2)
			#plt.plot(x1,y1,'o',x2,y2,'-')
			#plt.show()
		return y2


	def remove_jitter(self, threshold):
		# for kinefly L wing and R wing data, detect any outliers (jitter), and replace them by a interpolated value.
		
		for key in ['left','right','head']:
			for j in range(0, self.number_of_speeds):
				print (key + str(j))

				list_to_filter = self.kinefly_data_separated_on_speed[key][j]
				#smooth 1 - prev,next
				for idx in range(1,len(list_to_filter) - 1):
					if abs(list_to_filter[idx] - list_to_filter[idx - 1]) > threshold and abs(list_to_filter[idx] - list_to_filter[idx + 1]) > threshold:
						self.counter += 1
						list_to_filter[idx] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 1]))  

				print ('smooth 1:',self.counter)
				self.counter = 0
				# smooth 2 - prev,next+1 and replace two indices
				for idx in range(1,len(list_to_filter) - 2):
					if abs(list_to_filter[idx] - list_to_filter[idx - 1]) > threshold and abs(list_to_filter[idx] - list_to_filter[idx + 2]) > threshold:
						self.counter += 1
						list_to_filter[idx] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 2]))  
						list_to_filter[idx + 1] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 2]))  
				print ('smooth 2:',self.counter)
				self.counter = 0


				# smooth 3 - prev, next+2
				for idx in range(1,len(list_to_filter) - 3):
					if abs(list_to_filter[idx] - list_to_filter[idx - 1]) > threshold and abs(list_to_filter[idx] - list_to_filter[idx + 3]) > threshold:
						self.counter += 1
						list_to_filter[idx] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 3]))  
						list_to_filter[idx + 1] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 3]))
						list_to_filter[idx + 2] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 3]))  

				print ('smooth 3:',self.counter)
				self.counter = 0
				self.kinefly_data_separated_on_speed[key][j] = np.array(list_to_filter)


	def remove_jitter_baseline(self, threshold):
		# exactly the same as remove_jitter but just for baseline kinefly values
		
		for key in ['left','right','head']:
			for j in range(0, self.number_of_speeds):
				print ('baseline: ' + key + str(j))

				list_to_filter = self.kinefly_data_separated_on_speed_baseline[key][j]
				#smooth 1 - prev,next
				for idx in range(1,len(list_to_filter) - 1):
					if abs(list_to_filter[idx] - list_to_filter[idx - 1]) > threshold and abs(list_to_filter[idx] - list_to_filter[idx + 1]) > threshold:
						self.counter += 1
						list_to_filter[idx] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 1]))  

				print ('smooth 1:',self.counter)
				self.counter = 0
				# smooth 2 - prev,next+1 and replace two indices
				for idx in range(1,len(list_to_filter) - 2):
					if abs(list_to_filter[idx] - list_to_filter[idx - 1]) > threshold and abs(list_to_filter[idx] - list_to_filter[idx + 2]) > threshold:
						self.counter += 1
						list_to_filter[idx] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 2]))  
						list_to_filter[idx + 1] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 2]))  
				print ('smooth 2:',self.counter)
				self.counter = 0


				# smooth 3 - prev, next+2
				for idx in range(1,len(list_to_filter) - 3):
					if abs(list_to_filter[idx] - list_to_filter[idx - 1]) > threshold and abs(list_to_filter[idx] - list_to_filter[idx + 3]) > threshold:
						self.counter += 1
						list_to_filter[idx] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 3]))  
						list_to_filter[idx + 1] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 3]))
						list_to_filter[idx + 2] = np.mean((list_to_filter[idx - 1],list_to_filter[idx + 3]))  

				print ('smooth 3:',self.counter)
				self.counter = 0
				self.kinefly_data_separated_on_speed_baseline[key][j] = np.array(list_to_filter)



	def intersection(self, a1, a2):
		a3 = [value for value in a1 if value in a2]
		return np.array(a3)

	def segregate_based_on_command_timings(self):
		# use the command and command timings to segregate the kinefly and autostep data into different chunks of different periods and baselines
		### remember : autostep only publishes during sinusoid motion ...not while moving to 0 not while moving to sinusoid start..only during periodic motion
		
		for i in range(0,self.number_of_speeds*2, 2):
			# get kinefly and autostep data in the period ranges
			print (self.trajectory_speed_segregation_indices[i],self.trajectory_speed_segregation_indices[i+1])
						
			self.trajectory_data_separated_on_speed['position'].append(self.IMU_data['gyro_z'][self.trajectory_speed_segregation_indices[i]:self.trajectory_speed_segregation_indices[i+1]])
			self.trajectory_data_separated_on_speed['timestamp'].append(self.IMU_data['timestamp'][self.trajectory_speed_segregation_indices[i]:self.trajectory_speed_segregation_indices[i+1]])

			a1 = np.where(self.kinefly_data['timestamp'] >= self.IMU_data['timestamp'][self.trajectory_speed_segregation_indices[i]])[0]
			a2 = np.where(self.kinefly_data['timestamp'] < self.IMU_data['timestamp'][self.trajectory_speed_segregation_indices[i+1]])[0]
			a3 = self.intersection(a1, a2) # these are the timestamps for the first period for autostep and kinefly data..this is the data to use to segregate the kinefly data
			self.a3.append(a3)

			self.kinefly_data_separated_on_speed['left'].append(self.kinefly_data['left'][a3])
			self.kinefly_data_separated_on_speed['right'].append(self.kinefly_data['right'][a3])
			self.kinefly_data_separated_on_speed['head'].append(self.kinefly_data['head'][a3])
			self.kinefly_data_separated_on_speed['aux'].append(self.kinefly_data['aux'][a3])
			self.kinefly_data_separated_on_speed['timestamp'].append(self.kinefly_data['timestamp'][a3])

			if self.head_data:
				h1 = np.where(self.top_camera_data['timestamp'] >= self.IMU_data['timestamp'][self.trajectory_speed_segregation_indices[i]])[0]
				h2 = np.where(self.top_camera_data['timestamp'] <= self.IMU_data['timestamp'][self.trajectory_speed_segregation_indices[i+1]])[0]
				h3 = self.intersection(h1, h2)
				#import ipdb; ipdb.set_trace()
				self.h3.append(h3)
				#import ipdb;ipdb.set_trace()
				self.top_camera_data_separated_on_speed['roll'].append(self.top_camera_data['roll'][h3])
				self.top_camera_data_separated_on_speed['pitch'].append(self.top_camera_data['pitch'][h3])
				
				self.top_camera_data_separated_on_speed['timestamp'].append(self.top_camera_data['timestamp'][h3])

	def segregate_get_kinefly_baselines(self):
		# to get the left and right kinefly baseline values during the 15 second sleep intervals between speeds in order to use this for baseline subtraction
		idx = 0
		for i in range(0, self.number_of_speeds):
			a1 = np.where(self.kinefly_data['timestamp'] >= self.command_data['timestamp'][idx])[0]
			a2 = np.where(self.kinefly_data['timestamp'] < self.command_data['timestamp'][idx + 1])[0]
			a3 = self.intersection(a1, a2) # these are the timestamp indices for the first period for autostep and kinefly data..this is the data to use to segregate the kinefly data

			self.kinefly_data_separated_on_speed_baseline['left'].append(self.kinefly_data['left'][a3])
			self.kinefly_data_separated_on_speed_baseline['right'].append(self.kinefly_data['right'][a3])
			self.kinefly_data_separated_on_speed_baseline['head'].append(self.kinefly_data['head'][a3])
			self.kinefly_data_separated_on_speed_baseline['aux'].append(self.kinefly_data['aux'][a3])
			self.kinefly_data_separated_on_speed_baseline['timestamp'].append(self.kinefly_data['timestamp'][a3])

			a1 = np.where(self.IMU_data['timestamp'] >= self.command_data['timestamp'][idx])[0]
			a2 = np.where(self.IMU_data['timestamp'] < self.command_data['timestamp'][idx + 1])[0]
			a3 = self.intersection(a1, a2) ### indices of timestamps in the IMU data


			self.trajectory_data_separated_on_speed_baseline['position'].append(self.IMU_data['gyro_z'][a3])
			self.trajectory_data_separated_on_speed_baseline['timestamp'].append(self.IMU_data['timestamp'][a3])

			## resample the IMU baseline data and kinefly to a common evenly spaced time vector
			t_baseline = np.arange(max(self.trajectory_data_separated_on_speed_baseline['timestamp'][-1][0], self.kinefly_data_separated_on_speed_baseline['timestamp'][-1][0]), 
				min(self.trajectory_data_separated_on_speed_baseline['timestamp'][-1][-1], self.kinefly_data_separated_on_speed_baseline['timestamp'][-1][-1]), 1/72)

			self.kinefly_data_separated_on_speed_baseline['left'][-1] = self.resample(self.kinefly_data_separated_on_speed_baseline['timestamp'][-1], 
				self.kinefly_data_separated_on_speed_baseline['left'][-1], t_baseline, kind='linear')

			self.kinefly_data_separated_on_speed_baseline['right'][-1] = self.resample(self.kinefly_data_separated_on_speed_baseline['timestamp'][-1], 
				self.kinefly_data_separated_on_speed_baseline['right'][-1], t_baseline, kind='linear')
			self.kinefly_data_separated_on_speed_baseline['head'][-1] = self.resample(self.kinefly_data_separated_on_speed_baseline['timestamp'][-1], 
				self.kinefly_data_separated_on_speed_baseline['head'][-1], t_baseline, kind='linear')

			self.kinefly_data_separated_on_speed_baseline['timestamp'][-1] = t_baseline

			self.trajectory_data_separated_on_speed_baseline['position'][-1] = self.resample(self.trajectory_data_separated_on_speed_baseline['timestamp'][-1], 
				self.trajectory_data_separated_on_speed_baseline['position'][-1], t_baseline, kind='linear')

			self.trajectory_data_separated_on_speed_baseline['timestamp'][-1] = t_baseline




			#if i==0:
			#	plt.clf()
			#	plt.title('baseline L in blue R in green')
			#	plt.plot(range(0,len(self.kinefly_data['left'][a3])), self.kinefly_data['left'][a3], c='b')
			#	plt.plot(range(0,len(self.kinefly_data['right'][a3])), self.kinefly_data['right'][a3], c='g')
			#
			#	plt.show()

			if self.head_data:
				h1 = np.where(self.top_camera_data['timestamp'] >= self.command_data['timestamp'][idx])[0]
				h2 = np.where(self.top_camera_data['timestamp'] <= self.command_data['timestamp'][idx + 1])[0]
				h3 = self.intersection(h1, h2)

				self.top_camera_data_separated_on_speed_baseline['roll'].append(self.top_camera_data['roll'][h3])
				self.top_camera_data_separated_on_speed_baseline['pitch'].append(self.top_camera_data['pitch'][h3])
				
				self.top_camera_data_separated_on_speed_baseline['timestamp'].append(self.top_camera_data['timestamp'][h3])
			
			# because commands are three i.e 'sleeping', period, 'returning'
			idx += 3



	def nan_interp(self, baseline=False):
		# interpolate to remove nan values..I am doing this per speed separately because otherwise there is baselines and speeds and im not sure if the interp will work well
		# only kinefly might have nans and not trajectory data
		for i in range(0, self.number_of_speeds):
			for key in ['left','right','head']:
				y = self.kinefly_data_separated_on_speed[key][i]
				nans, x = self._nan_helper(y)

				if len(nans) == 1:
					continue
				# interpolate the nan values using the non nan values
				y[nans] = np.interp(x(nans), x(~nans), y[~nans])
				self.kinefly_data_separated_on_speed[key][i] = y
				
				# repeat for baseline values
				if baseline:
					y = self.kinefly_data_separated_on_speed_baseline[key][i]
					nans, x = self._nan_helper(y)
					y[nans] = np.interp(x(nans), x(~nans), y[~nans])
					self.kinefly_data_separated_on_speed_baseline[key][i] = y
		        

	def _nan_helper(self, y):
		# helper function
		if sum(np.isnan(y)) > 0:
			print ('NANS are present')
			import ipdb;ipdb.set_trace()
			return np.isnan(y), lambda z: z.nonzero()[0]
		else:
			return [-1],[-1]

	def butterworth_filter(self, order=5.0, cutoff=12.0):
		# low pass filter to remove high freq vibrations
		# for each speed get kf_dt, assume order and cutoff and see how results look

		for i in range(0, self.number_of_speeds):
			# this is after resampling so should be constant
			kf_dt = np.mean(np.diff(self.kinefly_data_separated_on_speed['timestamp'][i])) 
			nyquist = 1./(2 * kf_dt)
			bfilt = sps.butter(order, cutoff/nyquist, btype='lowpass')
			self.kinefly_data_separated_on_speed['left'][i] = sps.filtfilt(bfilt[0], bfilt[1], self.kinefly_data_separated_on_speed['left'][i])
			self.kinefly_data_separated_on_speed['right'][i] = sps.filtfilt(bfilt[0], bfilt[1], self.kinefly_data_separated_on_speed['right'][i])
			self.kinefly_data_separated_on_speed['head'][i] = sps.filtfilt(bfilt[0], bfilt[1], self.kinefly_data_separated_on_speed['head'][i])

			if self.head_data:
				self.top_camera_data_separated_on_speed['roll'][i] = sps.filtfilt(bfilt[0], bfilt[1], self.top_camera_data_separated_on_speed['roll'][i])
				self.top_camera_data_separated_on_speed['pitch'][i] = sps.filtfilt(bfilt[0], bfilt[1], self.top_camera_data_separated_on_speed['pitch'][i])
				

	def fourier_transform_to_get_cutoff(self, signal, t):
		N = signal.shape[0]
		dt = t[2] - t[1]
		fft = 1.0/N * np.fft.fft(signal)
		fft = fft[:N//2]
		fftfreq = np.fft.fftfreq(N, dt)[:N//2] 
		idx = 1 + np.argmax(np.abs(fft[1:]))
		amp = np.abs(fft[idx]) * 2
		phase =  np.angle(fft[idx]) * 180/np.pi
		freq = np.abs(fftfreq[idx])
		## add offset np.mean(signal) just for visualizing the fits better
		signal_reconstucted =  amp*np.cos(2*np.pi*freq*t + phase*np.pi/180) + np.mean(signal)
		plt.clf()
		return signal_reconstucted, amp, phase, signal
		
		

	def subtract_baseline(self):
		# subtract the baselines for kinefly data so as to in a way normalize
		for i in range(0, self.number_of_speeds):
			for key in ['left','right','head']:
				#import ipdb;ipdb.set_trace()
				self.kinefly_data_separated_on_speed[key][i] -= np.mean(self.kinefly_data_separated_on_speed_baseline[key][i])
			if self.head_data:
				self.top_camera_data_separated_on_speed['roll'][i] -= np.mean(self.top_camera_data_separated_on_speed_baseline['roll'][i])
				self.top_camera_data_separated_on_speed['pitch'][i] -= np.mean(self.top_camera_data_separated_on_speed_baseline['pitch'][i])


	def fit_sine_curve(self):
		# use the segregated based on period data...try reconstructing the data 
		self.fourier_transform_to_get_cutoff()
		# highest value in fft
		amp = yf[0]
		#freq = xf[0]
		#phase = 90
		a1 = np.zeros_like(yf)
		a1[0] = amp
		x_ = scipy.fftpack.ifft(yf)
		#import ipdb;ipdb.set_trace()
		
		t = self.kinefly_data_separated_on_speed['timestamp'][-1]
		#x_ = amp * np.sin(2*np.pi*freq*t + phase)
		
		if show_plots:
			plt.plot(t, x_.real)
			plt.show()


	def check_no_fly(self):
		# use the kinefly aux to ignore data in a sinusoid cycle when the fly is not flying
		pass



def run_functions(rka_object):
	rka_object.find_trajectory_speed_change_points()
	rka_object.segregate_based_on_command_timings()
	rka_object.segregate_get_kinefly_baselines()


	#rka_object.visualize_autostep_motion_and_kinefly_data()

	# i am interpolating the nans for each period separately to be more accurate
	rka_object.nan_interp(baseline=True)

	### TODO: commenting out remove_jitter for now but I had run it with it on initially
	rka_object.remove_jitter(5.0)
	rka_object.remove_jitter_baseline(5.0)

	# apply the low pass filter to remove high freq vibration effects
	# for this i need to first fourier transform the signal to determine cutoff freq
	#rka_object.fourier_transform_to_get_cutoff()
	
	### TODO: commenting out butterworth_filter for now but I had run it with it on initially
	rka_object.butterworth_filter(order = 5.0, cutoff = 6)
	
	#rka_object.subtract_baseline()

	## before resampling removal of no fly times using aux should be done
	## resample wings and head roll to trajectory timestamp
	rka_object.resample_kinefly_and_autostep_data()

	#rka_object.visualize_autostep_motion_and_kinefly_data()
	
	#rka_object.fit_sine_curve()
	#rka_object.fourier_transform_to_get_cutoff()




if __name__ == '__main__':
	# calling all the essential functions for 3 separate trials of the same fly in order to average (because this is open loop)
	show_plots = False
	genotype = 'UXJ00yawneg/'
	subdirs = glob.glob('../hdf5/' + genotype +'*')
	
	#subdirs = ['/home/tarun/catkin_ws/src/trajectories_autostep_ros/bagfiles/hdf5/KirJ00yawpos/Nov-1-fly9']
	#subdirs = ['../hdf5/KirJ00yawneg/Nov-3-fly20']
	
	for s in subdirs:
		print ('#######################################################################################')
		print (s)
		print ('#######################################################################################')
		subdir = genotype + s.split('/')[-1] + '/'

		#subdir = 'UXJ96yawneg/Jan-13-fly18/'
		bagfile_names = glob.glob('../hdf5/'+subdir+'*.hdf5')
		#bagfile_names = ['/home/tarun/catkin_ws/src/trajectories_autostep_ros/bagfiles/hdf5/KirJ00yawpos/Nov-1-fly9/final_tests_fly1_2022-11-01-16-26-49_compressed.hdf5']
		

		for i in range(0,len(bagfile_names)):
			bagfile_names[i] = bagfile_names[i].split('/')[-1]

		fly_objects = []

		for i in range(0, len(bagfile_names)):
			fly_trial = RkfAnalysis()
			fly_trial.read_hdf5("../hdf5/"+subdir+ bagfile_names[i])
			run_functions(fly_trial)
			fly_objects.append(fly_trial)

		fly_trial1 = fly_objects[0]



		#### visualize all three trials separately
		# for k,fly in enumerate(fly_objects):
			
		# 	for i in range(0,fly.number_of_speeds):
		# 		figure, axes = plt.subplots()
		# 		axes.set_title('trial ' + str(k+1))
		# 		axes.plot(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.trajectory_data_separated_on_speed['position'][i], 'b')
		# 		axes2 = axes.twinx()
		# 		axes2.plot(fly.kinefly_data_separated_on_speed['timestamp'][i] - fly.kinefly_data_separated_on_speed['timestamp'][i][0], fly.kinefly_data_separated_on_speed['left'][i] - fly.kinefly_data_separated_on_speed['right'][i], c='g', label='left - right')
		# 		#axes2.plot(fly.kinefly_data_separated_on_speed['timestamp'][i] - fly.kinefly_data_separated_on_speed['timestamp'][i][0], fly.kinefly_data_separated_on_speed['head'][i], c='r', label='left - right')
				
		# 		### collapse over the 5 cycles


		# 		optimized_fitted_sin_curve, amp, w_phase, _ = fly_trial1.fourier_transform_to_get_cutoff(fly.kinefly_data_separated_on_speed['left'][i] - fly.kinefly_data_separated_on_speed['right'][i],fly.kinefly_data_separated_on_speed['timestamp'][i] - fly.kinefly_data_separated_on_speed['timestamp'][i][0])
		# 		print ('wing phase:', w_phase)	
		# 		#axes2.plot(fly.kinefly_data_separated_on_speed['timestamp'][i] - fly.kinefly_data_separated_on_speed['timestamp'][i][0], optimized_fitted_sin_curve, 'k')
		# 		#optimized_fitted_sin_curve, amp, h_phase, _ = fly.fourier_transform_to_get_cutoff(fly.kinefly_data_separated_on_speed['head'][i],fly.kinefly_data_separated_on_speed['timestamp'][i] - fly.kinefly_data_separated_on_speed['timestamp'][i][0])
		# 		#axes2.plot(fly.kinefly_data_separated_on_speed['timestamp'][i] - fly.kinefly_data_separated_on_speed['timestamp'][i][0], optimized_fitted_sin_curve, 'y')	
		# 		#print ('head phase:', h_phase)	
		# 		plt.show()

		# try plotting the trajectories of all three flies to see if there is complete and perfect overlap - I did this and there is perfect overlap in trajectories because they are just sinusoid functions

		
		#### average across the three trials...have to resample first..resampling based on trajectory position does not make sense This is because the trajectory is not a linearly increasing thing like time...it is 5 sinusoid cycles so this will not work
		### resampling based on time still makes sense because the trajectories are perfectly aligned I checked.
		dict_for_pickle = {}
		
		
		### resample all to the shortest time list among flies
		for i in range(0, fly_trial1.number_of_speeds):

			#resampling_dt = np.mean(np.diff(fly_trial1.trajectory_data_separated_on_speed['timestamp'][i]))
			#t = np.arange(len(fly_trial1.trajectory_data_separated_on_speed['timestamp'][i])) * resampling_dt
			
			## across the 3 trials, take the shortest time and resample everything to that at 72Hz.

			t = np.arange(0, min([fly_objects[x].trajectory_data_separated_on_speed['timestamp'][i][-1] for x in range(len(fly_objects))]), 1/72)

			average_left_minus_right_list = np.zeros((len(fly_objects),t.shape[0]))
			average_left_list = np.zeros_like(average_left_minus_right_list)
			average_right_list = np.zeros_like(average_left_minus_right_list)
			average_yaw = np.zeros_like(average_left_minus_right_list)
			average_trajectory = np.zeros_like(average_left_minus_right_list)

			############ Here I averaging the baselines of the trials. I am doing this solely for the purpose of making a plot of baseline WBA vs stabilization response #####################
			############ The baseline is just a single value per speed as I average across the baseline period #
			
			#resampling_dt_baseline = np.mean(np.diff(fly_trial1.kinefly_data_separated_on_speed_baseline['timestamp'][i]))
			resampling_dt_baseline = 1/72
			t_baseline = np.arange(0, min([fly_objects[x].kinefly_data_separated_on_speed_baseline['timestamp'][i][-1] - fly_objects[x].kinefly_data_separated_on_speed_baseline['timestamp'][i][0] for x in range(len(fly_objects))]) 
				, resampling_dt_baseline)
			
			average_left_baseline = np.zeros((len(fly_objects),t_baseline.shape[0]))
			average_right_baseline = np.zeros_like(average_left_baseline)
			average_head_baseline = np.zeros_like(average_left_baseline)

			average_trajectory_baseline = np.zeros_like(average_left_baseline)
			
			for k,fly in enumerate(fly_objects):
				#average_left_baseline.append(fly.kinefly_data_separated_on_speed_baseline['left'][i][0:495])
				#average_right_baseline.append(fly.kinefly_data_separated_on_speed_baseline['right'][i][0:495])
				average_left_baseline[k] = fly.resample(fly.kinefly_data_separated_on_speed_baseline['timestamp'][i] - fly.kinefly_data_separated_on_speed_baseline['timestamp'][i][0], fly.kinefly_data_separated_on_speed_baseline['left'][i], t_baseline, kind='linear')
				average_right_baseline[k] = fly.resample(fly.kinefly_data_separated_on_speed_baseline['timestamp'][i] - fly.kinefly_data_separated_on_speed_baseline['timestamp'][i][0], fly.kinefly_data_separated_on_speed_baseline['right'][i], t_baseline, kind='linear')
				average_head_baseline[k] = fly.resample(fly.kinefly_data_separated_on_speed_baseline['timestamp'][i] - fly.kinefly_data_separated_on_speed_baseline['timestamp'][i][0], fly.kinefly_data_separated_on_speed_baseline['head'][i], t_baseline, kind='linear')
				average_trajectory_baseline[k] = fly.resample(fly.trajectory_data_separated_on_speed_baseline['timestamp'][i] - fly.kinefly_data_separated_on_speed_baseline['timestamp'][i][0], fly.trajectory_data_separated_on_speed_baseline['position'][i], t_baseline, kind='linear')
				
			mean_baseline_left_trace = np.mean(average_left_baseline, axis=0)
			mean_baseline_right_trace = np.mean(average_right_baseline, axis=0)
			mean_baseline_head_trace = np.mean(average_head_baseline, axis=0)
			mean_baseline_trajectory = np.mean(average_trajectory_baseline, axis=0)

			if show_plots:
				plt.clf()
				plt.plot(t_baseline, mean_baseline_left_trace, 'b')
				plt.plot(t_baseline, mean_baseline_right_trace, 'g')
				plt.show()
			
			mean_baseline_left = np.mean(mean_baseline_left_trace, axis=0)
			mean_baseline_right = np.mean(mean_baseline_right_trace, axis=0)

			##########################################################################################################################################################################################################

			if fly_trial1.head_data:
				average_roll = np.zeros_like(average_left_minus_right_list)
				average_pitch = np.zeros_like(average_left_minus_right_list)

			for k,fly in enumerate(fly_objects):
				#### for some reason last value is nan..removing last value
				fly.kinefly_data_separated_on_speed['left'][i][-1] = fly.kinefly_data_separated_on_speed['left'][i][-2]
				fly.kinefly_data_separated_on_speed['right'][i][-1] = fly.kinefly_data_separated_on_speed['right'][i][-2]
				fly.kinefly_data_separated_on_speed['head'][i][-1] = fly.kinefly_data_separated_on_speed['head'][i][-2]

				left = fly.resample(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.kinefly_data_separated_on_speed['left'][i], t, kind='linear')
				right = fly.resample(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.kinefly_data_separated_on_speed['right'][i], t, kind='linear')
				head = fly.resample(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.kinefly_data_separated_on_speed['head'][i], t, kind='linear')
				
				trajectory = fly.resample(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.trajectory_data_separated_on_speed['position'][i], t, kind='linear')

				average_left_minus_right_list[k] = left - right
				average_left_list[k] = left
				average_right_list[k] = right
				average_trajectory[k] = trajectory
				
				average_yaw[k] = head
				if fly.head_data:
					roll = fly.resample(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.top_camera_data_separated_on_speed['roll'][i], t, kind='linear')
					pitch = fly.resample(fly.trajectory_data_separated_on_speed['timestamp'][i] - fly.trajectory_data_separated_on_speed['timestamp'][i][0], fly.top_camera_data_separated_on_speed['pitch'][i], t, kind='linear')
					average_roll[k] = roll
					average_pitch[k] = pitch

			#std = np.std(average_left_minus_right_list, axis = 0)
			mean_left_minus_right = np.mean(average_left_minus_right_list, axis=0)
			mean_left = np.mean(average_left_list, axis=0)
			mean_right = np.mean(average_right_list, axis=0)
			mean_yaw = np.mean(average_yaw, axis=0)
			mean_trajectory = np.mean(average_trajectory, axis=0)


			# plt.plot(t, mean_left_minus_right, 'b')
			# plt.plot(t, mean_yaw, 'r')
			# plt.show()

			if fly_trial1.head_data:
				mean_roll = np.mean(average_roll, axis=0)
				mean_pitch = np.mean(average_pitch, axis=0)
				
				#plt.clf()
				#plt.plot(t, mean_roll,c='g')
				#plt.show()

			## get sin fit of mean trace(using fourier transform)
			#fly_trial1.fourier_transform_to_get_cutoff(mean_left_minus_right, t, i)
			#plt.clf()
			#plt.close()


			# plot the average
			if show_plots:
				figure, axes = plt.subplots()
				#import ipdb;ipdb.set_trace()
				axes.set_title('left WBA of ' +genotype[:-1] + ' ' + (subdir[:-1].split('/')[-1]).split('-')[-1] + ' trials max speed = ' + str(round(fly_trial1.list_of_peak_velocities_used[i],2)))
				
				axes.set_ylim(0,80)
				axes.set_ylabel('left WBA')
				axes.set_xlabel('Time')
				colors = ['b', 'g', 'r']

				for k in range(0, len(fly_objects)):
					axes.plot(t, average_left_minus_right_list[k], c=colors[k], alpha=0.5)
					#axes.plot(t, average_left_list[k], c=colors[k], alpha=0.5)

				axes.plot(t, mean_left_minus_right, c='k', label='left - right (avg)')
				#axes.plot(t, mean_left, c='k', label='left (avg)', alpha=0.8)
				
				#axes.fill_between(t, mean-std, mean+std,color='green',alpha=0.2)
				
				## plot each of the trials as well in alpha
				axes2 = axes.twinx()
				axes2.plot(t, fly_trial1.trajectory_data_separated_on_speed['position'][i])
				axes2.set_ylabel('angular velocity')
				#plt.savefig('/home/tarun/Desktop/individual_WBA_plots/'+ genotype[:-1] + ' '+ (subdir[:-1].split('/')[-1]).split('-')[-1] + '_left_speed_'+str(i+1)+'.png')
				plt.clf()
				plt.close(figure)
				#plt.show()

			## along with the plot also save the values to a pickle file so that it can be used to average across flies..also save trajectory to pickle file for average plot
			if fly_trial1.head_data:
				dict_for_pickle[i] = [t, 1.0*mean_left_minus_right, 1.0*mean_roll, 1.0*mean_pitch, 1.0*mean_yaw]
			else:
				dict_for_pickle[i] = [t, 1.0*mean_left_minus_right, 1.0*mean_yaw]


			### trajectory position and timestamp and exactly same for all trials, I checked. Hence only saving it from the first trial.
			dict_for_pickle['trajectory_'+str(i)] = [1.0*mean_trajectory]
			dict_for_pickle['trajectory_'+str(i)+'_time'] = [t]

			dict_for_pickle['mean_baseline_left_'+str(i)] = mean_baseline_left
			dict_for_pickle['mean_baseline_right_'+str(i)] = mean_baseline_right

			dict_for_pickle['mean_baseline_left_trace_'+str(i)] = mean_baseline_left_trace
			dict_for_pickle['mean_baseline_right_trace_'+str(i)] = mean_baseline_right_trace
			dict_for_pickle['mean_baseline_head_trace_'+str(i)] = mean_baseline_head_trace
			dict_for_pickle['baseline_time_'+str(i)] = t_baseline
			dict_for_pickle['mean_baseline_trajectory_' + str(i)] = mean_baseline_trajectory


			dict_for_pickle['mean_left_'+str(i)] = 1.0*mean_left
			dict_for_pickle['mean_right_'+str(i)] = 1.0*mean_right



		pickle.dump(dict_for_pickle, open('../hdf5/'+subdir+'avg_left_minus_right_'+subdir.split('/')[-2] +'.p', 'wb'))
		