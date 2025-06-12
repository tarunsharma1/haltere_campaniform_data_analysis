import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import itertools
import scipy.stats as st
import random


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

#pickle_files = glob.glob('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_baseline_correlation_plot/*.p')

# list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ88yawbothdirections', 'UX28C05yawbothdirections',  'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']

# pickle_files = ['./data_for_baseline_correlation_plot_new/' + f + '-speed-2.p' for f in list_of_genotypes]

# master_baseline_left = []
# master_baseline_right = []
# master_amplitudes = []



# for pickle_file in pickle_files:
	
# 	temp_amplitudes, temp_baseline_left, temp_baseline_right = pickle.load(open(pickle_file, "rb"))
# 	print (pickle_file)
# 	print (temp_amplitudes)
# 	for i in range(0, len(temp_amplitudes)):
# 		master_amplitudes.append(temp_amplitudes[i])
# 		master_baseline_left.append(temp_baseline_left[i])
# 		master_baseline_right.append(temp_baseline_right[i])
# 		plt.scatter(temp_baseline_left[i] + temp_baseline_right[i],temp_amplitudes[i], c=color)

# #plt.plot(master_baseline_left, master_amplitudes, 'o')
# #plt.xlim(0,150)
# #plt.ylim(-10,50)
# plt.xlabel('Baseline L+R WBA')
# plt.ylabel('Amplitude of sinusoid fits to wing stabilization responses')
# plt.title('Baseline L+R WBA vs stabilization response for all flies')
# plt.show()



genotype_names = ['Control','R14B04', 'R17G01', 'R28C05', 'R60B12','R58F02', 'R31A09','R22E04','R86H12']
list_of_genotypes = ['KirJ00yawbothdirections', 'KirJ70yawbothdirections', 'KirJ73yawbothdirections', 'Kir28C05yawbothdirections','KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections',  'KirJ75yawbothdirections', 'KirJ96yawbothdirections']


# list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ88yawbothdirections', 'UX28C05yawbothdirections',  'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']
# genotype_names = [ 'Control', 'R22E04','R31A09','R60B12','R28C05','R74B09','R86H12']



for speed in ['2']:
	x_list = []
	y_list = []
	x_list_left = []
	y_list_left = []
	x_list_right = []
	y_list_right = []
	
	x_labels = []
	y_conf_min = []
	y_conf_max = []
	for k,genotype in enumerate(list_of_genotypes):

		temp_amplitudes, temp_baseline_left, temp_baseline_right = pickle.load(open('./data_for_baseline_correlation_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
		

		y_list_temp = [temp_baseline_left[i] + temp_baseline_right[i] for i in range(0,len(temp_baseline_left))]
		bootstrapped_data = bootstrap(y_list_temp, 10000)
		conf_min, conf_max = confidence_interval(bootstrapped_data)
		y_conf_min.append(conf_min)
		y_conf_max.append(conf_max)

		x_list.append([k]*len(y_list_temp))
		y_list.append(y_list_temp)
		x_labels.append(genotype_names[k])

		# for j,b in enumerate(['L', 'R']):
		# 	if b=='L':		
		# 		y_list_temp = [temp_baseline_left[i] for i in range(0,len(temp_baseline_left))]
		# 		#y_list_temp = [(temp_amplitudes[i]*temp_baseline[i]) for i in range(0,len(temp_baseline))]
		# 		#print (genotype)
		# 		#print ([round(i,2) for i in y_list_temp])
		# 		bootstrapped_data = bootstrap(y_list_temp, 10000)
		# 		conf_min, conf_max = confidence_interval(bootstrapped_data)
		# 		y_conf_min.append(conf_min)
		# 		y_conf_max.append(conf_max)

		# 		x_list_left.append([(k+j/2.0)*2.0]*len(y_list_temp))
		# 		y_list_left.append(y_list_temp)
		# 		x_labels.append(genotype_names[k]+'_'+b)
		# 	else:
		# 		y_list_temp = [temp_baseline_right[i] for i in range(0,len(temp_baseline_right))]
		# 		#y_list_temp = [(temp_amplitudes[i]*temp_baseline[i]) for i in range(0,len(temp_baseline))]
		# 		#print (genotype)
		# 		#print ([round(i,2) for i in y_list_temp])
		# 		bootstrapped_data = bootstrap(y_list_temp, 10000)
		# 		conf_min, conf_max = confidence_interval(bootstrapped_data)
		# 		y_conf_min.append(conf_min)
		# 		y_conf_max.append(conf_max)

		# 		x_list_right.append([(k+j/2.0)*2.0]*len(y_list_temp))
		# 		y_list_right.append(y_list_temp)
		# 		x_labels.append(genotype_names[k]+'_'+b)
			

	
	### plot (L-R)baseline on X axis and (L-R) signal vertical offset on y axis
	#L_minus_R_baseline = [y_list_left[0][i] - y_list_right[0][i] for i in range(0,len(y_list_left[0]))]

	#plt.plot(L_minus_R_baseline, temp_amplitudes,'o')
	#plt.xlabel('Mean of left baseline - mean of right baseline')
	#plt.ylabel('Vertical offset in L-R stabilization response')
	#plt.show()

	plt.scatter(list(itertools.chain.from_iterable(x_list)), list(itertools.chain.from_iterable(y_list)))
	## convert the confidence interval values into mean - value so as to use matplotlib errbar function
	y_err_lower = []
	y_err_upper = []
	for i,y in enumerate(y_list):
		mean = np.mean(y)
		y_err_lower.append(mean - y_conf_min[i])
		y_err_upper.append(y_conf_max[i] - mean)
	
	plt.errorbar(range(0,len(list_of_genotypes)), [np.mean(y) for y in y_list], yerr=[y_err_lower, y_err_upper], fmt="o", c='r')

			
	plt.xticks(range(0,len(list_of_genotypes)), x_labels, rotation=25)


	
	#plt.scatter(list(itertools.chain.from_iterable(x_list_left)),list(itertools.chain.from_iterable(y_list_left)),c='b')
	#plt.scatter(list(itertools.chain.from_iterable(x_list_right)),list(itertools.chain.from_iterable(y_list_right)),c='g')

	## draw line connecting each L and R point
	#for x in range(0, len(y_list_left)):
	#	plt.plot([x_list_left[x], x_list_right[x]], [y_list_left[x], y_list_right[x]], linewidth=1, color='grey')

	
	## plot the means and conf intervals
	#plt.scatter(range(0,len(list_of_genotypes)), [np.mean(y) for y in y_list],c='r')
	
	## convert the confidence interval values into mean - value so as to use matplotlib errbar function
	for y_list in [y_list_left, y_list_right]:
		y_err_lower = []
		y_err_upper = []
		for i,y in enumerate(y_list_left):
			mean = np.mean(y)
			y_err_lower.append(mean - y_conf_min[i])
			y_err_upper.append(y_conf_max[i] - mean)
		
		#plt.errorbar(range(0,len(list_of_genotypes)), [np.mean(y) for y in y_list], yerr=[y_err_lower, y_err_upper], fmt="o", c='r')

			
	#plt.xticks(range(0,2*len(list_of_genotypes)), x_labels, rotation=25)
	#plt.ylim(0,40)
	plt.ylabel('Baseline L and R WBA')
	plt.xlabel('genotypes')
	plt.title(' Baseline L and R WBA per genotype ')
	#plt.ylim(0,80)
	plt.show()
	plt.clf()