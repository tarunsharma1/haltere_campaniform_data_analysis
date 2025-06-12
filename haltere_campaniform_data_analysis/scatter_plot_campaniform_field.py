## make scatter plot of flies across genotypes
import pickle
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy.stats as st

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


parameters = {'axes.labelsize':15, 'axes.titlesize':18, 'xtick.labelsize':12}
plt.rcParams.update(parameters)



############# cells silenced from cartoons ###############
# campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':1, 'vscab':0, 'vped':7, 'lcho':0} ,'60B12':{'dscab':0, 'dped':5, 'vscab':3, 'vped':5, 'lcho':7},
#  '31A09':{'dscab':6, 'dped':16, 'vscab':0, 'vped':11, 'lcho':7}, '74B09':{'dscab':16, 'dped':9, 'vscab':0, 'vped':9, 'lcho':0}, '58F02':{'dscab':16, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':12, 'dped':10, 'vscab':0, 'vped':7, 'lcho':0}, '14B04':{'dscab':4, 'dped':2, 'vscab':2, 'vped':0, 'lcho':0},
#  '86H12':{'dscab':7, 'dped':6, 'vscab':5, 'vped':13, 'lcho':0}, '17G01':{'dscab':3, 'dped':7, 'vscab':1, 'vped':3, 'lcho':0}}

## cells remaining
#campaniform_fields = {'UXS00yawbothdirections':{'dscab':42, 'dped':43, 'vscab':5, 'vped':46, 'lcho':7}, 'UX28C05yawbothdirections':{'dscab':40, 'dped':42, 'vscab':5, 'vped':39, 'lcho':7},
#'UXJ88yawbothdirections':{'dscab':42, 'dped':38, 'vscab':2, 'vped':41, 'lcho':0}, 'UXJ79yawredobothdirections':{'dscab':36, 'dped':27, 'vscab':5, 'vped':35, 'lcho':0}, 
#'UXJ90yawredobothdirections': {'dscab':26, 'dped':34, 'vscab':5, 'vped':37, 'lcho':7}}



#### cells silenced average of mine and annes counts ####
campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}

## cells remaining
# campaniform_fields = {'J00':{'dscab':42, 'dped':43, 'vscab':5, 'vped':46, 'lcho':12},'28C05':{'dscab':41, 'dped':42, 'vscab':5, 'vped':41, 'lcho':12} ,'60B12':{'dscab':42, 'dped':42, 'vscab':5, 'vped':40, 'lcho':0},
#  '31A09':{'dscab':38, 'dped':27, 'vscab':5, 'vped':35, 'lcho':3}, '74B09':{'dscab':31, 'dped':31, 'vscab':5, 'vped':36, 'lcho':12}, '58F02':{'dscab':27, 'dped':36, 'vscab':5, 'vped':38, 'lcho':12}, '22E04':{'dscab':32, 'dped':32, 'vscab':5, 'vped':36, 'lcho':12}, '14B04':{'dscab':40, 'dped':42, 'vscab':4, 'vped':46, 'lcho':12},
#  '86H12':{'dscab':35, 'dped':39, 'vscab':0, 'vped':37, 'lcho':12}, '17G01':{'dscab':42, 'dped':34, 'vscab':5, 'vped':45, 'lcho':6}}



#genotype_mapping = {'UXJ00yawbothdirections':'J00', 'UX28C05yawbothdirections':'28C05', 'UXJ88yawbothdirections':'60B12', 'UXJ79yawredobothdirections':'31A09', 'UXJ90yawredobothdirections':'74B09', 'UXJ75yawbothdirections':'22E04', 'UXJ96yawbothdirections':'86H12'}
genotype_mapping = {'KirJ00yawbothdirections':'J00','KirJ73yawbothdirections': '17G01' ,'KirJ86yawbothdirections': '58F02', 'KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ70yawbothdirections':'14B04', 'Kir28C05yawbothdirections':'28C05', 'KirJ88yawbothdirections':'60B12', 'KirJ75yawbothdirections':'22E04'}

#genotype_mapping = {'KirJ00yawbothdirections':'J00','KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ75yawbothdirections':'22E04'}

field = 'vscab'

for speed in ['1','2']:
	x_list = []
	y_list = []
	x_labels = []
	y_conf_min = []
	y_conf_max = []
	for k,genotype in enumerate(genotype_mapping.keys()):

		y_list_temp = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
		#print (genotype)
		#print (y_list_temp)
		bootstrapped_data = bootstrap(y_list_temp, 10000)
		conf_min, conf_max = confidence_interval(bootstrapped_data)
		y_conf_min.append(conf_min)
		y_conf_max.append(conf_max)

		## number of cells in a field of interest
		x_list.append([campaniform_fields[genotype_mapping[genotype]][field]] * len(y_list_temp))
		y_list.append(y_list_temp)
		
		x_labels.append(genotype_mapping[genotype])
		#import ipdb;ipdb.set_trace()

	#plt.xlim(-5,20)
	for i in range(0,len(y_list)):
		 
		plt.plot(x_list[i], y_list[i], 'o', mfc='none')
		plt.plot(x_list[i][0], np.mean(y_list[i]),'o')

	plt.xlabel('Number of cells silenced in ' + field)
	plt.ylabel('Amplitude of fit to L-R WBA')

	plt.ylim(0, 40)
	#plt.show()
	#import ipdb;ipdb.set_trace()

	x = np.array([item for sublist in x_list for item in sublist])
	y = np.array([item for sublist in y_list for item in sublist])
	a,b = np.polyfit(x,y,1)
	plt.plot(x, a*x + b)
	corr, pvalue = st.pearsonr(x,y)
	plt.title('Pearson R:'+ str(round(corr,2)) + ' p value:' + str(round(pvalue,5)))
	plt.show()

	print (st.pearsonr(x,y))
	


	####################### IDEA ###################
	## using the stabilization magnitude as Y values and the number of cells in each field as features, do a linear regression fit to see what fields are significant for wing and for head ##




	#################################################


	## plot the means and conf intervals
	#plt.scatter(range(0,len(list_of_genotypes)), [np.mean(y) for y in y_list],c='r')
	
	## convert the confidence interval values into mean - value so as to use matplotlib errbar function
	y_err_lower = []
	y_err_upper = []
	for i,y in enumerate(y_list):
		mean = np.mean(y)
		y_err_lower.append(mean - y_conf_min[i])
		y_err_upper.append(y_conf_max[i] - mean)
	
	plt.errorbar(range(0,len(genotype_mapping.keys())), [np.mean(y) for y in y_list], yerr=[y_err_lower, y_err_upper], fmt="o", c='r')

			
	plt.xticks(range(0,len(genotype_mapping.keys())), x_labels, rotation=25)
	plt.ylim(0,40)
	plt.ylabel('Amplitude of sinusoid fits to average stabilization response')
	plt.xlabel('genotypes')
	plt.title(' UX speed '+ speed)
	#plt.show()
	plt.clf()