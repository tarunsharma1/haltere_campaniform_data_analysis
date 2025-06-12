import pickle
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy.stats as st
import matplotlib as mpl
from scipy import stats
import matplotlib as mpl
import statsmodels.api as sm
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
    if 'top' in spines or 'bottom' in spines:
        ax.set_xticks(xticks)
    
    for line in ax.get_xticklines() + ax.get_yticklines():
        #line.set_markersize(6)
        line.set_markeredgewidth(linewidth)




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




parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8}
plt.rcParams.update(parameters)

#genotype_names = [ 'Control', '22E04','31A09','60B12','28C05','74B09','86H12']
#list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ88yawbothdirections', 'UX28C05yawbothdirections',  'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']
#genotype_mapping = {'UXJ00yawbothdirections':'Control', 'UX28C05yawbothdirections':'28C05', 'UXJ88yawbothdirections':'60B12', 'UXJ79yawredobothdirections':'31A09', 'UXJ90yawredobothdirections':'74B09', 'UXJ75yawbothdirections':'22E04', 'UXJ96yawbothdirections':'86H12'}

list_of_genotypes = ['KirJ00yawbothdirections','KirJ70yawbothdirections','KirJ73yawbothdirections', 'Kir28C05yawbothdirections', 'KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections', 'KirJ75yawbothdirections', 'KirJ96yawbothdirections']
genotype_names = ['KirJ00','Kir17G01', 'Kir58F02', 'Kir31A09', 'Kir86H12', 'Kir14B04', 'Kir28C05', 'Kir60B12', 'Kir22E04']
genotype_mapping = {'KirJ00yawbothdirections':'Control','KirJ73yawbothdirections': '17G01' ,'KirJ86yawbothdirections': '58F02', 'KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ70yawbothdirections':'14B04', 'Kir28C05yawbothdirections':'28C05', 'KirJ88yawbothdirections':'60B12', 'KirJ75yawbothdirections':'22E04'}


campaniform_fields = {'Control':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}



for speed in ['2']:
	x_list = []
	y_list = []
	x_labels = []
	y_conf_min = []
	y_conf_max = []
	bootstrapped_means = []
	for k,genotype in enumerate(list(genotype_mapping.keys())):

		y_list_temp_wing = pickle.load(open('./data_for_scatter_plot/'+ genotype + '-speed-' + speed +'.p', 'rb'))
		y_list_temp_head = pickle.load(open('./data_for_scatter_plot/'+ genotype + '-speed-' + speed +'-head-yaw.p', 'rb'))
		
		y_list_temp = np.array(y_list_temp_wing)
		#y_list_temp = np.array(y_list_temp_head)


		#### BASELINE 
		#y_list_baseline = pickle.load(open('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_baseline_correlation_plot/'+ genotype + '-speed-' + speed +'.p', 'rb'))
		#y_list_temp = np.array(y_list_baseline[1]) + np.array(y_list_baseline[2]) ## L + R
		
		print (genotype)
		#print ([round(i,2) for i in y_list_temp])
		bootstrapped_data = bootstrap(y_list_temp, 10000)
		conf_min, conf_max = confidence_interval(bootstrapped_data)
		y_conf_min.append(conf_min)
		y_conf_max.append(conf_max)
		number_of_campaniform_silenced = campaniform_fields[genotype_mapping[genotype]]['dscab'] + campaniform_fields[genotype_mapping[genotype]]['vped'] + campaniform_fields[genotype_mapping[genotype]]['dped'] + campaniform_fields[genotype_mapping[genotype]]['lcho'] + campaniform_fields[genotype_mapping[genotype]]['vscab']
		x_list.extend([number_of_campaniform_silenced]*len(y_list_temp))
		y_list.extend(y_list_temp)
		bootstrapped_means.append(np.mean(bootstrapped_data))

		x_labels.append(genotype_names[k])
		

	plt.plot(x_list, y_list, 'o')
	plt.show()
	slope, intercept, r, p, se = stats.linregress(x_list, y_list)
	slope = round(slope, 3)
	intercept = round(intercept, 3)
	r = round(r, 3)
	#p = round(p, 3)
	se = round(se, 3)
	print ('linregress: ', slope,intercept, r, p, se)
	continue



	#print (y_error_min, y_error_max)
	
	fig, axes = plt.subplots()
	fig.set_size_inches(2.2,1.5)
	#fig.set_size_inches(0.2,0.5)
	
	x_axis = list(itertools.chain.from_iterable(x_list))
	## add some jitter for visaulization purposes
	
	#x_axis = x_axis + np.random.uniform(low=-0.11, high=0.11, size=len(x_axis))
	
	axes.scatter(x_axis, list(itertools.chain.from_iterable(y_list)), marker='.', c=[[0.5,0.5,0.5]],s=7)

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
	adjust_spines(axes,['bottom', 'left'],xticks=range(0,len(list_of_genotypes)), yticks=[0,5], linewidth=0.6,spineColor='k')
	
	#axes.set_xticks(range(0,len(list_of_genotypes)), x_labels)
	axes.set_xticklabels(x_labels, rotation=90)
	#axes.set_xlim(-0.12,0.12)
	#axes.set_ylim(0,70)
	#plt.ylim(0,4)
	axes.set_ylabel('Wing \n response \n (degs)')
	#plt.ylabel('Ratio of sinusoid fits to wing and head yaw')
	#plt.xlabel('genotypes',fontsize=8, font='Arial')
	#plt.title('  Baseline left + right of wba UX  ')
	plt.show()
	plt.clf()