## make scatter plot of flies across genotypes
import pickle
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy.stats as st
import matplotlib as mpl

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




parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8, 'svg.fonttype':'none'}
plt.rcParams.update(parameters)







### sorted according to WBA ####
genotype_names = [ 'Control', 'R22E04','R31A09','R60B12','R28C05','R74B09','R86H12']
list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ88yawbothdirections', 'UX28C05yawbothdirections',  'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']


### sorted according to WBA ##
# genotype_names = ['Control','R14B04', 'R17G01', 'R28C05', 'R60B12','R58F02', 'R31A09','R22E04','R86H12']
# list_of_genotypes = ['KirJ00yawbothdirections', 'KirJ70yawbothdirections', 'KirJ73yawbothdirections', 'Kir28C05yawbothdirections','KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections',  'KirJ75yawbothdirections', 'KirJ96yawbothdirections']


#genotype_names = ['Wing response']
#list_of_genotypes = ['KirJ00yawbothdirections']


if __name__ == '__main__':

	for speed in ['2']:
		x_list = []
		y_list = []
		x_labels = []
		y_conf_min = []
		y_conf_max = []
		bootstrapped_means = []
		for k,genotype in enumerate(list_of_genotypes):

			y_list_temp_wing = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
			y_list_temp_head = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'-head-yaw.p', 'rb'))
			
			y_list_temp = np.array(y_list_temp_head)
			#y_list_temp = np.array(y_list_temp_wing) - np.array(y_list_temp_head)


			#### BASELINE 
			#y_list_baseline = pickle.load(open('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_baseline_correlation_plot/'+ genotype + '-speed-' + speed +'.p', 'rb'))
			#y_list_temp = np.array(y_list_baseline[1]) + np.array(y_list_baseline[2]) ## L + R
			
			print (genotype)
			#print ([round(i,2) for i in y_list_temp])
			bootstrapped_data = bootstrap(y_list_temp, 10000)
			conf_min, conf_max = confidence_interval(bootstrapped_data)
			y_conf_min.append(conf_min)
			y_conf_max.append(conf_max)

			x_list.append([k]*len(y_list_temp))
			y_list.append(y_list_temp)
			bootstrapped_means.append(np.mean(bootstrapped_data))

			x_labels.append(genotype_names[k])
			
		print (list_of_genotypes)
		print (bootstrapped_means)
		print (y_conf_min)
		print (y_conf_max)
		#print (y_error_min, y_error_max)
		
		fig, axes = plt.subplots()
		#fig.set_size_inches(2.2,1.5)
		fig.set_size_inches(2.0,1.3)
		
		#fig.set_size_inches(0.2,0.5)
		
		x_axis = list(itertools.chain.from_iterable(x_list))
		## add some jitter for visaulization purposes
		
		x_axis = x_axis + np.random.uniform(low=-0.11, high=0.11, size=len(x_axis))
		print ('max response value:' + str(max(list(itertools.chain.from_iterable(y_list))) ))
		
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

				
		adjust_spines(axes,['bottom', 'left'],xticks=range(0,len(list_of_genotypes)), yticks=[8,14], linewidth=0.6,spineColor='k')
		#adjust_spines(axes,['left'], yticks=[20,30], linewidth=0.6,spineColor='k')
		#adjust_spines(axes,['bottom', 'left'],xticks=range(0,len(list_of_genotypes)), yticks=[8, 12], linewidth=0.6,spineColor='k')
		#adjust_spines(axes,['bottom', 'left'],xticks=range(0,len(list_of_genotypes)), yticks=[20, 30], linewidth=0.6,spineColor='k')
		
		#axes.set_xticks(range(0,len(list_of_genotypes)), x_labels)
		axes.set_xticklabels(x_labels, rotation=90)
		#axes.set_xlim(-0.12,0.12)
		axes.set_ylim(0,18)
		#plt.ylim(0,4)
		axes.set_ylabel('Head response \n (degs)')
		#plt.ylabel('Ratio of sinusoid fits to wing and head yaw')
		#plt.xlabel('genotypes',fontsize=8, font='Arial')
		#plt.title('  Baseline left + right of wba UX  ')
		plt.show()
		plt.clf()