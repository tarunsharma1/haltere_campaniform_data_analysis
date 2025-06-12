import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib as mpl
import statsmodels.api as sm
def f(x, A, B): # this is your 'straight line' y=f(x)
	return A*x + B


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


parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'ytick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8}
plt.rcParams.update(parameters)




all_slopes = []

wba_all_genotypes = []
head_all_genotypes = []
campaniforms_silenced = []

wba_controls = []
head_controls = []

## per genotype, make plots of wing vs head as scatter plot dots
#list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UX28C05yawbothdirections', 'UXJ88yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']
#genotype_names = [ 'UXJ00', 'UX22E04', 'UX28C05' , 'UX60B12', 'UX31A09', 'UX74B09', 'UX86H12']


#list_of_genotypes = ['KirJ00yawbothdirections','KirJ73yawbothdirections' ,'KirJ86yawbothdirections', 'KirJ79yawbothdirections', 'KirJ96yawbothdirections', 'KirJ70yawbothdirections', 'Kir28C05yawbothdirections', 'KirJ88yawbothdirections', 'KirJ75yawbothdirections']
#genotype_names = ['KirJ00','Kir17G01', 'Kir58F02', 'Kir31A09', 'Kir86H12', 'Kir14B04', 'Kir28C05', 'Kir60B12', 'Kir22E04']

genotype_mapping = {'KirJ00yawbothdirections':'J00','KirJ73yawbothdirections': '17G01' ,'KirJ86yawbothdirections': '58F02', 'KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ70yawbothdirections':'14B04', 'Kir28C05yawbothdirections':'28C05', 'KirJ88yawbothdirections':'60B12', 'KirJ75yawbothdirections':'22E04'}
#genotype_mapping = {'UXJ00yawbothdirections':'J00', 'UX28C05yawbothdirections':'28C05', 'UXJ88yawbothdirections':'60B12', 'UXJ79yawredobothdirections':'31A09', 'UXJ90yawredobothdirections':'74B09', 'UXJ75yawbothdirections':'22E04', 'UXJ96yawbothdirections':'86H12'}

campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}



#list_of_genotypes = [ 'KirJ00yawbothdirections']
list_of_genotypes = ['KirJ00yawbothdirections', 'KirJ73yawbothdirections' ,'KirJ86yawbothdirections', 'KirJ79yawbothdirections', 'KirJ96yawbothdirections', 'KirJ70yawbothdirections', 'Kir28C05yawbothdirections', 'KirJ88yawbothdirections', 'KirJ75yawbothdirections']

#list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UX28C05yawbothdirections', 'UXJ88yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']

for genotype in list_of_genotypes:
    speeds = ['2']
    for speed in speeds:
        wba = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
        head = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'-head-yaw.p', 'rb'))
        if genotype == 'KirJ00yawbothdirections':
            wba_controls.extend(wba)
            head_controls.extend(head)
        else:
            wba_all_genotypes.extend(wba)
            head_all_genotypes.extend(head)
            #import ipdb;ipdb.set_trace()
            campaniforms_silenced.extend([campaniform_fields[genotype_mapping[genotype]]['dscab'] + campaniform_fields[genotype_mapping[genotype]]['vped'] + campaniform_fields[genotype_mapping[genotype]]['dped'] + campaniform_fields[genotype_mapping[genotype]]['lcho'] + campaniform_fields[genotype_mapping[genotype]]['vscab']]*len(wba))

    
    slope, intercept, r, p, se = stats.linregress(head, wba)
    all_slopes.append(slope)
    slope = round(slope, 3)
    intercept = round(intercept, 3)
    r = round(r, 3)
    #p = round(p, 3)
    se = round(se, 3)
    print ('linregress: ', slope,intercept, r, p, se, genotype_mapping[genotype])


head_all_genotypes = np.array(head_all_genotypes)
wba_all_genotypes = np.array(wba_all_genotypes)

head_controls = np.array(head_controls)
wba_controls = np.array(wba_controls)

print ('pearson r:', stats.pearsonr(head_all_genotypes, wba_all_genotypes))

A2 = sm.add_constant(head_all_genotypes)
#A2 = head_all_genotypes
model = sm.OLS(wba_all_genotypes,A2, hasconst=True)
multi_model = model.fit()

print (multi_model.summary())
print (multi_model.conf_int(alpha=0.05, cols=None))


fig, axes = plt.subplots()
fig.set_size_inches(2.5,1.5)

adjust_spines(axes,['bottom','left'], xticks=[0,15], yticks=[10,20,30,40], linewidth=0.6,spineColor='k')

#plt.scatter(head_all_genotypes, wba_all_genotypes, c=np.array(campaniforms_silenced), cmap='plasma', s=7)
#plt.colorbar()

plt.scatter(head_all_genotypes, wba_all_genotypes, c='gray', s=7)



axes.set_xlim(0,17)
axes.set_ylim(0,45)
axes.set_xlabel('Head response \n (degs)')
axes.set_ylabel('Wing response \n (degs)')



##extend the line a little bit
#head_all_genotypes = np.concatenate((head_all_genotypes, np.array([2,4,5,6,7,15])))
axes.plot(head_all_genotypes, head_all_genotypes*multi_model.params[1] + multi_model.params[0], c=[0.5,0.5,0.5])

pred = multi_model.get_prediction(A2).summary_frame()
sorted_indices = np.argsort(head_all_genotypes)
axes.fill_between(head_all_genotypes[sorted_indices], pred['mean_ci_lower'][sorted_indices], pred['mean_ci_upper'][sorted_indices], color='gray', alpha=0.3)


####### only controls on the same plot in a different color ###
plt.scatter(head_controls, wba_controls, c='black', s=7)

print('### CONTROLS ###')
A2 = sm.add_constant(head_controls)

model_controls = sm.OLS(wba_controls,A2, hasconst=True)
multi_model_controls = model_controls.fit()

print (multi_model_controls.summary())
print ('conf interval:')
print (multi_model_controls.conf_int(alpha=0.05, cols=None))

axes.plot(head_controls, head_controls*multi_model_controls.params[1] + multi_model_controls.params[0], c=[0,0,0])

pred = multi_model_controls.get_prediction(A2).summary_frame()
sorted_indices = np.argsort(head_controls)
axes.fill_between(head_controls[sorted_indices], pred['mean_ci_lower'][sorted_indices], pred['mean_ci_upper'][sorted_indices], color='black', alpha=0.3)

print ('pearson r:', stats.pearsonr(head_controls, wba_controls))

plt.show()

