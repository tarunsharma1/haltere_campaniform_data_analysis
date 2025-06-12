import numpy as np
import pickle
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.compat import lzip
import statsmodels.stats.api as sms
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats
import matplotlib as mpl
import scipy.stats

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



def bootstrap_slopes(X,Y, n=10000):
    bootstrapped_data = np.zeros(n)
    indices = list(range(0,len(X)-1))
    for i in range(0,n):
        sample = np.random.choice(indices, size=len(indices))
        slope, intercept, r, p, se = stats.linregress(num_campaniforms_silenced[sample], response[sample])
        bootstrapped_data[i] = slope
    print ('bootstrapped slope')
    print (np.mean(bootstrapped_data))
    conf_interval = np.percentile(bootstrapped_data,[2.5,97.5])
    print (conf_interval)




parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'ytick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8}
plt.rcParams.update(parameters)


genotype_mapping = {'KirJ00yawbothdirections':'J00','KirJ73yawbothdirections': '17G01' ,'KirJ86yawbothdirections': '58F02', 'KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ70yawbothdirections':'14B04', 'Kir28C05yawbothdirections':'28C05', 'KirJ88yawbothdirections':'60B12', 'KirJ75yawbothdirections':'22E04'}
#genotype_mapping = {'UXJ00yawbothdirections':'J00', 'UX28C05yawbothdirections':'28C05', 'UXJ88yawbothdirections':'60B12', 'UXJ79yawredobothdirections':'31A09', 'UXJ90yawredobothdirections':'74B09', 'UXJ75yawbothdirections':'22E04', 'UXJ96yawbothdirections':'86H12'}



### cells silenced annes counts ###
# campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':1,'dped':1, 'vscab':0, 'vped':5, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':0, 'vped':6, 'lcho':12},
# '31A09':{'dscab':4, 'dped':16, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':11, 'dped':12, 'vscab':0, 'vped':10, 'lcho':0}, '58F02':{'dscab': 15, 'dped':7, 'vscab':0, 'vped':8, 'lcho':0}, '22E04':{'dscab':10, 'dped':11, 'vscab':0, 'vped':10, 'lcho':0}, '14B04':{'dscab':2, 'dped':1, 'vscab':1, 'vped':0, 'lcho':0},
# '86H12':{'dscab':7, 'dped':4, 'vscab':5, 'vped':9, 'lcho':0}, '17G01':{'dscab':0, 'dped':9, 'vscab':0, 'vped':1, 'lcho':6}}

#### cells silenced average of mine and annes counts ####
campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}


speed = '2'


num_campaniforms_silenced = []
response = []

bootstrapped_mean_response = []
unique_num_campaniforms_silenced = []

fig, axes = plt.subplots()
fig.set_size_inches(2.5,1.5)


for k,genotype in enumerate(list(genotype_mapping.keys())):
    print (genotype_mapping[genotype])
	
    y_list_wing = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'-phases.p', 'rb'))
    y_list_head = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed + '-head-yaw-phases.p', 'rb'))
    #y_list_head = [x for x in y_list_head if x <50]
    y_list_temp = np.array(y_list_wing)
    print (f'number of flies for {genotype_mapping[genotype]} is {len(y_list_temp)}')

    count = campaniform_fields[genotype_mapping[genotype]]['dscab'] + campaniform_fields[genotype_mapping[genotype]]['vped'] + campaniform_fields[genotype_mapping[genotype]]['dped'] + campaniform_fields[genotype_mapping[genotype]]['vscab']# + campaniform_fields[genotype_mapping[genotype]]['lcho']
    num_campaniforms_silenced.extend([count]*len(y_list_temp))
    unique_num_campaniforms_silenced.append(count)

    bootstrapped_data = bootstrap(y_list_temp, 10000)
    conf_min, conf_max = confidence_interval(bootstrapped_data)

    #B[k] = np.mean(bootstrapped_data)
    response.extend(y_list_temp)
    bootstrapped_mean_response.append(np.mean(bootstrapped_data))
    print (f'bootstrapped mean {np.mean(bootstrapped_data)}')

    #axes.plot(A[k], B[k], 'k.')
    #yerr = np.array([np.mean(bootstrapped_data)- conf_min, conf_max - np.mean(bootstrapped_data)]).reshape(2,-1)
    #axes.errorbar(A[k], B[k], xerr=0.2, yerr=yerr, color='k', elinewidth=1)
    

## center the variables (subtract mean)
#A = A - np.mean(A, 0)

axes.plot(num_campaniforms_silenced, response, 'o', markersize=1, c='gray')
axes.plot(unique_num_campaniforms_silenced, bootstrapped_mean_response, 'o', markersize=4, c='black')

num_campaniforms_silenced = np.array(num_campaniforms_silenced)
response = np.array(response)

slope, intercept, r, p, se = stats.linregress(num_campaniforms_silenced, response)
print ('t test, slope')
print (p, slope)

bootstrap_slopes(num_campaniforms_silenced, response)

A2 = sm.add_constant(num_campaniforms_silenced)
model = sm.OLS(response,A2, hasconst=True)
multi_model = model.fit()

predictions = multi_model.predict(A2)
residuals = predictions - response

pred = multi_model.get_prediction(A2).summary_frame()
sorted_indices = np.argsort(num_campaniforms_silenced)
axes.fill_between(num_campaniforms_silenced[sorted_indices], pred['mean_ci_lower'][sorted_indices], pred['mean_ci_upper'][sorted_indices], color='gray', alpha=0.3)

print (multi_model.summary())
print (multi_model.rsquared)
print ('conf', multi_model.conf_int(alpha=0.05, cols=None))
Ftest = np.identity(len(multi_model.params))
Ftest = Ftest[1:,:]
print(multi_model.f_test(Ftest))

print ('pearson r:', stats.pearsonr(num_campaniforms_silenced, response))


  
# creating regression plots
# fig = sm.graphics.plot_regress_exog(multi_model, 1, fig=fig)

# plt.show()

# plt.plot(B, residuals, 'o')
# plt.show()

# import ipdb;ipdb.set_trace()

####### to test for heteroskedascity ###########
test_result = sms.het_breuschpagan(multi_model.resid, multi_model.model.exog)
names = ['Lagrange multiplier statistic', 'p-value',
         'f-value', 'f p-value']
print (lzip(names, test_result))


######## to test for multicollinearity ###########
for i in range(0, A2.shape[1]):
	print (variance_inflation_factor(A2, i))



### plot prediction ###
# for i in range(0, A.shape[0]):
# 	#axes.plot(A[i][0], B[i],'bo')
# 	axes.plot(A[i][3], multi_model.predict(A2[i]), 'gx')

axes.plot(num_campaniforms_silenced, multi_model.params[0] + multi_model.params[1]*num_campaniforms_silenced, c=[0,0,0])
# slope, intercept, r, p, se = stats.linregress(A, multi_model.predict(A2))
# print(f"R-squared: {r**2:.6f}")

# slope = round(slope, 3)
# intercept = round(intercept, 3)
# r = round(r, 3)
# #p = round(p, 3)
# se = round(se, 3)
# print ('linregress: ', slope,intercept, r, p, se)
adjust_spines(axes,['bottom', 'left'],xticks=[0,10,20,30,40], yticks=[5,10,15], linewidth=0.6,spineColor='k')


axes.set_xlabel('Number of campaniforms targeted')
axes.set_ylabel('Head responses \n (degs)')
#axes.set_ylim(0,20)
axes.set_xlim(-5,45)
#plt.legend(loc="upper right")
#plt.title('Kir speed 1; Adj R squared: 0.91. P value of x4 =0.012')
plt.show()