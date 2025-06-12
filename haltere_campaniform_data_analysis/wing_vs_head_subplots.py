import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib as mpl
import statsmodels.api as sm
import seaborn as sns
import pandas as pd

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


def bootstrap_slope(data, x_col, y_col, n_reps=1000):
    """Calculates the bootstrapped slope and confidence intervals."""

    slopes = []
    for _ in range(n_reps):
        # Resample the data with replacement
        sample = data.sample(n=len(data), replace=True)

        # Fit a linear regression model
        #slope, _, _, _, _ = stats.linregress(sample[x_col], sample[y_col])
        #slopes.append(slope)

        bs_model = sm.OLS(sample[y_col],sample[x_col], hasconst=False)
        bs_multi_model = bs_model.fit()
        ci = bs_multi_model.conf_int(alpha=0.05, cols=None)
        slopes.append(bs_multi_model.params[0])


    # Calculate the confidence interval
    lower = np.percentile(slopes, 2.5)
    upper = np.percentile(slopes, 97.5)

    return np.mean(slopes), lower, upper, slopes











parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8}
plt.rcParams.update(parameters)




all_slopes = []

means = []
lower_ci = []
upper_ci = []

wba_genotype_for_subplot = []
head_genotype_for_subplot = []

## per genotype, make plots of wing vs head as scatter plot dots
#list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UX28C05yawbothdirections', 'UXJ88yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']
#genotype_names = [ 'UXJ00', 'UX22E04', 'UX28C05' , 'UX60B12', 'UX31A09', 'UX74B09', 'UX86H12']


list_of_genotypes = ['KirJ00yawbothdirections','KirJ70yawbothdirections','KirJ73yawbothdirections', 'Kir28C05yawbothdirections', 'KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections', 'KirJ75yawbothdirections', 'KirJ96yawbothdirections']
genotype_names = ['KirJ00','Kir17G01', 'Kir58F02', 'Kir31A09', 'Kir86H12', 'Kir14B04', 'Kir28C05', 'Kir60B12', 'Kir22E04']

genotype_mapping = {'KirJ00yawbothdirections':'Control','KirJ73yawbothdirections': '17G01' ,'KirJ86yawbothdirections': '58F02', 'KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ70yawbothdirections':'14B04', 'Kir28C05yawbothdirections':'28C05', 'KirJ88yawbothdirections':'60B12', 'KirJ75yawbothdirections':'22E04'}
#genotype_mapping = {'UXJ00yawbothdirections':'Control', 'UX28C05yawbothdirections':'28C05', 'UXJ88yawbothdirections':'60B12', 'UXJ79yawredobothdirections':'31A09', 'UXJ90yawredobothdirections':'74B09', 'UXJ75yawbothdirections':'22E04', 'UXJ96yawbothdirections':'86H12'}

campaniform_fields = {'Control':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}



fig, axes = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(1,4/6)

df1 = pd.DataFrame()
for i,genotype in enumerate(list_of_genotypes):    
    wba_both_speeds = []
    head_both_speeds = []
    campaniforms_silenced = []
    speeds = ['2']
    for speed in speeds:
        wba = pickle.load(open('./data_for_scatter_plot/'+ genotype + '-speed-' + speed +'.p', 'rb'))
        head = pickle.load(open('./data_for_scatter_plot/'+ genotype + '-speed-' + speed +'-head-yaw.p', 'rb'))
        wba_both_speeds.extend(wba)
        head_both_speeds.extend(head)
        #import ipdb;ipdb.set_trace()
        campaniforms_silenced.extend([campaniform_fields[genotype_mapping[genotype]]['dscab'] + campaniform_fields[genotype_mapping[genotype]]['vped'] + campaniform_fields[genotype_mapping[genotype]]['dped'] + campaniform_fields[genotype_mapping[genotype]]['lcho'] + campaniform_fields[genotype_mapping[genotype]]['vscab']]*len(wba))

    wba = wba_both_speeds
    head = head_both_speeds
    slope, intercept, r, p, se = stats.linregress(head, wba)

    
    slope = round(slope, 3)
    intercept = round(intercept, 3)
    r = round(r, 3)
    #p = round(p, 3)
    se = round(se, 3)
    #print ('linregress: ', slope,intercept, r, p, se, genotype_mapping[genotype])


    #### make a scatter plot of wing vs head only for one genotype, the genotype_for_subplot

    head = np.array(head)
    wba = np.array(wba)
    A2 = head
    #A2 = sm.add_constant(head)
    model = sm.OLS(wba,A2, hasconst=False)
    multi_model = model.fit()
    ci = multi_model.conf_int(alpha=0.05, cols=None)
    #print (ci)
    
    print (multi_model.summary())

    adjust_spines(axes,['bottom','left'], xticks=[5,15], yticks=[20,30], linewidth=0.6,spineColor='k')

    #import ipdb;ipdb.set_trace()
    axes.scatter(head, wba, c='k', s=7)
    #plt.colorbar()
    axes.set_xlim(3,18)
    axes.set_ylim(15,40)
    #axes[i//3, i%3].set_xlabel('Head response \n (degs)')
    #axes[i//3, i%3].set_ylabel('Wing response \n (degs)')
    #axes[i//3, i%3].set_title(genotype_mapping[genotype] + ' pearson_coeff:' + str(round(stats.pearsonr(head, wba)[0], 2)) + ' p_value:' + str(round(stats.pearsonr(head, wba)[1], 3)))
    #plt.ylim(5,45)
    #plt.xlim(3,20)

    ##extend the line a little bit
    #head = np.concatenate((head, np.array([2,4,5,6,7,15])))
    axes.plot(head, head*multi_model.params[0], c=[0.5,0.5,0.5])

    pred = multi_model.get_prediction(A2).summary_frame()
    sorted_indices = np.argsort(head)
    axes.fill_between(head[sorted_indices], pred['mean_ci_lower'][sorted_indices], pred['mean_ci_upper'][sorted_indices], color='gray', alpha=0.3)


    df2 = pd.DataFrame({'head':head, 'wing':wba, 'genotype':genotype_mapping[genotype]})
    #df1 = df1.append(df2)

    bootstrapped_slope, lower, upper, all_slopes_bootstrapped = bootstrap_slope(df2, 'head', 'wing', 10000)
    means.append(bootstrapped_slope)
    lower_ci.append(lower)
    upper_ci.append(upper)
    all_slopes.append(all_slopes_bootstrapped)

    print (genotype_mapping[genotype])
    print ('bootstrapped slope and CI')
    print (bootstrapped_slope, lower, upper)
    print ('corr coef:')
    print (stats.pearsonr(head, wba) )
    break
#g = sns.lmplot(x='head', y='wing', data=df1, fit_reg=True, ci=95, n_boot=10000, col='genotype', height=3)
#g.fig.suptitle()

axes.set_xlabel('Head response \n (degs)')
axes.set_ylabel('Wing response \n (degs)')



# axes[2, 0].set_xlabel('Head response \n (degs)')
# axes[2, 1].set_xlabel('Head response \n (degs)')
# axes[2, 2].set_xlabel('Head response \n (degs)')
# axes[0, 0].set_ylabel('Wing response \n (degs)')
# axes[1, 0].set_ylabel('Wing response \n (degs)')
# axes[2, 0].set_ylabel('Wing response \n (degs)')

#plt.title('Wing vs head per genotype using Kir')

plt.show()






plt.bar([genotype_mapping[x] for x in list_of_genotypes], means, yerr=[np.array(means) - np.array(lower_ci), np.array(upper_ci) - np.array(means)], capsize=5)
plt.ylim(0,4)
plt.show()

## plot differences between bootstrapped slopes of controls and all other genotypes
for i in range(1, len(list_of_genotypes)):
    print (genotype_mapping[list_of_genotypes[i]])
    diff = np.array(all_slopes[i]) - np.array(all_slopes[0])
    lower = np.percentile(diff, 2.5)
    upper = np.percentile(diff, 97.5)
    print ('lower - upper', lower, upper)



