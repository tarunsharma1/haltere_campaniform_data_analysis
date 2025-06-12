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
		interp = spi.interp1d(x1, y1, kind=kind, bounds_error=False, fill_value=fill_value)
		y2 = interp(x2)
	return y2


parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8, 'svg.fonttype':'none'}
plt.rcParams.update(parameters)

time1, WBA1, head1, trajectory1 = pickle.load(open('KirJ00yawbothdirections_data_for_lineplot.p', 'rb'))
time2, WBA2, head2, trajectory2 = pickle.load(open('KirJ96yawbothdirections_data_for_lineplot.p', 'rb'))

# WBA2_resampled = []
# for w2 in WBA2:
# 	WBA2_resampled.append(resample(time2 , w2, time1, kind='linear'))

# figure, axes = plt.subplots()
# #figure.set_size_inches(1.3,0.9)
# figure.set_size_inches(1.4,0.9)

# sns.lineplot(list(time1)*len(WBA1), list(np.array(WBA1).flatten()),ax = axes, sort=False,color='k', linewidth=0.6)
# axes.set_xlim([0,time1[-1]])
# axes.set_ylim([-40,40])


# axes2 = axes.twinx()
# sns.lineplot(list(time1)*len(WBA2_resampled), list(np.array(WBA2_resampled).flatten()),ax = axes2, sort=False,color='k', linewidth=0.6)
# axes2.set_xlim([0,time1[-1]])
# axes2.set_ylim([-40,40])


# axes3 = axes.twinx()

# axes3.plot(time1, trajectory1, 'b', linewidth=0.6)

# adjust_spines(axes,['bottom', 'left'],xticks=[0,5],yticks=[-20, 20], linewidth=0.6,spineColor='k')
# adjust_spines(axes3,['bottom', 'right'],xticks=[0,5], yticks=[-200,200], linewidth=0.6,spineColor='k')

# adjust_spines(axes2,['bottom'],xticks=[0,5], linewidth=0.6,spineColor='k')


# axes.set_xlabel('5 (s)')
# axes.set_ylabel('Wing response \n (degs)')
# axes3.set_ylabel('Angular velocity \n' + r'(deg $s^{-1}$)')
# plt.show()



# # plt.clf()
head2_resampled = []
for w2 in head2:
	head2_resampled.append(resample(time2 , w2, time1, kind='linear'))

figure, axes = plt.subplots()
figure.set_size_inches(1.4,0.9)

sns.lineplot(list(time1)*len(head1), list(np.array(head1).flatten()),ax = axes, sort=False,color='k', linewidth=0.6)
axes.set_xlim([0,time1[-1]])
axes.set_ylim([-15,15])


axes2 = axes.twinx()
sns.lineplot(list(time1)*len(head2_resampled), list(np.array(head2_resampled).flatten()),ax = axes2, sort=False,color='k', linewidth=0.6)
axes2.set_xlim([0,time1[-1]])
axes2.set_ylim([-15,15])


axes3 = axes.twinx()

axes3.plot(time1, trajectory1, 'b', linewidth=0.6)

adjust_spines(axes,['bottom', 'left'],xticks=[0,5],yticks=[-10, 10], linewidth=0.6,spineColor='k')
adjust_spines(axes3,['bottom', 'right'],xticks=[0,5], yticks=[-200,200], linewidth=0.6,spineColor='k')

adjust_spines(axes2,['bottom'],xticks=[0,5], linewidth=0.6,spineColor='k')


axes.set_xlabel('5 (s)')
axes.set_ylabel('Head response \n (degs)')
axes3.set_ylabel('Angular velocity \n' + r'(deg $s^{-1}$)')
plt.show()