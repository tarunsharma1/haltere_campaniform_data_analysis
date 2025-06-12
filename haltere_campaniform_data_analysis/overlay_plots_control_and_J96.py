import pickle
import matplotlib.pyplot as plt
import numpy as np
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

average_left_minus_right1, average_head_yaw1, resampled_trajectories1, t1 = pickle.load(open('./data_for_wing_head_phase_new/KirJ00yawbothdirections-speed-2.p', 'rb'))
average_left_minus_right2, average_head_yaw2, resampled_trajectories2, t2 = pickle.load(open('./data_for_wing_head_phase_new/KirJ96yawbothdirections-speed-2.p', 'rb'))

t = np.arange(0,min(t1[-1], t2[-1]), 1/72)



for k in range(0,len(average_left_minus_right1)):
	average_left_minus_right1[k] = resample(t1, average_left_minus_right1[k], t, kind='linear')
	average_head_yaw1[k] = resample(t1, average_head_yaw1[k], t, kind='linear')
	resampled_trajectories1[k] = resample(t1, resampled_trajectories1[k], t, kind='linear')

for k in range(0,len(average_left_minus_right2)):
	average_left_minus_right2[k] = resample(t2, average_left_minus_right2[k], t, kind='linear')
	average_head_yaw2[k] = resample(t2, average_head_yaw2[k], t, kind='linear')
	resampled_trajectories2[k] = resample(t2, resampled_trajectories2[k], t, kind='linear')



mean_trajectory = np.mean(resampled_trajectories1,axis=0)

fig, ax1 = plt.subplots(figsize=(2.4,1))
	
#ax1.set_ylim(-60,60)
ax1.set_ylim(-20,20)

adjust_spines(ax1,['bottom', 'left'],xticks=[0,20],  yticks=[-10,10], linewidth=0.6,spineColor='k')


for k in range(0,len(average_left_minus_right1)):
	average_head_yaw1[k] = average_head_yaw1[k] - np.mean(average_head_yaw1[k])
	ax1.plot(t, average_head_yaw1[k], c='r', alpha=0.07, linewidth=0.5)


ax2 = ax1.twinx()

for k in range(0,len(average_left_minus_right2)):
	average_head_yaw2[k] = average_head_yaw2[k] - np.mean(average_head_yaw2[k])
	ax2.plot(t, average_head_yaw2[k] - np.mean(average_head_yaw2[k]),c='g', alpha=0.07, linewidth=0.5)


mean1 = np.mean(average_left_minus_right1, axis=0)
mean_yaw1 = np.mean(average_head_yaw1, axis=0)
mean2 = np.mean(average_left_minus_right2, axis=0)
mean_yaw2 = np.mean(average_head_yaw2, axis=0)

ax1.plot(t, mean_yaw1, c='r', alpha=1, linewidth=0.5)
ax2.plot(t, mean_yaw2, c='g', alpha=1, linewidth=0.5)
#ax2.set_ylim(-60,60)
ax2.set_ylim(-20,20)

adjust_spines(ax2,['bottom'],xticks=[0,20], linewidth=0.6,spineColor='k')


ax4 = ax1.twinx()
ax4.plot(t, mean_trajectory, 'b', linewidth=0.5)
ax4.set_ylim(-530,530)
adjust_spines(ax4,['bottom', 'right'],xticks=[0,20],  yticks=[-200,200], linewidth=0.6,spineColor='k')


plt.show()
