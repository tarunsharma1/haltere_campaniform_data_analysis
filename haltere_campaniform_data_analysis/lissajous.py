from numpy import sin,pi,linspace,cos
from pylab import plot,show,subplot
import matplotlib.pyplot as plt
import matplotlib as mpl

parameters = {'axes.labelsize':8,'axes.titlesize':8, 'xtick.labelsize':8, 'font.family':"sans-serif", 'font.sans-serif':['Arial'], 'font.size':8}
plt.rcParams.update(parameters)

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

fig, axes = plt.subplots()
fig.set_size_inches(1,1)

a = 0.055 # plotting the curves for
b = 0.055 # different values of a/b
delta = 96.28*pi/180
delta2 = 117.13*pi/180
t = linspace(0,6*pi,300)

wings = 27.08*cos(2*pi*a*t + delta)
head = 9.46*cos(2*pi*b*t + delta2)

axes.plot(head,wings, c='k', linewidth=0.6)
adjust_spines(axes,['bottom', 'left'],xticks=[-10,10], yticks=[-20,20], linewidth=0.6,spineColor='k')

axes.set_xlabel('Head response')
axes.set_ylabel('Wing response')

#plt.scatter(head[0], wings[0], color='red', label='Start Point (t=0)', zorder=5)
#plt.scatter(head[-1], wings[-1], color='blue', label='End Point (t=6π)', zorder=5)

# Plot arrows showing direction (optional)
#plt.arrow(head[0], wings[0], head[1] - head[0], wings[1] - wings[0], 
#          head_width=1, head_length=2, fc='black', ec='black')


# for i in range(0,4):
#     if i==0:
#         x = 27.08*cos(2*pi*a*t + delta)
#         y = 9.46*cos(2*pi*b*t + delta2)
#         plt.subplot(2,2,i+1)
#         plt.plot(x,y)
#     else:
#         x = 24.68*sin(2*pi*a*t + delta)
#         y = 8.57*sin(2*pi*b*t + delta2)
#         plt.subplot(2,2,i+1)
#         plt.plot(x,y)
plt.show()

