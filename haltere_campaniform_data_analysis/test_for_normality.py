import numpy as np
import scipy.stats as stats
import pickle
from numpy.random import seed
from numpy.random import randn
from scipy.stats import shapiro
from scipy.stats import normaltest
from scipy.stats import anderson
import matplotlib.pyplot as plt

list_of_genotypes = ['KirJ00yawbothdirections','KirJ73yawbothdirections' ,'KirJ86yawbothdirections', 'KirJ79yawbothdirections', 'KirJ96yawbothdirections', 'KirJ287yawbothdirections', 'KirJ70yawbothdirections', 'Kir28C05yawbothdirections', 'KirJ88yawbothdirections', 'KirJ75yawbothdirections']
#list_of_genotypes = [ 'UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UX28C05yawbothdirections', 'UXJ88yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']

all_data = []
for genotype in list_of_genotypes:
	print (genotype)
	#genotype = 'KirJ00yawbothdirections'

	speeds = ['1']
	for speed in speeds:
		print ('Speed: ' + speed)
		data = pickle.load(open('/home/tarun/catkin_ws/src/trajectories_autostep_ros/src/data_for_scatter_plot/'+ genotype + '-speed-' + speed +'-head-yaw.p', 'rb'))
		all_data.extend(data)
data = all_data
plt.hist(data, edgecolor='black', bins=20)
plt.show()
stat, p = shapiro(data)
print('Statistics=%.3f, p=%.3f' % (stat, p))
alpha = 0.05
if p > alpha:
	print('Shapiro : Sample looks Gaussian (fail to reject H0)')
else:
	print('############# Shapiro : Sample does not look Gaussian (reject H0) ##########################')


## normality test
stat, p = normaltest(data)
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('K squared : Sample looks Gaussian (fail to reject H0)')
else:
	print(' #############K squared : Sample does not look Gaussian (reject H0)############################')

# anderson test
result = anderson(data)
print('Statistic: %.3f' % result.statistic)
p = 0
for i in range(len(result.critical_values)):
	sl, cv = result.significance_level[i], result.critical_values[i]
	if result.statistic < result.critical_values[i]:
		print('%.3f: %.3f, Andersen data looks normal (fail to reject H0)' % (sl, cv))
	else:
		print('%.3f: %.3f, ############################ Andersen data does not look normal (reject H0) ############################' % (sl, cv))