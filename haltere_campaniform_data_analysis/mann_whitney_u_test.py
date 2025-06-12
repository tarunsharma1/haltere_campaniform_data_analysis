import numpy as np
import scipy.stats as stats
import pickle
from scipy.stats import mannwhitneyu

controls_genotype = 'UXJ00yawbothdirections'
#controls_genotype = 'KirJ00yawbothdirections'

list_of_genotypes = ['UXJ75yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ88yawbothdirections', 'UX28C05yawbothdirections',  'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']
#list_of_genotypes = ['KirJ70yawbothdirections', 'KirJ73yawbothdirections', 'Kir28C05yawbothdirections','KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections',  'KirJ75yawbothdirections', 'KirJ96yawbothdirections']


alpha = 0.05
alpha_corrected = alpha/(1.0*len(list_of_genotypes))

for genotype in list_of_genotypes:
	wba_both_speeds = []
	head_both_speeds = []

	speeds = ['2']
	for speed in speeds:
		#print ('Speed: ' + speed)
		# wba_controls = pickle.load(open('./data_for_scatter_plot_new/'+ controls_genotype + '-speed-' + speed +'-phases.p', 'rb'))
		# head_controls = pickle.load(open('./data_for_scatter_plot_new/'+ controls_genotype + '-speed-' + speed +'-head-yaw-phases.p', 'rb'))
		
		# wba = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'-phases.p', 'rb'))
		# head = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'-head-yaw-phases.p', 'rb'))


		### baselines
		wba_controls_baseline = pickle.load(open('./data_for_baseline_correlation_plot_new/'+ controls_genotype + '-speed-' + speed +'.p', 'rb'))
		#wba_controls = np.array(wba_controls_baseline[1]) + np.array(wba_controls_baseline[2]) ## L + R
		
		
		### baselines
		wba_baseline = pickle.load(open('./data_for_baseline_correlation_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
		#wba = np.array(wba_baseline[1]) + np.array(wba_baseline[2]) ## L + R
		

		
		wba_controls = np.array(wba_controls)
		#head_controls = np.array(head_controls)

		wba = np.array(wba)
		#head = np.array(head)
		

		U1, p = mannwhitneyu(wba,wba_controls,alternative='two-sided')
		if p<= alpha_corrected:
			print ('Wing: ' + genotype + ' is significant with p value:' + str(p))


		# U1, p = mannwhitneyu(head,head_controls,alternative='two-sided')
		# if p<= alpha_corrected:
		# 	print ('Head: ' + genotype + ' is significant with p value:' + str(p))


		#print (U1, p)
		#nx, ny = len(wba_controls), len(wba)
		#print ('counts:',nx, ny)
		#U2 = nx*ny - U1
		#print ('U2', U2)
		

		# if nx < ny:
		# 	print ('test statistic is:', U1)
		# elif nx > ny:
		# 	print ('test statistic is:', U2)
		# else:
		# 	print ('equal', U1, U2)


		# print ('data :')
		# print (wba_controls)
		# print (wba)

		# print ('Wings :')
		# print(stats.ttest_ind(wba_controls, wba, equal_var = False))
		# print (' Head : ')
		# print(stats.ttest_ind(head_controls, head, equal_var = False))
