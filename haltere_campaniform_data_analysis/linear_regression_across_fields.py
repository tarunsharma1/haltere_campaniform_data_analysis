import numpy as np
from sklearn.linear_model import LinearRegression
import pickle
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor

def bootstrap(data, n):
	## randomly resample n times and take the mean each time
	bootstrapped_data = np.zeros(n)
	for i in range(0,n):
		sample = np.random.choice(data, size=len(data))
		bootstrapped_data[i] = np.mean(np.array(sample))
	return bootstrapped_data





#### standard linear regression - only 9 rows for Kir and 5 predictors so statistics can vary a lot with small changes #####################
def standard_regression():
	Y =[]
	campaniform_fields_present = {}

	for k,genotype in enumerate(list(genotype_mapping.keys())):
	    campaniform_fields_present[genotype_mapping[genotype]] = campaniform_fields[genotype_mapping[genotype]]
		
	    y_list_wing = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
	    y_list_head = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed + '-head-yaw.p', 'rb'))
	    

	    bootstrapped_data = bootstrap(np.array(y_list_wing), 10000)

	    Y.append(np.mean(bootstrapped_data)) 


	df = pd.DataFrame(campaniform_fields_present)
	df = df.transpose()


	df['mean_responses'] = Y

	X = df[['dscab', 'dped', 'vped', 'vscab']]
	Y = df['mean_responses']



	X = sm.add_constant(X)  # adds intercept term
	model = sm.OLS(Y, X).fit()
	print(model.summary())

	correlation_matrix = X.corr(method='pearson')
	print ('correlation matrix:')
	print (correlation_matrix)

	print ([variance_inflation_factor(X.values, i)
                   for i in range(X.shape[1])])


############# Mixed effects model with genotype as random effect and includes all flies as rows so about 9 x 20 rows for Kir ##############
def mixed_effects_regression():

	rows = []

	for k,genotype in enumerate(list(genotype_mapping.keys())):

		y_list_wing = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed +'.p', 'rb'))
		y_list_head = pickle.load(open('./data_for_scatter_plot_new/'+ genotype + '-speed-' + speed + '-head-yaw.p', 'rb'))
		for i in range(len(y_list_wing)):
			rows.append([campaniform_fields[genotype_mapping[genotype]]['dscab'], campaniform_fields[genotype_mapping[genotype]]['dped'], campaniform_fields[genotype_mapping[genotype]]['vscab'], 
				campaniform_fields[genotype_mapping[genotype]]['vped'], campaniform_fields[genotype_mapping[genotype]]['lcho'], y_list_head[i], genotype_mapping[genotype]])

	df = pd.DataFrame(rows, columns=['dscab', 'dped', 'vscab', 'vped', 'lcho', 'response', 'genotype'])

	model = smf.mixedlm("response ~ dscab + dped + vscab + vped + lcho",
                    data=df,
                    groups="genotype")
	result = model.fit()
	print(result.summary())



if __name__ == "__main__":
	#### cells silenced average of mine and annes counts ####
	campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
	'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
	'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}

	speed = '2'

	genotype_mapping = {'KirJ00yawbothdirections':'J00','KirJ73yawbothdirections': '17G01' ,'KirJ86yawbothdirections': '58F02', 'KirJ79yawbothdirections': '31A09', 'KirJ96yawbothdirections': '86H12', 'KirJ70yawbothdirections':'14B04', 'Kir28C05yawbothdirections':'28C05', 'KirJ88yawbothdirections':'60B12', 'KirJ75yawbothdirections':'22E04'}
	#genotype_mapping = {'UXJ00yawbothdirections':'J00', 'UX28C05yawbothdirections':'28C05', 'UXJ88yawbothdirections':'60B12', 'UXJ79yawredobothdirections':'31A09', 'UXJ90yawredobothdirections':'74B09', 'UXJ75yawbothdirections':'22E04', 'UXJ96yawbothdirections':'86H12'}

	#mixed_effects_regression()
	standard_regression()