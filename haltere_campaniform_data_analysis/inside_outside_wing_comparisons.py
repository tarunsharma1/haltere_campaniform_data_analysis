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


# genotype_mapping = {'KirJ00yawneg':'J00', 'KirJ00yawpos':'J00', 'KirJ96yawneg':'86H12', 'KirJ96yawpos':'86H12', 'KirJ73yawneg': '17G01', 'KirJ73yawpos': '17G01', 'KirJ86yawneg': '58F02',
# 'KirJ86yawpos': '58F02', 'KirJ79yawneg': '31A09', 'KirJ79yawpos': '31A09', 'KirJ70yawneg':'14B04', 'KirJ70yawpos':'14B04', 'Kir28C05yawneg':'28C05', 'Kir28C05yawpos':'28C05',
# 'KirJ88yawneg':'60B12', 'KirJ88yawpos':'60B12', 'KirJ75yawneg':'22E04', 'KirJ75yawpos':'22E04'}


genotype_mapping = {'UXJ00yawneg':'J00','UXJ00yawpos':'J00', 'UX28C05yawneg':'28C05', 'UX28C05yawpos':'28C05', 'UXJ88yawneg':'60B12', 'UXJ88yawpos':'60B12',
'UXJ79yawredoneg':'31A09', 'UXJ79yawredopos':'31A09', 'UXJ90yawredoneg':'74B09', 'UXJ90yawredopos':'74B09', 'UXJ75yawneg':'22E04', 'UXJ75yawpos':'22E04',
'UXJ96yawneg':'86H12', 'UXJ96yawpos':'86H12'}



#### cells silenced average of mine and annes counts ####
campaniform_fields = {'J00':{'dscab':0, 'dped':0, 'vscab':0, 'vped':0, 'lcho':0},'28C05':{'dscab':2,'dped':2, 'vscab':0, 'vped':6, 'lcho':0} ,'60B12':{'dscab':0, 'dped':1, 'vscab':2, 'vped':6, 'lcho':10},
'31A09':{'dscab':4, 'dped':14, 'vscab':0, 'vped':11, 'lcho':9}, '74B09':{'dscab':14, 'dped':12, 'vscab':0, 'vped':11, 'lcho':0}, '58F02':{'dscab': 14, 'dped':8, 'vscab':0, 'vped':9, 'lcho':0}, '22E04':{'dscab':11, 'dped':11, 'vscab':0, 'vped':12, 'lcho':0}, '14B04':{'dscab':3, 'dped':2, 'vscab':1, 'vped':1, 'lcho':0},
'86H12':{'dscab':7, 'dped':8, 'vscab':4, 'vped':10, 'lcho':0}, '17G01':{'dscab':2, 'dped':9, 'vscab':0, 'vped':3, 'lcho':6}}

speed = '2'


num_campaniforms_silenced = []

response = []

unique_num_campaniforms_silenced = []
fig, axes = plt.subplots()
fig.set_size_inches(2.5,1.5)


for k,genotype in enumerate(list(genotype_mapping.keys())):
    print (genotype_mapping[genotype])

    inside, outside = pickle.load(open('./data_for_ipsi_contra_new/' + genotype + '-inside-outside-speed-' + speed + '.p', 'rb'))
    
    y_list_temp = np.abs(np.array(inside))

    count = campaniform_fields[genotype_mapping[genotype]]['dscab'] + campaniform_fields[genotype_mapping[genotype]]['vped'] + campaniform_fields[genotype_mapping[genotype]]['dped'] + campaniform_fields[genotype_mapping[genotype]]['vscab']# + campaniform_fields[genotype_mapping[genotype]]['lcho']
    num_campaniforms_silenced.extend([count]*len(y_list_temp))
    unique_num_campaniforms_silenced.append(count)

    response.extend(y_list_temp)
    



axes.plot(num_campaniforms_silenced, response, 'o', markersize=1, c='gray')


num_campaniforms_silenced = np.array(num_campaniforms_silenced)
response = np.array(response)

slope, intercept, r, p, se = stats.linregress(num_campaniforms_silenced, response)
print ('t test, slope')
print (p, slope)

A2 = sm.add_constant(num_campaniforms_silenced)
model = sm.OLS(response,A2, hasconst=True)
multi_model = model.fit()
predictions = multi_model.predict(A2)
residuals = predictions - response

pred = multi_model.get_prediction(A2).summary_frame()
sorted_indices = np.argsort(num_campaniforms_silenced)
axes.fill_between(num_campaniforms_silenced[sorted_indices], pred['mean_ci_lower'][sorted_indices], pred['mean_ci_upper'][sorted_indices], color='gray', alpha=0.3)
axes.plot(num_campaniforms_silenced, multi_model.params[0] + multi_model.params[1]*num_campaniforms_silenced, c=[0,0,0])

print (multi_model.summary())
plt.show()