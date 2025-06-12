import numpy as np
import scipy.stats as stats
from scipy.stats import mannwhitneyu


def calculate_ratio_confidence_interval(A_point, A_lower, A_upper, B_point, B_lower, B_upper,covariance, confidence_level=0.95):
    # Calculate the standard errors from the confidence intervals
    z_score = stats.norm.ppf(1 - (1 - confidence_level) / 2)  # z for the given confidence level
    
    SE_A = (A_upper - A_lower) / (2 * z_score)
    SE_B = (B_upper - B_lower) / (2 * z_score)

    var_A = SE_A**2
    var_B = SE_B**2
    
    # Calculate the point estimate of the ratio
    R = A_point / B_point

    term1 = var_A / B_point**2
    term2 = (A_point**2 * var_B) / B_point**4
    term3 = -2 * A_point * covariance / (B_point**3)

    var_R = term1 + term2 + term3
    #print (var_R)
    SE_R = np.sqrt(abs(var_R))
    
    # Propagate the uncertainties (standard errors) to the ratio
    #SE_R = np.sqrt((1 / B_point**2) * SE_A**2 + (A_point / B_point**2)**2 * SE_B**2)
    
    # Calculate the confidence interval for the ratio
    margin_of_error = z_score * SE_R
    CI_lower = R - margin_of_error
    CI_upper = R + margin_of_error
    
    return R, CI_lower, CI_upper

# Example usage:
### Kir
genotypes = ['KirJ00yawbothdirections', 'KirJ70yawbothdirections', 'KirJ73yawbothdirections', 'Kir28C05yawbothdirections', 'KirJ88yawbothdirections', 'KirJ86yawbothdirections', 'KirJ79yawbothdirections', 'KirJ75yawbothdirections', 'KirJ96yawbothdirections']
A_point = [27.229012708680898, 30.05242997292106, 29.203215921520123, 24.82584208851877, 23.958673467724807, 22.881234038549376, 20.952477067973888, 20.788509786488454, 18.08183549732626]
A_lower = [25.56402648267131, 28.157472606232044, 26.48451511260686, 22.822114171922404, 21.85087178873235, 20.597191795669836, 19.375125511328317, 18.626376863484776, 15.850941299562027]
A_upper = [28.964376506307715, 31.95720404824618, 31.7790265496718, 26.69589299153346, 26.022334732925103, 25.177126281683876, 22.416949149929128, 22.95833887957429, 20.409336477931504]

B_point = [9.481838949615813, 9.6906003935501, 10.300777372787325, 8.591444949519328, 10.540544849568775, 10.20634020279452, 7.6051101762672975, 7.102728317367063, 6.373236746639519]
B_lower = [8.689069082473965, 8.87300115841861, 9.264143484030484, 7.758442598406445, 9.823045788524807, 9.539421912412758, 6.956614633724738, 6.659385493180754, 5.885524920868067]
B_upper = [10.333802537767005, 10.540020971303308, 11.36834272129144, 9.394164326458572, 11.243652312468104, 10.885791914698949, 8.241714253775122, 7.565876066913203, 6.93600153429438]

### UX
# genotypes = ['UXJ00yawbothdirections', 'UXJ75yawbothdirections', 'UXJ79yawredobothdirections', 'UXJ88yawbothdirections', 'UX28C05yawbothdirections', 'UXJ90yawredobothdirections', 'UXJ96yawbothdirections']
# A_point = [30.13242977703665, 29.063201053843116, 28.726879433164378, 28.549452727352705, 27.412542815001082, 21.61870633023243, 16.716367325708497]
# A_lower = [27.474187734792267, 26.53771409006447, 25.711021585756495, 26.372789794245193, 24.591387808660887, 18.965989582854235, 13.519832153181829]
# A_upper = [32.624163284501776, 31.48135237566113, 31.530512559241874, 30.572879218056954, 30.12906298152973, 24.11303301183178, 20.06607787320253]

# B_point = [12.671004475145333, 11.180571191768554, 11.096417982521091, 12.571377957197948, 13.136827903156068, 9.478794823423852, 6.731050884495325]
# B_lower = [11.861688516124266, 10.103281375998185, 9.75123151986668, 11.576179278538726, 11.978136782599183, 8.413689435306088, 5.07031338613447]
# B_upper = [13.486107223639802, 12.223095092400083, 12.398427023783524, 13.572693501231289, 14.159907294734335, 10.484597470097835, 8.315730268946322]


confidence_level = 0.95  # 95% confidence level

covariance = np.cov(A_point, B_point)[0,1]
print ('covariance is: ',covariance)

# Calculate the ratio and its confidence interval
for i in range(0,len(A_point)):
    print (genotypes[i])
    
    R, CI_lower, CI_upper = calculate_ratio_confidence_interval(A_point[i], A_lower[i], A_upper[i], B_point[i], B_lower[i], B_upper[i], covariance, confidence_level)
    print (R, CI_lower, CI_upper)

    U1, p = mannwhitneyu(wba,wba_controls,alternative='two-sided')