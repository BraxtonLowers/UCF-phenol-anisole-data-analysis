#! python3
# The purpose of this script is to calculate reference corrections for CDCl3 samples using TMS standardized data
# Author: Braxton Lowers

import pandas
from sklearn import linear_model
import json

# Load data
filename = r'C:\Users\Braxton Lowers\Desktop\Raw spectral data\TMS standardized.csv'
data = pandas.read_csv(filename, header=0)
# Initialize lists
analyte_list = ['anisole', 'phenol', 'thioanisole', 'thiophenol']
peak_list = ['molality', 'ipso', 'meta', 'para', 'ortho', 'solvent peak']
json_data = {}
# Iterate over all analytes
for analyte in analyte_list:
    # Prepare data for linear regression
    solvent_shift = data[data['analyte'] == analyte]['solvent peak'].values.reshape(-1, 1)
    solvent_concentration = data[data['analyte'] == analyte]['molality'].values.reshape(-1, 1)
    # Perform linear regression
    linear_regression = linear_model.LinearRegression()
    linear_regression.fit(solvent_concentration, solvent_shift)
    print(77.23 - linear_regression.intercept_[0])
    # Calculate intercept correction
    interceptCorrection = 77.23 - linear_regression.intercept_[0]
    if analyte == 'phenol':
        isLinear = False
        # TODO insert log corrections here in an if-else chain
    else:
        isLinear = True
    # Add log factors manually
    logFactor = -0.16671387
    logIntercept = 155.11514
    # Create corrections dictionary
    json_data[analyte] = {'slope':float(linear_regression.coef_[0][0]), 'intercept': float(interceptCorrection),
                          'isLinear': isLinear, 'logFactor': logFactor, 'logIntercept': logIntercept}
# Convert corrections dictionary to json format for saving
outputJson = json.dumps(json_data)
# Save calculated corrections as output file
with open(r'.\corrections', 'w') as correctionsFile:
    correctionsFile.write(outputJson)
