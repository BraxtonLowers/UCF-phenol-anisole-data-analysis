#! python3
# The purpose of this script is to utilize data from NMRPeaks.csv and analyte_info.json to produce plots of
# concentration vs corrected chemical shift, perform regression, and output the results in a table
#
# author: Braxton Lowers

import numpy as np
import pandas as pd
from matplotlib import pyplot
import seaborn as sns
from sklearn import linear_model
from sklearn.metrics import r2_score
import json



def create_plot(x_corrected, y_corrected, x_ground_truth, y_ground_truth, analyte, peak,
                ground_truth_uncertainty, corrected_uncertainty, is_linear, regression):
    """
Creates a plot with corrected datapoints, ground truth datapoints, a best fit line, and a fully formatted legend
    :param numpy.Series x_corrected: Series with train x values
    :param numpy.Series y_corrected: Series with train y values
    :param numpy.Series x_ground_truth: Series with test x values
    :param numpy.Series y_ground_truth: Series with test y values
    :param str analyte: name of analyte
    :param str peak: name of peak plotted
    :param numpy.Series ground_truth_uncertainty: Series with x uncertainty for ground truth
    :param numpy.Series corrected_uncertainty: Series with x uncertainty for corrected values
    :param bool is_linear: is linear or logarithmic
    :param linear_model.LinearRegression regression: LinearRegression fitted to train values
    """
    best_fit_line_x = np.arange(x_corrected.iloc[0], x_corrected.iloc[-1], 0.1).reshape(-1,1)
    if is_linear is False:
        best_fit_line_y = regression.predict(np.log(best_fit_line_x))
        r2_train = r2_score(y_corrected, np.log(x_corrected) * corrections[analyte_to_view]['slope'] +
                            corrections[analyte_to_view]['intercept'])
        r2_test = r2_score(y_ground_truth, np.log(x_corrected) * corrections[analyte_to_view]['slope'] +
                           corrections[analyte_to_view]['intercept'])
        equation_of_line = 'best fit line: f(x) = ' + str(regression.coef_[0][0]) + '*ln(x) + ' + \
                           str(regression.intercept_[0])
    else:
        best_fit_line_y = regression.predict(best_fit_line_x)
        r2_train = r2_score(np.array(y_corrected).reshape(-1,1),
                            regression.predict(np.array(x_corrected).reshape(-1, 1)))
        r2_test = r2_score(np.array(y_ground_truth).reshape(-1,1),
                           regression.predict(np.array(x_ground_truth).reshape(-1, 1)))
        equation_of_line = 'best fit line: f(x) = ' + str(regression.coef_[0][0]) + 'x + ' + str(regression.intercept_[0])
    best_fit_line = pyplot.plot(best_fit_line_x, best_fit_line_y, ':', color='black',
                                label=equation_of_line)
    corrected_axis_label = 'corrected assuming linear CDCl$_3$ concentration dependence\nR$^2$: ' + str(r2_train)
    pyplot.scatter(x=x_corrected, y=y_corrected, c='blue', marker='v', label=corrected_axis_label)
    ground_truth_axis_label = 'referenced to TMS\nR$^2$: ' + str(r2_test)
    pyplot.scatter(x=x_ground_truth, y=y_ground_truth, c='orange', marker='^', label=ground_truth_axis_label)
    title_str = analyte + ' in CDCl$_3$: ' + peak + ' peak. Molality vs chemical shift'
    pyplot.title(title_str, fontdict={'fontsize': 24})
    pyplot.xlabel('Molality (mol/kg)', fontdict={'fontsize': 16})
    pyplot.ylabel('Chemical shift (ppm)', fontdict={'fontsize': 16})
    pyplot.legend()
    pyplot.errorbar(x=x_ground_truth, y=y_ground_truth, xerr=ground_truth_uncertainty, capsize=3,
                    color='black', linestyle='None')
    pyplot.errorbar(x=x_corrected, y=y_corrected, xerr=corrected_uncertainty, capsize=3, color='black',
                    linestyle='None')
    pyplot.show()


# Loads corrections to reference data
################################################## commented out for dummy test purposes
# correctionsFile = '.\path to corrections output'
# with open(correctionsFile) as json_data:
#     corrections = json.load(json_data)
json_data = '{"phenol":{"slope":0.02323674138873531, "intercept":0.21047156916347376, "isLinear":"False"}}'
corrections = json.loads(json_data)
print(corrections)
# Open uncorrected dataset as a pandas dataframe and prepare an averaged dataset
filename = r'C:\Users\Braxton Lowers\Desktop\Raw spectral data\allNMRpeaksWithUncertainty.csv'
dataset = pd.read_csv(filename, header=0)
# TODO needs a better strategy manage NaN values
dataset = dataset.dropna()
grouped_dataset = dataset.groupby(['analyte', 'solvent', 'molality'])
averaged_dataset = grouped_dataset.mean().reset_index()
# Open ground truth dataset as a pandas dataframe
TMSfilename = r'C:\Users\Braxton Lowers\Desktop\Raw spectral data\TMS standardized - uncertainties.csv'
TMS_referenced_data = pd.read_csv(TMSfilename, header=0)
# Iterates over peaks and analytes and perform linear correction on each
analyteList = ['phenol', 'anisole', 'thiophenol', 'thioanisole']
for analyte_to_view in analyteList:
    peakList = ['ipso', 'ortho', 'meta', 'para']
    for peak_to_view in peakList:
        averaged_dataset['corrected'] = averaged_dataset[peak_to_view] + (averaged_dataset.molality *
                                                                     corrections[analyte_to_view]['slope'] -
                                                                     corrections[analyte_to_view]['intercept'])
        # Prepare data for plotting
        correctedX = averaged_dataset[(averaged_dataset['analyte'] == analyte_to_view) &
                                       (averaged_dataset['solvent'] == 'cdcl3')]['molality']
        correctedY = averaged_dataset[(averaged_dataset['analyte'] == analyte_to_view) &
                                      (averaged_dataset['solvent'] == 'cdcl3')]['corrected']
        groundX = TMS_referenced_data[TMS_referenced_data['analyte'] == analyte_to_view]['molality']
        groundY = TMS_referenced_data[TMS_referenced_data['analyte'] == analyte_to_view][peak_to_view]
        groundUncertain = TMS_referenced_data[TMS_referenced_data['analyte'] == analyte_to_view]['uncertain']
        correctedUncertain = averaged_dataset[(averaged_dataset['analyte'] == analyte_to_view) &
                                              (averaged_dataset['solvent'] == 'cdcl3')]['molaluncertainty']
        # Prepare a linear regression model and fit it to the corrected data
        regression = linear_model.LinearRegression()
        regression.fit(np.array(correctedX).reshape(-1,1), np.array(correctedY).reshape(-1,1))
        is_linear = corrections[analyte_to_view]['isLinear']
        # Plot data
        # Contains dummy data
        create_plot(x_corrected=correctedX, y_corrected=correctedY, x_ground_truth=groundX, y_ground_truth=groundY,
                    analyte=analyte_to_view, peak=peak_to_view, ground_truth_uncertainty=groundUncertain,
                    corrected_uncertainty=correctedUncertain, is_linear=True, regression=regression)
