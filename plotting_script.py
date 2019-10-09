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



def create_plot(x_corrected, y_corrected, x_ground_truth, y_ground_truth, r2_train, r2_test, analyte, peak,
                ground_truth_uncertainty, corrected_uncertainty, is_linear, regression):
    """
Creates a plot with corrected datapoints, ground truth datapoints, a best fit line, and a fully formatted legend
    :param numpy.Series x_corrected: Series with train x values
    :param numpy.Series y_corrected: Series with train y values
    :param numpy.Series x_ground_truth: Series with test x values
    :param numpy.Series y_ground_truth: Series with test y values
    :param float r2_train: r squared score for train values with best fit line
    :param float r2_test: r squared score for test values with best fit line
    :param str analyte: name of analyte
    :param str peak: name of peak plotted
    :param numpy.Series ground_truth_uncertainty: Series with x uncertainty for ground truth
    :param numpy.Series corrected_uncertainty: Series with x uncertainty for corrected values
    :param bool is_linear: is linear or logarithmic
    :param linear_model.LinearRegression regression: LinearRegression fitted to train values
    """
    best_fit_line_x = np.arange(x_corrected[0], x_corrected[-1], 0.1).reshape(-1,1)
    if is_linear is False:
        best_fit_line_y = regression.predict(np.log(best_fit_line_x))
        equation_of_line = 'best fit line: f(x) = ' + str(regression.coef_[0][0]) + '*ln(x) + ' + str(regression.intercept_[0])
    else:
        best_fit_line_y = regression.predict(best_fit_line_x)
        equation_of_line = 'best fit line: f(x) = ' + str(regression.coef_[0][0]) + 'x + ' + str(regression.intercept_[0])
    best_fit_line = pyplot.plot(best_fit_line_x, best_fit_line_y, ':', color='black',
                                label=equation_of_line)
    corrected_axis_label = 'corrected assuming linear CDCl$_3$ concentration dependence\nR$^2$: ' + r2_train
    pyplot.scatter(x=x_corrected, y=y_corrected, c='blue', marker='v', label=corrected_axis_label)
    ground_truth_axis_label = 'referenced to TMS\nR$^2$: ' + r2_test
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


# Open dataset as a pandas dataframe
filename = r'C:\Users\Braxton Lowers\Desktop\Raw spectral data\allNMRpeaksWithUncertainty.csv'
dataset = pd.read_csv(filename, header=0)
grouped_dataset = dataset.groupby(['analyte', 'solvent', 'molality'])
averaged_dataset = grouped_dataset.mean().reset_index()