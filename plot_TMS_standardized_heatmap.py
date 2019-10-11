#! python3
# The purpose of this script is to calculate and plot the difference between the low and high concentration chemical
# shifts
# Author: Braxton Lowers

import pandas
from matplotlib import pyplot
import seaborn

# Loads data as dataframe
filename = r'C:\Users\Braxton Lowers\Desktop\Raw spectral data\TMS standardized.csv'
data = pandas.read_csv(filename, header=0)
# Sorts and groups data in preparation for taking difference
data.sort_values(['analyte', 'molality'])
grouped = data.groupby(by='analyte')
# Removes NaN
difference = grouped.diff().dropna()
# Creates a new column that maps the indicies to the analyte column
difference['index'] = data.loc[difference.index]['analyte']
# Turns analyte column into index
difference = difference.set_index(['index'])
# Fixes sign to be intuitive
difference = difference.apply(lambda x: x*(-1))
# Plot heatmap
plot = seaborn.heatmap(difference.loc[:,['ipso', 'meta', 'para', 'ortho', 'solvent peak']], annot=True, center=0,
                       cmap='coolwarm_r', robust=True)
plot.set_title('Range of chemical shift in ppm. Note that concentrations are NOT controlled for!')
pyplot.xlabel('ppm', fontdict={'size':16})
pyplot.ylabel('Analyte', fontdict={'size':16})
pyplot.show()
# TODO perform molal control calculations (df.apply?) and replot with these settings
plot = seaborn.heatmap(difference, annot=True, center=0, cmap='coolwarm_r', robust=True, vmin=-0.04)
plot.set_title('Change of chemical shift/molal in CDCl$_3$ upon dilution (ppm/molal)', fontdict={'size':16})
pyplot.xlabel('ppm/molal', fontdict={'size':16})
pyplot.ylabel('Analyte', fontdict={'size':16})
pyplot.xlabel('ppm')
pyplot.show()
print(difference)
