import sys # Standard Packages
from copy import deepcopy
sys.path.insert(0, '../')

import numpy as np # 3rd Party Packages

from mmm_package import variables,constants # Local Packages

# Converts a variable in CDF format into the format needed for mmm
# This interpolates variables onto a new grid and converts some units
# Also applies optional smoothing and removal of outliers
def convert_variable(cdf_var, num_points):
    # deepcopy needed to create a new variable, else var is just a reference to cdf_var
    var = deepcopy(cdf_var)

    # Convert units to those needed by MMM
    units = var.units
    if units == 'CM':
        var.values /= 100
        var.units = 'M'
    elif units == 'CM/SEC':
        var.values /= 100
        var.units = 'M/SEC'
    elif units == 'N/CM**3':
        var.values = var.values*(10**6)
        var.units = 'N/M**3'
    elif units == 'EV':
        var.values /= 1000
        var.units = 'kEV'
    elif units == 'CM**2/SEC':
        var.values /= 10**4
        var.units = 'M**2/SEC'

    return var

# Calculates input variables for the MMM script from CDF data
def calc_inputs(cdf_vars, num_points=200):
    # Input variables for MMM will be stored in new vars object
    vars = variables.Variables()

    # Copy independent variables
    vars.time = deepcopy(cdf_vars.time)
    vars.x = deepcopy(cdf_vars.x)
    vars.xb = deepcopy(cdf_vars.xb)

    # Add the origin to the boundary grid (xb.values.shape[0] == x.values.shape[0] + 1)
    vars.xb.values = np.concatenate((np.zeros((1, cdf_vars.get_ntimes())), cdf_vars.xb.values), axis=0)

    # Store single column arrays for future calculations
    x = vars.x.values[:, 0]
    xb = vars.xb.values[:, 0]

    # Check that interpolation points is not smaller than the number of boundary points
    num_points = max(num_points, xb.size)

    # Get list of CDF variables to convert to MMM variables
    cdf_var_list = cdf_vars.get_cdf_variables()
    already_copied = ['time', 'x', 'xb']

    # Remove independent variables already copied from CDF variable list
    for var in already_copied:
        if var in cdf_var_list:
            cdf_var_list.remove(var)

    # Convert remaining CDF variables into MMM variables
    for var in cdf_var_list:
        converted_var = convert_variable(getattr(cdf_vars, var), num_points)
        # Can't assign getattr() to the function call, so need to set members individually
        getattr(vars, var).desc = converted_var.desc
        getattr(vars, var).units = converted_var.units
        getattr(vars, var).values = converted_var.values

if __name__ == '__main__':
    pass
