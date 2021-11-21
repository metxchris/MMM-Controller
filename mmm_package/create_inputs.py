import sys # Standard Packages
from copy import deepcopy
sys.path.insert(0, '../')

import numpy as np # 3rd Party Packages
from scipy.interpolate import interp1d

from mmm_package import variables,constants # Local Packages

# Store sizes of the number of points for different dimensions
class NumPoints(object):
    def __init__(self, interpolation_points, boundary_points, time_points):
        self.interpolation = interpolation_points # TODO: unused
        self.boundary = boundary_points
        self.time = time_points 

# Converts a variable in CDF format into the format needed for mmm
# This interpolates variables onto a new grid and converts some units
# Also applies optional smoothing and removal of outliers
# Assumes values of cdf_var are not None
def convert_variable(cdf_var, num_points, x, xb):
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
        var.values *= 10**6
        var.units = 'N/M**3'
    elif units == 'EV':
        var.values /= 1000
        var.units = 'kEV'
    elif units == 'CM**2/SEC':
        var.values /= 10**4
        var.units = 'M**2/SEC'

    # Reshape non-scalar variables so that their shape matches (XB, TIME)
    if var.values.ndim == 0:
        pass
    # Tile 1-dim time arrays into 2-dim arrays, in the format of (XB, TIME)
    # This also adds the origin to the X-axis
    elif var.values.shape[0] == num_points.time:
        var.values = np.tile(var.values, (num_points.boundary, 1))
    # Some variables (i.e. VPOL) are mirrored around the X-axis, so take non-zero X values
    elif var.values.shape[0] == 2 * num_points.boundary - 1:
        var.values = var.values[num_points.boundary - 1:, :]
    # Interpolate remaining variables from X to XB (adds origin to the X-axis)
    else:
        set_interp = interp1d(x, var.values.T, kind='cubic', fill_value="extrapolate")
        var.values = set_interp(xb).T 

    # TODO: Apply smoothing using moving average

    return var

# Calculates input variables for the MMM script from CDF data
def create_inputs(cdf_vars, num_interp_points=200):
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

    # Store sizes and check that interpolation points is not smaller than the number of boundary points
    num_points = NumPoints(max(num_interp_points, xb.size), xb.size, vars.get_ntimes())

    # Get list of CDF variables to convert to MMM variables
    cdf_var_list = cdf_vars.get_cdf_variables()
    already_copied = ['time', 'x', 'xb']

    # Remove independent variables already copied from CDF variable list
    for var in already_copied:
        if var in cdf_var_list:
            cdf_var_list.remove(var)

    # Convert remaining CDF variables into MMM variables
    for var in cdf_var_list:
        cdf_var = getattr(cdf_vars, var)
        # Variables previously not found in the CDF will not have values
        if cdf_var.values is not None:
            setattr(vars, var, convert_variable(cdf_var, num_points, x, xb))

    return vars

if __name__ == '__main__':
    pass
