# Standard Packages
import sys
from copy import deepcopy
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np 
from scipy.interpolate import interp1d
import scipy.ndimage
# Local Packages
from mmm_package import variables 

# Store sizes of the number of points for different dimensions
class NumPoints(object):
    def __init__(self, interpolation_points, boundary_points, time_points):
        self.interpolation = interpolation_points # TODO: unused
        self.boundary = boundary_points
        self.time = time_points

# Stores single dimension arrays of the values of X, XB, and XB + origin (xbo)
class XValues(object):
    def __init__(self, x, xb, xbo):
        self.x = x
        self.xb = xb
        self.xbo = xbo

# Converts a variable in CDF format into the format needed for mmm
# This interpolates variables onto a new grid and converts some units
# Also applies optional smoothing and removal of outliers
# Assumes values of cdf_var are not None
def convert_variable(cdf_var, num_points, xvals):
    # deepcopy needed to create a new variable instead of a reference
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
    elif units == 'AMPS':
        var.values /= 10**6
        var.units = 'MAMPS'

    # Reshape all non-scalar variables so that their shape matches (XBo, TIME)
    # This allows us to vectorize our calculations later, making them much faster
    xdim = var.get_xdim()
    # 0-dimensional variables are not reshaped
    if xdim is None or var.values.ndim == 0:
        pass
    # Tile 1-dim time arrays into 2-dim arrays, in the format of (XBo, TIME)
    # This also adds the origin to the X-axis  
    elif xdim in ['TIME', 'TIME3']:
        var.values = np.tile(var.values, (num_points.boundary, 1))
        var.set_dims(['XBO', xdim])
    # Some variables (i.e. VPOL) are mirrored around the X-axis, so take non-negative XB values
    # TODO: Handle this case better
    elif xdim in ['RMAJM']:
        var.values = var.values[num_points.boundary - 1:, :]
        var.set_xdim('XBO')
    # Interpolate/Extrapolate variable from X or XB to XBo using a cubic spline
    elif xdim in ['X', 'XB']:
        set_interp = interp1d(getattr(xvals, xdim.lower()), var.values.T, kind='cubic', fill_value="extrapolate")
        var.values = set_interp(xvals.xbo).T
        var.set_xdim('XBO')
    else:
        print('[create_inputs] *** Warning: Unsupported interpolation xdim type for variable', var.name, xdim)

    # Variable smoothing using a Gaussian filter
    var.apply_smoothing()

    return var

# Calculates input variables for the MMM script from CDF data
def convert_inputs(cdf_vars, num_interp_points=200):
    # Input variables for MMM will be stored in new vars object
    vars = variables.Variables()

    # Copy independent variables
    vars.time = deepcopy(cdf_vars.time)
    vars.x = deepcopy(cdf_vars.x)
    vars.xb = deepcopy(cdf_vars.xb)

    # Cache single column arrays; xbo is xb with the origin tacked on
    xvals = XValues(x=vars.x.values[:, 0], 
                    xb=vars.xb.values[:, 0], 
                    xbo=np.append([0], vars.xb.values[:, 0]))

    # Add the origin to the boundary grid
    vars.xb.values = np.concatenate((np.zeros((1, cdf_vars.get_ntimes())), cdf_vars.xb.values), axis=0)
    vars.xb.set_xdim('XBo')

    # Cache sizes and check that interpolation points is not smaller than the number of boundary points
    num_points = NumPoints(max(num_interp_points, xvals.xbo.size), xvals.xbo.size, vars.get_ntimes())

    # Get list of CDF variables to convert to the format needed for MMM
    cdf_var_list = cdf_vars.get_cdf_variables()

    # Remove independent variables that were already copied
    for var in ['time', 'x', 'xb']:
        cdf_var_list.remove(var)

    # Convert remaining CDF variables into the format needed for MMM
    # TODO: Add option to use TIPRO, TEPRO in place of TI, TE
    for var in cdf_var_list:
        cdf_var = getattr(cdf_vars, var)
        # Variables previously not found in the CDF will not have values
        if cdf_var.values is not None:
            setattr(vars, var, convert_variable(cdf_var, num_points, xvals))

    return vars

if __name__ == '__main__':
    pass
