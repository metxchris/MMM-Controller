# Standard Packages
import sys
from copy import deepcopy
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np 
from scipy.interpolate import interp1d
# Local Packages
from main import variables
from main.options import Options
import settings

class XValues:
    '''Stores single dimension arrays of the values of X, XB, and XB + origin (xbo)'''
    def __init__(self, xvar, xbvar):
        self.x = xvar.values[:, 0]
        self.xb = xbvar.values[:, 0] # implicitly used in the interpolation step
        self.xbo = np.append([0], self.xb)

def convert_variable(cdf_var, xvals):
    '''
    Converts variable units from CDF format to MMM format, then interpolates data onto the new XBO grid

    A new XB grid which also contains the origin (XBO) is created, since the origin is needed when later 
    calculating gradients. All variables are mapped from their current positional dimension to this
    new XBO grid, which makes plotting and variable comparisons much easier. All variables of a single
    dimension have their values copied over to their missing dimension, which allows for fast vectorized
    calculations when later calculating input variables.

    Parameters:
    * cdf_var (Variable): A single variable object containing data from the CDF
    * xvals (XValues): Cached x-dimension values needed for the interpolation process

    Returns:
    * input_var (Variable): A single variable object containing the converted data from cdf_var
    '''

    # deepcopy needed to create a new variable instead of a reference
    input_var = deepcopy(cdf_var)

    # Convert units to those needed by MMM
    units = input_var.units
    if units == 'CM':
        input_var.set_variable(input_var.values / 100, 'M')
    elif units == 'CM/SEC':
        input_var.set_variable(input_var.values / 100, 'M/SEC')
    elif units == 'N/CM**3':
        input_var.set_variable(input_var.values * 10**6, 'N/M**3')
    elif units == 'EV':
        input_var.set_variable(input_var.values / 1000, 'kEV')
    elif units == 'CM**2/SEC':
        input_var.set_variable(input_var.values / 10**4, 'M**2/SEC')
    elif units == 'AMPS':
        input_var.set_variable(input_var.values / 10**6, 'MAMPS')

    # Reshape all non-scalar variables so that their shape matches (XBo, TIME)
    xdim = input_var.get_xdim()
    # 0-dimensional variables are not reshaped
    if xdim is None or input_var.values.ndim == 0:
        pass
    # Tile 1-dim time arrays into 2-dim arrays, in the format of (XBO, TIME) 
    elif xdim in ['TIME', 'TIME3']:
        input_var.set_variable(np.tile(input_var.values, (xvals.xbo.size, 1)))
        input_var.dimensions = ['XBO', xdim]
    # Some variables (i.e. VPOL) are mirrored around the X-axis, so take non-negative XB values
    # TODO: Handle this case better
    elif xdim in ['RMAJM']:
        input_var.set_variable(input_var.values[xvals.xbo.size - 1:, :])
        input_var.set_xdim('XBO')
    # Interpolate/Extrapolate variable from X or XB to XBO
    elif xdim in ['X', 'XB']:
        set_interp = interp1d(getattr(xvals, xdim.lower()), input_var.values, kind='cubic', fill_value="extrapolate", axis=0)
        input_var.set_variable(set_interp(xvals.xbo))
        input_var.set_xdim('XBO')
    else:
        print('[initial_conversion] *** Warning: Unsupported interpolation xdim type for variable', input_var.name, xdim)
    
    # Apply smoothing, then verify minimum values (fixes errors due to interpolation)
    input_var.apply_smoothing()
    input_var.set_minvalue()

    return input_var

def initial_conversion(cdf_vars):
    '''
    Initializes the process of converting variables from CDF format to MMM format

    The values of cdf_vars are copied to input_vars, then various values are stored before the conversion
    process begins.  The measurement time and index of this time are also obtained from the input time, where
    the time closest to the input time is taken as the measurement time.

    Parameters:
    * cdf_vars (InputVariables): Raw data from the CDF of all InputVariables with a specifed cdf_var value

    Returns:
    * input_vars (InputVariables): Data from the CDF converted into a format needed to run MMM
    '''

    input_options = Options.instance
    input_vars = variables.InputVariables()

    # Directly copy independent variables, which don't need to be converted
    input_vars.time = deepcopy(cdf_vars.time)
    input_vars.x = deepcopy(cdf_vars.x)
    input_vars.xb = deepcopy(cdf_vars.xb)

    # Cache single column arrays of x-values
    xvals = XValues(input_vars.x, input_vars.xb)

    # Add the origin to the boundary grid
    input_vars.xb.set_variable(np.concatenate((np.zeros((1, cdf_vars.get_ntimes())), cdf_vars.xb.values), axis=0))

    # Set the array index and measurement time value corresponding to the input time
    input_options.set_measurement_time(input_vars.time)

    # Update the number of input points, if needed
    if input_options.input_points is None:
        input_options.input_points = input_vars.get_nboundaries()

    # Get list of CDF variables to convert to the format needed for MMM, and remove variables already copied
    cdf_var_list = cdf_vars.get_cdf_variables()
    for var in ['time', 'x', 'xb']:
        cdf_var_list.remove(var)

    # Convert remaining CDF variables into the format needed for MMM
    for var in cdf_var_list:
        cdf_var = getattr(cdf_vars, var)
        if cdf_var.values is not None:
            setattr(input_vars, var, convert_variable(cdf_var, xvals))

    # Use TEPRO, TIPRO in place of TE, TI
    if settings.USE_TEMPERATURE_PROFILES:
        input_vars.use_temperature_profiles()

    return input_vars

def interp_to_input_points(input_vars):
    '''
    Interpolates from XB to a grid determined by input_points, if the input_point grid differs from XB

    Since the XB grid is constant over time, the interpolation process can be done in one step for all timeslices
    for each input variable.

    Parameters:
    * input_vars (InputVariables): Contains all variables from CDF + extra calculated variables

    Returns:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file + extra calculated variables
    '''

    input_options = Options.instance
    mmm_vars = deepcopy(input_vars)

    # Interpolation only needed if input_points != xb points
    if input_options.input_points != input_vars.get_nboundaries():
        # Single column arrays for interpolation
        xb = mmm_vars.xb.values[:, 0]
        xb_mmm = np.arange(input_options.input_points) / (input_options.input_points - 1)

        # Get list of CDF variables to convert to the format needed for MMM
        full_var_list = mmm_vars.get_nonzero_variables()

        # Remove independent variables
        for var in ['time', 'x']:
            full_var_list.remove(var)

        # Interpolate variables onto grid specified by input_options.input_points
        for var in full_var_list:
            mmm_var = getattr(mmm_vars, var)
            if mmm_var.values is not None:
                set_interp = interp1d(xb, mmm_var.values, kind='cubic', fill_value="extrapolate", axis=0)
                mmm_var.set_variable(set_interp(xb_mmm))
                mmm_var.set_minvalue()
            else:
                print(f'ERROR: Trying to interpolate variable {var} with values equal to None')

    return mmm_vars

def interp_to_uniform_rho(input_vars):
    '''
    Interpolates each timeslice onto an evenly spaced rho, determined by input_points

    Since the value of rmin varies over time, the interpolation onto a grid of evenly spaced rho values
    requires that interpolation be carried out over each individual timeslice of variable data, which can 
    increase the time needed to interpolate the data by a factor of 10.

    Parameters:
    * input_vars (InputVariables): Contains all variables from CDF + extra calculated variables

    Returns:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file + extra calculated variables
    '''

    input_options = Options.instance
    mmm_vars = deepcopy(input_vars)

    # Single column arrays for interpolation
    old_rho = mmm_vars.rho.values
    new_rho = np.arange(input_options.input_points) / (input_options.input_points - 1)

    # Get list of CDF variables to convert to the format needed for MMM
    full_var_list = mmm_vars.get_nonzero_variables()

    # Remove independent variables
    for var in ['time', 'x']:
        full_var_list.remove(var)

    # Interpolate variables onto grid specified by input_options.input_points
    for var in full_var_list:
        mmm_var = getattr(mmm_vars, var)
        interp_values = np.empty((len(new_rho), mmm_var.values.shape[1]))
        if mmm_var.values is not None:
            for time_idx in range(mmm_var.values.shape[1]):
                set_interp = interp1d(old_rho[:, time_idx], mmm_var.values[:, time_idx], 
                    kind='cubic', fill_value="extrapolate", axis=0)
                interp_values[:, time_idx] = set_interp(new_rho)
            mmm_var.set_variable(interp_values)
            mmm_var.set_minvalue()
        else:
            print(f'ERROR: Trying to interpolate variable {var} with values equal to None')

    return mmm_vars

def final_interpolation(input_vars):
    '''
    Interpolates to either a new grid of input_points, or a grid of uniform rho values

    Parameters:
    * input_vars (InputVariables): Contains all variables from CDF + extra calculated variables

    Returns:
    * (InputVariables): Contains all variables needed to write MMM input file + extra calculated variables
    '''

    uniform_rho = Options.instance.uniform_rho
    return interp_to_uniform_rho(input_vars) if uniform_rho else interp_to_input_points(input_vars)

if __name__ == '__main__':
    pass
