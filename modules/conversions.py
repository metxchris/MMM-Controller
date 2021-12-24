# Standard Packages
from copy import deepcopy

# 3rd Party Packages
import numpy as np
from scipy.interpolate import interp1d

# Local Packages
import modules.options as options


class XValues:
    '''Stores single dimension arrays of the values of X, XB, and XB + origin (xbo)'''

    def __init__(self, xvar, xbvar):
        self.x = xvar.values[:, 0]
        self.xb = xbvar.values[:, 0]  # implicitly used in the interpolation step
        self.xbo = np.append([0], self.xb)


def convert_units(input_var):
    '''
    Converts variable units from CDF format to MMM format

    Parameters:
    * input_var (Variable): A reference to a variable object in input_vars, containing data from the CDF
    '''

    units = input_var.units
    if units == 'CM':
        input_var.set(values=input_var.values / 100, units='m')
    elif units == 'CM/SEC':
        input_var.set(values=input_var.values / 100, units='m/s')
    elif units == 'N/CM**3':
        input_var.set(values=input_var.values * 10**6, units='m^-3')
    elif units == 'EV':
        input_var.set(values=input_var.values / 1000, units='keV')
    elif units == 'CM**2/SEC':
        input_var.set(values=input_var.values / 10**4, units='m^2/s')
    elif units == 'AMPS':
        input_var.set(values=input_var.values / 10**6, units='MA')
    elif units == 'TESLA*CM':
        input_var.set(values=input_var.values / 100, units='T*m')
    elif units == 'SEC**-1':
        input_var.set(units='s^-1')
    elif units == 'RAD/SEC':
        input_var.set(units='rad/s')
    elif units == 'PASCALS':
        input_var.set(units='Pa')
    elif units == 'SECONDS':
        input_var.set(units='s')
    elif units == 'TESLA':
        input_var.set(units='T')


def interp_to_boundarygrid(input_var, xvals):
    '''
    Reshape/Interpolate all non-scalar variables so that their shape matches (XBO, TIME), where XBO = XB + origin

    A new XB grid which also contains the origin (XBO) is created, since the origin is needed when later
    calculating gradients. All variables are mapped from their current positional dimension to this
    new XBO grid, which makes plotting and variable comparisons much easier. All variables of a single
    dimension have their values copied over to their missing dimension, which allows for fast vectorized
    calculations when later calculating input variables.

    Parameters:
    * input_var (Variable): A reference to a variable object in input_vars, containing data from the CDF
    * xvals (XValues): Cached x-dimension values needed for the interpolation process
    '''

    xdim = input_var.get_xdim()
    # 0-dimensional variables are not reshaped
    if xdim is None or input_var.values.ndim == 0:
        pass
    # Tile 1-dim time arrays into 2-dim arrays, in the format of (XBO, TIME)
    elif xdim in ['TIME', 'TIME3']:
        input_var.set(values=np.tile(input_var.values, (xvals.xbo.size, 1)))
        input_var.dimensions = ['XBO', xdim]
    # Some variables (i.e. VPOL) are mirrored around the X-axis, so take non-negative XB values
    # TODO: Handle this case better
    elif xdim in ['RMAJM']:
        input_var.set(values=input_var.values[xvals.xbo.size - 1:, :])
        input_var.set_xdim('XBO')
    # Interpolate/Extrapolate variable from X or XB to XBO
    elif xdim in ['X', 'XB']:
        set_interp = interp1d(getattr(xvals, xdim.lower()), input_var.values, kind='cubic', fill_value="extrapolate", axis=0)
        input_var.set(values=set_interp(xvals.xbo))
        input_var.set_xdim('XBO')
    else:
        print('[initial_conversion] *** Warning: Unsupported interpolation xdim type for variable', input_var.name, xdim)


def interp_to_input_points(input_vars):
    '''
    Interpolates from XB to a grid determined by input_points, if the input_point grid differs from XB

    Since the XB grid is constant over time, the interpolation process can be done in one step for all time slices
    for each input variable.

    Parameters:
    * input_vars (InputVariables): Contains all variables from CDF + extra calculated variables

    Returns:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file + extra calculated variables
    '''
    opts = options.instance
    mmm_vars = deepcopy(input_vars)

    # Interpolation only needed if input_points != xb points
    if opts.input_points != input_vars.get_nboundaries():
        # Single column arrays for interpolation
        xb = mmm_vars.xb.values[:, 0]
        xb_new = np.arange(opts.input_points) / (opts.input_points - 1)

        # Get list of CDF variables to convert to the format needed for MMM
        full_var_list = mmm_vars.get_nonzero_variables()

        # Remove independent variables
        for var in ['time', 'x']:
            full_var_list.remove(var)

        # Interpolate variables onto grid specified by opts.input_points
        for var in full_var_list:
            mmm_var = getattr(mmm_vars, var)
            if mmm_var.values is None:
                raise ValueError(f'Trying to interpolate variable {var} with values equal to None')

            if mmm_var.values.size > 1:
                set_interp = interp1d(xb, mmm_var.values, kind='cubic', fill_value="extrapolate", axis=0)
                mmm_var.set(values=set_interp(xb_new))

    return mmm_vars


def interp_to_uniform_rho(input_vars):
    '''
    Interpolates each time slice onto an evenly spaced rho, determined by input_points

    Since the value of rmin varies over time, the interpolation onto a grid of evenly spaced rho values
    requires that interpolation be carried out over each individual time slice of variable data, which can
    increase the time needed to interpolate the data by a factor of 10.

    Parameters:
    * input_vars (InputVariables): Contains all variables from CDF + extra calculated variables

    Returns:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file + extra calculated variables
    '''
    opts = options.instance
    mmm_vars = deepcopy(input_vars)

    # Single column arrays for interpolation
    rho_old = mmm_vars.rho.values
    rho_new = np.arange(opts.input_points) / (opts.input_points - 1)

    # Get list of CDF variables to convert to the format needed for MMM
    full_var_list = mmm_vars.get_nonzero_variables()

    # Remove independent variables
    for var in ['time', 'x']:
        full_var_list.remove(var)

    # Interpolate variables onto grid specified by opts.input_points
    for var in full_var_list:
        mmm_var = getattr(mmm_vars, var)
        interp_values = np.empty((len(rho_new), mmm_var.values.shape[1]))
        if mmm_var.values is None:
            raise ValueError(f'Trying to interpolate variable {var} with values equal to None')

        if mmm_var.values.size > 1:
            for time_idx in range(mmm_var.values.shape[1]):
                set_interp = interp1d(rho_old[:, time_idx], mmm_var.values[:, time_idx],
                                      kind='cubic', fill_value="extrapolate", axis=0)
                interp_values[:, time_idx] = set_interp(rho_new)

            mmm_var.set(values=interp_values)

    return mmm_vars


def initial_conversion(cdf_vars):
    '''
    Initializes the process of converting variables from CDF format to MMM format

    The values of cdf_vars are copied to input_vars, then various values are stored before the conversion
    process begins.  The measurement time and index of this time are also obtained from the input time, where
    the time closest to the input time is taken as the measurement time.

    Parameters:
    * cdf_vars (InputVariables): Raw data from the CDF of all InputVariables with a specified cdf_var value

    Returns:
    * input_vars (InputVariables): Data from the CDF converted into a format needed to run MMM
    '''
    opts = options.instance
    input_vars = deepcopy(cdf_vars)

    # Cache single column arrays of x-values
    xvals = XValues(input_vars.x, input_vars.xb)

    # Add the origin to the boundary grid
    input_vars.xb.set(values=np.concatenate((np.zeros((1, cdf_vars.get_ntimes())), cdf_vars.xb.values), axis=0))

    # Set the array index and measurement time value corresponding to the input time
    opts.set_measurement_time(input_vars.time)

    # Update the number of input points, if needed
    if opts.input_points is None:
        opts.input_points = input_vars.get_nboundaries()

    # Get list of CDF variables to convert to the format needed for MMM
    # Independent variables listed below don't need to be converted
    cdf_var_list = cdf_vars.get_cdf_variables()
    for var_name in ['time', 'x', 'xb']:
        cdf_var_list.remove(var_name)

    # Convert remaining CDF variables into the format needed for MMM
    for var_name in cdf_var_list:
        input_var = getattr(input_vars, var_name)
        if input_var.values is not None:
            convert_units(input_var)
            interp_to_boundarygrid(input_var, xvals)

    # Use TEPRO, TIPRO in place of TE, TI
    if opts.temperature_profiles:
        input_vars.use_temperature_profiles()

    # Set value of rmin at origin to 0
    input_vars.rmin.values[0, :] = 0

    return input_vars


def convert_variables(cdf_vars):
    '''
    Initializes the process of converting variables from CDF format to MMM format

    The values of cdf_vars are copied to input_vars, then various values are stored before the conversion
    process begins.  The measurement time and index of this time are also obtained from the input time, where
    the time closest to the input time is taken as the measurement time.

    Parameters:
    * cdf_vars (InputVariables): Raw data from the CDF of all InputVariables with a specified cdf_var value

    Returns:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file + extra calculated variables
    '''
    opts = options.instance
    input_vars = initial_conversion(cdf_vars)
    input_vars.set_rho_values()

    uniform_rho = opts.uniform_rho
    mmm_vars = interp_to_uniform_rho(input_vars) if uniform_rho else interp_to_input_points(input_vars)

    full_var_list = mmm_vars.get_nonzero_variables()
    # Apply smoothing, then verify minimum values (fixes errors due to interpolation)
    for var_name in full_var_list:
        mmm_var = getattr(mmm_vars, var_name)
        if opts.apply_smoothing:
            mmm_var.apply_smoothing(opts.input_points)
        # Since interpolation can create multiple nonphysical values, no exceptions are raised for nonphysical values
        mmm_var.set_minvalue(raise_exception=False)

    # Update x from xb (x is the grid between xb, and has one fewer point than xb)
    mmm_vars.x.values = (mmm_vars.xb.values[0:-1, :] + mmm_vars.xb.values[1:, :]) / 2
    mmm_vars.set_rho_values()

    return mmm_vars
