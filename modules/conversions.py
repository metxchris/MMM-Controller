"""Handles the conversion of variables from TRANSP format to MMM format

The conversion process includes converting units used in TRANSP to those used
by MMM, followed by one or two interpolation steps.  The first interpolation
step is used to put all TRANSP variables on the same radial grid, since some
variables are defined on X (r/a centers), where as other variables are
defined on XB (r/a boundaries).  The new radial grid is these grids are
interpolated to is just the same XB grid, but with the origin added. In doing
so, we create some unavoidable error at the end-points of our calculations
(rho=0 and rho=1), however this error can be reduced by increasing the amount
of radial points used in a TRANSP run.

The second interpolation step is optional, and allows the user to create more
data points to send to the MMM driver.  In particular, the amount of TRANSP
data points should be increased (roughly doubled) when the uniform rho option
is specified, to minimize errors from creating uniformly spaced radial
points.

All variables are interpolated onto their final radial grid before any
calculations are performed.  In our experience, interpolating variables after
performing calculations will lead to inconsistencies in values if those
variables are later recalculated during the variable adjustment process.  For
example, if variable A was interpolated onto a larger grid first and then
calculated, and variable B was calculated first and then interpolated onto a
larger grid, comparing the values of A and B will show a difference greater
than that attributed to floating point errors.
"""

# Standard Packages
from copy import deepcopy

# 3rd Party Packages
import numpy as np
from scipy.interpolate import interp1d

# Local Packages
import modules.options as options


class _XValues:
    '''Stores single dimension arrays of the values of X, XB, and XB + origin (xbo)'''

    def __init__(self, xvar, xbvar):
        self.x = xvar.values[:, 0]
        self.xb = xbvar.values[:, 0]  # implicitly used in the interpolation step
        self.xbo = np.append([0], self.xb)


def _convert_units(input_var):
    '''
    Converts variable units from CDF format to MMM format

    Parameters:
    * input_var (Variable): Object containing variable data
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


def _interp_to_boundarygrid(input_var, xvals):
    '''
    Reshape/Interpolate all non-scalar variables so that their shape matches
    (XBO, TIME), where XBO = XB + origin

    A new XB grid which also contains the origin (XBO) is created, since the
    origin is needed when later calculating gradients. All variables are
    mapped from their current positional dimension to this new XBO grid,
    which makes plotting and variable comparisons much easier. All variables
    of a single dimension have their values copied over to their missing
    dimension, which allows for fast vectorized calculations when later
    calculating new variables.

    Parameters:
    * input_var (Variable): Object containing variable data
    * xvals (_XValues): Cached x-dimension values needed for the interpolation process

    Raises:
    * NotImplementedError: If no interpolation is defined for the xdim of the variable
    '''

    xdim = input_var.get_xdim()
    # 0-dimensional variables are not reshaped
    if xdim is None or input_var.values.ndim == 0:
        pass

    # Tile 1-dim time arrays into 2-dim arrays, in the format of (XBO, TIME)
    elif xdim in ['TIME', 'TIME3']:
        input_var.set(values=np.tile(input_var.values, (xvals.xbo.size, 1)))
        input_var.dimensions = ['XBO', xdim]

    # Some variables (i.e. VPOL) exist on a full diameter, rather than the radius of rmin
    elif xdim in ['RMAJM']:
        input_var.set(values=input_var.values[xvals.xbo.size - 1:, :])
        input_var.set_xdim('XBO')

    # Interpolate/Extrapolate variable from X or XB to XBO
    elif xdim in ['X', 'XB']:
        set_interp = interp1d(getattr(xvals, xdim.lower()), input_var.values, kind='cubic', fill_value="extrapolate", axis=0)
        input_var.set(values=set_interp(xvals.xbo))
        input_var.set_xdim('XBO')

    else:
        raise NotImplementedError(f'Unsupported interpolation xdim {xdim} for variable {input_var.name}')


def _interp_to_input_points(input_vars):
    '''
    Interpolates from XB to a grid of size determined by input points

    Since the XB grid is constant over time, the interpolation process can be
    done in one step for all time slices for each input variable.  This
    function only executes if the specified input points differs from the
    number of boundaries already used by the data.

    Parameters:
    * input_vars (InputVariables): Object containing variable data

    Returns:
    * mmm_vars (InputVariables): Object containing interpolated variable data

    Raises:
    * ValueError: If variable to interpolate is None
    '''

    input_points = options.instance.input_points
    mmm_vars = deepcopy(input_vars)

    # Interpolation only needed if input_points != xb points
    if input_points != input_vars.get_nboundaries():
        # Single column arrays for interpolation
        xb = mmm_vars.xb.values[:, 0]
        xb_new = np.arange(input_points) / (input_points - 1)

        full_var_list = mmm_vars.get_nonzero_variables()
        for var in ['time', 'x']:  # Remove independent variables
            full_var_list.remove(var)

        # Interpolate variables onto grid specified by input_points
        for var in full_var_list:
            mmm_var = getattr(mmm_vars, var)
            if mmm_var.values is None:
                raise ValueError(f'Trying to interpolate variable {var} with values equal to None')

            if mmm_var.values.size > 1:
                set_interp = interp1d(xb, mmm_var.values, kind='cubic', fill_value="extrapolate", axis=0)
                mmm_var.set(values=set_interp(xb_new))

    return mmm_vars


def _interp_to_uniform_rho(input_vars):
    '''
    Interpolates each time slice onto uniformly spaced grid determined by
    input points

    Since the value of rmin varies over time, the interpolation onto a grid of
    evenly spaced rho values requires that interpolation be carried out over
    each individual time slice of variable data, which can increase the time
    needed to interpolate the data by a factor of 10.

    Parameters:
    * input_vars (InputVariables): Object containing variable data

    Returns:
    * mmm_vars (InputVariables): Object containing interpolated variable data

    Raises:
    * ValueError: If variable to interpolate is None
    '''

    input_points = options.instance.input_points
    mmm_vars = deepcopy(input_vars)

    # Single column arrays for interpolation
    rho_old = mmm_vars.rho.values
    rho_new = np.arange(input_points) / (input_points - 1)

    full_var_list = mmm_vars.get_nonzero_variables()
    for var in ['time', 'x']:  # Remove independent variables
        full_var_list.remove(var)

    # Interpolate variables onto uniform grid specified by input_points
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


def _initial_conversion(cdf_vars):
    '''
    Initializes the process of converting variables from CDF format to MMM format

    The values of cdf_vars are copied to input_vars, then various values are
    stored before the conversion process begins.  The measurement time and
    index of this time are also obtained from the input time, where the time
    closest to the input time is taken as the measurement time.

    Parameters:
    * cdf_vars (InputVariables): Object containing TRANSP formatted variable data

    Returns:
    * input_vars (InputVariables): Object containing MMM formatted variable data
    '''

    opts = options.instance
    input_vars = deepcopy(cdf_vars)

    # Cache single column arrays of x-values
    xvals = _XValues(input_vars.x, input_vars.xb)

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
            _convert_units(input_var)
            _interp_to_boundarygrid(input_var, xvals)

    # Use TEPRO, TIPRO in place of TE, TI
    if opts.temperature_profiles:
        input_vars.use_temperature_profiles()

    input_vars.rmin.values[0, :] = 0  # Set value of rmin at origin to 0

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

    input_vars = _initial_conversion(cdf_vars)
    input_vars.set_rho_values()

    uniform_rho = options.instance.uniform_rho
    mmm_vars = _interp_to_uniform_rho(input_vars) if uniform_rho else _interp_to_input_points(input_vars)

    full_var_list = mmm_vars.get_nonzero_variables()
    # Apply smoothing, then verify minimum values (fixes errors due to interpolation)
    for var_name in full_var_list:
        mmm_var = getattr(mmm_vars, var_name)
        if options.instance.apply_smoothing:
            mmm_var.apply_smoothing(options.instance.input_points)

        # Since interpolation can create multiple expected nonphysical values,
        # no exceptions are raised for fixing these issues
        mmm_var.set_minvalue(raise_exception=False)

    mmm_vars.set_x_values()
    mmm_vars.set_rho_values()

    return mmm_vars
