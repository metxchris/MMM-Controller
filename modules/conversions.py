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

# 3rd Party Packages
import numpy as np
from scipy.interpolate import interp1d, Akima1DInterpolator

# Local Packages
import settings
import modules.datahelper as datahelper
import modules.constants as constants


class _XValues:
    '''Stores single dimension arrays of the values of X, XB, and XB + origin (xbo)'''

    def __init__(self, xvar, xbvar):
        self.x = xvar.values[:, 0]
        self.xb = xbvar.values[:, 0]  # implicitly used in the interpolation step
        self.xbo = np.append([0], self.xb)
        self.time = xvar.values[0, :]


def _choose_variables(input_vars):
    '''
    Chooses a CDF variable to use in cases where CDF variables don't always
    contain data

    Parameters:
    * input_vars (InputVariables): Object containing variable data
    '''

    # For some reason, SREXBSMOD values appear to be shifted by an entire radial point
    input_vars.wexbsmod.values[:-1, :] = input_vars.wexbsmod.values[1:, :]

    # Choose wexb (CDF values are actually 1/s even if they say rad/s)
    wexbsv2 = input_vars.wexbsv2
    wexbsmod = input_vars.wexbsmod
    wexbsa = input_vars.wexbsa
    tolerance = 1.001

    if wexbsmod.values is not None and not (wexbsmod.values <= tolerance * wexbsv2.default_values).all():
        input_vars.wexb.values = wexbsmod.values
    elif wexbsv2.values is not None and not (wexbsv2.values <= tolerance * wexbsv2.default_values).all():
        input_vars.wexb.values = wexbsv2.values
    elif wexbsa.values is not None and not (wexbsa.values <= tolerance * wexbsv2.default_values).all():
        input_vars.wexb.values = wexbsa.values


def convert_units(input_var):
    '''
    Converts variable units from CDF format to MMM format

    Parameters:
    * input_var (Variable): Object containing variable data
    '''

    # TODO: Implement lowercase checks to remove redundancies
    units = input_var.units.upper()
    if units == 'CM':
        input_var.set(values=input_var.values / 100, units='m')
    elif units == 'CM**-1':
        input_var.set(values=input_var.values * 100, units='m^-1')
    elif units == 'CM**-2':
        input_var.set(values=input_var.values * 10**4, units='m^-2')
    elif units == 'CM**2' or units == 'CM2':
        input_var.set(values=input_var.values / 10**4, units='m^2')
    elif units == 'CM**3' or units == 'CM3':
        input_var.set(values=input_var.values / 10**6, units='m^3')
    elif units == 'CM/SEC':
        input_var.set(values=input_var.values / 100, units='m/s')
    elif units == 'N/CM**3' or units == '#/CM**3':
        input_var.set(values=input_var.values * 10**6, units='m^-3')
    elif units == 'EV' or units =='eV':
        input_var.set(values=input_var.values / 1000, units='keV')
    elif units == 'CM**2/SEC':
        input_var.set(values=input_var.values / 10**4, units='m^2/s')
    elif units == 'AMPS' or units == 'A':
        input_var.set(values=input_var.values / 10**6, units='MA')
    elif units == 'TESLA*CM':
        input_var.set(values=input_var.values / 100, units='T*m')
    elif units == 'V/CM':
        input_var.set(values=input_var.values * 100, units='V/m')
    elif units == 'AMPS/CM2' or units == 'AMPS/CM**2':
        input_var.set(values=input_var.values / 100, units='MA/m^2')
    elif units == 'OHM*CM':
        input_var.set(values=input_var.values / 100, units='ohm*m')
    elif units == 'JLES/CM3':
        input_var.set(values=input_var.values * 10**6 / constants.ZCKB, units='keV/m^3')
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
    elif units == 'WEBERS':
        input_var.set(units='T*m^2')
    elif units == 'HOURS':
        input_var.set(units='h')


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
    # Check if variable was found
    if xdim is None:
        if input_var.default_values is not None:
            input_var.set(values=np.tile(input_var.values, (xvals.xbo.size, xvals.time.size)))
            input_var.dimensions = ['XBO', 'TIME']
        else:
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
        # set_interp = interp1d(getattr(xvals, xdim.lower()), input_var.values,
        #                       kind=settings.INTERPOLATION_METHOD, fill_value="extrapolate", axis=0)
        # input_var.set(values=set_interp(xvals.xbo))
        set_interp = Akima1DInterpolator(getattr(xvals, xdim.lower()), input_var.values, axis=0)
        input_var.set(values=set_interp(xvals.xbo, extrapolate=True))
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

    mmm_vars = datahelper.deepcopy_data(input_vars)
    input_points = mmm_vars.options.input_points

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

            if isinstance(mmm_var.values, np.ndarray) and mmm_var.values.size > 1:
                set_interp = interp1d(xb, mmm_var.values, kind=settings.INTERPOLATION_METHOD,
                                      fill_value="extrapolate", axis=0)
                mmm_var.set(values=set_interp(xb_new))

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

    input_vars = datahelper.deepcopy_data(cdf_vars)

    # Cache single column arrays of x-values
    xvals = _XValues(input_vars.x, input_vars.xb)

    # Add the origin to the boundary grid
    input_vars.xb.set(values=np.concatenate((np.zeros((1, cdf_vars.get_ntimes())), cdf_vars.xb.values), axis=0))

    # Update the number of input points, if needed
    if not input_vars.options.input_points:
        input_vars.options.input_points = input_vars.get_nboundaries()

    # Get list of CDF variables to convert to the format needed for MMM
    # Independent variables listed below don't need to be converted
    cdf_var_list = cdf_vars.get_cdf_variables()

    for var_name in cdf_var_list:
        input_var = getattr(input_vars, var_name)
        convert_units(input_var)
        if input_var.values is not None and var_name not in ['time', 'x', 'xb']:
            _interp_to_boundarygrid(input_var, xvals)

    # Use TEPRO, TIPRO in place of TE, TI
    if input_vars.options.use_experimental_profiles:
        input_vars.use_experimental_profiles()

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
    input_vars.set_radius_values()
    mmm_vars = _interp_to_input_points(input_vars)
    _choose_variables(mmm_vars)

    full_var_list = mmm_vars.get_nonzero_variables()
    # Apply smoothing, then verify minimum values (fixes errors due to interpolation)
    for var_name in full_var_list:
        mmm_var = getattr(mmm_vars, var_name)
        if mmm_vars.options.apply_smoothing:
            mmm_var.apply_smoothing()

        # Since interpolation can create multiple expected nonphysical values,
        # no exceptions are raised for fixing these issues
        mmm_var.set_minvalue(ignore_exceptions=True)

    mmm_vars.set_x_values()
    mmm_vars.set_radius_values()

    return mmm_vars
