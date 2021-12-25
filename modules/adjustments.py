"""Adjusts a variable during a variable scan.

Variable Classifications:
* Base variable: A non-gradient variable that is either used as input to MMM,
  or is needed to calculate an MMM input variable.  There are no base
  variables that depend on gradient values.
* Gradient variable: A gradient variable that is used as input to MMM.
* Additional variable: A variable that is not needed for MMM input, or input
  calculations. Additional variables can depend on gradient values.

Adjustment Methodology:
* To keep this process as physical as possible, we recalculate all variables
  as necessary following an adjustment. However, a full recalculation will
  overwrite any gradient adjustments made, so we only recalculate additional
  variables in this case.  As such, gradient adjustments are somewhat
  nonphysical, but the error from this will be minimal.  For example, we
  would recalculate all variables following an adjustment to electron
  temperature (te), but would only recalculate additional variables following
  an adjustment to the electron temperature gradient (gte).
* We note that recalculating every variable, even those that would have
  clearly not changed from an adjustment, does increase the time needed to
  conduct the adjustment process.  However, it also eliminates human error
  associated with choosing only specific variables that need to be
  recalculated.  Furthermore, recalculating everything also helps with error
  checking, in the event that some variable values change that were supposed
  to remain constant.
* When adjusting a derived quantity, adjustments are made directly to the
  variables that determine that derived quantity.  When multiple variables
  need to be adjusted in this manner, each variable is adjusted by the same
  base amount (either multiplied or divided by the same number).  For
  example, if we wanted to double the value of the temperature ratio (tau =
  te / ti), the adjustment factor would be sqrt(2); te would be multiplied by
  this factor, and ti would be divided by it.  If any variables have powers
  attached to them, then the adjustment factor is updated as needed, but all
  variables are still adjusted by the same factor.

Example Usage:
  * new_vars = adjust_scanned_variable(original_vars, 'tau', 2.5)
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
from math import log10
from copy import deepcopy

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.calculations as calculations
import modules.options as options


EQUALITY_ERROR_TOLERANCE = 1e-8
SCAN_FACTOR_TOLERANCE = 1e-4
PRINT_FMT_STR = '>7.2f'


def _print_factors(scan_factor, adjusted_factor):
    '''Prints the scan factor and adjusted factor'''
    print(f'{scan_factor:{PRINT_FMT_STR}}', f'{adjusted_factor:{PRINT_FMT_STR}}')


def _get_nonzero_idx(values):
    '''
    Gets the indices of a nonzero value from values at the previously specified measurement time

    Parameters:
    * values (np.ndarray): An array of variable data

    Returns:
    * nonzero_values[0] (int): The index along the radial dimension of first nonzero variable value
    * time_idx (int): The index of the measurement time

    Raises:
    * ValueError: If there are no nonzero values
    '''

    time_idx = options.instance.time_idx
    nonzero_values = np.where(values[:, time_idx] != 0)[0]
    if not len(nonzero_values):
        raise ValueError('Cannot adjust variable that is equal to 0 at all radial points')

    return nonzero_values[0], time_idx


def _check_adjusted_factor(scan_factor, base_var, adjusted_var):
    '''
    Checks that the adjusted factor equals scan_factor, within an allowable tolerance

    Parameters:
    * scan_factor (float): The intended factor that a variable was adjusted by
    * base_var (Variable): The unmodified variable object
    * adjusted_var (Variable): The adjusted variable object

    Raises:
    * ValueError: If the adjusted factor does not equal the scan factor within allowable tolerance
    '''

    r, t = _get_nonzero_idx(base_var.values)
    adjusted_factor = adjusted_var.values[r, t] / base_var.values[r, t]
    if abs(adjusted_factor / scan_factor - 1) > SCAN_FACTOR_TOLERANCE:
        var_to_scan = options.instance.var_to_scan
        fmt_length = abs(int(log10(SCAN_FACTOR_TOLERANCE)))
        raise ValueError(
            f'{var_to_scan} did not change within the allowable tolerance level\n'
            f'    adjusted_factor: {adjusted_factor}\n'
            f'    scan_factor:     {scan_factor}\n'
            f'    tolerance:       {SCAN_FACTOR_TOLERANCE:.{fmt_length}f}\n'
        )

    if __name__ == '__main__':  # For testing purposes
        _print_factors(scan_factor, adjusted_factor)


def _check_equality(base_var, adjusted_var):
    '''
    Checks if two variables are equal in value

    This equality check assumes that a linear adjustment has been made, so
    only one point of each variable is checked to be equal.  This point is
    specifically chosen to be nonzero in value, since adjustments to 0 would
    be trivial, and would also create divide by 0 errors in the error check.

    Parameters:
    * base_var (Variable): The unmodified variable
    * adjusted_var (Variable): The adjusted variable

    Raises:
    * ValueError: If the two variables are not equal within allowable tolerance
    '''

    r, t = _get_nonzero_idx(base_var.values)
    variable_error = np.absolute(adjusted_var.values[r, t] / base_var.values[r, t] - 1)
    if variable_error > EQUALITY_ERROR_TOLERANCE:
        fmt_length = abs(int(log10(EQUALITY_ERROR_TOLERANCE)))
        raise ValueError(
            f'{base_var.name} did not remain constant\n'
            f'    Variable Error: {variable_error}\n'
            f'    Allowed Error:  {EQUALITY_ERROR_TOLERANCE:.{fmt_length}f}'
        )


def _adjust_nuei(mmm_vars, scan_factor):
    '''
    Adjusts the Collision Frequency

    Adjustment Details:
    * nuei is adjusted indirectly by adjusting the values of ne and te by
      adjustment_total. The final value of adjustment_total is obtained by
      running a loop until nuei is adjusted by the target scan factor.
    * As part of this scan, we also require that pressure (p) and tau remain
      constant as nuei is adjusted. Consequently, we also adjust ti, nz, nd,
      nf, nh, and ni accordingly, and then check that both p and tau have
      remained constant afterwards.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file

    Raises:
    * ValueError: If the total adjustment could not be found in the specified max attempts
    '''

    max_adjustment_attempts = 10
    r, t = _get_nonzero_idx(mmm_vars.nuei.values)

    adjusted_vars = deepcopy(mmm_vars)
    adjustment_total = scan_factor**(-2 / 5)  # based on the formula for nuei
    adjustment_step = adjustment_total

    # Make adjustments to ne and te until the adjustment_total is found
    # The adjustment finding loop typically finds the target within 6 attempts for a tolerance of 1e-4
    for i in range(max_adjustment_attempts):

        adjusted_vars.ne.values /= adjustment_step
        adjusted_vars.te.values *= adjustment_step

        # loge must be recalculated to properly recalculate nuei
        calculations.loge(adjusted_vars)
        calculations.nuei(adjusted_vars)

        # Check if the current_factor is within our allowable tolerance for the scan_factor
        current_factor = adjusted_vars.nuei.values[r, t] / mmm_vars.nuei.values[r, t]
        if abs(current_factor / scan_factor - 1) < SCAN_FACTOR_TOLERANCE:
            if __name__ == '__main__':  # For testing purposes
                print(i + 1, end=' ')
            break

        if i + 1 == max_adjustment_attempts:
            print(i + 1, end=' ')
            _print_factors(scan_factor, adjustment_total)
            raise ValueError(f'nuei factor could not be found after {i + 1} attempts')

        # Make an adjustment_step based on the difference in the current_factor and our desired scan_factor
        adjustment_step = 1 + (current_factor - scan_factor) / (2 * scan_factor)
        adjustment_total *= adjustment_step

    # Additional adjustments to keep tau and pressure constant
    adjusted_vars.ti.values *= adjustment_total
    adjusted_vars.nz.values /= adjustment_total
    adjusted_vars.nd.values /= adjustment_total
    adjusted_vars.nf.values /= adjustment_total

    adjusted_vars = calculations.calculate_inputs(adjusted_vars)
    _check_adjusted_factor(scan_factor, mmm_vars.nuei, adjusted_vars.nuei)
    _check_equality(mmm_vars.p, adjusted_vars.p)
    _check_equality(mmm_vars.tau, adjusted_vars.tau)

    return adjusted_vars


def _adjust_tau(mmm_vars, scan_factor):
    '''
    Adjusts the temperature ratio

    Adjustment Details:
    * Tau is adjusted indirectly by adjusting te and ti by adjustment_total

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify tau by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    adjustment_total = scan_factor**(1 / 2)  # based on the formula for tau
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.ti.values /= adjustment_total
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)
    _check_adjusted_factor(scan_factor, mmm_vars.tau, adjusted_vars.tau)

    return adjusted_vars


def _adjust_zeff(mmm_vars, scan_factor):
    '''
    Adjusts the Effective Charge

    Adjustment Details:
    * zeff is adjusted indirectly by adjusting the values of nz by
      scan_factor. Since ne depends on nz, ne is updated by the change in nz
      values (ne += z*delta(nz)).  However, the values of nh and nd are
      intentionally held constant, and are not updated from the adjustment of ne.
    * Unlike other adjustments, the amount that zeff changes is different at
      each point of rho (nonlinear). Therefore, the adjusted factor is
      checked for the change in nz after all recalculations, instead of for
      the change in zeff.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify nz by (not zeff)

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    adjusted_vars.nz.values *= scan_factor
    adjusted_vars.ne.values += adjusted_vars.zimp.values * (adjusted_vars.nz.values - mmm_vars.nz.values)
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    if __name__ == '__main__':  # For testing purposes
        r, t = _get_nonzero_idx(mmm_vars.zeff.values)
        adjusted_factor = adjusted_vars.zeff.values[r, t] / mmm_vars.zeff.values[r, t]
        print(f'{adjusted_factor:{PRINT_FMT_STR}}', end=' ')

    _check_adjusted_factor(scan_factor, mmm_vars.nz, adjusted_vars.nz)

    return adjusted_vars


def _adjust_etae(mmm_vars, scan_factor):
    '''
    Adjust etae = gte/gne

    Adjustment Details:
    * etae is adjusted indirectly by adjusting gte and gne by adjustment_total

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify etae by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    adjustment_total = scan_factor**(1 / 2)  # based on the formula for etae
    adjusted_vars.gte.values *= adjustment_total
    adjusted_vars.gne.values /= adjustment_total

    # Only additional variables are recalculated, so that the gte, gne adjustments aren't overwritten
    calculations.calculate_additional_variables(adjusted_vars)
    _check_adjusted_factor(scan_factor, mmm_vars.etae, adjusted_vars.etae)

    return adjusted_vars


def _adjust_shear(mmm_vars, scan_factor):
    '''
    Adjust Magnetic Shear

    Adjustment Details:
    * shear is adjusted indirectly by adjusting gq by scan_factor

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify shear by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    adjusted_vars.gq.values *= scan_factor

    # Only additional variables are recalculated, so that the gq adjustment isn't overwritten
    calculations.calculate_additional_variables(adjusted_vars)
    _check_adjusted_factor(scan_factor, mmm_vars.shear, adjusted_vars.shear)

    return adjusted_vars


def _adjust_btor(mmm_vars, scan_factor):
    '''
    Adjust Toroidal Magnetic Field

    Adjustment Details:
    * btor is adjusted indirectly by adjusting bz by scan_factor

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify btor by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    adjusted_vars.bz.values *= scan_factor
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)
    _check_adjusted_factor(scan_factor, mmm_vars.btor, adjusted_vars.btor)

    return adjusted_vars


def _adjust_betae(mmm_vars, scan_factor):
    '''
    Adjust Electron Pressure Ratio

    Adjustment Details:
    * betae is adjusted indirectly by adjusting ne, te, and bz by
      adjustment_total. The values of nd, nz, and nf depend on changes in ne,
      so these are adjusted as well

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify betae by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    adjustment_total = scan_factor**(1 / 4)  # based on the formula for betae
    adjusted_vars.ne.values *= adjustment_total
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.bz.values /= adjustment_total
    adjusted_vars.nd.values *= adjustment_total
    adjusted_vars.nz.values *= adjustment_total
    adjusted_vars.nf.values *= adjustment_total
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)
    _check_adjusted_factor(scan_factor, mmm_vars.betae, adjusted_vars.betae)

    return adjusted_vars


def adjust_scanned_variable(mmm_vars, var_to_scan, scan_factor):
    '''
    Adjusts the variable being scanned, as well as any necessary dependencies

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * var_to_scan (str): The name of the variable being adjusted
    * scan_factor (float): The factor to modify var_to_scan by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    if var_to_scan == 'nuei':
        adjusted_vars = _adjust_nuei(mmm_vars, scan_factor)

    elif var_to_scan == 'zeff':
        adjusted_vars = _adjust_zeff(mmm_vars, scan_factor)

    elif var_to_scan == 'tau':
        adjusted_vars = _adjust_tau(mmm_vars, scan_factor)

    elif var_to_scan == 'etae':
        adjusted_vars = _adjust_etae(mmm_vars, scan_factor)

    elif var_to_scan == 'shear' or var_to_scan == 'gq':
        adjusted_vars = _adjust_shear(mmm_vars, scan_factor)

    elif var_to_scan == 'btor' or var_to_scan == 'bz':
        adjusted_vars = _adjust_btor(mmm_vars, scan_factor)

    elif var_to_scan == 'betae':
        adjusted_vars = _adjust_betae(mmm_vars, scan_factor)

    else:
        # Simple Scan (no advanced logic needed)
        adjusted_vars = deepcopy(mmm_vars)
        base_var = getattr(mmm_vars, var_to_scan)
        scanned_var = getattr(adjusted_vars, var_to_scan)
        scanned_var.values = scan_factor * base_var.values

        if var_to_scan in ['gte', 'gti', 'gne', 'gnh', 'gni', 'gnz', 'gvpar', 'gvpol', 'gvtor']:
            # Only recalculate additional variables when adjusting gradients (see module docstring)
            calculations.calculate_additional_variables(adjusted_vars)

        else:
            adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    return adjusted_vars


if __name__ == '__main__':  # For Testing Purposes
    from modules.utils import initialize_variables
    options.instance.set(
        runid='138536A01',
        input_points=51,
        apply_smoothing=True,
        uniform_rho=True,
        input_time=.63,
    )

    mmm_vars, __, __ = initialize_variables()

    # Check that all scan_factors can be found in scan_range (failures will raise a ValueError)
    scan_range = np.hstack(
        (np.arange(0.05, 1, 0.05),
         np.arange(1, 5, 0.2),
         np.arange(5, 20, 1),
         np.arange(20, 105, 5)),
    )

    advanced_scans = ['nuei', 'zeff', 'tau', 'etae', 'shear', 'btor', 'betae']
    for var_name in advanced_scans:
        options.instance.set(var_to_scan=var_name)
        print(f'\n{var_name} scan factors:')
        for scan_factor in scan_range:
            adjust_scanned_variable(mmm_vars, var_name, scan_factor)
