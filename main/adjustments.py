# Standard Packages
import sys; sys.path.insert(0, '../')
from math import log10
from copy import deepcopy

# 3rd Party Packages
import numpy as np

# Local Packages
import main.calculations as calculations
import main.options as options


VARIABLE_ERROR_TOLERANCE = 1e-8
SCAN_FACTOR_TOLERANCE = 1e-4
PRINT_FMT_STR = '>7.2f'


def print_factors(scan_factor, adjusted_factor):
    '''Prints the scan and adjusted factors'''
    print(f'{scan_factor:{PRINT_FMT_STR}}', f'{adjusted_factor:{PRINT_FMT_STR}}')


def get_nonzero_idx(values):
    '''
    Gets the indices of a nonzero value from values at the specified measurement time

    Parameters:
    * values (np.ndarray): An array of variable data

    Returns:
    * nonzero_values[0] (int): The index along the radial dimension of first nonzero variable value
    * time_idx (int): The index of the measurement time
    '''

    time_idx = options.instance.time_idx
    nonzero_values = np.where(values[:, time_idx] != 0)[0]
    if not len(nonzero_values):
        raise ValueError('Cannot adjust variable that is equal to 0 at all radial points')

    return nonzero_values[0], time_idx


def check_adjusted_factor(scan_factor, base_vals, new_vals):
    '''
    Checks that the adjusted factor equals scan_factor, within an allowable tolerance.
    An exception is raised if the check fails.

    Parameters:
    * adjusted_factor: The actual factor that a variable was adjusted by
    * scan_factor: The intended factor to adjust a variable by
    '''

    r, t = get_nonzero_idx(base_vals)
    adjusted_factor = new_vals[r, t] / base_vals[r, t]
    if abs(adjusted_factor / scan_factor - 1) > SCAN_FACTOR_TOLERANCE:
        var_to_scan = options.instance.var_to_scan
        fmt_length = abs(int(log10(SCAN_FACTOR_TOLERANCE)))
        raise ValueError(
            f'{var_to_scan} did not change within the allowable tolerance level\n'
            f'    adjusted_factor: {adjusted_factor}\n'
            f'    scan_factor:     {scan_factor}\n'
            f'    tolerance:       {SCAN_FACTOR_TOLERANCE:.{fmt_length}f}\n'
        )

    # For testing purposes (only executes when running this file directly)
    if __name__ == '__main__':
        print_factors(scan_factor, adjusted_factor)


def check_remained_constant(var_name, base_vars, new_vars):
    '''
    Raises an exception if the specified variable did not remain constant

    Parameters:
    * var_name (str): The name of the variable
    * base_vars (Variables): The base variables object
    * new_vars (Variables): The new variables object
    * r (int): The index along the radial dimension to check
    * t (int): The index along the time dimension to check
    '''

    base_var = getattr(base_vars, var_name)
    new_var = getattr(new_vars, var_name)
    r, t = get_nonzero_idx(base_var.values)
    variable_error = np.absolute(new_var.values[r, t] / base_var.values[r, t] - 1)
    if variable_error > VARIABLE_ERROR_TOLERANCE:
        fmt_length = abs(int(log10(VARIABLE_ERROR_TOLERANCE)))
        raise ValueError(
            f'{var_name} did not remain constant\n'
            f'    Variable Error: {variable_error}\n'
            f'    Allowed Error:  {VARIABLE_ERROR_TOLERANCE:.{fmt_length}f}'
        )


def adjust_nuei(mmm_vars, scan_factor):
    '''
    Collision Frequency Scan:
    * Since nuei is not an input variable, we must adjust ne and te to adjust nuei
    * Dividing ne and multiplying te by a adjustment_total will adjust nuei by the scan_factor
    * Our goal is to adjust by the scan_factor nuei, while keeping tau and pressure constant
    * Keeping tau constant also keeps alphamhd and gave constant
    * In order to keep tau and pressure constant, we must also adjust ti, nz, nd, nf, nh, and ni
    * These density adjustments keep average masses and zeff constant as well
    * Note that normalized gradients remain unchanged when multiplying their base variable by a constant factor
    * The adjustment finding loop typically finds the target adjustment within 6 attempts for a tolerance of 1e-4

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): Contains adjusted mmm_vars needed to write MMM input file
    '''

    max_adjustment_attempts = 10
    adjusted_vars = deepcopy(mmm_vars)
    r, t = get_nonzero_idx(mmm_vars.nuei.values)

    # The initial guess for the adjustment_total is based on the formula for nuei
    adjustment_total = scan_factor**(-2 / 5)
    adjustment_step = adjustment_total

    # Make adjustments to ne and te until the adjustment_total is found
    for i in range(max_adjustment_attempts):

        adjusted_vars.ne.values /= adjustment_step
        adjusted_vars.te.values *= adjustment_step

        # loge must be recalculated to properly calculate nuei
        calculations.loge(adjusted_vars)
        calculations.nuei(adjusted_vars)

        # Check if the current_factor is within our allowable tolerance for the scan_factor
        current_factor = adjusted_vars.nuei.values[r, t] / mmm_vars.nuei.values[r, t]
        if abs(current_factor / scan_factor - 1) < SCAN_FACTOR_TOLERANCE:
            # For testing purposes (only executes when running this file directly)
            if __name__ == '__main__':
                print(i + 1, end=' ')
            break

        if i + 1 == max_adjustment_attempts:
            print(i + 1, end=' ')
            print_factors(scan_factor, adjustment_total)
            raise ValueError(f'nuei factor could not be found after {i + 1} attempts')

        # Make an adjustment_step based on the difference in the current_factor and our desired scan_factor
        adjustment_step = 1 + (current_factor - scan_factor) / (2 * scan_factor)
        adjustment_total *= adjustment_step

    # Adjust ti to keep tau constant
    adjusted_vars.ti.values *= adjustment_total

    # Adjust densities needed to recalculate ni (to keep pressure constant)
    adjusted_vars.nz.values /= adjustment_total
    adjusted_vars.nd.values /= adjustment_total
    adjusted_vars.nf.values /= adjustment_total

    # Recalculate all variables
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    check_adjusted_factor(scan_factor, mmm_vars.nuei.values, adjusted_vars.nuei.values)
    check_remained_constant('p', mmm_vars, adjusted_vars)
    check_remained_constant('tau', mmm_vars, adjusted_vars)

    return adjusted_vars


def adjust_tau(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    # adjustment_total is based on the formula for tau
    adjustment_total = scan_factor**(1 / 2)

    # Adjust te and ti by the adjustment total to adjust tau
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.ti.values /= adjustment_total

    # Recalculate all variables
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    check_adjusted_factor(scan_factor, mmm_vars.tau.values, adjusted_vars.tau.values)

    return adjusted_vars


def adjust_zeff(mmm_vars, scan_factor):
    '''
    Effective Charge Scan:
    * Since zeff is a derived quantity, we adjust zeff by adjusting nz, and then update ne accordingly
    * This adjustment is nonlinear, so we don't have control over how much zeff changes at each point of rho
    * Consequently, all input variables are being recalculated after nz and ne are updated
    * These recalculations are our error checks, because exceptions will be thrown if nonphysical values emerge

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * scan_factor (float): The factor to modify nz by

    Returns:
    * adjusted_vars (InputVariables): Contains adjusted mmm_vars needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)

    # Adjust nz scan_factor, then increase ne by the change in nz
    adjusted_vars.nz.values *= scan_factor
    adjusted_vars.ne.values += adjusted_vars.zimp.values * (adjusted_vars.nz.values - mmm_vars.nz.values)

    # Recalculate all variables
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    # For testing purposes (only executes when running this file directly)
    if __name__ == '__main__':
        r, t = get_nonzero_idx(mmm_vars.zeff.values)
        adjusted_factor = adjusted_vars.zeff.values[r, t] / mmm_vars.zeff.values[r, t]
        print(f'{adjusted_factor:{PRINT_FMT_STR}}', end=' ')

    # nz is checked instead of zeff, since the change in zeff is different at each radial point
    check_adjusted_factor(scan_factor, mmm_vars.nz.values, adjusted_vars.nz.values)

    return adjusted_vars


def adjust_etae(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    # adjustment_total is based on the formula for etae
    adjustment_total = scan_factor**(1 / 2)
    adjusted_vars.gte.values *= adjustment_total
    adjusted_vars.gne.values /= adjustment_total

    # Only recalculate additional variables, since recalculating all would undo the gradient adjustments
    calculations.calculate_additional_variables(adjusted_vars)

    check_adjusted_factor(scan_factor, mmm_vars.etae.values, adjusted_vars.etae.values)

    return adjusted_vars


def adjust_shear(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    # Adjust gq by the scan_factor to adjust shear
    adjusted_vars.gq.values *= scan_factor

    # Only recalculate additional variables, since recalculating all would undo the gradient adjustment
    calculations.calculate_additional_variables(adjusted_vars)

    check_adjusted_factor(scan_factor, mmm_vars.shear.values, adjusted_vars.shear.values)

    return adjusted_vars


def adjust_btor(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    # Adjust bz by the scan_factor to adjust btor
    adjusted_vars.bz.values *= scan_factor

    # Recalculate all variables
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    check_adjusted_factor(scan_factor, mmm_vars.btor.values, adjusted_vars.btor.values)

    return adjusted_vars


def adjust_betae(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    # Our adjustment_total is based on the formula for betae
    adjustment_total = scan_factor**(1 / 4)

    adjusted_vars.ne.values *= adjustment_total
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.bz.values /= adjustment_total

    # adjust dependencies of ne
    adjusted_vars.nd.values *= adjustment_total
    adjusted_vars.nz.values *= adjustment_total
    adjusted_vars.nf.values *= adjustment_total

    # Recalculate all variables
    adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    check_adjusted_factor(scan_factor, mmm_vars.betae.values, adjusted_vars.betae.values)

    return adjusted_vars


def adjust_scanned_variable(mmm_vars, var_to_scan, scan_factor):
    '''
    Adjusts the variable being scanned, as well as any necessary dependencies

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * var_to_scan (str): The name of the variable being adjusted
    * scan_factor (float): The factor to modify var_to_scan by

    Returns:
    * adjusted_vars (InputVariables): All adjusted and unaltered mmm_vars needed to write MMM input file
    '''

    if var_to_scan == 'nuei':
        '''Collision Frequency Scan'''
        adjusted_vars = adjust_nuei(mmm_vars, scan_factor)

    elif var_to_scan == 'zeff':
        '''Effective Charge Scan'''
        adjusted_vars = adjust_zeff(mmm_vars, scan_factor)

    elif var_to_scan == 'tau':
        '''Temperature Ratio Scan'''
        adjusted_vars = adjust_tau(mmm_vars, scan_factor)

    elif var_to_scan == 'etae':
        '''etae = gte/gne scan'''
        adjusted_vars = adjust_etae(mmm_vars, scan_factor)

    elif var_to_scan == 'shear' or var_to_scan == 'gq':
        '''Shear Scan (gq scan)'''
        adjusted_vars = adjust_shear(mmm_vars, scan_factor)

    elif var_to_scan == 'btor' or var_to_scan == 'bz':
        '''Toroidal Magnetic Field Scan (bz scan)'''
        adjusted_vars = adjust_btor(mmm_vars, scan_factor)

    elif var_to_scan == 'betae':
        '''Betae Scan'''
        adjusted_vars = adjust_betae(mmm_vars, scan_factor)

    else:
        '''Simple Scan (no advanced logic needed)'''
        adjusted_vars = deepcopy(mmm_vars)
        base_var = getattr(mmm_vars, var_to_scan)
        scanned_var = getattr(adjusted_vars, var_to_scan)
        scanned_var.values = scan_factor * base_var.values

        '''
        Simple Gradient Adjustments:
        * Only recalculate additional variables, since recalculating all would undo the gradient adjustment
        * Non-gradient input variables do not depend on gradient values
        '''
        if var_to_scan in ['gte', 'gti', 'gne', 'gnh', 'gni', 'gnz', 'gvpar', 'gvpol', 'gvtor']:
            calculations.calculate_additional_variables(adjusted_vars)

        else:
            adjusted_vars = calculations.calculate_inputs(adjusted_vars)

    return adjusted_vars


# For Testing Purposes
if __name__ == '__main__':
    from utils import initialize_variables
    opts = options.instance
    opts.set(
        runid='138536A01',
        input_points=51,
        apply_smoothing=True,
        uniform_rho=False,
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
        opts.set(var_to_scan=var_name)
        print(f'\n{var_name} scan factors:')
        for scan_factor in scan_range:
            adjust_scanned_variable(mmm_vars, opts.var_to_scan, scan_factor)
