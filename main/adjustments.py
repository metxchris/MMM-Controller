# Standard Packages
import sys; sys.path.insert(0, '../')
from copy import deepcopy

# 3rd Party Packages
import numpy as np

# Local Packages
import main.calculations as calculations
import main.options as options


VARIABLE_ERROR_TOLERANCE = 1e-8
SCAN_FACTOR_TOLERANCE = 1e-4


def check_adjusted_factor(adjusted_factor, scan_factor):
    if abs(adjusted_factor / scan_factor - 1) > SCAN_FACTOR_TOLERANCE:
        print(scan_factor, adjusted_factor)
        var_to_scan = options.Options.instance
        raise ValueError(f'{var_to_scan} did not change within the allowable tolerance level')
    # For testing purposes (only executes when running this file directly)
    elif __name__ == '__main__':
        print(round(scan_factor, 4), round(adjusted_factor, 4))


def get_nonzero_idx(var):
    nonzero_values = np.where(var.values[:, 0] != 0)[0]
    if not len(nonzero_values):
        raise ValueError('Cannot adjust variable that is equal to 0 everywhere')

    return nonzero_values[0]


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
    * The adjustment finding loop typically finds the target adjustment within 5 attempts for a tolerance of 1e-4

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): Contains adjusted mmm_vars needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)
    max_adjustment_attempts = 10

    nonzero_idx = get_nonzero_idx(mmm_vars.nuei)

    # Store initial variable values needed for calculating the total adjustment_step, and error checking
    nuei0 = mmm_vars.nuei.values[nonzero_idx, 0]
    p0 = mmm_vars.p.values[nonzero_idx, 0]
    tau0 = mmm_vars.tau.values[nonzero_idx, 0]

    # Our initial guess for the adjustment_total is based on the formula for nuei
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
        current_factor = adjusted_vars.nuei.values[nonzero_idx, 0] / nuei0
        if abs(current_factor / scan_factor - 1) < SCAN_FACTOR_TOLERANCE:
            # For testing purposes
            if __name__ == '__main__':
                print(i, round(scan_factor, 4), round(current_factor, 4))
            break

        if i + 1 == max_adjustment_attempts:
            print(i, scan_factor, adjustment_total)
            raise ValueError(f'nuei factor could not be found after {i + 1} attempts.')

        # Make an adjustment_step based on the difference in the current_factor and our desired scan_factor
        adjustment_step = 1 + (current_factor - scan_factor) / (2 * scan_factor)
        adjustment_total *= adjustment_step

    # Adjust ti to keep tau constant
    adjusted_vars.ti.values *= adjustment_total

    # Adjust densities needed to recalculate ni
    adjusted_vars.nz.values /= adjustment_total
    adjusted_vars.nd.values /= adjustment_total
    adjusted_vars.nf.values /= adjustment_total

    # Recalculate ni to keep pressure constant
    calculations.calculate_variable(calculations.nh0, adjusted_vars)
    calculations.calculate_variable(calculations.nh, adjusted_vars)
    calculations.calculate_variable(calculations.ni, adjusted_vars)

    # Recalculate constant variables for error checking
    calculations.calculate_variable(calculations.p, adjusted_vars)
    calculations.calculate_variable(calculations.tau, adjusted_vars)

    # Take ratios of recalculated constant variables and check for errors
    variable_ratio = np.array((
        adjusted_vars.p.values[nonzero_idx, 0] / p0,
        adjusted_vars.tau.values[nonzero_idx, 0] / tau0,
    ))
    error_check = np.absolute(variable_ratio - 1) > VARIABLE_ERROR_TOLERANCE
    if error_check.any():
        if error_check[0]:
            raise ValueError(f'Pressure did not remain constant for scan factor {scan_factor}: {variable_ratio[0]}')
        if error_check[1]:
            raise ValueError(f'Tau did not remain constant for scan factor {scan_factor}: {variable_ratio[1]}')

    return adjusted_vars


def adjust_tau(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    nonzero_idx = get_nonzero_idx(mmm_vars.tau)
    tau0 = mmm_vars.tau.values[nonzero_idx, 0]

    # Our adjustment_total is based on the formula for tau
    adjustment_total = scan_factor**(1 / 2)

    # Adjust te and ti by the adjustment total to adjust tau
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.ti.values /= adjustment_total

    calculations.calculate_variable(calculations.tau, adjusted_vars)

    # Check that the adjustment is correct
    adjusted_factor = adjusted_vars.tau.values[nonzero_idx, 0] / tau0
    check_adjusted_factor(adjusted_factor, scan_factor)

    return adjusted_vars


def adjust_etae(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    nonzero_idx = get_nonzero_idx(mmm_vars.etae)
    etae0 = mmm_vars.etae.values[nonzero_idx, 0]

    # Our adjustment_total is based on the formula for etae
    adjustment_total = scan_factor**(1 / 2)

    # Adjust te and ti by the adjustment total to adjust etae
    adjusted_vars.gte.values *= adjustment_total
    adjusted_vars.gne.values /= adjustment_total

    calculations.calculate_variable(calculations.etae, adjusted_vars)

    # Check that the adjustment is correct
    adjusted_factor = adjusted_vars.etae.values[nonzero_idx, 0] / etae0
    check_adjusted_factor(adjusted_factor, scan_factor)

    return adjusted_vars


def adjust_shear(mmm_vars, scan_factor):
    adjusted_vars = deepcopy(mmm_vars)

    nonzero_idx = get_nonzero_idx(mmm_vars.shear)
    shear0 = mmm_vars.shear.values[nonzero_idx, 0]

    # Adjust gq by the scan_factor to adjust shear (no need for an adjustment total)
    adjusted_vars.gq.values *= scan_factor
    calculations.calculate_variable(calculations.shear, adjusted_vars)

    # Check that the adjustment is correct
    adjusted_factor = adjusted_vars.shear.values[nonzero_idx, 0] / shear0
    check_adjusted_factor(adjusted_factor, scan_factor)

    return adjusted_vars


def adjust_betae(mmm_vars, scan_factor):

    adjusted_vars = deepcopy(mmm_vars)

    nonzero_idx = get_nonzero_idx(mmm_vars.betae)
    betae0 = mmm_vars.betae.values[nonzero_idx, 0]

    # Our adjustment_total is based on the formula for betae
    adjustment_total = scan_factor**(1 / 4)

    # Adjust te and ti by the adjustment total to adjust etae
    adjusted_vars.ne.values *= adjustment_total
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.btor.values /= adjustment_total

    calculations.calculate_variable(calculations.betae, adjusted_vars)

    # Check that the adjustment is correct
    adjusted_factor = adjusted_vars.betae.values[nonzero_idx, 0] / betae0
    check_adjusted_factor(adjusted_factor, scan_factor)

    return adjusted_vars


def recalculate_dependencies(adjusted_vars, var_to_scan):
    '''
    Recalculates variables dependent on previously made adjustments as per the scan

    Parameters:
    * adjusted_vars (InputVariables): Contains adjusted mmm_vars needed to write MMM input file
    * var_to_scan (str): The variable being scanned
    '''

    if var_to_scan == 'nuei':
        '''Collision Frequency Scan'''
        calculations.calculate_variable(calculations.nuei2, adjusted_vars)
        calculations.calculate_variable(calculations.vthe, adjusted_vars)
        calculations.calculate_variable(calculations.vthi, adjusted_vars)
        calculations.calculate_variable(calculations.nuste, adjusted_vars)
        calculations.calculate_variable(calculations.nusti, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
    elif var_to_scan == 'gte':
        '''Electron Temperature Gradient Scan'''
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
        calculations.calculate_variable(calculations.etae, adjusted_vars)
    elif var_to_scan == 'gti':
        '''Ion Temperature Gradient Scan'''
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
        calculations.calculate_variable(calculations.etai, adjusted_vars)
    elif var_to_scan == 'q':
        '''
        Safety Factor Scan
        * vpar is not recalculated since it's considered as an independent input variable
        '''
        calculations.calculate_variable(calculations.bpol, adjusted_vars)
        calculations.calculate_variable(calculations.nuste, adjusted_vars)
        calculations.calculate_variable(calculations.nusti, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
    elif var_to_scan == 'zeff':
        '''Effective Charge Scan'''
        calculations.calculate_variable(calculations.loge, adjusted_vars)
        calculations.calculate_variable(calculations.nuei, adjusted_vars)
        calculations.calculate_variable(calculations.nuei2, adjusted_vars)
        calculations.calculate_variable(calculations.nuste, adjusted_vars)
        calculations.calculate_variable(calculations.nusti, adjusted_vars)
    elif var_to_scan == 'gnh':
        '''Hydrogenic Density Gradient Scan'''
        ...
    elif var_to_scan == 'gnz':
        '''Impurity Density Gradient Scan (No dependencies)'''
        ...
    elif var_to_scan == 'tau':
        calculations.calculate_variable(calculations.p, adjusted_vars)
        calculations.calculate_variable(calculations.betae, adjusted_vars)
        calculations.calculate_variable(calculations.beta, adjusted_vars)
        calculations.calculate_variable(calculations.loge, adjusted_vars)
        calculations.calculate_variable(calculations.nuei, adjusted_vars)
        calculations.calculate_variable(calculations.nuei2, adjusted_vars)
        calculations.calculate_variable(calculations.vthe, adjusted_vars)
        calculations.calculate_variable(calculations.vthi, adjusted_vars)
        calculations.calculate_variable(calculations.nuste, adjusted_vars)
        calculations.calculate_variable(calculations.nusti, adjusted_vars)
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
    elif var_to_scan == 'btor':
        '''
        Toroidal Magnetic Field Scan
        * vpar is not recalculated since it's considered as an independent input variable
        '''
        calculations.calculate_variable(calculations.bpol, adjusted_vars)
        calculations.calculate_variable(calculations.betae, adjusted_vars)
        calculations.calculate_variable(calculations.beta, adjusted_vars)
        calculations.calculate_variable(calculations.gyrfi, adjusted_vars)
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
    elif var_to_scan == 'etae':
        '''etae = gte/gne scan'''
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
    elif var_to_scan == 'gnh':
        '''Hydrogenic Ion Density Gradient Scan'''
        ...
    elif var_to_scan == 'q':
        '''
        Safety Factor Scan
        * vpar is not recalculated since it's considered as an independent input variable
        '''
        calculations.calculate_variable(calculations.bpol, adjusted_vars)
        calculations.calculate_variable(calculations.nuste, adjusted_vars)
        calculations.calculate_variable(calculations.nusti, adjusted_vars)
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
    elif var_to_scan == 'shear':
        '''Shear Scan'''
        calculations.calculate_variable(calculations.shat, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)
    elif var_to_scan == 'betae':
        '''
        Beta_e Scan
        * zeff, vpar are intentionally not recalculated
        '''
        calculations.calculate_variable(calculations.tau, adjusted_vars)
        calculations.calculate_variable(calculations.bpol, adjusted_vars)
        calculations.calculate_variable(calculations.p, adjusted_vars)
        calculations.calculate_variable(calculations.beta, adjusted_vars)
        calculations.calculate_variable(calculations.loge, adjusted_vars)
        calculations.calculate_variable(calculations.nuei, adjusted_vars)
        calculations.calculate_variable(calculations.nuei2, adjusted_vars)
        calculations.calculate_variable(calculations.vthe, adjusted_vars)
        calculations.calculate_variable(calculations.gyrfi, adjusted_vars)
        calculations.calculate_variable(calculations.alphamhd, adjusted_vars)
        calculations.calculate_variable(calculations.gmax, adjusted_vars)
        calculations.calculate_variable(calculations.gave, adjusted_vars)


def adjust_scanned_variable(mmm_vars, var_to_scan, scan_factor):
    '''
    Adjusts the variable being scanned, as well as any necessary dependencies

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * var_to_scan (str): The name of the variable being adjusted
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): All adjusted and unaltered mmm_vars needed to write MMM input file
    '''

    if var_to_scan == 'nuei':
        '''Collision Frequency Scan'''
        adjusted_vars = adjust_nuei(mmm_vars, scan_factor)
    elif var_to_scan == 'tau':
        '''Temperature Ratio Frequency Scan'''
        adjusted_vars = adjust_tau(mmm_vars, scan_factor)
    elif var_to_scan == 'etae':
        '''etae = gte/gne scan'''
        adjusted_vars = adjust_etae(mmm_vars, scan_factor)
    elif var_to_scan == 'shear' or var_to_scan == 'gq':
        '''Shear Scan (gq scan)'''
        adjusted_vars = adjust_shear(mmm_vars, scan_factor)
    elif var_to_scan == 'betae':
        '''Betae Scan'''
        adjusted_vars = adjust_betae(mmm_vars, scan_factor)
    else:
        '''Simple Scan (no advanced logic needed)'''
        adjusted_vars = deepcopy(mmm_vars)
        base_var = getattr(mmm_vars, var_to_scan)
        scanned_var = getattr(adjusted_vars, var_to_scan)
        scanned_var.values = scan_factor * base_var.values

    recalculate_dependencies(adjusted_vars, var_to_scan)

    return adjusted_vars


# For Testing Purposes
if __name__ == '__main__':
    from utils import initialize_variables
    opts = options.Options.instance
    opts.set(
        var_to_scan='betae',
        runid='120982A09',
        input_points=51,
        apply_smoothing=True,
        uniform_rho=False,
        input_time=0,  # Value isn't used here, but needs to be some number
    )

    mmm_vars, __, __ = initialize_variables()

    # Check that all scan_factors can be found in scan_range (failures will raise a ValueError)
    scan_range = np.hstack((np.arange(1e-6, 5, 0.05), np.arange(5, 25, 1), np.arange(25, 105, 5)))
    for scan_factor in scan_range:
        adjust_scanned_variable(mmm_vars, opts.var_to_scan, scan_factor)
