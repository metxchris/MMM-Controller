# Standard Packages
from copy import deepcopy
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
from main import calculations


def adjust_nuei(mmm_vars, scan_factor):
    '''
    Collision Frequency Scan:
    * Since nuei is not an input variable, we must adjust ne and te to adjust nuei
    * Dividing ne and multiplying te by a adjustment_total will adjust nuei by the scan_factor
    * Our goal is to adjust by the scan_factor nuei, while keeping tau and pressure constant
    * Keeping tau constant also keeps alphamhd and gave constant
    * In order to keep tau and pressure constant, we must also adjust ti, nz, nd, nf, nh, and ni
    * These density adjustments keep average masses and zeff constant as well
    * Note that normalized gradients are constant when multiplying their base by a constant factor
    * The adjustment finding loop typically finds the target adjustment within 4 attemps for a tolerance of 1e-3

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): All adjusted and unaltered mmm_vars needed to write MMM input file
    '''

    adjusted_vars = deepcopy(mmm_vars)

    max_adjustment_attempts = 10
    scan_factor_tolerance = 1e-3
    variable_error_tolerance = 1e-8

    # Store initial variable values needed for calculating the total adjustment_step, and error checking
    nuei0 = mmm_vars.nuei.values[0, 0]
    p0 = mmm_vars.p.values[0, 0]
    tau0 = mmm_vars.tau.values[0, 0]

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
        current_factor = adjusted_vars.nuei.values[0, 0] / nuei0
        if abs(current_factor / scan_factor - 1) < scan_factor_tolerance:
            # For testing purposes
            if __name__ == '__main__':
                print(i, round(current_factor, 4), round(scan_factor, 4))
            break

        if i + 1 == max_adjustment_attempts:
            print(i, scan_factor, adjustment_total)
            raise ValueError(f'nuei factor could not be found after {i + 1} attempts.')

        # Make an adjustment_step based on the difference in the current_factor and our desired scan_factor
        adjustment_step = 1 + (current_factor - scan_factor) / (2 * scan_factor)
        adjustment_total *= adjustment_step

    # Adjust ti to keep tau constant
    adjusted_vars.ti.values *= adjustment_total

    # Adjust densities to adjust ni, which is needed to keep pressure constant
    adjusted_vars.nz.values /= adjustment_total
    adjusted_vars.nd.values /= adjustment_total
    adjusted_vars.nf.values /= adjustment_total

    # Recalculate dependent variables
    calculations.calculate_variable(calculations.nh0, adjusted_vars)
    calculations.calculate_variable(calculations.nh, adjusted_vars)
    calculations.calculate_variable(calculations.ni, adjusted_vars)

    # Recalculate variables that are dependent on previous changes
    calculations.calculate_variable(calculations.vthe, adjusted_vars)
    calculations.calculate_variable(calculations.vthi, adjusted_vars)
    calculations.calculate_variable(calculations.nuste, adjusted_vars)
    calculations.calculate_variable(calculations.nusti, adjusted_vars)

    # Recalculate constant variables for error checking
    calculations.calculate_variable(calculations.p, adjusted_vars)
    calculations.calculate_variable(calculations.tau, adjusted_vars)

    # Take ratios of recalculated constant variables and check for errors
    variable_ratio = np.array((
        adjusted_vars.p.values[0, 0] / p0,
        adjusted_vars.tau.values[0, 0] / tau0,
    ))
    error_check = np.absolute(variable_ratio - 1) > variable_error_tolerance
    if error_check.any():
        if error_check[0]:
            raise ValueError(f'Pressure did not remain constant for scan factor {scan_factor}: {variable_ratio[0]}')
        if error_check[1]:
            raise ValueError(f'Tau did not remain constant for scan factor {scan_factor}: {variable_ratio[1]}')

    return adjusted_vars


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
    else:
        '''Simple Scan (no dependent recalculations needed)'''
        adjusted_vars = deepcopy(mmm_vars)
        base_var = getattr(mmm_vars, var_to_scan)
        scanned_var = getattr(adjusted_vars, var_to_scan)
        scanned_var.values = scan_factor * base_var.values

    return adjusted_vars


# For Testing Purposes
if __name__ == '__main__':
    from mmm_controller import initialize_variables
    from main.options import Options

    Options.instance.set(
        var_to_scan='nuei',
        runid='120982A09',
        input_points=51,
        apply_smoothing=True,
        uniform_rho=False,
        input_time=0,  # Value isn't used here, but needs to be some number
    )

    mmm_vars, _, _ = initialize_variables()

    # Check that all scan_factors can be found in scan_range (failures will raise a ValueError)
    scan_range = np.hstack((np.arange(1e-6, 5, 0.01), np.arange(5, 25, 1), np.arange(25, 105, 5)))
    for scan_factor in scan_range:
        adjust_scanned_variable(mmm_vars, Options.instance.var_to_scan, scan_factor)
