"""Adjusts a variable during a variable scan.

Variable Classifications:
* Base variable: A non-gradient variable that is either used as input to MMM,
  or is needed to calculate an MMM input variable.  Base variables do not
  depend on values of gradient or additional variables.
* Gradient variable: A variable that is a normalized gradient, and may be used
  as an input to MMM.  These variables depend on values of base variables,
  but not on values of additional variables.
* Additional variable: A variable that is not needed for MMM input nor input
  calculations. Additional variables can depend on both base and gradient
  values, so they are always recalculated when an adjustment to any variable
  type is made.

Adjustment Methodology:
* Our goal is to keep adjustments to variables as physical as possible, within
  reason.  As such, we use the following variable recalculation rules when
  adjusting different variable types:
  - Adjusting base variables: If the adjustment is linear with respect to rho,
    then normalized gradient values will remain constant, by definition.  In
    this case, only base variables and additional variables are recalculated,
    which saves computational time from also having to recalculate gradients.
    Alternatively, if the adjustment is nonlinear with respect to rho(such as
    for zeff), then all variables are recalculated, since gradient values
    will change.
  - Adjusting gradient variables: We only recalculate additional variables
    following a gradient adjustment, as recalculating the gradients would
    overwrite the gradient adjustment.  Consequently, our adjustment process
    for handling gradients is somewhat nonphysical, but the overall nature of
    scan results will remain the same.  Note that base variables do not
    depend on gradient values, so their recalculations are not needed.
  - Adjusting additional variables: Additional variables are not used as input
    for MMM, so they are always adjusted indirectly by adjusting the base or
    gradient variables that determine the value of the additional variable.

* Our recalculation process aims to find a balance between computational
  performance and minimizing human error.  For example, if we were to adjust
  the value of the electron temperature (te), then all base variables would
  still be recalculated, even though there are no base variables that depend
  on the value of te.  This approach eliminates the need for the user to
  correctly determine the dependencies of each variable being adjusted.
  However, we would not to recalculate gradient values if the adjustment to
  te is linear, as the value of the electron temperature gradient will remain
  constant, by definition of normalized gradients.
* When adjusting a derived quantity, adjustments are made directly to the
  variables that determine that derived quantity.  When multiple variables
  need to be adjusted in this manner, each variable is adjusted
  (either multiplied or divided) by the same base amount, which is the
  adjustment factor.  For example, if we wanted to double the value of the
  temperature ratio (te_ti = te / ti), the scan factor for te_ti would be 2, and
  the adjustment factor for te and ti would be sqrt(2). By multiplying te and
  dividing ti by this adjustment factor, we indirectly multiply te_ti by the
  value of the scan factor.
* If any variables have powers attached to them, then the adjustment factor is
  updated as needed, but all variables are still adjusted by the same
  adjustment factor.  For example, if we wanted to double the value of betae,
  where betae is proportional to (ne * te / bzxr**2), then our scan factor
  would be 2, and our adjustment factor would be 2**(1 / 4).  In this case,
  by multiplying both ne and te and dividing bzxr by the adjustment factor, we
  indirectly multiply betae by the value of the scan factor.

Example Usage:
* new_vars = adjust_scanned_variable(original_vars, 'te_ti', 2.5)
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
from math import log10

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.calculations as calculations
import modules.datahelper as datahelper
from modules.enums import SaveType


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
    * values (np.ndarray): A one-dimensional array of variable data

    Returns:
    * nonzero_values[0] (int): The index along the radial dimension of first nonzero variable value

    Raises:
    * ValueError: If there are no nonzero values
    '''

    nonzero_values = np.where(values != 0)[0]
    if not len(nonzero_values):
        raise ValueError('Cannot adjust variable that is equal to 0 at all radial points')

    return nonzero_values[0]


def _check_adjusted_factor(scan_factor, base_var, adjusted_var, t):
    '''
    Checks that the adjusted factor equals scan_factor, within an allowable tolerance

    Parameters:
    * scan_factor (float): The intended factor that a variable was adjusted by
    * base_var (Variable): The unmodified variable object
    * adjusted_var (Variable): The adjusted variable object
    * t (int): The time index

    Raises:
    * ValueError: If the adjusted factor does not equal the scan factor within allowable tolerance
    '''

    r = _get_nonzero_idx(base_var.values[:, t])
    adjusted_factor = adjusted_var.values[r, t] / base_var.values[r, t]
    if abs(adjusted_factor / scan_factor - 1) > SCAN_FACTOR_TOLERANCE:
        fmt_length = abs(int(log10(SCAN_FACTOR_TOLERANCE)))
        raise ValueError(
            f'Scanned variable did not change within the allowable tolerance level\n'
            f'    adjusted_factor: {adjusted_factor}\n'
            f'    scan_factor:     {scan_factor}\n'
            f'    tolerance:       {SCAN_FACTOR_TOLERANCE:.{fmt_length}f}\n'
        )

    if __name__ == '__main__':  # For testing purposes
        _print_factors(scan_factor, adjusted_factor)


def _check_equality(base_var, adjusted_var, t):
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

    r = _get_nonzero_idx(base_var.values)
    variable_error = np.absolute(adjusted_var.values[r, t] / base_var.values[r, t] - 1)
    if variable_error > EQUALITY_ERROR_TOLERANCE:
        fmt_length = abs(int(log10(EQUALITY_ERROR_TOLERANCE)))
        raise ValueError(
            f'{base_var.name} did not remain constant\n'
            f'    Variable Error: {variable_error}\n'
            f'    Allowed Error:  {EQUALITY_ERROR_TOLERANCE:.{fmt_length}f}'
        )


def _adjust_ah(mmm_vars, scan_factor):
    '''
    Adjusts the hydrogenic atomic mass

    The ion atomic mass is also updated.  Full base variables are not
    recalculated because the ah adjustment would be overwritten.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify ah by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    print((adjusted_vars.ah.values[:, :].min()))
    adjusted_vars.ah.values *= scan_factor

    calculations.ai(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    # Check that ah increased by the scan_factor
    _check_adjusted_factor(scan_factor, mmm_vars.ah, adjusted_vars.ah, t)

    return adjusted_vars

def _adjust_ai(mmm_vars, scan_factor):
    '''
    Adjusts the hydrogenic atomic mass

    The ion atomic mass is also updated.  Full base variables are not
    recalculated because the ah adjustment would be overwritten.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify ah by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjusted_vars.ai.values *= scan_factor

    calculations.calculate_additional_variables(adjusted_vars)

    # Check that ah increased by the scan_factor
    _check_adjusted_factor(scan_factor, mmm_vars.ai, adjusted_vars.ai, t)

    return adjusted_vars

def _adjust_ne(mmm_vars, scan_factor):
    '''
    Adjusts the electron density

    All densities dependent on electron density are also adjusted so as to
    remain physically consistent

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify te_ti by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjusted_vars.ne.values *= scan_factor
    adjusted_vars.nz.values *= scan_factor
    adjusted_vars.nd.values *= scan_factor
    adjusted_vars.nf.values *= scan_factor
    adjusted_vars.nh0.values *= scan_factor

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    # Check that nh increased by the scan_factor
    _check_adjusted_factor(scan_factor, mmm_vars.nh, adjusted_vars.nh, t)

    return adjusted_vars


def _adjust_nuei_alphaconst(mmm_vars, scan_factor):
    '''
    Adjusts the Collision Frequency

    nuei is adjusted indirectly by adjusting the values of ne and te. The
    final value of the adjustment total is obtained by running a loop until
    nuei is adjusted by the target scan factor. As part of this scan, we also
    require that alphamhd remain constant as nuei is adjusted, so additional
    adjustments are made to ti and density variables afterwards.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file

    Raises:
    * ValueError: If the total adjustment could not be found in the specified max attempts
    '''

    max_adjustment_attempts = 10
    t = mmm_vars.options.time_idx
    r = _get_nonzero_idx(mmm_vars.nuei.values[:, t])

    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjustment_total = scan_factor**(-2 / 5)  # based on the formula for nuei
    adjustment_step = adjustment_total

    # Make adjustments to ne and te until the adjustment_total is found The
    # adjustment finding loop typically finds the target within 6 attempts
    # for a tolerance of 1e-4

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

        # Make an adjustment_step based on the difference in the current_factor and scan_factor
        adjustment_step = 1 + (current_factor - scan_factor) / (2 * scan_factor)
        adjustment_total *= adjustment_step

    # Additional adjustments to keep alphamhd constant
    adjusted_vars.ti.values *= adjustment_total
    adjusted_vars.nz.values /= adjustment_total
    adjusted_vars.nd.values /= adjustment_total
    adjusted_vars.nf.values /= adjustment_total
    adjusted_vars.nh0.values /= adjustment_total

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.nuei, adjusted_vars.nuei, t)
    _check_equality(mmm_vars.alphamhd, adjusted_vars.alphamhd, t)

    return adjusted_vars


def _adjust_nuei_lareunitconst(mmm_vars, scan_factor):
    '''
    Adjusts the Collision Frequency

    nuei is adjusted indirectly by adjusting the value te. The final value of
    the adjustment total is obtained by running a loop until nuei is adjusted
    by the target scan factor. As part of this scan, the Electron Gyroradius
    (unit) will be held constant by ensuring that te / bu**2 is constant.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify nuei by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file

    Raises:
    * ValueError: If the total adjustment could not be found in the specified max attempts
    '''

    max_adjustment_attempts = 10
    t = mmm_vars.options.time_idx
    r = _get_nonzero_idx(mmm_vars.nuei.values[:, t])

    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjustment_total = scan_factor**(-2 / 3)  # based on the formula for nuei
    adjustment_step = adjustment_total

    # Make adjustments to ne and te until the adjustment_total is found The
    # adjustment finding loop typically finds the target within 6 attempts
    # for a tolerance of 1e-4

    for i in range(max_adjustment_attempts):

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

        # Make an adjustment_step based on the difference in the current_factor and scan_factor
        adjustment_step = 1 + (current_factor - scan_factor) / (2 * scan_factor)
        adjustment_total *= adjustment_step

    # Additional adjustments to electron gyroradius (unit) constant
    adjusted_vars.ti.values *= adjustment_total
    adjusted_vars.bftor.values *= adjustment_total**(1 / 2)

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.nuei, adjusted_vars.nuei, t)
    _check_equality(mmm_vars.alphamhdu, adjusted_vars.alphamhdu, t)
    _check_equality(mmm_vars.lareu, adjusted_vars.lareu, t)

    return adjusted_vars


def _adjust_tau(mmm_vars, scan_factor):
    '''
    Adjusts the temperature ratio

    te_ti is adjusted indirectly by adjusting te and ti

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify te_ti by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjustment_total = scan_factor**(1 / 2)  # based on the formula for te_ti
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.ti.values /= adjustment_total

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.te_ti, adjusted_vars.te_ti, t)

    return adjusted_vars


def _adjust_zeff(mmm_vars, scan_factor):
    '''
    Adjusts the Effective Charge

    zeff is adjusted indirectly by adjusting the values of nz. Since ne
    depends on nz, ne is updated by the change in nz values(ne += z*delta
    (nz)).  However, the values of nh and nd are intentionally held constant,
    and are not updated from the adjustment of ne. Unlike other adjustments,
    the amount that zeff changes is different at each point of rho, meaning
    the adjustment is nonlinear. Therefore, the adjusted factor is checked
    against the change in nz after all recalculations, instead of against the
    change in zeff.  This means that zeff *DOES NOT* change by the value of
    the scan factor, which is a departure for how adjustments generally
    work.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify nz by (not zeff)

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx

    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjusted_vars.nz.values *= scan_factor
    adjusted_vars.ne.values += adjusted_vars.zz.values * (adjusted_vars.nz.values - mmm_vars.nz.values)

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_gradient_variables(adjusted_vars)  # needed as ne change was nonlinear
    calculations.calculate_additional_variables(adjusted_vars)

    if __name__ == '__main__':  # For testing purposes
        r = _get_nonzero_idx(mmm_vars.zeff.values)
        adjusted_factor = adjusted_vars.zeff.values[r, t] / mmm_vars.zeff.values[r, t]
        print(f'{adjusted_factor:{PRINT_FMT_STR}}', end=' ')

    _check_adjusted_factor(scan_factor, mmm_vars.nz, adjusted_vars.nz, t)

    return adjusted_vars


def _adjust_etae(mmm_vars, scan_factor):
    '''
    Adjust etae = gte/gne

    etae is adjusted indirectly by adjusting gte and gne

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify etae by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)

    # adjustment_total = scan_factor**(1 / 2)  # based on the formula for etae
    # adjusted_vars.gte.values *= adjustment_total
    # adjusted_vars.gne.values /= adjustment_total

    # adjusted_vars.gte.values *= scan_factor**(3 / 2)
    # adjusted_vars.gne.values *= scan_factor**(1 / 2)

    # adjusted_vars.gte.values /= scan_factor**(1 / 2)
    # adjusted_vars.gne.values /= scan_factor**(3 / 2)

    adjusted_vars.gte.values *= scan_factor
    adjusted_vars.gne.values *= scan_factor

    calculations.calculate_additional_variables(adjusted_vars)

    # _check_adjusted_factor(scan_factor, mmm_vars.etae, adjusted_vars.etae, t)

    return adjusted_vars


def _adjust_shear(mmm_vars, scan_factor):
    '''
    Adjust Magnetic Shear

    shear is adjusted indirectly by adjusting gq

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify shear by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjusted_vars.gq.values *= scan_factor

    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.shear, adjusted_vars.shear, t)

    return adjusted_vars


def _adjust_btor(mmm_vars, scan_factor):
    '''
    Adjust Toroidal Magnetic Field

    btor is adjusted indirectly by adjusting bftor.  This will not adjust
    bu, since the adjustment to btor cancels out in the bu equation.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify btor by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjusted_vars.bzxr.values *= scan_factor

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.btor, adjusted_vars.btor, t)

    return adjusted_vars


def _adjust_bunit(mmm_vars, scan_factor):
    '''
    Adjust Toroidal Magnetic Field (unit)

    bu is adjusted indirectly by adjusting bftor

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify btor by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjusted_vars.bftor.values *= scan_factor

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.bu, adjusted_vars.bu, t)

    return adjusted_vars


def _adjust_betae(mmm_vars, scan_factor):
    '''
    Adjust Electron Pressure Ratio

    betae is adjusted indirectly by adjusting ne, te, and bzxr. The values of
    nd, nz, and nf depend on changes in ne, so these variables are updated
    accordingly.  In addition, alphamhd is held constant.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify betae by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjustment_total = scan_factor**(1 / 4)  # based on the formula for betae
    adjusted_vars.ne.values *= adjustment_total
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.bzxr.values /= adjustment_total
    adjusted_vars.nd.values *= adjustment_total
    adjusted_vars.nz.values *= adjustment_total
    adjusted_vars.nf.values *= adjustment_total

    # NEW: Keep alphamhd constant
    adjusted_vars.ti.values *= adjustment_total
    adjusted_vars.gte.values /= scan_factor
    adjusted_vars.gti.values /= scan_factor
    adjusted_vars.gne.values /= scan_factor
    adjusted_vars.gni.values /= scan_factor

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.betae, adjusted_vars.betae, t)
    _check_equality(mmm_vars.alphamhd, adjusted_vars.alphamhd, t)

    return adjusted_vars


def _adjust_betaeu(mmm_vars, scan_factor):
    '''
    Adjust Electron Pressure Ratio

    betae is adjusted indirectly by adjusting ne, te, and bzxr. The values of
    nd, nz, and nf depend on changes in ne, so these variables are updated
    accordingly.

    Parameters:
    * mmm_vars (InputVariables): Contains unmodified variables
    * scan_factor (float): The factor to modify betae by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file
    '''

    t = mmm_vars.options.time_idx
    adjusted_vars = datahelper.deepcopy_data(mmm_vars)
    adjustment_total = scan_factor**(1 / 4)  # based on the formula for betae
    adjusted_vars.ne.values *= adjustment_total
    adjusted_vars.te.values *= adjustment_total
    adjusted_vars.bftor.values /= adjustment_total
    adjusted_vars.nd.values *= adjustment_total
    adjusted_vars.nz.values *= adjustment_total
    adjusted_vars.nf.values *= adjustment_total

    calculations.calculate_base_variables(adjusted_vars)
    calculations.calculate_additional_variables(adjusted_vars)

    _check_adjusted_factor(scan_factor, mmm_vars.betaeu, adjusted_vars.betaeu, t)

    return adjusted_vars


def adjust_scanned_variable(mmm_vars, scan_factor):
    '''
    Adjusts the variable being scanned, as well as any necessary dependencies

    If the adjustment requires any advanced logic, then the corresponding
    private function is called.  These functions handle all adjustments,
    recalculations, and error checks.  Otherwise, simple adjustments are
    handled in a generalized manner, with no error checking needed.  In
    either case, adjustments are made using a deepcopy of the base variable,
    so that adjustments do not unintentionally alter base variable values.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * scan_factor (float): The factor to modify var_to_scan by

    Returns:
    * adjusted_vars (InputVariables): Adjusted variables needed to write MMM input file

    Raises:
    * ValueError: If scanned variable is not of SaveType.Input and does not
      have an advanced scan defined
    '''

    adjustment_name = mmm_vars.options.adjustment_name

    # adjustment names that don't correspond to variable names need to be
    # specified in _adjustment_name_to_var_dict within the options module
    if adjustment_name == 'nuei_alphaconst':
        adjusted_vars = _adjust_nuei_alphaconst(mmm_vars, scan_factor)

    elif adjustment_name == 'nuei_lareunitconst':
        adjusted_vars = _adjust_nuei_lareunitconst(mmm_vars, scan_factor)

    elif adjustment_name == 'ah':
        adjusted_vars = _adjust_ah(mmm_vars, scan_factor)

    elif adjustment_name == 'ai':
        adjusted_vars = _adjust_ai(mmm_vars, scan_factor)

    elif adjustment_name == 'ne':
        adjusted_vars = _adjust_ne(mmm_vars, scan_factor)

    elif adjustment_name == 'zeff':
        adjusted_vars = _adjust_zeff(mmm_vars, scan_factor)

    elif adjustment_name == 'te_ti':
        adjusted_vars = _adjust_tau(mmm_vars, scan_factor)

    elif adjustment_name == 'etae':
        adjusted_vars = _adjust_etae(mmm_vars, scan_factor)

    elif adjustment_name == 'shear' or adjustment_name == 'gq' or adjustment_name == 'shat_gxi':
        adjusted_vars = _adjust_shear(mmm_vars, scan_factor)

    elif adjustment_name == 'btor' or adjustment_name == 'bzxr':
        adjusted_vars = _adjust_btor(mmm_vars, scan_factor)

    elif adjustment_name == 'bu' or adjustment_name == 'bftor':
        adjusted_vars = _adjust_bunit(mmm_vars, scan_factor)

    elif adjustment_name == 'betae':
        adjusted_vars = _adjust_betae(mmm_vars, scan_factor)

    elif adjustment_name == 'betaeu':
        adjusted_vars = _adjust_betaeu(mmm_vars, scan_factor)

    else:

        # Simple Scan (no advanced logic needed)
        var_to_scan = mmm_vars.options.var_to_scan
        adjusted_vars = datahelper.deepcopy_data(mmm_vars)
        base_var = getattr(mmm_vars, var_to_scan)
        scanned_var = getattr(adjusted_vars, var_to_scan)
        scanned_var.values = scan_factor * base_var.values

        if var_to_scan in ['gte', 'gti', 'gne', 'gnh', 'gni', 'gnz', 'gvpar', 'gvpol', 'gvtor']:
            calculations.calculate_additional_variables(adjusted_vars)

        else:
            calculations.calculate_base_variables(adjusted_vars)
            calculations.calculate_additional_variables(adjusted_vars)

    return adjusted_vars


if __name__ == '__main__':  # For Testing Purposes
    from modules.options import Options

    options = Options()
    options.set(
        runid='138536A01',
        input_points=51,
        apply_smoothing=True,
        uniform_rho=False,
        input_time=.63,
    )

    mmm_vars = datahelper.initialize_variables(options, 1, 1)

    # Check that all scan_factors can be found in scan_range (failures will raise a ValueError)
    scan_range = np.hstack(
        (np.arange(0.05, 1, 0.05),
         np.arange(1, 5, 0.2),
         np.arange(5, 20, 1),
         np.arange(20, 105, 5)),
    )

    advanced_scans = [
        'ne', 'nuei_alphaconst', 'nuei_lareunitconst', 'zeff',
        'te_ti', 'etae', 'shear', 'btor',
        'bu', 'betae', 'betaeu',
    ]

    for var_name in advanced_scans:
        mmm_vars.options.set(adjustment_name=var_name)
        print(f'\n{var_name} scan, {mmm_vars.options.var_to_scan} factors:')
        for scan_factor in scan_range:
            adjust_scanned_variable(mmm_vars, scan_factor)
