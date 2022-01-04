"""Acts as an interface between the various classes that contain data

This module provides the primary function to load variable data from a CDF and
into Variables class objects.  Additionally, it allows for interfacing
between data classes of different types.

Example Usage:
* mmm_vars, cdf_vars, raw_cdf_vars = initialize_variables()
"""

from copy import deepcopy

import modules.variables as variables
import modules.controls as controls
import modules.calculations as calculations
import modules.conversions as conversions
import modules.cdfreader as cdfreader
import modules.utils as utils
from modules.enums import SaveType, ScanType


def initialize_variables(options):
    '''
    Initializes all input variables needed to run the MMM Driver and plot
    variable profiles

    Parameters:
    * options (Options): Contains user specified options

    Returns:
    * mmm_vars (InputVariables): All calculated variables
    * cdf_vars (InputVariables): All interpolated CDF variables
    * raw_cdf_vars (InputVariables): All unedited CDF variables
    '''

    raw_cdf_vars = cdfreader.extract_data(options)
    cdf_vars = conversions.convert_variables(raw_cdf_vars)
    mmm_vars = calculations.calculate_new_variables(cdf_vars)

    return mmm_vars, cdf_vars, raw_cdf_vars


def deepcopy_data(obj):
    '''
    Creates a deepcopy of the given object and reference between their
    options

    A deepcopy is needed to avoid creating a reference between two object
    classes.  However, we do want to create a reference between the options
    stored in each class, so that only one options object exists between all
    copied objects.

    Parameters:
    * obj (InputControls | InputVariables | OutputVariables): The obj to deepcopy

    Returns
    * new_obj (InputControls | InputVariables | OutputVariables): deepcopy of obj
    '''

    new_obj = deepcopy(obj)
    new_obj.options = obj.options

    return new_obj


def get_all_rho_data(options):
    '''
    Creates dictionaries that map rho values to InputVariables and
    OutputVariables objects

    Data is loaded from CSVs stored in the rho folder of the runid, scan_num,
    and var_to_scan, which are supplied via the options parameter. A list of
    rho values for the scan is created from the filenames of the CSVs.

    Parameters:
    * options (Options): Object containing user options

    Returns:
    * input_vars_dict (dict): Dictionary mapping rho values (str) to InputVariables input data
    * output_vars_dict (dict): Dictionary mapping rho values (str) to OutputVariables data
    * input_controls (InputControls or None): InputControls object with np.ndarray for values
    '''

    input_vars_dict, output_vars_dict = {}, {}
    rho_values = utils.get_rho_strings(options, SaveType.OUTPUT)

    # Stores InputVariables and OutputVariables data objects for each rho_value
    for rho in rho_values:
        input_vars = variables.InputVariables(options)
        output_vars = variables.OutputVariables(options)

        input_vars.load_from_csv(SaveType.INPUT, rho_value=rho)
        input_vars.load_from_csv(SaveType.ADDITIONAL, rho_value=rho)
        output_vars.load_from_csv(SaveType.OUTPUT, rho_value=rho)

        input_vars_dict[rho] = input_vars
        output_vars_dict[rho] = output_vars

    # Get control_file from rho folder (there's at most one control file, as controls are independent of rho values)
    input_controls = controls.InputControls(options)
    input_controls.load_from_csv(use_rho=True)

    return input_vars_dict, output_vars_dict, input_controls


def get_data_objects(options, scan_factor=None, rho_value=None):
    '''
    Get InputVariables, OutputVariables, and InputControls data objects

    Parameters:
    * options (Options): Object containing user options
    * scan_factor (float): The scan factor to load (Optional)
    * rho_value (str): The rho value to load (Optional)

    Returns:
    * input_vars (InputVariables): Object containing base input variable data
    * output_vars (OutputVariables): Object containing base output variable data
    * input_controls (InputControls): Object containing base input control data
    '''

    input_vars = variables.InputVariables(options)
    output_vars = variables.OutputVariables(options)
    input_controls = controls.InputControls(options)

    input_vars.load_from_csv(SaveType.INPUT, scan_factor, rho_value)
    input_vars.load_from_csv(SaveType.ADDITIONAL, scan_factor, rho_value)
    output_vars.load_from_csv(SaveType.OUTPUT, scan_factor, rho_value)

    use_rho = True if rho_value is not None else False
    input_controls.load_from_csv(scan_factor=scan_factor, use_rho=use_rho)

    return input_vars, output_vars, input_controls


def get_scan_type(var_to_scan):
    '''
    Gets the scan type from the variable being scanned

    Parameters:
    * var_to_scan (str): The variable being scanned
    * options (Options): Object containing user options

    Raises:
    * TypeError: If var_to_scan is not a member of InputVariables or InputControls
    '''

    scan_type = ScanType.NONE
    if var_to_scan is not None:
        if hasattr(variables.InputVariables(), var_to_scan):
            scan_type = ScanType.VARIABLE
        elif hasattr(controls.InputControls(), var_to_scan):
            scan_type = ScanType.CONTROL
        else:
            raise TypeError(f'Variable {var_to_scan} is not defined under InputVariables or InputControls')

    return scan_type
