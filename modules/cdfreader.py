"""Reads variable data from a CDF produced by TRANSP

This module has been tested to work with how data is formatted for both DIII-D
and NSTX TRANSP discharge types, and may not work for other types of
discharges.  The CDF's must be saved in the "cdfs" folder of the top level
directory in order for this module to find them.

Example Usage:
* cdf_vars = extract_data('120968A02.CDF')
* print_variables('120968A02.CDF')
* print_dimensions('120968A02.CDF')
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import os.path
import logging

# 3rd Party Packages
from netCDF4 import Dataset
import numpy as np

# Local Packages
import modules.variables as variables
import modules.utils as utils


_log = logging.getLogger(__name__)


def extract_data(options, print_warnings=False):
    '''
    Extracts variable data from a CDF and stores it in a variables object

    Parameters:
    * options (Options): Object containing user options
    * print_warnings (bool): Prints warning messages

    Returns:
    * cdf_vars (InputVariables): Object containing extracted variable data from the CDF

    Raises:
    * FileNotFoundError: If the CDF file was not found
    '''

    cdf_file = utils.get_cdf_path(options.runid, options.shot_type)
    if not os.path.exists(cdf_file):
        raise FileNotFoundError(
            f'CDF {options.runid} could not be found in the cdf folder'
            f'\n\tPath: {cdf_file}'
        )

    cdf = Dataset(cdf_file)

    # Runid from CDF should match input runid, else CDF file might be named incorrectly
    if options.runid != cdf.Runid.strip() and options.runid != 'TEST':
        _log.warning(f'\n\tThe CDF Runid {cdf.Runid.strip()} does not match runid {options.runid}\n')

    cdf_vars = variables.InputVariables(options)
    cdf_vars_to_get = cdf_vars.get_cdf_variables()

    # Get values for all specified CDF variables
    for var_name in cdf_vars_to_get:
        var = getattr(cdf_vars, var_name)
        cdfvar = var.cdfvar.upper()

        if cdfvar in cdf.variables:
            # Transpose to put the values in the format needed for calculations: (X, T)
            values = np.array(cdf.variables[cdfvar][:].T)

            # Not all variable values in the CDF are arrays
            var.values = values[:] if values.size > 1 else values
            var.units = (cdf.variables[cdfvar].units).strip()
            var.desc = (cdf.variables[cdfvar].long_name).strip()

            # Store variable dimensions in reverse order, since we transposed the values above
            cdf_dimensions = cdf.variables[cdfvar].get_dims()
            var.dimensions = [dim.name for dim in cdf_dimensions]
            var.dimensions.reverse()

        elif print_warnings and var.default_values is None:
            _log.error(f'\n\t{var.cdfvar} not found in CDF and no default values were set\n')

    cdf_vars.options.set_measurement_time(cdf_vars.time.values)

    return cdf_vars


def print_variables(runid):
    '''
    Print names, descriptions, units, and dimensions of all variables in the CDF

    Parameters:
    * runid (str): The file name of the CDF (without the path)
    '''

    cdf = Dataset(utils.get_cdf_path(runid))
    cdf_cdf_vars = sorted(cdf.variables.keys())

    for var_name in cdf_cdf_vars:
        var = cdf.variables[var_name]
        var_dims = [dim.name for dim in var.get_dims()]
        print(f'{var.name}, {var.long_name.strip()}, {var.units.strip()}, {str(var_dims)}')


def print_dimensions(runid):
    '''
    Print names and sizes of dimensions in the CDF

    Parameters:
    * runid (str): The file name of the CDF (without the path)
    '''

    cdf = Dataset(utils.get_cdf_path(runid))
    cdf_dims = sorted(cdf.dimensions.keys())

    for dim_name in cdf_dims:
        dim = cdf.dimensions[dim_name]
        print(f'{dim.name}, {dim.size}')


if __name__ == '__main__':
    # For testing purposes
    import modules.options
    utils.init_logging()
    options = modules.options.Options(runid='129016Z41')
    cdf_cdf_vars = extract_data(options, print_warnings=True)
    print_dimensions(options.runid)
    print_variables(options.runid)
