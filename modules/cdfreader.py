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

    if hasattr(cdf, 'Runid') and options.runid != cdf.Runid.strip() and options.runid != 'TEST':
        _log.warning(f'\n\tThe CDF Runid {cdf.Runid.strip()} does not match runid {options.runid}\n')

    cdf_vars = variables.InputVariables(options)
    cdf_vars_to_get = cdf_vars.get_cdf_variables()
    if not options.load_all_vars:
        cdf_vars_to_get = [c for c in cdf_vars_to_get if getattr(cdf_vars, c).required]

    # Find indices of specified time values, so all times don't load
    time_values = np.array(cdf.variables[cdf_vars.time.cdfvar.upper()][:].T)
    cdf_vars.options.set_time_ranges(time_values)
    cdf_vars.options.set_measurement_time(time_values)
    time_idxs = cdf_vars.options.scan_range_idxs or options.time_idx
    if time_idxs is not None and not isinstance(time_idxs, list):
        time_idxs = [time_idxs]

    # Get values for all specified CDF variables
    for var_name in cdf_vars_to_get:
        var = getattr(cdf_vars, var_name)
        cdfvar = var.cdfvar.upper()

        if cdfvar in cdf.variables:
            # Transpose to put the values in the format needed for calculations: (X, T)
            values = np.array(cdf.variables[cdfvar][:].T)

            # Not all variable values in the CDF are arrays
            var.values = values[:] if values.size > 1 else values
            if time_idxs is not None:
                var.values = values[:, time_idxs] if values.ndim == 2 else values[time_idxs]
            if hasattr(cdf.variables[cdfvar], 'units'):
                var.units = (cdf.variables[cdfvar].units).strip()
            if hasattr(cdf.variables[cdfvar], 'long_name'):
                var.desc = (cdf.variables[cdfvar].long_name).strip()

            # Store variable dimensions in reverse order, since we transposed the values above
            cdf_dimensions = cdf.variables[cdfvar].get_dims()
            var.dimensions = [dim.name for dim in cdf_dimensions]
            var.dimensions.reverse()

        elif print_warnings and var.default_values is None:
            _log.error(f'\n\t{var.cdfvar} not found in CDF and no default values were set\n')

    # Update measurement time using loaded time values
    cdf_vars.options.set_measurement_time(cdf_vars.time.values)

    return cdf_vars


def print_variables(runid, cdf):
    '''
    Print names, descriptions, units, and dimensions of all variables in the CDF

    Parameters:
    * cdf (Dataset): The cdf data
    '''

    cdf_cdf_vars = sorted(cdf.variables.keys())
    if 'PH' not in runid:
        for var_name in cdf_cdf_vars:
            var = cdf.variables[var_name]
            var_dims = [dim.name for dim in var.get_dims()]
            print(f'{var.name:<16}{var.long_name.strip():<60}{var.units.strip():<20}{var_dims}')
            # print(f'{var.name:<16}')
    else:
        for var_name in cdf_cdf_vars:
            var = cdf.variables[var_name]
            var_dims = [dim.name for dim in var.get_dims()]
            print(f'{var.name:<16}{var_dims}')

def print_dimensions(cdf):
    '''
    Print names and sizes of dimensions in the CDF

    Parameters:
    * cdf (Dataset): The cdf data
    '''

    cdf_dims = sorted(cdf.dimensions.keys())
    for dim_name in cdf_dims:
        dim = cdf.dimensions[dim_name]
        print(f'{dim.name:<16}{dim.size}')


if __name__ == '__main__':
    # For testing purposes
    runid = 'ps_tb.debug'
    runid = '16325A00'
    cdf = Dataset(utils.get_cdf_path(runid))
    print_variables(runid, cdf)
    print_dimensions(cdf)

    # for v in cdf.variables['omegat']:
    #     print(v)
