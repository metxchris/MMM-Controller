"""Reads variable data from a CDF produced by TRANSP

This module has been tested to work with how data is formatted for both DIII-D
and NSTX TRANSP discharge types, and may not work for other types of
discharges.  The CDF's must be saved in the "cdfs" folder of the top level
directory in order for this module to find them.

Example Usage:
* cdf_vars = read_cdf('120968A02.CDF')
* print_cdf_variables('120968A02.CDF')
* print_cdf_dimensions('120968A02.CDF')
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import os.path

# 3rd Party Packages
from netCDF4 import Dataset
import numpy as np

# Local Packages
import modules.options as options
import modules.variables as variables
import modules.utils as utils


def read_cdf(print_warnings=False):
    '''
    Reads variable data from a CDF and stores it in a variables object

    Parameters:
    * print_warnings (bool): Prints warning messages

    Returns:
    * cdf_vars (InputVariables): Object containing extracted variable data from the CDF

    Raises:
    * FileNotFoundError: If the CDF file was not found
    '''

    runid = options.instance.runid

    cdf_file = utils.get_cdf_path(runid)
    if not os.path.exists(cdf_file):
        raise FileNotFoundError(f'CDF {runid} could not be found in the cdf folder')

    cdf = Dataset(cdf_file)

    # Runid from CDF should match input runid, else CDF file might be named incorrectly
    if runid != cdf.Runid.strip() and runid != 'TEST':
        print(f'Warning: The CDF Runid {cdf.Runid.strip()} does not match runid {runid}')

    cdf_vars = variables.InputVariables()
    cdf_vars_to_get = cdf_vars.get_cdf_variables()

    # Get values for all specified CDF variables
    for var_name in cdf_vars_to_get:
        if getattr(cdf_vars, var_name).cdfvar in cdf.variables:
            # Transpose to put the values in the format needed for calculations: (X, T)
            values = np.array(cdf.variables[getattr(cdf_vars, var_name).cdfvar][:].T)

            # Not all variable values in the CDF are arrays
            getattr(cdf_vars, var_name).values = values[:] if values.size > 1 else values
            getattr(cdf_vars, var_name).units = (cdf.variables[getattr(cdf_vars, var_name).cdfvar].units).strip()
            getattr(cdf_vars, var_name).desc = (cdf.variables[getattr(cdf_vars, var_name).cdfvar].long_name).strip()

            # Store variable dimensions in reverse order, since we transposed the values above
            cdf_dimensions = cdf.variables[getattr(cdf_vars, var_name).cdfvar].get_dims()
            getattr(cdf_vars, var_name).dimensions = [dim.name for dim in cdf_dimensions]
            getattr(cdf_vars, var_name).dimensions.reverse()

        elif print_warnings:
            # Not all variables will be found in the CDF, which is expected
            print(f'*** [read_cdf] WARNING: {getattr(cdf_vars, var_name).cdfvar} not found in CDF')

    return cdf_vars


def print_cdf_variables(cdf_name):
    '''
    Print names, descriptions, units, and dimensions of all variables in the CDF

    Parameters:
    * cdf_name (str): The file name of the CDF (without the path)
    '''

    cdf = Dataset(utils.get_cdf_path(cdf_name))
    cdf_cdf_vars = sorted(cdf.variables.keys())

    for var_name in cdf_cdf_vars:
        var = cdf.variables[var_name]
        var_dims = [dim.name for dim in var.get_dims()]
        print(f'{var.name}, {var.long_name.strip()}, {var.units.strip()}, {str(var_dims)}')


def print_cdf_dimensions(cdf_name):
    '''
    Print names and sizes of dimensions in the CDF

    Parameters:
    * cdf_name (str): The file name of the CDF (without the path)
    '''

    cdf = Dataset(utils.get_cdf_path(cdf_name))
    cdf_dims = sorted(cdf.dimensions.keys())

    for dim_name in cdf_dims:
        dim = cdf.dimensions[dim_name]
        print(f'{dim.name}, {dim.size}')


if __name__ == '__main__':
    # For testing purposes
    opts = options.instance
    opts.runid = '132017T01'
    cdf_cdf_vars = read_cdf(True)
    print_cdf_dimensions(opts.runid)
    print_cdf_variables(opts.runid)
