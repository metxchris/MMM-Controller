# Standard Packages
import sys; sys.path.insert(0, '../')
from os.path import exists

# 3rd Party Packages
from netCDF4 import Dataset
import numpy as np

# Local Packages
import modules.options as options
import modules.variables as variables
import modules.utils as utils


# Reads CDF variables specified by Variables().cdfname and a Variables() object
def read_cdf(print_warnings=False):
    opts = options.instance

    cdf_file = utils.get_cdf_path(opts.runid)

    # Check if file exists
    if not exists(cdf_file):
        raise FileNotFoundError(f'CDF {opts.runid} could not be found in the cdf folder')

    # Load CDF into memory
    cdf = Dataset(cdf_file)

    # Runid from CDF should match input runid
    if cdf.Runid.strip() != opts.runid and opts.runid != 'TEST':
        # TODO: Save all warnings strings and print at the end of code execution
        print(f'Warning: cdf.Runid {cdf.Runid.strip()} does not match opts.runid {opts.runid}')

    # Variables object to store CDF values
    cdf_vars = variables.InputVariables()

    # List all variables that have a specified CDF variable name in the Variables class
    cdf_vars_to_get = cdf_vars.get_cdf_variables()

    # Get values for all specified CDF variables
    for var_name in cdf_vars_to_get:
        if getattr(cdf_vars, var_name).cdfvar in cdf.variables:
            # Transpose to put the values in the format needed for calculations: (X, T)
            values = np.array(cdf.variables[getattr(cdf_vars, var_name).cdfvar][:].T)

            # Not all values in the CDF are arrays
            getattr(cdf_vars, var_name).values = values[:] if values.size > 1 else values

            # Store units of values and strip extra white space
            getattr(cdf_vars, var_name).units = (cdf.variables[getattr(cdf_vars, var_name).cdfvar].units).strip()

            # Store long name of values and strip extra white space
            getattr(cdf_vars, var_name).desc = (cdf.variables[getattr(cdf_vars, var_name).cdfvar].long_name).strip()

            # Store variable dimensions in reverse order, since we transposed the values above
            cdf_dimensions = cdf.variables[getattr(cdf_vars, var_name).cdfvar].get_dims()
            getattr(cdf_vars, var_name).dimensions = [dim.name for dim in cdf_dimensions]
            getattr(cdf_vars, var_name).dimensions.reverse()

        elif print_warnings:
            print(f'*** [read_cdf] WARNING: {getattr(cdf_vars, var_name).cdfvar} not found in CDF')

    if len(cdf_vars.get_nonzero_variables()) == 0:
        print('*** [read_cdf] ERROR: no variables were saved from ' + cdf_name)

    return cdf_vars

# Print all variable names, descriptions, units, and dimensions in the CDF
def print_cdf_variables(cdf_name):
    cdf = Dataset(utils.get_cdf_path(cdf_name))
    cdf_cdf_vars = sorted(cdf.variables.keys())

    for var_name in cdf_cdf_vars:
        var = cdf.variables[var_name]
        var_dims = [dim.name for dim in var.get_dims()]
        print(f'{var.name}, {var.long_name.strip()}, {var.units.strip()}, {str(var_dims)}')

# Print all dimension names and sizes in the CDF
def print_cdf_dimensions(cdf_name):
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
