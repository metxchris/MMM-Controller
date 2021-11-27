# Standard Packages
from os.path import exists, dirname
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
from netCDF4 import Dataset
import numpy as np
# Local Packages
from main import variables
import cdfs

# Returns the path to the CDF folder
def get_cdf_path(cdf_name):
    return "{0}/{1}.CDF".format(dirname(cdfs.__file__), cdf_name)

# Reads CDF variables specified by Variables().cdfname and a Variables() object
def read_cdf(input_options, print_warnings=False):
    cdf_file = get_cdf_path(input_options.cdf_name)

    # Check if file exists
    if not exists(cdf_file):
        raise FileNotFoundError("CDF {0} could not be found in the cdf folder".format(input_options.cdf_name))

    # Load CDF into memory
    cdf = Dataset(cdf_file)

    # Set runid from CDF (should match cdf_name)
    input_options.runid = cdf.Runid

    # Variables object to store CDF values
    vars = variables.Variables()

    # List all variables that have a specfied CDF variable name in the Variables class
    vars_to_get = vars.get_cdf_variables()

    # Get values for all specified CDF variables
    for var in vars_to_get:
        if getattr(vars, var).cdfvar in cdf.variables:
            # Transpose to put the values in the format needed for calculations: (X, T)
            values = np.array(cdf.variables[getattr(vars, var).cdfvar][:].T)

            # Not all values in the CDF are arrays
            getattr(vars, var).values = values[:]

            # Store units of values and strip extra whitespace
            getattr(vars, var).units = (cdf.variables[getattr(vars, var).cdfvar].units).strip()

            # Store long name of values and strip extra whitespace
            getattr(vars, var).desc = (cdf.variables[getattr(vars, var).cdfvar].long_name).strip()

            # Store variable dimensions in reverse order, since we transposed the values above
            cdf_dimensions = cdf.variables[getattr(vars, var).cdfvar].get_dims()
            getattr(vars, var).dimensions = [dim.name for dim in cdf_dimensions]
            getattr(vars, var).dimensions.reverse()

        elif print_warnings:
            print('*** [read_cdf] WARNING: {0} not found in CDF'.format(getattr(vars, var).cdfvar))

    if len(vars.get_nonzero_variables()) == 0:
        print('*** [read_cdf] ERROR: no variables were saved from ' + cdf_name)

    return vars

# Print all variable names, descriptions, units, and dimensions in the CDF
def print_cdf_variables(cdf_name):
    cdf = Dataset(get_cdf_path(cdf_name))
    cdf_vars = sorted(cdf.variables.keys())

    for var_name in cdf_vars:
        var = cdf.variables[var_name]
        var_dims = [dim.name for dim in var.get_dims()]
        print("{0}, {1}, {2}, {3}".format(var.name, var.long_name.strip(), var.units.strip(), str(var_dims)))

# Print all dimension names and sizes in the CDF
def print_cdf_dimensions(cdf_name):
    cdf = Dataset(get_cdf_path(cdf_name))
    cdf_dims = sorted(cdf.dimensions.keys())

    for dim_name in cdf_dims:
        dim = cdf.dimensions[dim_name]
        print("{0}, {1}".format(dim.name, str(dim.size)))

if __name__ == '__main__':
    # For testing purposes
    cdf_name = '132017T01'
    read_cdf(variables.InputOptions(cdf_name), True)
    print_cdf_dimensions(cdf_name)
    print_cdf_variables(cdf_name)
