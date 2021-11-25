# Standard Packages
from os.path import exists
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
from netCDF4 import Dataset
import numpy as np
# Local Packages
from mmm_package import variables

def read_cdf(cdfname, print_warnings=False):
    # Check if file exists
    if not exists(cdfname):
        raise FileNotFoundError("CDF " + cdfname + " could not be found in the cdf folder")

    # Load CDF into memory
    cdf = Dataset(cdfname)

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
            print('[read_cdf] *** WARNING:', getattr(vars, var).cdfvar, 'not found in CDF')

    if len(vars.get_nonzero_variables()) == 0:
        print('[read_cdf] *** ERROR: no variables were saved from ' + cdfname)

    return vars

if __name__ == '__main__':
    # For testing purposes
    read_cdf('../cdf/132017T01.CDF', True)

# Note: All variables in CDF can be viewed using
# for dimobj in cdf.variables.values():
    #     print(dimobj)
