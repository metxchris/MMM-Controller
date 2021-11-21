from os.path import exists # Standard Packages
import sys
sys.path.insert(0, '../')

from netCDF4 import Dataset # 3rd Party Packages

from mmm_package import variables # Local Packages

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
            # Transpose to put the values in the format needed for calculations: (X, Time)
            values = cdf.variables[getattr(vars, var).cdfvar][:].T

            # Not all values in the CDF are arrays
            getattr(vars, var).values = values if values.ndim == 0 else values[:]

            # Store units of values and strip extra whitespace
            getattr(vars, var).units = (cdf.variables[getattr(vars, var).cdfvar].units).strip()

            # Store long name of values and strip extra whitespace
            getattr(vars, var).desc = (cdf.variables[getattr(vars, var).cdfvar].long_name).strip()

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