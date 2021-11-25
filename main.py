# Standard Packages
import time
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
# Local Packages
from mmm_package import read_cdf, convert_inputs, calculate_inputs, variables, constants
from plots import plot2d

def main(cdfname, t):
    # Read variables from specified CDF
    cdf_vars = read_cdf.read_cdf(cdfname)

    # Convert variables from CDF format to MMM format
    input_vars = convert_inputs.convert_inputs(cdf_vars)

    # Calculate new variables from CDF variables
    calculate_inputs.calculate_inputs(input_vars)

    input_vars.print_nonzero_variables()

    plot2d.plot_input_profiles(input_vars, t)
    # plot2d.plot2d(input_vars.xb.values, input_vars.shat.values)

if __name__ == '__main__':
    t = 2.1
    main('cdf/132017T01.CDF', t)
