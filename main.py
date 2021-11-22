# Standard Packages
import time
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
# Local Packages
from mmm_package import read_cdf, create_inputs, variables, constants
from plots import plot2d

def main(cdfname):
    cdf_vars = read_cdf.read_cdf(cdfname)
    input_vars = create_inputs.create_inputs(cdf_vars)
    input_vars.print_nonzero_variables()
    # plot2d.plot2d(input_vars.xb.values, input_vars.ti.values, input_vars.xb.values, input_vars.te.values)

if __name__ == '__main__':
    main('cdf/132017T01.CDF')
