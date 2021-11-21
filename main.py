import time # Standard Packages

import numpy as np # 3rd Party Packages
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

from mmm_package import read_cdf, create_inputs, variables, constants # Local Packages

def main(cdfname):
    vars = read_cdf.read_cdf(cdfname)
    vars = create_inputs.create_inputs(vars)
    vars.print_nonzero_variables()

if __name__ == '__main__':
    main('cdf/132017T01.CDF')
