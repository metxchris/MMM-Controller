import time # Standard Packages

import numpy as np # 3rd Party Packages

from mmm_package import read_cdf, calc_inputs, variables, constants # Local Packages

def main(cdfname):
    
    vars = read_cdf.read_cdf(cdfname)
    # vars.print_variable_descriptions()

    calc_inputs.calc_inputs(vars)



if __name__ == '__main__':
    main('cdf/132017T01.CDF')
