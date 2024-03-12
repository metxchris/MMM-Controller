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


def write_cdf(mmm_vars):

    mmm_vars.options.set_time_ranges(mmm_vars.time.values)
    print(mmm_vars.options.scan_range_idxs)
    print(mmm_vars.time.values[mmm_vars.options.scan_range_idxs])
    # print(mmm_vars.te.values[:, mmm_vars.options.scan_range_idxs])

    try: ncfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    ncfile = Dataset('cdfs/new.CDF',mode='w',format='NETCDF4_CLASSIC') 
    all_vars = mmm_vars.get_variables()

    xb_dim = ncfile.createDimension('xb', getattr(mmm_vars, 'xb').values.shape[0])    # longitude axis
    time_dim = ncfile.createDimension('time', len(mmm_vars.options.scan_range_idxs)) # unlimited axis (can be appended to).

    xb = ncfile.createVariable('xb', np.float64, ('xb', 'time'))
    xb[:] = getattr(mmm_vars, 'xb').values[:, mmm_vars.options.scan_range_idxs]
    time = ncfile.createVariable('time', np.float64, ('time'))
    time[:] = getattr(mmm_vars, 'time').values[mmm_vars.options.scan_range_idxs]

    for dim in ncfile.dimensions.items():
        print(dim)



    print('--- Dimensions Complete ---')
    for name in all_vars:
        dims = getattr(mmm_vars, name).dimensions or ['XB', 'TIME']
        dims = ['xb' if d == 'XB0' else d for d in dims]
        dims = ['xb' if d == 'XBO' else d for d in dims]
        dims = ['time' if d == 'TIME3' else d for d in dims]
        dims = tuple([d.lower() for d in dims])
        # print(name, dims)
        if name == 'time' or name == 'xb' or name == 'x':
            continue
        if isinstance(getattr(mmm_vars, name).values, np.ndarray):
            var = ncfile.createVariable(name.lower(), np.float64, dims)
            var[:] = getattr(mmm_vars, name).values[:, mmm_vars.options.scan_range_idxs]
    # print(ncfile)

    mmm_vars.options.set_time_ranges(mmm_vars.time.values[mmm_vars.options.scan_range_idxs])
    print(mmm_vars.time.values)
    print(mmm_vars.options.scan_range_idxs)

if __name__ == '__main__':

    try: ncfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    ncfile = Dataset('new.CDF',mode='w',format='NETCDF4_CLASSIC') 
    print(ncfile)