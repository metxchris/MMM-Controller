# 3rd Party Packages
import numpy as np
import scipy.ndimage
from multipledispatch import dispatch

class Variable(object):
    def __init__(self, name, cdfvar=None, mmmvar=None, smooth=None, label=None, values=None, desc=None, units=None, dimensions=None):
        self.name = name
        self.cdfvar = cdfvar # Name of variable as used in CDF's
        self.mmmvar = mmmvar # Name of variable as used in MMM
        self.smooth = smooth # None to disable smoothing, or n = 1, 2, 3, ...  
        self.values = values
        self.label = label if label is not None else '' # LaTeX Format
        self.desc = desc if desc is not None else ''
        self.units = units if units is not None else ''
        self.dimensions = dimensions if dimensions is not None else ['','']

    def __str__(self):
        return str(self.name)

    def get_xdim(self):
        return self.dimensions[0] if self.dimensions is not None and len(self.dimensions) > 0 else None

    def set_xdim(self, xdim):
        if self.dimensions is not None and len(self.dimensions) > 0:
            self.dimensions[0] = xdim
        else:
            print('[variables] *** Error: Failed to set xdim on variable', self.name)

    def get_dims(self):
        return self.dimensions

    def get_units(self):
        return self.units

    @dispatch(list)
    def set_dims(self, dimensions):
        self.dimensions = dimensions

    @dispatch(str)
    def set_units(self, units):
        self.units = units

    @dispatch(np.ndarray)
    def set_variable(self, values):
        self.set_variable(values, self.units, self.dimensions)

    @dispatch(np.ndarray, str)
    def set_variable(self, values, units):
        self.set_variable(values, units, self.dimensions)

    # Set variable values, units, dimensions
    @dispatch(np.ndarray, str, list)
    def set_variable(self, values, units, dimensions):
        self.values = values
        self.units = units
        self.dimensions = dimensions

    # Variable smoothing using a Gaussian filter
    def apply_smoothing(self):
        if self.smooth is not None:
            self.values = scipy.ndimage.gaussian_filter(self.values, sigma=self.smooth)

class Variables(object):
    def __init__(self):
        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME') # TODO: What is TIME3 in CDF?
        self.x = Variable('X', cdfvar='X')
        self.xb = Variable('XB', cdfvar='XB')

        # CDF Variables
        self.aimp = Variable('AIMP', cdfvar='AIMP', smooth=1)
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT', smooth=1)
        self.bz = Variable('BZ', cdfvar='BZ', smooth=1)
        self.elong = Variable('Elongation', cdfvar='ELONG', smooth=1)
        self.omega = Variable('OMEGA', cdfvar='OMEGA', smooth=1)
        self.ne = Variable('Electron Density', cdfvar='NE', smooth=1)
        self.nf = Variable('NF', cdfvar='BDENS', smooth=1)
        self.nh = Variable('Hydrogenic Ion Density', cdfvar='NH', smooth=1)
        self.nd = Variable('ND', cdfvar='ND', smooth=1)
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', smooth=1)
        self.nz = Variable('NZ', cdfvar='NIMP', smooth=1)
        self.pcur = Variable('PCUR', cdfvar='PCUR', smooth=1)
        self.q = Variable('Safety Factor', cdfvar='Q', smooth=1)
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP', smooth=None)
        self.te = Variable('Electron Temperature', cdfvar='TE', smooth=1)
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI', smooth=1)
        self.triang = Variable('TRIANG', cdfvar='TRIANG', smooth=1)
        self.vpold = Variable('VPOL', cdfvar='VPOLD_NC', smooth=3)
        self.vpolh = Variable('VPOL', cdfvar='VPOLH_NC', smooth=3)
        self.wexbs = Variable('ExB Shear Rate', cdfvar='SREXBA', smooth=1)
        self.zimp = Variable('ZIMP', cdfvar='XZIMP', smooth=1)

        # Calculated Variables (some are also in the CDF)
        # TODO: Check that calculated values match CDF values
        self.aimass = Variable('AIMASS')
        self.alphamhd = Variable('Alpha_MHD')
        self.beta = Variable('Beta')
        self.betae = Variable('Electron Beta') # cdfvar='BETAE'
        self.btor = Variable('Toroidal Magnetic Field')
        self.eps = Variable('Inverse Aspect Ratio')
        self.nuei = Variable('Collision Frequency')
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality') # cdfvar='NUSTE'
        self.nusti = Variable('Ion Collisionality') # cdfvar='NUSTI'
        self.p = Variable('Plasma Pressure') # cdfvar='PPLAS'
        self.raxis = Variable('RAXIS')
        self.rmin = Variable('Minor Radius')
        self.shat = Variable('Effective Magnetic Shear') # cdfvar='SHAT'
        self.shear = Variable('Magnetic Shear')
        self.vpar = Variable('VPAR')
        self.vpol = Variable('VPOL')
        self.vtor = Variable('VTOR')
        self.tau = Variable('Temperature Ratio')
        self.zeff = Variable('Effective Charge') # cdfvar='ZEFF'
        self.zgmax = Variable('ZGMAX')
        self.zgyrfi = Variable('Ion Gyrofrequency')
        self.zlog = Variable('Coulomb Logarithm')
        self.zvthe = Variable('Electron Thermal Velocity')
        self.zvthi = Variable('Ion Thermal Velocity')

        # Calculated Gradients
        self.gne = Variable('Electron Density Gradient')
        self.gnh = Variable('Hydrogenic Ion Density Gradient')
        self.gni = Variable('Thermal Ion Density Gradient')
        self.gnz = Variable('NZ Gradient')
        self.gq = Variable('Safety Factor Gradient')
        self.gte = Variable('Electron Temperature Gradient')
        self.gti = Variable('Thermal Ion Temperature Gradient')
        self.gvpar = Variable('VPAR Gradient')
        self.gvpol = Variable('VPOL Gradient')
        self.gvtor = Variable('VTOR Gradient')

    def get_variables(self):
        return [var for var in dir(self) if not callable(getattr(self, var)) and not var.startswith("__")]
        
    def get_nonzero_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).values is not None]

    def get_cdf_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).cdfvar is not None]

    def print_nonzero_variables(self):
        vars = self.get_nonzero_variables()
        for var in vars:
            print(var + ", "
                  + str(getattr(self, var).name) + ", "
                  + str(getattr(self, var).desc) + ", " 
                  + str(getattr(self, var).units) + ", "
                  + str(getattr(self, var).values.shape) + ", "
                  + str(getattr(self, var).dimensions))

    def __str__(self):
        return str(self.get_nonzero_variables())

    def get_nzones(self):
        return self.x.values.shape[0] if self.xb.values is not None and self.xb.values.ndim > 0 else 0

    def get_ntimes(self):
        return self.x.values.shape[1] if self.xb.values is not None and self.xb.values.ndim > 1 else 0
