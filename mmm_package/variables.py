# Standard Packages
from copy import deepcopy
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np
import scipy.ndimage
from multipledispatch import dispatch
# Local Packages
import settings

class Variables:
    def __init__(self):
        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME') # TODO: What is TIME3 in CDF?
        self.x = Variable('X', cdfvar='X')
        self.xb = Variable('XB', cdfvar='XB')

        # CDF Variables
        self.aimp = Variable('AIMP', cdfvar='AIMP', mmmvar='aimp', smooth=1)
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT', smooth=1)
        self.bz = Variable('BZ', cdfvar='BZ', smooth=1)
        self.elong = Variable('Elongation', cdfvar='ELONG', mmmvar='elong', label=r'$\kappa$', smooth=1)
        self.omega = Variable('OMEGA', cdfvar='OMEGA', smooth=1)
        self.ne = Variable('Electron Density', cdfvar='NE', mmmvar='ne', label=r'$n_\mathrm{e}$', smooth=2)
        self.nf = Variable('Fast Ion Density', cdfvar='BDENS', mmmvar='nf', label=r'$n_\mathrm{f}$', smooth=1)
        self.nd = Variable('ND', cdfvar='ND', label=r'$n_d$', smooth=1)
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', label=r'$n_\mathrm{i}$', smooth=2)
        self.nz = Variable('Impurity Density', cdfvar='NIMP', mmmvar='nz', label=r'$n_z$', smooth=2)
        self.pcur = Variable('PCUR', cdfvar='PCUR', smooth=1)
        self.q = Variable('Safety Factor', cdfvar='Q', mmmvar='q', label=r'$q$', smooth=2)
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP', mmmvar='rmaj', smooth=None)
        self.te = Variable('Electron Temperature', cdfvar='TE', mmmvar='te', label=r'$T_\mathrm{e}$', smooth=2)
        self.tepro = Variable('Electron Temperature', cdfvar='TEPRO', mmmvar='te', label=r'$T_\mathrm{e}$', smooth=1)
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI', mmmvar='ti', label=r'$T_\mathrm{i}$', smooth=2)
        self.tipro = Variable('Thermal Ion Temperature', cdfvar='TIPRO', mmmvar='ti', label=r'$T_\mathrm{i}$', smooth=1)
        self.triang = Variable('TRIANG', cdfvar='TRIANG', smooth=1)
        self.vpold = Variable('VPOL', cdfvar='VPOLD_NC', smooth=1)
        self.vpolh = Variable('VPOL', cdfvar='VPOLH_NC', smooth=1)
        self.wexbs = Variable('ExB Shear Rate', cdfvar='SREXBA', mmmvar='wexbs', label=r'$\omega_{E \times B}$', smooth=1)
        self.zimp = Variable('Mean Charge of Impurities', cdfvar='XZIMP', mmmvar='zimp', smooth=1)

        # Calculated Variables (some are also in the CDF)
        # TODO: Check that calculated values match CDF values
        self.aimass = Variable('Thermal Ion Mean Atomic Mass', mmmvar='aimass')
        self.ahyd = Variable('Hydrogenic Ion Mean Atomic Mass', mmmvar='ahyd')
        self.alphamhd = Variable('Alpha_MHD', label=r'$\alpha_\mathrm{MHD}$')
        self.beta = Variable('Pressure Ratio', label=r'$\beta$')
        self.betae = Variable('Electron Pressure Ratio', label=r'$\beta_\mathrm{\,e}$') # cdfvar='BETAE'
        self.btor = Variable('Toroidal Magnetic Field', mmmvar='btor')
        self.eps = Variable('Inverse Aspect Ratio')
        self.etae = Variable('Electron Gradient Ratio', label=r'$\eta_\mathrm{\,e}$')
        self.etai = Variable('Ion Gradient Ratio', label=r'$\eta_\mathrm{\,i}$')
        self.nh = Variable('Hydrogenic Ion Density', mmmvar='nh', label=r'$n_\mathrm{h}$',smooth=1) # cdfvar='NH'
        self.nuei = Variable('Collision Frequency')
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality', label=r'$\nu^{*}_\mathrm{e}$') # cdfvar='NUSTE'
        self.nusti = Variable('Ion Collisionality', label=r'$\nu^{*}_\mathrm{i}$') # cdfvar='NUSTI'
        self.p = Variable('Plasma Pressure') # cdfvar='PPLAS'
        self.raxis = Variable('RAXIS')
        self.rho = Variable('Radius', label=r'$\rho$')
        self.rmin = Variable('Minor Radius', mmmvar='rmin')
        self.shat = Variable('Effective Magnetic Shear', label=r'$\hat{s}$') # cdfvar='SHAT'
        self.shear = Variable('Magnetic Shear', label=r'$s$')
        self.vpar = Variable('Parallel Velocity', mmmvar='vpar')
        self.vpol = Variable('Poloidal Velocity', mmmvar='vpol', label=r'$v_\theta$')
        self.vtor = Variable('Toroidal Velocity', mmmvar='vtor', label=r'$v_\phi$')
        self.tau = Variable('Temperature Ratio', label=r'$\tau$')
        self.zeff = Variable('Effective Charge', mmmvar='zeff') # cdfvar='ZEFF'
        self.zgmax = Variable('ZGMAX')
        self.zgyrfi = Variable('Ion Gyrofrequency')
        self.zlog = Variable('Coulomb Logarithm')
        self.zvthe = Variable('Electron Thermal Velocity')
        self.zvthi = Variable('Ion Thermal Velocity')

        # Calculated Gradients
        self.gne = Variable('Electron Density Gradient', mmmvar='gne', label=r'$g_{n_\mathrm{e}}$')
        self.gnh = Variable('Hydrogenic Ion Density Gradient', mmmvar='gnh', label=r'$g_{n_\mathrm{h}}$')
        self.gni = Variable('Thermal Ion Density Gradient', mmmvar='gni', smooth=0, label=r'$g_{n_\mathrm{i}}$')
        self.gnz = Variable('Impurity Density Gradient', mmmvar='gnz', label=r'$g_{n_\mathrm{z}}$')
        self.gq = Variable('Safety Factor Gradient', mmmvar='gq', label=r'$g_{q}$')
        self.gte = Variable('Electron Temperature Gradient', mmmvar='gte', label=r'$g_{T_\mathrm{e}}$')
        self.gti = Variable('Thermal Ion Temperature Gradient', mmmvar='gti', smooth=0, label=r'$g_{T_\mathrm{i}}$')
        self.gvpar = Variable('Parallel Velocity Gradient', mmmvar='gvpar')
        self.gvpol = Variable('Poloidal Velocity Gradient', mmmvar='gvpol', label=r'$g_{\nu_\theta}$')
        self.gvtor = Variable('Toroidal Velocity Gradient', mmmvar='gvtor', label=r'$g_{\nu_\phi}$')

    def get_variables(self):
        return [var for var in dir(self) if not callable(getattr(self, var)) and not var.startswith("__")]
        
    def get_nonzero_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).values is not None]

    def get_cdf_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).cdfvar is not None]

    def get_mmm_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).mmmvar is not None]

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

    def get_nboundaries(self):
        return self.xb.values.shape[0] if self.xb.values is not None and self.xb.values.ndim > 0 else 0

    def get_ntimes(self):
        return self.x.values.shape[1] if self.xb.values is not None and self.xb.values.ndim > 1 else 0

    def use_temperature_profiles(self):
        if self.tepro.values is not None:
            self.te = deepcopy(self.tepro)
        else:
            raise ValueError('Failed to set TEPRO since TEPRO is None')

        if self.tipro.values is not None:
            self.ti = deepcopy(self.tipro)
        else:
            raise ValueError('Failed to set TIPRO since TIPRO is None')

class Variable:
    def __init__(self, name, cdfvar=None, mmmvar=None, smooth=None, label=None, desc=None, units=None, dimensions=None, values=None):
        # Public
        self.name = name
        self.cdfvar = cdfvar # Name of variable as used in CDF's
        self.mmmvar = mmmvar # Name of variable as used in MMM
        self.smooth = smooth # None to disable smoothing, or n = 1, 2, 3, ...  
        self.label = label if label is not None else '' # LaTeX Format
        self.desc = desc if desc is not None else ''
        # Private
        self._units = units if units is not None else ''
        self._dimensions = dimensions if dimensions is not None else ['','']
        self._values = values

    def __str__(self):
        return str(self.name)

    def get_xdim(self):
        return self.dimensions[0] if self.dimensions is not None and len(self.dimensions) > 0 else None

    def set_xdim(self, xdim):
        if self.dimensions is not None and len(self.dimensions) > 0:
            self.dimensions[0] = xdim
        else:
            raise ValueError('Failed to set xdim on variable {0}'.format(self.name))

    @property
    def dimensions(self):
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions):
        if type(dimensions) == list:
            self._dimensions = dimensions
        else:
            raise ValueError('Variable dimensions must be {0} and not {1}'.format(list, type(dimensions)))
    
    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        self._units = units

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        if type(values) == np.ndarray:
            self._values = values
        else:
            raise ValueError('Variable values must be {0} and not {1}'.format(np.ndarray, type(values)))
    
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
        if self.smooth is not None and settings.APPLY_SMOOTHING:
            self.values = scipy.ndimage.gaussian_filter(self.values, sigma=self.smooth)

    # Clamps values between -value and value, and sets origin value to apprximately 0
    def clamp_gradient(self, value):
        self.values[0, :] = 1e-6
        self.values[self.values > value] = value
        self.values[self.values < -value] = -value

    # Removes values outside of m standard deviations
    # TODO: Currently not ideal since removed values are replaced with None,
    # which turns everything into nan after smoothing or intepolating again
    def reject_outliers(self, m=4):
        if settings.REMOVE_OUTLIERS:
            self.values[(np.abs(self.values - np.mean(self.values)) > m * np.std(self.values))] = None

class InputOptions:
    def __init__(self, cdf_name, shot_type=None, input_time=None, input_points=None):
        # Public
        self.cdf_name = cdf_name
        self.shot_type = shot_type
        self.input_time = input_time
        self.input_points = input_points
        # Private
        self._runid = None
        self._interp_points = None
        self._time = None
        self._time_idx = None

    @property
    def runid(self):
        return self._runid

    @runid.setter
    def runid(self, runid):
        self._runid = runid.strip()
        if self._runid != self.cdf_name:
            print('*** WARNING: runid {0} does not match cdf_name {1}'.format(self.runid, self.cdf_name))

    @property
    def interp_points(self):
        return self._interp_points

    @interp_points.setter
    def interp_points(self, points):
        self._interp_points = points
        if self._interp_points < self.input_points:
            print('*** WARNING: possible interpolation points ({0}) is less than specified input points ({1})'
                .format(self.interp_points, self.input_points))

    @property
    def time_idx(self):
        return self._time_idx

    @property
    def time(self):
        return self._time

    # Find the index of the measurement time closest to the input_time and index and the value
    def set_measurement_time(self, tvar):
        self._time_idx = np.argmin(np.abs(tvar.values - self.input_time))
        self._time = "{:.3f}".format(tvar.values[self.time_idx])