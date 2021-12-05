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

# Parent class for input and output variables
class Variables:
    def __str__(self):
        return str(self.get_nonzero_variables())

    def get_variables(self):
        return [var for var in dir(self) if not callable(getattr(self, var)) and not var.startswith("__")]
        
    def get_nonzero_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).values is not None]

    def print_nonzero_variables(self):
        vars = self.get_nonzero_variables()
        for var in vars:
            print(var + ", "
                  + str(getattr(self, var).name) + ", "
                  + str(getattr(self, var).desc) + ", " 
                  + str(getattr(self, var).units) + ", "
                  + str(getattr(self, var).values.shape) + ", "
                  + str(getattr(self, var).dimensions))

# Variables obtained from a CDF
class InputVariables(Variables):
    def __init__(self):
        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME') # TODO: What is TIME3 in CDF?
        self.x = Variable('X', cdfvar='X', label=r'$x$')
        self.xb = Variable('XB', cdfvar='XB', label=r'$x_\mathrm{B}$')

        # CDF Variables needed for calculations
        self.aimp = Variable('Impurity Mean Mass', cdfvar='AIMP', label=r'$\overline{M}_\mathrm{imp}$', minvalue=1e-6, smooth=2)
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT', smooth=1)
        self.bz = Variable('BZ', cdfvar='BZ', smooth=1)
        self.elong = Variable('Elongation', cdfvar='ELONG', label=r'$\kappa$', smooth=1)
        self.omega = Variable('Toroidal Angular Velocity', cdfvar='OMEGA', smooth=2)
        self.ne = Variable('Electron Density', cdfvar='NE', label=r'$n_\mathrm{e}$', minvalue=1e-6, smooth=2)
        self.nf = Variable('Fast Ion Density', cdfvar='BDENS', label=r'$n_\mathrm{f}$', minvalue=1e-6, smooth=2)
        self.nd = Variable('Deuterium Ion Density', cdfvar='ND', label=r'$n_d$', minvalue=1e-6, smooth=2)
        self.nz = Variable('Impurity Density', cdfvar='NIMP', label=r'$n_z$', minvalue=1e-6, smooth=2)
        self.q = Variable('Safety Factor', cdfvar='Q', label=r'$q$', minvalue=1e-6, smooth=2)
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP', label=r'$R$', smooth=None)
        self.te = Variable('Electron Temperature', cdfvar='TE', label=r'$T_\mathrm{e}$', minvalue=1e-6, smooth=2)
        self.tepro = Variable('Electron Temperature', cdfvar='TEPRO', label=r'$T_\mathrm{e}$', minvalue=1e-6, smooth=2)
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI',label=r'$T_\mathrm{i}$', minvalue=1e-6, smooth=2)
        self.tipro = Variable('Thermal Ion Temperature', cdfvar='TIPRO', label=r'$T_\mathrm{i}$', minvalue=1e-6, smooth=2)
        self.vpold = Variable('VPOL', cdfvar='VPOLD_NC', smooth=3)
        self.vpolh = Variable('VPOL', cdfvar='VPOLH_NC', smooth=3)
        self.wexbs = Variable(r'ExB Shear Rate', cdfvar='SREXBA', label=r'$\omega_{E \times B}$', smooth=2)
        self.zimp = Variable('Mean Charge of Impurities', cdfvar='XZIMP', label=r'$\overline{Z}_\mathrm{imp}$', smooth=2)

        # Additional CDF variables for comparisons
        self.betat = Variable('BETAT', cdfvar='BETAT', smooth=1)

        # Calculated Variables (some are also in the CDF)
        self.aimass = Variable('Thermal Ion Mean Mass', label=r'$\overline{M}_\mathrm{i}$')
        self.ahyd = Variable('Hydrogenic Ion Mean Mass', label=r'$\overline{M}_\mathrm{h}$')
        self.alphamhd = Variable('Alpha MHD', label=r'$\alpha_\mathrm{MHD}$')
        self.beta = Variable('Pressure Ratio', cdfvar='BETAT', label=r'$\beta$')
        self.betae = Variable('Electron Pressure Ratio', cdfvar='BETAE', label=r'$\beta_\mathrm{\,e}$') 
        self.btor = Variable('Toroidal Magnetic Field', cdfvar='BTTOT', label=r'$B_\mathrm{T}$')
        self.eps = Variable('Inverse Aspect Ratio', label=r'$\epsilon$')
        self.etae = Variable('Electron Gradient Ratio', cdfvar='ETAE', label=r'$\eta_\mathrm{\,e}$')
        self.etai = Variable('Ion Gradient Ratio', label=r'$\eta_\mathrm{\,i}$') # ETAI in CDF is not gTI/gNI
        self.etaih = Variable('Hydrogenic Gradient Ratio', cdfvar='ETAIH', label=r'$\eta_\mathrm{\,ih}$')
        self.etaid = Variable('ETAID', label=r'$\eta_\mathrm{\,id}$')
        self.etaie = Variable('ETAIE', cdfvar='ETAIE', label=r'$\eta_\mathrm{\,ie}$')
        self.gave = Variable('Avg Curvature, Magnetic Field', label=r'$g_\mathrm{ave}$')
        self.gmax = Variable('Max Gradient', label=r'$g_\mathrm{max}$')
        self.gyrfi = Variable('Ion Gyrofrequency', label=r'$\omega_\mathrm{ci}$')
        self.loge = Variable('Coulomb Logarithm', cdfvar='CLOGE', label=r'$\ln\, \Lambda_\mathrm{e}$')
        self.nh = Variable('Hydrogenic Ion Density', cdfvar='NH', label=r'$n_\mathrm{h}$', smooth=2)
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', label=r'$n_\mathrm{i}$', smooth=1)
        self.nuei = Variable('Electron Collision Frequency', label=r'$\nu_\mathrm{ei}$')
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality', cdfvar='NUSTE', label=r'$\nu^{*}_\mathrm{e}$')
        self.nusti = Variable('Ion Collisionality', cdfvar='NUSTI', label=r'$\nu^{*}_\mathrm{i}$')
        self.p = Variable('Plasma Pressure', cdfvar='PPLAS', label=r'$p$') 
        self.raxis = Variable('RAXIS')
        self.rho = Variable('Radius', label=r'$\rho$')
        self.rmin = Variable('Minor Radius', label=r'$r$')
        self.shat = Variable('Effective Magnetic Shear', cdfvar='SHAT', label=r'$\hat{s}$')
        self.shear = Variable('Magnetic Shear', label=r'$s$')
        self.tau = Variable('Temperature Ratio', label=r'$\tau$')
        self.vpar = Variable('Parallel Velocity', label=r'$v_\mathrm{par}$', absminvalue=1e-6, smooth=2)
        self.vpol = Variable('Poloidal Velocity', label=r'$v_\theta$', absminvalue=1e-6, smooth=2)
        self.vtor = Variable('Toroidal Velocity', label=r'$v_\phi$', absminvalue=1e-6, smooth=2)
        self.vthe = Variable('Electron Thermal Velocity', label=r'$v_{T_\mathrm{e}}$')
        self.vthi = Variable('Ion Thermal Velocity', label=r'$v_{T_\mathrm{i}}$')
        self.zeff = Variable('Effective Charge', cdfvar='ZEFF', label=r'$Z_\mathrm{eff}$')

        # Calculated Gradients
        self.gne = Variable('Electron Density Gradient', label=r'$g_{n_\mathrm{e}}$')
        self.gnh = Variable('Hydrogenic Ion Density Gradient', label=r'$g_{n_\mathrm{h}}$')
        self.gni = Variable('Thermal Ion Density Gradient', smooth=0, label=r'$g_{n_\mathrm{i}}$')
        self.gnz = Variable('Impurity Density Gradient', label=r'$g_{n_\mathrm{z}}$')
        self.gnd = Variable('Deuterium Ion Density Gradient', label=r'$g_{n_\mathrm{d}}$')
        self.gq = Variable('Safety Factor Gradient', label=r'$g_{q}$')
        self.gte = Variable('Electron Temperature Gradient', label=r'$g_{T_\mathrm{e}}$')
        self.gti = Variable('Thermal Ion Temperature Gradient', smooth=0, label=r'$g_{T_\mathrm{i}}$')
        self.gvpar = Variable('Parallel Velocity Gradient', label=r'$g_{v_\mathrm{par}}$')
        self.gvpol = Variable('Poloidal Velocity Gradient', label=r'$g_{v_\theta}$')
        self.gvtor = Variable('Toroidal Velocity Gradient', label=r'$g_{v_\phi}$')

        # Test Variables
        self.test = Variable('Test Variable')
        self.test2 = Variable('Test Variable 2')
        self.gtest = Variable('Test Variable Gradient')

    def get_cdf_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).cdfvar is not None]

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

# Variables obtained from MMM output
class OutputVariables(Variables):
    def __init__(self):
        self.rho = Variable('rho', label=r'$\rho$')
        self.rmin = Variable('rmin', label=r'$r_\mathrm{min}$')
        self.xti = Variable('xti', label='xti')
        self.xdi = Variable('xdi', label='xdi')
        self.xte = Variable('xte', label='xte')
        self.xdz = Variable('xdz', label='xdz')
        self.xvt = Variable('xvt', label='xvt')
        self.xvp = Variable('xvp', label='xvp')
        self.xtiW20 = Variable('xtiW20', label='xtiW20')
        self.xdiW20 = Variable('xdiW20', label='xdiW20')
        self.xteW20 = Variable('xteW20', label='xteW20')
        self.xtiDBM = Variable('xtiDBM', label='xtiDBM')
        self.xdiDBM = Variable('xdiDBM', label='xdiDBM')
        self.xteDBM = Variable('xteDBM', label='xteDBM')
        self.xteETG = Variable('xteETG', label='xteETG')
        self.xteMTM = Variable('xteMTM', label='xteMTM')
        self.xteETGM = Variable('xteETGM', label='xteETGM')
        self.xdiETGM = Variable('xdiETGM', label='xdiETGM')
        self.gmaW20ii = Variable('gmaW20ii', label='gmaW20ii')
        self.omgW20ii = Variable('omgW20ii', label='omgW20ii')
        self.gmaW20ie = Variable('gmaW20ie', label='gmaW20ie')
        self.omgW20ie = Variable('omgW20ie', label='omgW20ie')
        self.gmaW20ei = Variable('gmaW20ei', label='gmaW20ei')
        self.omgW20ei = Variable('omgW20ei', label='omgW20ei')
        self.gmaW20ee = Variable('gmaW20ee', label='gmaW20ee')
        self.omgW20ee = Variable('omgW20ee', label='omgW20ee')
        self.gmaDBM = Variable('gmaDBM', label='gmaDBM')
        self.omgDBM = Variable('omgDBM', label='omgDBM')
        self.gmaMTM = Variable('gmaMTM', label='gmaMTM')
        self.omgMTM = Variable('omgMTM', label='omgMTM')
        self.gmaETGM = Variable('gmaETGM', label='gmaETGM')
        self.omgETGM = Variable('omgETGM', label='omgETGM')
        self.dbsqprf = Variable('dbsqprf', label='dbsqprf')

class Variable:
    def __init__(self, name, cdfvar=None, smooth=None, label='', desc='', minvalue=None, absminvalue=None, units='', dimensions=None, values=None):
        # Public
        self.name = name
        self.cdfvar = cdfvar # Name of variable as used in CDF's
        self.smooth = smooth # None to disable smoothing, or n = 1, 2, 3, ...  
        self.label = label # Plot label in LaTeX Format
        self.desc = desc
        self.minvalue = minvalue
        self.absminvalue = absminvalue
        # Private
        self._units_label = ''
        self._units = units
        self._dimensions = dimensions if dimensions is not None else ['','']
        self._values = values

    def __str__(self):
        return str(self.name)

    # Minimum values are used to handle variables that cannot take values below a minimum amount (such as negative Temperature)
    # Absolute minimum values are used to handle variables that can be negative, but get too close to zero
    def set_minvalue(self):
        if self.minvalue is not None:
            self.values[self.values < self.minvalue] = self.minvalue
        if self.absminvalue is not None:
            too_small = np.absolute(self.values) < self.absminvalue
            if too_small.any():
                # np.sign(0) = 0, so set these to +1
                value_signs = np.sign(self.values[too_small])
                value_signs[value_signs == 0] = 1
                self.values[too_small] = self.absminvalue * value_signs

    def get_xdim(self):
        return self.dimensions[0] if self.dimensions is not None and len(self.dimensions) > 0 else None

    def set_xdim(self, xdim):
        if self.dimensions is not None and len(self.dimensions) > 0:
            self.dimensions[0] = xdim
        else:
            raise ValueError(f'Failed to set xdim on variable {self.name}')

    @property
    def units_label(self):
        return self._units_label if self._units_label is not None else self._units

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        self._units = units

        # Set units_label in LaTeX format
        if units != '':
            # First item is the search string, second item is the replacement string
            unit_strs = [
                ['N/M**3', r'$\mathrm{m}^{-3}$'],
                ['M**2/SEC', r'$\mathrm{m}^{2}/s$'],
                ['M/SEC', r'm/s'],
                ['M', r'm'],
                ['SEC**-1', r's$^{-1}$'],
                ['MAMPS', r'MA'],
                ['RAD/SEC', r'rad/s'],
                ['PASCALS', r'Pa'],
                ['SECONDS', r's'],
                ['TESLA', r'T'],
                ['EV', r'eV'],
                ['kEV', r'keV'],
                ['m/s^2', r'm/s$^2$'],
                ['m^2/s', r'm$^2$/s'],
                ['s^-1', r's$^{-1}$']]

            for unit_str in unit_strs:
                if (unit_str[0] == self._units):
                    self._units_label = unit_str[1]
                    break

    @property
    def dimensions(self):
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions):
        if type(dimensions) == list:
            self._dimensions = dimensions
        else:
            raise ValueError(f'Variable dimensions must be {list} and not {type(dimensions)}')

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        if type(values) == np.ndarray:
            self._values = values
        else:
            raise ValueError(f'Variable values must be {np.ndarray} and not {type(values)}')
    
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
            self.values = scipy.ndimage.gaussian_filter(self.values, sigma=(self.smooth, 0))

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

    def remove_nan(self):
        if np.isnan(self.values).any():
            print('nan values found for var ' + self.name)
            self.values[np.isnan(self.values)] = 0
            self.set_minvalue()
