# Standard Packages
from copy import deepcopy
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import scipy.ndimage
from multipledispatch import dispatch

# Local Packages
from main import constants


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
        self.aimp = Variable('Mean Mass of Impurities', cdfvar='AIMP', label=r'$\overline{M}_\mathrm{imp}$', minvalue=1e-6, smooth=2)
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
        self.vpolavg = Variable('VPOL', cdfvar='VPOL_AVG', smooth=3)
        self.vpold = Variable('VPOL', cdfvar='VPOLD_NC', smooth=3)
        self.vpolh = Variable('VPOL', cdfvar='VPOLH_NC', smooth=3)
        self.wexbs = Variable(r'ExB Shear Rate', cdfvar='SREXBA', label=r'$\omega_{E \times B}$', smooth=2)
        self.zimp = Variable('Mean Charge of Impurities', cdfvar='XZIMP', label=r'$\overline{Z}_\mathrm{imp}$', smooth=2)

        # Additional CDF variables for comparisons
        self.betat = Variable('BETAT', cdfvar='BETAT', smooth=1)
        self.bzxr = Variable('BZXR', cdfvar='BZXR')

        # Calculated Variables (some are also in the CDF)
        self.aimass = Variable('Mean Mass of Thermal Ions', label=r'$\overline{M}_\mathrm{i}$')
        self.ahyd = Variable('Mean Mass of Hydrogenic Ions', label=r'$\overline{M}_\mathrm{h}$')
        self.alphamhd = Variable('Alpha MHD', label=r'$\alpha_\mathrm{MHD}$')
        self.beta = Variable('Pressure Ratio', cdfvar='BTPL', label=r'$\beta$')
        self.betae = Variable('Electron Pressure Ratio', cdfvar='BTE', label=r'$\beta_\mathrm{\,e}$') # cdfvar='BETAE' is a scalar
        self.bpol = Variable('Poloidal Magnetic Field', cdfvar='BPOL', label=r'$B_\theta$')
        self.btor = Variable('Toroidal Magnetic Field', cdfvar='', label=r'$B_\phi$')
        self.eps = Variable('Inverse Aspect Ratio', label=r'$\epsilon$')
        self.etae = Variable('Electron Gradient Ratio', cdfvar='ETAE', label=r'$\eta_\mathrm{\,e}$')
        self.etai = Variable('Ion Gradient Ratio', label=r'$\eta_\mathrm{\,i}$') # cdfvar='ETAI' in CDF is not gTI/gNI
        self.etaih = Variable('Hydrogenic Gradient Ratio', cdfvar='ETAIH', label=r'$\eta_\mathrm{\,ih}$')
        self.etaid = Variable('ETAID', label=r'$\eta_\mathrm{\,id}$')
        self.etaie = Variable('ETAIE', cdfvar='ETAIE', label=r'$\eta_\mathrm{\,ie}$')
        self.gave = Variable('Avg Curvature of Magnetic Field', label=r'$G_\mathrm{ave}$')
        self.gmax = Variable('Max Gradient', label=r'$g_\mathrm{max}$')
        self.gyrfi = Variable('Ion Gyrofrequency', label=r'$\omega_\mathrm{ci}$')
        self.loge = Variable('Electron Coulomb Logarithm', cdfvar='CLOGE', label=r'$\lambda_\mathrm{e}$')
        self.logi = Variable('Ion Coulomb Logarithm', cdfvar='CLOGI', label=r'$\lambda_\mathrm{i}$')
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', label=r'$n_\mathrm{i}$', smooth=1)
        self.nh0 = Variable('Hydrogen Ion Density', cdfvar='NH', label=r'$n_\mathrm{h}$', smooth=2)
        self.nh = Variable('Total Hydrogenic Ion Density', label=r'$n_\mathrm{h,T}$')
        self.nuei = Variable('Electron Collision Frequency', label=r'$\nu_\mathrm{ei}$')
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality', cdfvar='NUSTE', label=r'$\nu^{*}_\mathrm{e}$')
        self.nusti = Variable('Ion Collisionality', cdfvar='NUSTI', label=r'$\nu^{*}_\mathrm{i}$')
        self.p = Variable('Plasma Pressure', cdfvar='PPLAS', label=r'$p$') 
        self.raxis = Variable('RAXIS')
        self.rho = Variable('Radius', label=r'$\rho$')
        self.rmin = Variable('Minor Radius', cdfvar='RMNMP', label=r'$r$')
        self.shat = Variable('Effective Magnetic Shear', cdfvar='SHAT', label=r'$\hat{s}$') # MMM uses a different definition of shat than what cdfvar='SHAT' uses
        self.shear = Variable('Magnetic Shear', label=r'$s$')
        self.tau = Variable('Temperature Ratio', label=r'$\tau$')
        self.vpar = Variable('Parallel Velocity', label=r'$v_\parallel$', absminvalue=1e-6, smooth=2)
        self.vpol = Variable('Poloidal Velocity', label=r'$v_\theta$', absminvalue=1e-6, smooth=2)
        self.vtor = Variable('Toroidal Velocity', cdfvar='VTOR_AVG', label=r'$v_\phi$', absminvalue=1e-6, smooth=2) # cdfvar='VTOR_AVG' is a slightly different VTOR than what we are using
        self.vthe = Variable('Electron Thermal Velocity', label=r'$v_{T_\mathrm{e}}$')
        self.vthi = Variable('Ion Thermal Velocity', label=r'$v_{T_\mathrm{i}}$')
        self.zeff = Variable('Effective Charge', cdfvar='ZEFFP', label=r'$Z_\mathrm{eff}$')

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
        self.rho = Variable('rho', units='', label=r'$\rho$')
        self.rmin = Variable('rmin', units='m', label=r'$r_\mathrm{min}$')
        self.xti = Variable('xti', units='m^2/s', label='xti')
        self.xdi = Variable('xdi', units='m^2/s', label='xdi')
        self.xte = Variable('xte', units='m^2/s', label='xte')
        self.xdz = Variable('xdz', units='m^2/s', label=r'$D_\mathrm{nw}$')
        self.xvt = Variable('xvt', units='m^2/s', label=r'$\chi_{v_\phi}$')
        self.xvp = Variable('xvp', units='m^2/s', label=r'$\chi_{v_\theta}$')
        self.xtiW20 = Variable('xtiW20', units='m^2/s', label=r'$\chi_\mathrm{iw}$')
        self.xdiW20 = Variable('xdiW20', units='m^2/s', label='xdiW20')
        self.xteW20 = Variable('xteW20', units='m^2/s', label='xteW20')
        self.xtiDBM = Variable('xtiDBM', units='m^2/s', label='xtiDBM')
        self.xdiDBM = Variable('xdiDBM', units='m^2/s', label='xdiDBM')
        self.xteDBM = Variable('xteDBM', units='m^2/s', label='xteDBM')
        self.xteETG = Variable('xteETG', units='m^2/s', label='xteETG')
        self.xteMTM = Variable('xteMTM', units='m^2/s', label='xteMTM')
        self.xteETGM = Variable('xteETGM', units='m^2/s', label='xteETGM')
        self.xdiETGM = Variable('xdiETGM', units='m^2/s', label='xdiETGM')
        self.gmaW20ii = Variable('gmaW20ii', units='s^-1', label='gmaW20ii')
        self.omgW20ii = Variable('omgW20ii', units='s^-1', label='omgW20ii')
        self.gmaW20ie = Variable('gmaW20ie', units='s^-1', label='gmaW20ie')
        self.omgW20ie = Variable('omgW20ie', units='s^-1', label='omgW20ie')
        self.gmaW20ei = Variable('gmaW20ei', units='s^-1', label='gmaW20ei')
        self.omgW20ei = Variable('omgW20ei', units='s^-1', label='omgW20ei')
        self.gmaW20ee = Variable('gmaW20ee', units='s^-1', label='gmaW20ee')
        self.omgW20ee = Variable('omgW20ee', units='s^-1', label='omgW20ee')
        self.gmaDBM = Variable('gmaDBM', units='s^-1', label='gmaDBM')
        self.omgDBM = Variable('omgDBM', units='s^-1', label='omgDBM')
        self.gmaMTM = Variable('gmaMTM', units='s^-1', label='gmaMTM')
        self.omgMTM = Variable('omgMTM', units='s^-1', label='omgMTM')
        self.gmaETGM = Variable('gmaETGM', units='s^-1', label='gmaETGM')
        self.omgETGM = Variable('omgETGM', units='s^-1', label='omgETGM')
        self.dbsqprf = Variable('dbsqprf', units='', label='dbsqprf')

    def get_vars_to_plot(self):
        vars_to_plot = self.get_variables()
        vars_to_plot.remove('rho')
        vars_to_plot.remove('rmin')
        return vars_to_plot


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
        self._units = ''
        self._dimensions = dimensions if dimensions is not None else ['','']
        self._values = values

        # Call units setter to also set units_label
        self.units = units

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
            for unit_str in constants.UNIT_STRINGS:
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
        if self.smooth is not None:
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
        self.values[(np.abs(self.values - np.mean(self.values)) > m * np.std(self.values))] = None

    def remove_nan(self):
        if np.isnan(self.values).any():
            print('nan values found for var ' + self.name)
            self.values[np.isnan(self.values)] = 0
            self.set_minvalue()
