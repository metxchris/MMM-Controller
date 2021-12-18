# Standard Packages
import sys; sys.path.insert(0, '../')
from copy import deepcopy

# 3rd Party Packages
import numpy as np
import scipy.ndimage
from multipledispatch import dispatch

# Local Packages
import main.constants as constants
import main.utils as utils
from main.enums import SaveType


# Parent class for input and output variables
class Variables:
    def __init__(self):
        self.rho = None
        self.rmin = None

    def __str__(self):
        return str(self.get_nonzero_variables())

    def get_variables(self):
        return [var for var in dir(self) if not callable(getattr(self, var)) and not var.startswith("__")]

    def get_nonzero_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).values is not None]

    def print_nonzero_variables(self):
        var_names = self.get_nonzero_variables()
        for v in var_names:
            print(f'{v}, '
                  f'{getattr(self, v).name}, '
                  f'{getattr(self, v).desc}, '
                  f'{getattr(self, v).units}, '
                  f'{getattr(self, v).values.shape}, '
                  f'{getattr(self, v).dimensions}')

    def set_rho_values(self):
        if self.rmin.values.ndim == 2:
            self.rho.values = self.rmin.values / self.rmin.values[-1, :]
        elif self.rmin.values.ndim == 1:
            # This is expected when loading data from rho files with rho = 0
            if self.rmin.values[-1] == 0:
                self.rho.values = np.zeros_like(self.rmin.values)
            else:
                self.rho.values = self.rmin.values / self.rmin.values[-1]

    def get_data_as_array(self, var_list, time_idx=None):
        '''
        Gets data from requested variables in array format

        Parameters:
        * var_list (list): List of names (str) of variables to get data for
        * time_idx (int): Index of the measurement_time, if applicable

        Returns:
        * data (np.ndarray): The requested data in a 2-dimensional array
        * header (str): The string of variable names corresponding to data in the array
        '''

        num_points = self.rmin.values.shape[0]
        num_vars = len(var_list)
        data = np.zeros((num_points, num_vars), dtype=float)
        header = ','.join(var_list)

        # InputVariables data uses time_idx
        if time_idx is not None:
            for i, var_name in enumerate(var_list):
                data[:, i] = getattr(self, var_name).values[:, time_idx]

        # OutputVariables data does not use time_idx
        else:
            for i, var_name in enumerate(var_list):
                data[:, i] = getattr(self, var_name).values

        return data, header

    def save_to_csv(self, data, header, save_type, options, scan_factor=None):
        '''
        Saves data in np.ndarray format to a CSV

        Parameters:
        * data (np.ndarray): The data to save
        * header (str): The header to be saved to the CSV
        * save_type (SaveType): The SaveType of the data being saved
        * options (OptionsData): Options.instance
        * scan_factor (float): The scan_factor, if doing a parameter scan
        '''

        # Save data to the sub folder of the parameter being scanned
        if scan_factor is not None:
            scan_factor_str = constants.SCAN_FACTOR_FMT_STR.format(scan_factor)
            save_dir = utils.get_var_to_scan_path(options.runid, options.scan_num, options.var_to_scan)
            file_name = (f'{save_dir}\\{save_type.name.capitalize()} {options.var_to_scan}'
                         f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}{scan_factor_str}.csv')

        # Save data to top level directory of scan folder
        else:
            save_dir = utils.get_scan_num_path(options.runid, options.scan_num)
            file_name = f'{save_dir}\\{options.runid} {save_type.name.capitalize()} Profiles.csv'

        utils.create_directory(save_dir)
        np.savetxt(file_name, data, header=header, fmt='%.6e', delimiter=',')

        print(f'{save_type.name.capitalize()} data saved to \n    {file_name}\n')

    def load_from_csv(self, save_type, runid, scan_num, var_to_scan=None, scan_factor=None, rho_value=None):
        '''
        Loads data from a CSV into the current Variables subclass object

        Parameters:
        * save_type (SaveType): The SaveType of the data being saved
        * runid (str): The runid of the CSV to use
        * scan_num (int): The scan number of the CSV to use
        * var_to_scan (str): The scanned variable of the CSV to use (optional)
        * scan_factor (float): The scan_factor, if doing a parameter scan (optional)
        * rho_value (str or float): The rho value of the CSV to use (optional)
        '''

        if rho_value is not None:
            rho_str = rho_value if type(rho_value) is str else constants.RHO_VALUE_FMT_STR.format(rho_value)
            dir_path = utils.get_rho_path(runid, scan_num, var_to_scan)
            file_path = (f'{dir_path}\\{save_type.name.capitalize()} '
                         f'rho{constants.RHO_VALUE_SEPARATOR}{rho_str}.csv')

        elif scan_factor is not None:
            scan_factor_str = constants.SCAN_FACTOR_FMT_STR.format(scan_factor)
            dir_path = utils.get_var_to_scan_path(runid, scan_num, var_to_scan)
            file_path = (f'{dir_path}\\{save_type.name.capitalize()} {var_to_scan}'
                         f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}{scan_factor_str}.csv')

        else:
            dir_path = utils.get_scan_num_path(runid, scan_num)
            file_path = f'{dir_path}\\{runid} {save_type.name.capitalize()} Profiles.csv'

        data_array = np.genfromtxt(file_path, delimiter=',', dtype=float, names=True)
        var_names = data_array.dtype.names

        if len(var_names) == 0:
            raise ValueError(f'No variable names were loaded from {file_path}')

        for var_name in var_names:
            getattr(self, var_name).values = data_array[var_name]

        if self.rmin.values is not None:
            self.set_rho_values()


# Variables obtained from a CDF
class InputVariables(Variables):
    def __init__(self):
        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME')  # TODO: What is TIME3 in CDF?
        self.x = Variable('X', cdfvar='X', label=r'$x$')
        self.xb = Variable('XB', cdfvar='XB', label=r'$x_\mathrm{B}$')

        # CDF Variables needed for calculations
        self.aimp = Variable('Mean Mass of Impurities',cdfvar='AIMP', label=r'$\overline{M}_\mathrm{imp}$', minvalue=1e-6, smooth=5, save_type=SaveType.INPUT)
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT', smooth=0)
        self.bz = Variable('BZ', cdfvar='BZ', smooth=None)
        self.elong = Variable('Elongation', cdfvar='ELONG', label=r'$\kappa$', smooth=5, save_type=SaveType.INPUT)
        self.omega = Variable('Toroidal Angular Velocity', cdfvar='OMEGA', smooth=8)
        self.ne = Variable('Electron Density', cdfvar='NE', label=r'$n_\mathrm{e}$', minvalue=1e-6, smooth=5, save_type=SaveType.INPUT)
        self.nf = Variable('Fast Ion Density', cdfvar='BDENS', label=r'$n_\mathrm{f}$', minvalue=1e-6, smooth=5, save_type=SaveType.INPUT)
        self.nd = Variable('Deuterium Ion Density', cdfvar='ND', label=r'$n_d$', minvalue=1e-6, smooth=5, save_type=SaveType.ADDITIONAL)
        self.nz = Variable('Impurity Density', cdfvar='NIMP', label=r'$n_z$', minvalue=1e-6, smooth=5, save_type=SaveType.INPUT)
        self.q = Variable('Safety Factor', cdfvar='Q', label=r'$q$', minvalue=1e-6, smooth=5, save_type=SaveType.INPUT)
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP', label=r'$R$', smooth=None, save_type=SaveType.INPUT)
        self.rmin = Variable('Minor Radius', cdfvar='RMNMP', label=r'$r$', smooth=None, save_type=SaveType.INPUT)
        self.te = Variable('Electron Temperature', cdfvar='TE', label=r'$T_\mathrm{e}$', minvalue=1e-6, smooth=8, save_type=SaveType.INPUT)
        self.tepro = Variable('Electron Temperature', cdfvar='TEPRO', label=r'$T_\mathrm{e}$', minvalue=1e-6, smooth=8)
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI', label=r'$T_\mathrm{i}$', minvalue=1e-6, smooth=8, save_type=SaveType.INPUT)
        self.tipro = Variable('Thermal Ion Temperature', cdfvar='TIPRO', label=r'$T_\mathrm{i}$', minvalue=1e-6, smooth=8)
        self.vpolavg = Variable('VPOL', cdfvar='VPOL_AVG', smooth=0)
        self.vpold = Variable('VPOL', cdfvar='VPOLD_NC', smooth=0)
        self.vpolh = Variable('VPOL', cdfvar='VPOLH_NC', smooth=0)
        self.wexbs = Variable(r'ExB Shear Rate', cdfvar='SREXBA', label=r'$\omega_{E \times B}$', smooth=5, save_type=SaveType.INPUT)
        self.zimp = Variable('Mean Charge of Impurities', cdfvar='XZIMP', label=r'$\overline{Z}_\mathrm{imp}$', smooth=5, save_type=SaveType.INPUT)

        # Additional CDF variables for comparisons
        self.betat = Variable('BETAT', cdfvar='BETAT')
        self.bzxr = Variable('BZXR', cdfvar='BZXR')

        # Calculated Variables (some are also in the CDF)
        self.aimass = Variable('Mean Mass of Thermal Ions', label=r'$\overline{M}_\mathrm{i}$', save_type=SaveType.INPUT)
        self.ahyd = Variable('Mean Mass of Hydrogenic Ions', label=r'$\overline{M}_\mathrm{h}$', save_type=SaveType.INPUT)
        self.alphamhd = Variable('Alpha MHD', label=r'$\alpha_\mathrm{MHD}$', save_type=SaveType.ADDITIONAL)
        self.beta = Variable('Pressure Ratio', cdfvar='BTPL', label=r'$\beta$', save_type=SaveType.ADDITIONAL)
        self.betae = Variable('Electron Pressure Ratio', cdfvar='BTE', label=r'$\beta_\mathrm{\,e}$', save_type=SaveType.ADDITIONAL)  # cdfvar='BETAE' is a scalar
        self.bpol = Variable('Poloidal Magnetic Field', cdfvar='BPOL', label=r'$B_\theta$', save_type=SaveType.ADDITIONAL)
        self.btor = Variable('Toroidal Magnetic Field', cdfvar='', label=r'$B_\phi$', save_type=SaveType.INPUT)
        self.eps = Variable('Inverse Aspect Ratio', label=r'$\epsilon$', save_type=SaveType.ADDITIONAL)
        self.etae = Variable('Electron Gradient Ratio', cdfvar='ETAE', label=r'$\eta_\mathrm{\,e}$', save_type=SaveType.ADDITIONAL)
        self.etai = Variable('Ion Gradient Ratio', label=r'$\eta_\mathrm{\,i}$', save_type=SaveType.ADDITIONAL)  # cdfvar='ETAI' in CDF is not gTI/gNI
        self.etaih = Variable('Hydrogenic Gradient Ratio', cdfvar='ETAIH', label=r'$\eta_\mathrm{\,ih}$')
        self.etaid = Variable('ETAID', label=r'$\eta_\mathrm{\,id}$')
        self.etaie = Variable('ETAIE', label=r'$\eta_\mathrm{\,ie}$')  # cdfvar='ETAIE' in CDF is not gTI/gNE
        self.gave = Variable('Avg Curvature of Magnetic Field', label=r'$G_\mathrm{ave}$', save_type=SaveType.ADDITIONAL)
        self.gmax = Variable('Max Gradient', label=r'$g_\mathrm{max}$', save_type=SaveType.ADDITIONAL)
        self.gyrfi = Variable('Ion Gyrofrequency', label=r'$\omega_\mathrm{ci}$', save_type=SaveType.ADDITIONAL)
        self.loge = Variable('Electron Coulomb Logarithm', cdfvar='CLOGE', label=r'$\lambda_\mathrm{e}$', save_type=SaveType.ADDITIONAL)
        self.logi = Variable('Ion Coulomb Logarithm', cdfvar='CLOGI', label=r'$\lambda_\mathrm{i}$')
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', label=r'$n_\mathrm{i}$', smooth=1, save_type=SaveType.ADDITIONAL)
        self.nh0 = Variable('Hydrogen Ion Density', cdfvar='NH', label=r'$n_\mathrm{h}$', smooth=5)
        self.nh = Variable('Total Hydrogenic Ion Density', label=r'$n_\mathrm{h,T}$', smooth=5, save_type=SaveType.INPUT)
        self.nuei = Variable('Electron Collision Frequency', label=r'$\nu_\mathrm{ei}$', save_type=SaveType.ADDITIONAL)
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality', cdfvar='NUSTE', label=r'$\nu^{*}_\mathrm{e}$', save_type=SaveType.ADDITIONAL)
        self.nusti = Variable('Ion Collisionality', cdfvar='NUSTI', label=r'$\nu^{*}_\mathrm{i}$', save_type=SaveType.ADDITIONAL)
        self.p = Variable('Plasma Pressure', cdfvar='PPLAS', label=r'$p$', save_type=SaveType.ADDITIONAL)
        self.rho = Variable('Normalized Radius', label=r'$\rho$')
        self.shat = Variable('Effective Magnetic Shear', cdfvar='SHAT', label=r'$\hat{s}$', save_type=SaveType.ADDITIONAL)  # MMM uses a different definition of shat than what cdfvar='SHAT' uses
        self.shear = Variable('Magnetic Shear', label=r'$s$', save_type=SaveType.ADDITIONAL)
        self.tau = Variable('Temperature Ratio', label=r'$\tau$', save_type=SaveType.ADDITIONAL)
        self.vpar = Variable('Parallel Velocity', label=r'$v_\parallel$', absminvalue=1e-1, smooth=None, save_type=SaveType.INPUT)
        self.vpol = Variable('Poloidal Velocity', label=r'$v_\theta$', absminvalue=1e-1, smooth=15, save_type=SaveType.INPUT)
        self.vtor = Variable('Toroidal Velocity', cdfvar='VTOR_AVG', label=r'$v_\phi$', absminvalue=1e-1, smooth=None, save_type=SaveType.INPUT)  # cdfvar='VTOR_AVG' is a slightly different VTOR than what we are using
        self.vthe = Variable('Electron Thermal Velocity', label=r'$v_{T_\mathrm{e}}$', save_type=SaveType.ADDITIONAL)
        self.vthi = Variable('Ion Thermal Velocity', label=r'$v_{T_\mathrm{i}}$', save_type=SaveType.ADDITIONAL)
        self.zeff = Variable('Effective Charge', cdfvar='ZEFFP', label=r'$Z_\mathrm{eff}$', save_type=SaveType.INPUT)

        # Calculated Gradients
        self.gne = Variable('Electron Density Gradient', label=r'$g_{n_\mathrm{e}}$', save_type=SaveType.INPUT)
        self.gnh = Variable('Hydrogenic Ion Density Gradient', label=r'$g_{n_\mathrm{h}}$', save_type=SaveType.INPUT)
        self.gni = Variable('Thermal Ion Density Gradient', smooth=0, label=r'$g_{n_\mathrm{i}}$', save_type=SaveType.INPUT)
        self.gnz = Variable('Impurity Density Gradient', label=r'$g_{n_\mathrm{z}}$', save_type=SaveType.INPUT)
        self.gnd = Variable('Deuterium Ion Density Gradient', label=r'$g_{n_\mathrm{d}}$')
        self.gq = Variable('Safety Factor Gradient', label=r'$g_{q}$', save_type=SaveType.INPUT)
        self.gte = Variable('Electron Temperature Gradient', label=r'$g_{T_\mathrm{e}}$', save_type=SaveType.INPUT)
        self.gti = Variable('Thermal Ion Temperature Gradient', smooth=0, label=r'$g_{T_\mathrm{i}}$', save_type=SaveType.INPUT)
        self.gvpar = Variable('Parallel Velocity Gradient', label=r'$g_{v_\mathrm{par}}$', save_type=SaveType.INPUT)
        self.gvpol = Variable('Poloidal Velocity Gradient', label=r'$g_{v_\theta}$', save_type=SaveType.INPUT)
        self.gvtor = Variable('Toroidal Velocity Gradient', label=r'$g_{v_\phi}$', save_type=SaveType.INPUT)

        # Test Variables
        self.test = Variable('Test Variable')
        self.test2 = Variable('Test Variable 2')
        self.gtest = Variable('Test Variable Gradient')

    def get_vars_of_type(self, save_type):
        ''' Returns (list): List of all variables with the specified save_type'''
        nonzero_variables = self.get_nonzero_variables()
        return [v for v in nonzero_variables if getattr(self, v).save_type == save_type]

    def get_cdf_variables(self):
        ''' Returns (list): List of all variables where cdfvar is not None'''
        all_variables = self.get_variables()
        return [v for v in all_variables if getattr(self, v).cdfvar is not None]

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

    def save_vars_of_type(self, save_type, options, scan_factor=None):
        # Put rmin at the front of the variable list
        var_list = self.get_vars_of_type(save_type)
        if 'rmin' not in var_list:
            var_list.insert(0, 'rmin')
        else:
            var_list.insert(0, var_list.pop(var_list.index('rmin')))

        data, header = self.get_data_as_array(var_list, options.time_idx)
        self.save_to_csv(data, header, save_type, options, scan_factor)

    def save_all_vars(self, options, scan_factor=None):
        self.save_vars_of_type(SaveType.INPUT, options, scan_factor)
        self.save_vars_of_type(SaveType.ADDITIONAL, options, scan_factor)


# Variables obtained from MMM output
class OutputVariables(Variables):
    def __init__(self):
        self.rho = Variable('rho', units='', label=r'$\rho$')
        self.rmin = Variable('rmin', units='m', label=r'$r_\mathrm{min}$')
        # Total Diffusivities
        self.xti = Variable('xti', units='m^2/s', label='xti')
        self.xdi = Variable('xdi', units='m^2/s', label='xdi')
        self.xte = Variable('xte', units='m^2/s', label='xte')
        self.xdz = Variable('xdz', units='m^2/s', label=r'$D_\mathrm{n, w}$')
        self.xvt = Variable('xvt', units='m^2/s', label=r'$\chi_{v_\phi}$')
        self.xvp = Variable('xvp', units='m^2/s', label=r'$\chi_{v_\theta}$')
        # Weiland Components
        self.xtiW20 = Variable('xtiW20', units='m^2/s', label=r'$\chi_\mathrm{i, w}$')
        self.xdiW20 = Variable('xdiW20', units='m^2/s', label='xdiW20')
        self.xteW20 = Variable('xteW20', units='m^2/s', label='xteW20')
        self.gmaW20ii = Variable('gmaW20ii', units='s^-1', label='gmaW20ii')
        self.omgW20ii = Variable('omgW20ii', units='s^-1', label='omgW20ii')
        self.gmaW20ie = Variable('gmaW20ie', units='s^-1', label='gmaW20ie')
        self.omgW20ie = Variable('omgW20ie', units='s^-1', label='omgW20ie')
        self.gmaW20ei = Variable('gmaW20ei', units='s^-1', label='gmaW20ei')
        self.omgW20ei = Variable('omgW20ei', units='s^-1', label='omgW20ei')
        self.gmaW20ee = Variable('gmaW20ee', units='s^-1', label='gmaW20ee')
        self.omgW20ee = Variable('omgW20ee', units='s^-1', label='omgW20ee')
        # DBM Components
        self.xtiDBM = Variable('xtiDBM', units='m^2/s', label='xtiDBM')
        self.xdiDBM = Variable('xdiDBM', units='m^2/s', label='xdiDBM')
        self.xteDBM = Variable('xteDBM', units='m^2/s', label='xteDBM')
        self.gmaDBM = Variable('gmaDBM', units='s^-1', label='gmaDBM')
        self.omgDBM = Variable('omgDBM', units='s^-1', label='omgDBM')
        # ETG Component
        self.xteETG = Variable('xteETG', units='m^2/s', label='xteETG')
        # MTM Components
        self.xteMTM = Variable('xteMTM', units='m^2/s', label='xteMTM')
        self.gmaMTM = Variable('gmaMTM', units='s^-1', label='gmaMTM')
        self.omgMTM = Variable('omgMTM', units='s^-1', label='omgMTM')
        # ETGM Components
        self.xteETGM = Variable('xteETGM', units='m^2/s', label=r'$\chi_\mathrm{e, etgm}$')
        self.xdiETGM = Variable('xdiETGM', units='m^2/s', label=r'$D_\mathrm{n, etgm}$')
        self.gmaETGM = Variable('gmaETGM', units='s^-1', label=r'$\gamma_\mathrm{etgm}$')
        self.omgETGM = Variable('omgETGM', units='s^-1', label=r'$\omega_\mathrm{etgm}$')

        self.dbsqprf = Variable('dbsqprf', units='', label=r'$|\delta B/B|^2$')

    def get_all_output_vars(self):
        all_vars = self.get_variables()
        all_vars.remove('rho')
        all_vars.remove('rmin')
        return all_vars

    def get_etgm_vars(self):
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'ETGM' in var]

    def get_mtm_vars(self):
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'MTM' in var]

    def get_dbm_vars(self):
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'DBM' in var]

    def get_etg_vars(self):
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'ETG' in var and 'ETGM' not in var]

    def get_weiland_vars(self):
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'W20' in var]

    def save_all_vars(self, options, scan_factor=None):
        # Put rmin at the front of the variable list
        var_list = self.get_variables()
        var_list.insert(0, var_list.pop(var_list.index('rmin')))
        var_list.remove('rho')

        data, header = self.get_data_as_array(var_list)
        self.save_to_csv(data, header, SaveType.OUTPUT, options, scan_factor)


class Variable:
    def __init__(self,
                 name,
                 cdfvar=None,
                 smooth=None,
                 label='',
                 desc='',
                 minvalue=None,
                 absminvalue=None,
                 save_type=None,
                 units='',
                 dimensions=None,
                 values=None):
        # Public
        self.name = name
        self.cdfvar = cdfvar  # Name of variable as used in CDF's
        self.smooth = smooth  # None to disable smoothing, or n = integer value (n in N)
        self.label = label  # Plot label in LaTeX Format
        self.desc = desc
        self.minvalue = minvalue
        self.absminvalue = absminvalue
        self.save_type = save_type
        # Private
        self._units_label = ''
        self._units = ''
        self._dimensions = dimensions if dimensions is not None else ['', '']
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
                value_signs = np.sign(self.values[too_small])
                value_signs[value_signs == 0] = 1  # np.sign(0) = 0, so set these to +1
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
        return self._units_label

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        self._units = units
        self._units_label = units

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

    def set(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def apply_smoothing(self, input_points):
        '''
        Variable smoothing using a Gaussian filter

        The value of sigma needs to increase nearly linearly as the amount of input_points increases, to maintain
        the same level of smoothing.  To achieve this, the self.smooth value is multiplied by input_points.

        Parameters:
        * input_points (int): Number of radial points each variable has
        '''

        if self.smooth is not None:
            sigma = int(input_points * self.smooth / 100)
            self.values = scipy.ndimage.gaussian_filter(self.values, sigma=(sigma, 0))

    # Clamps values between -value and value, and sets origin value to approximately 0
    def clamp_gradient(self, value):
        self.values[0, :] = 1e-6
        self.values[self.values > value] = value
        self.values[self.values < -value] = -value

    # Removes values outside of m standard deviations
    # TODO: Currently not ideal since removed values are replaced with None,
    # which turns everything into nan after smoothing or interpolating again
    def reject_outliers(self, m=4):
        self.values[(np.abs(self.values - np.mean(self.values)) > m * np.std(self.values))] = None

    def remove_nan(self):
        if np.isnan(self.values).any():
            print('nan values found for var ' + self.name)
            self.values[np.isnan(self.values)] = 0
            self.set_minvalue()


# For testing purposes
if __name__ == '__main__':
    from options import Options
    Options.instance.set(
        runid='TEST',
        scan_num=10,
        var_to_scan='gti',
        time_idx=0,
    )

    # Create InputVariables and OutputVariables, and populate with non-zero values
    ivars = InputVariables()
    input_var_names = ivars.get_variables()
    for i, var_name in enumerate(input_var_names):
        values = np.ones((5, 4), dtype=float) * i
        getattr(ivars, var_name).set(name='name', desc='desc', units='units', dimensions=['X', 'T'], values=values)

    ovars = OutputVariables()
    output_var_names = ovars.get_variables()
    for i, var_name in enumerate(output_var_names):
        values = np.ones((5), dtype=float) * i
        getattr(ovars, var_name).set(name='name', desc='desc', units='units', dimensions=['X', 'T'], values=values)

    # Save variable data to CSV
    # ivars.save_all_vars(Options.instance)
    # ovars.save_all_vars(Options.instance)

    ivars = InputVariables()
    ivars.load_from_csv(SaveType.INPUT, 'TEST', 25, 'tau', scan_factor=1.5)
    ivars.load_from_csv(SaveType.ADDITIONAL, 'TEST', 25, 'tau', scan_factor=1.5, rho_value=0.5)

    ivars.print_nonzero_variables()
