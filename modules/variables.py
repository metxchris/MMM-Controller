# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import scipy.ndimage

# Local Packages
import modules.constants as constants
import modules.utils as utils
from modules.enums import SaveType


# Used to create units labels to display on plots from units strings
UNITS_TO_UNITS_LABEL = {
    'T*m': r'T\,m',
    'm^-3': r'm$^{-3}$',
    'm/s^2': r'm/s$^2$',
    'm^2/s': r'm$^2$/s',
    's^-1': r's$^{-1}$',
}


# Parent class for input and output variables
class Variables:
    def __init__(self):
        self.rho = None
        self.rmin = None

    def __str__(self):
        return str(self.get_nonzero_variables())

    def get_variables(self):
        '''Returns (list of str): all variable names'''
        return [var for var in dir(self) if not callable(getattr(self, var)) and not var.startswith("__")]

    def get_nonzero_variables(self):
        '''Returns (list of str): variable names with nonzero values'''
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
        '''Sets rho from rmin'''
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

        # InputVariables data can use time_idx
        if time_idx is not None:
            for i, var_name in enumerate(var_list):
                data[:, i] = getattr(self, var_name).values[:, time_idx]

        # OutputVariables data does not use time_idx
        else:
            for i, var_name in enumerate(var_list):
                data[:, i] = getattr(self, var_name).values

        return data, header

    def save_to_csv(self, data, header, save_type, options, scan_factor=None, rho_value=None):
        '''
        Saves data in np.ndarray format to a CSV

        Parameters:
        * data (np.ndarray): The data to save
        * header (str): The header to be saved to the CSV
        * save_type (SaveType): The SaveType of the data being saved
        * options (Options): An instance of the Options class
        * scan_factor (float): The scan_factor, if doing a parameter scan
        '''

        args = (save_type, options.runid, options.scan_num, options.var_to_scan, scan_factor, rho_value)
        dir_path, file_path = self.get_csv_save_path(*args)
        utils.create_directory(dir_path)
        np.savetxt(file_path, data, header=header, fmt='%.6e', delimiter=',')

        print(f'{save_type.name.capitalize()} data saved to \n    {file_path}\n')

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

        __, file_path = self.get_csv_save_path(save_type, runid, scan_num, var_to_scan, scan_factor, rho_value)
        self.load_from_file_path(file_path)

    def load_from_file_path(self, file_path):
        '''
        Loads data from a file into the current Variables subclass object

        Parameters:
        * file_path (str): The path of the file to load
        '''

        # TODO: Add check if file exists
        data_array = np.genfromtxt(file_path, delimiter=',', dtype=float, names=True)
        var_names = data_array.dtype.names

        if len(var_names) == 0:
            raise ValueError(f'No variable names were loaded from {file_path}')

        for var_name in var_names:
            getattr(self, var_name).values = data_array[var_name]

        if self.rmin.values is not None:
            self.set_rho_values()

    def get_csv_save_path(self, save_type, runid, scan_num, var_to_scan=None, scan_factor=None, rho_value=None):
        '''
        Gets the path where a CSV of variable data will be saved, based on the input parameter values

        Parameters:
        * save_type (SaveType): The SaveType of the data being saved
        * runid (str): The runid of the CSV to use
        * scan_num (int): The scan number of the CSV to use
        * var_to_scan (str): The scanned variable of the CSV to use (optional)
        * scan_factor (float): The scan_factor, if doing a parameter scan (optional)
        * rho_value (str or float): The rho value of the CSV to use (optional)
        '''

        if rho_value is not None:
            rho_str = rho_value if isinstance(rho_value, str) else f'{rho_value:{constants.RHO_VALUE_FMT}}'
            dir_path = utils.get_rho_path(runid, scan_num, var_to_scan)
            file_path = (f'{dir_path}\\{save_type.name.capitalize()} '
                         f'rho{constants.RHO_VALUE_SEPARATOR}{rho_str}.csv')

        elif scan_factor is not None:
            scan_factor_str = f'{scan_factor:{constants.SCAN_FACTOR_FMT}}'
            dir_path = utils.get_var_to_scan_path(runid, scan_num, var_to_scan)
            file_path = (f'{dir_path}\\{save_type.name.capitalize()} {var_to_scan}'
                         f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}{scan_factor_str}.csv')

        else:
            dir_path = utils.get_scan_num_path(runid, scan_num)
            file_path = f'{dir_path}\\{runid} {save_type.name.capitalize()} Profiles.csv'

        return dir_path, file_path


class InputVariables(Variables):
    '''
    Input variables are defined as anything that isn't read as output from the MMM driver

    All members are defined using the Variable class.  See the Variable class definition for more info.
    '''

    def __init__(self):
        self.rho = Variable('Normalized Radius', label=r'$\rho$')
        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME')
        self.x = Variable('X', cdfvar='X', label=r'$x$')
        self.xb = Variable('XB', cdfvar='XB', label=r'$x_\mathrm{B}$')

        # CDF Variables needed for calculations
        self.aimp = Variable('Mean Mass of Impurities', cdfvar='AIMP', label=r'$\overline{M}_\mathrm{imp}$',
                             save_type=SaveType.INPUT, minvalue=1, smooth=1)
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT')
        self.bz = Variable('BZ', cdfvar='BZ')
        self.elong = Variable('Elongation', cdfvar='ELONG', label=r'$\kappa$', smooth=1,
                              save_type=SaveType.INPUT)
        self.ne = Variable('Electron Density', cdfvar='NE', label=r'$n_\mathrm{e}$', minvalue=1e-6, smooth=1,
                           save_type=SaveType.INPUT, units='m^-3')
        self.nf = Variable('Fast Ion Density', cdfvar='BDENS', label=r'$n_\mathrm{f}$', minvalue=1e-6, smooth=1,
                           save_type=SaveType.INPUT, units='m^-3')
        self.nd = Variable('Deuterium Ion Density', cdfvar='ND', label=r'$n_d$', minvalue=1e-6, smooth=1,
                           save_type=SaveType.ADDITIONAL, units='m^-3')
        self.nz = Variable('Impurity Density', cdfvar='NIMP', label=r'$n_z$', minvalue=1e-6, smooth=1,
                           save_type=SaveType.INPUT, units='m^-3')
        self.q = Variable('Safety Factor', cdfvar='Q', label=r'$q$', minvalue=1e-6, smooth=1,
                          save_type=SaveType.INPUT)
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP', label=r'$R$',
                             save_type=SaveType.INPUT, units='m', minvalue=0)
        self.rmin = Variable('Minor Radius', cdfvar='RMNMP', label=r'$r$',
                             save_type=SaveType.INPUT, units='m', minvalue=0)
        self.te = Variable('Electron Temperature', cdfvar='TE', label=r'$T_\mathrm{e}$', minvalue=1e-6, smooth=1,
                           save_type=SaveType.INPUT, units='keV')
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI', label=r'$T_\mathrm{i}$', minvalue=1e-6, smooth=1,
                           save_type=SaveType.INPUT, units='keV')
        self.vpol = Variable('Poloidal Velocity', cdfvar='VPOL_AVG', label=r'$v_\theta$', absminvalue=1e-6, smooth=3,
                             save_type=SaveType.INPUT, units='m/s')
        self.vtor = Variable('Toroidal Velocity', cdfvar='VTOR_AVG', label=r'$v_\phi$', absminvalue=1e-6, smooth=3,
                             save_type=SaveType.INPUT, units='m/s')
        self.wexbs = Variable(r'ExB Shear Rate', cdfvar='SREXBA', label=r'$\omega_{E \times B}$', smooth=1,
                              save_type=SaveType.INPUT, units='s^-1')
        self.zimp = Variable('Mean Charge of Impurities', cdfvar='XZIMP', label=r'$\overline{Z}_\mathrm{imp}$', smooth=1,
                             save_type=SaveType.INPUT)

        # Additional CDF variables
        self.betat = Variable('BETAT', cdfvar='BETAT')
        self.bzxr = Variable('BZXR', cdfvar='BZXR')
        self.tepro = Variable('Electron Temperature', cdfvar='TEPRO')
        self.tipro = Variable('Thermal Ion Temperature', cdfvar='TIPRO')

        # Calculated Variables (some are also in the CDF)
        self.aimass = Variable('Mean Mass of Thermal Ions', label=r'$\overline{M}_\mathrm{i}$',
                               save_type=SaveType.INPUT, minvalue=1)
        self.ahyd = Variable('Mean Mass of Hydrogenic Ions', label=r'$\overline{M}_\mathrm{h}$',
                             save_type=SaveType.INPUT, minvalue=1)
        self.alphamhd = Variable('Alpha MHD', label=r'$\alpha_\mathrm{MHD}$',
                                 save_type=SaveType.ADDITIONAL)
        self.beta = Variable('Pressure Ratio', cdfvar='BTPL', label=r'$\beta$',
                             save_type=SaveType.ADDITIONAL, minvalue=0)
        self.betae = Variable('Electron Pressure Ratio', cdfvar='BTE', label=r'$\beta_\mathrm{\,e}$',
                              save_type=SaveType.ADDITIONAL, minvalue=0)  # cdfvar='BETAE' is a scalar
        self.bpol = Variable('Poloidal Magnetic Field', cdfvar='BPOL', label=r'$B_\theta$',
                             save_type=SaveType.ADDITIONAL, units='T')
        self.btor = Variable('Toroidal Magnetic Field', cdfvar='', label=r'$B_\phi$',
                             save_type=SaveType.INPUT, units='T', absminvalue=1e-6)
        self.eps = Variable('Inverse Aspect Ratio', label=r'$\epsilon$',
                            save_type=SaveType.ADDITIONAL)
        self.etae = Variable('Electron Gradient Ratio', cdfvar='ETAE', label=r'$\eta_\mathrm{\,e}$',
                             save_type=SaveType.ADDITIONAL)
        self.etai = Variable('Ion Gradient Ratio', label=r'$\eta_\mathrm{\,i}$',
                             save_type=SaveType.ADDITIONAL)  # cdfvar='ETAI' in CDF is not gTI/gNI
        self.gave = Variable('Avg Curvature of Magnetic Field', label=r'$G_\mathrm{ave}$',
                             save_type=SaveType.ADDITIONAL)
        self.gmax = Variable('Max Gradient', label=r'$g_\mathrm{max}$',
                             save_type=SaveType.ADDITIONAL)
        self.gyrfi = Variable('Ion Gyrofrequency', label=r'$\omega_\mathrm{ci}$',
                              save_type=SaveType.ADDITIONAL, units='s^-1')
        self.loge = Variable('Electron Coulomb Logarithm', cdfvar='CLOGE', label=r'$\lambda_\mathrm{e}$',
                             save_type=SaveType.ADDITIONAL)
        self.logi = Variable('Ion Coulomb Logarithm', cdfvar='CLOGI', label=r'$\lambda_\mathrm{i}$')
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', label=r'$n_\mathrm{i}$', units='m^-3',
                           save_type=SaveType.ADDITIONAL, minvalue=1e-6)
        self.nh0 = Variable('Hydrogen Ion Density', cdfvar='NH', label=r'$n_\mathrm{h}$', units='m^-3',
                            save_type=SaveType.ADDITIONAL)
        self.nh = Variable('Total Hydrogenic Ion Density', label=r'$n_\mathrm{h,T}$',
                           save_type=SaveType.INPUT, units='m^-3', minvalue=1e-6)
        self.nuei = Variable('Electron Collision Frequency', label=r'$\nu_\mathrm{ei}$',
                             save_type=SaveType.ADDITIONAL)
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality', cdfvar='NUSTE', label=r'$\nu^{*}_\mathrm{e}$',
                              save_type=SaveType.ADDITIONAL)
        self.nusti = Variable('Ion Collisionality', cdfvar='NUSTI', label=r'$\nu^{*}_\mathrm{i}$',
                              save_type=SaveType.ADDITIONAL)
        self.p = Variable('Plasma Pressure', cdfvar='PPLAS', label=r'$p$',
                          save_type=SaveType.ADDITIONAL, minvalue=1e-6)
        self.shat = Variable('Effective Magnetic Shear', cdfvar='SHAT', label=r'$\hat{s}$',
                             save_type=SaveType.ADDITIONAL)  # MMM uses a different definition of shat than cdfvar='SHAT'
        self.shear = Variable('Magnetic Shear', label=r'$s$',
                              save_type=SaveType.ADDITIONAL)
        self.tau = Variable('Temperature Ratio', label=r'$\tau$',
                            save_type=SaveType.ADDITIONAL, minvalue=0)
        self.vpar = Variable('Parallel Velocity', label=r'$v_\parallel$', absminvalue=1e-6,
                             save_type=SaveType.INPUT, units='m/s')
        self.vthe = Variable('Electron Thermal Velocity', label=r'$v_{T_\mathrm{e}}$',
                             save_type=SaveType.ADDITIONAL, units='m/s')
        self.vthi = Variable('Ion Thermal Velocity', label=r'$v_{T_\mathrm{i}}$',
                             save_type=SaveType.ADDITIONAL, units='m/s')
        self.zeff = Variable('Effective Charge', cdfvar='ZEFFP', label=r'$Z_\mathrm{eff}$',
                             save_type=SaveType.INPUT, minvalue=1)

        # Calculated Gradients
        self.gne = Variable('Electron Density Gradient', label=r'$g_{n_\mathrm{e}}$',
                            save_type=SaveType.INPUT)
        self.gnh = Variable('Hydrogenic Ion Density Gradient', label=r'$g_{n_\mathrm{h, T}}$',
                            save_type=SaveType.INPUT)
        self.gni = Variable('Thermal Ion Density Gradient', label=r'$g_{n_\mathrm{i}}$',
                            save_type=SaveType.INPUT)
        self.gnz = Variable('Impurity Density Gradient', label=r'$g_{n_\mathrm{z}}$',
                            save_type=SaveType.INPUT)
        self.gq = Variable('Safety Factor Gradient', label=r'$g_{q}$',
                           save_type=SaveType.INPUT)
        self.gte = Variable('Electron Temperature Gradient', label=r'$g_{T_\mathrm{e}}$',
                            save_type=SaveType.INPUT)
        self.gti = Variable('Thermal Ion Temperature Gradient', label=r'$g_{T_\mathrm{i}}$',
                            save_type=SaveType.INPUT)
        self.gvpar = Variable('Parallel Velocity Gradient', label=r'$g_{v_\parallel}$',
                              save_type=SaveType.INPUT)
        self.gvpol = Variable('Poloidal Velocity Gradient', label=r'$g_{v_\theta}$',
                              save_type=SaveType.INPUT)
        self.gvtor = Variable('Toroidal Velocity Gradient', label=r'$g_{v_\phi}$',
                              save_type=SaveType.INPUT)

    def set_x_values(self):
        '''
        Sets x from xb

        x is the grid between xb, and has one fewer point than xb.
        '''
        self.x.values = (self.xb.values[0:-1, :] + self.xb.values[1:, :]) / 2

    def get_vars_of_type(self, save_type):
        '''Returns (list of str): List of all variables with the specified save_type'''
        nonzero_variables = self.get_nonzero_variables()
        return [v for v in nonzero_variables if getattr(self, v).save_type == save_type]

    def get_cdf_variables(self):
        '''Returns (list of str): List of all variables where cdfvar is not None'''
        all_variables = self.get_variables()
        return [v for v in all_variables if getattr(self, v).cdfvar is not None]

    def get_nboundaries(self):
        '''Returns (int): The number of boundary points in the radial dimension of xb'''
        return self.xb.values.shape[0] if self.xb.values is not None and self.xb.values.ndim > 0 else 0

    def get_ntimes(self):
        '''Returns (int): The number of points in the time dimension of xb'''
        return self.x.values.shape[1] if self.xb.values is not None and self.xb.values.ndim > 1 else 0

    def use_temperature_profiles(self):
        '''Attempts to use experimental temperature profiles in place of calculated profiles'''
        if self.tepro.values is not None:
            self.te.values = self.tepro.values
        else:
            raise ValueError('Failed to set TEPRO since TEPRO is None')

        if self.tipro.values is not None:
            self.ti.values = self.tipro.values
        else:
            raise ValueError('Failed to set TIPRO since TIPRO is None')

    def save_vars_of_type(self, save_type, options, scan_factor=None):
        '''
        Saves variable values of the specified save type to a CSV

        Parameters:
        * save_type (SaveType): The save type of the variables to save
        * options (Options): An instance of the Options class
        * scan_factor (scan_factor): The value of the scan factor (Optional)
        '''

        # Put rmin at the front of the variable list
        var_list = self.get_vars_of_type(save_type)
        if 'rmin' not in var_list:
            var_list.insert(0, 'rmin')
        else:
            var_list.insert(0, var_list.pop(var_list.index('rmin')))

        data, header = self.get_data_as_array(var_list, options.time_idx)
        self.save_to_csv(data, header, save_type, options, scan_factor)

    def save_all_vars(self, options, scan_factor=None):
        '''
        Saves variable values of all relevant save types to a CSV

        Parameters:
        * save_type (SaveType): The save type of the variables to save
        * options (Options): An instance of the Options class
        * scan_factor (scan_factor): The value of the scan factor (Optional)
        '''

        self.save_vars_of_type(SaveType.INPUT, options, scan_factor)
        self.save_vars_of_type(SaveType.ADDITIONAL, options, scan_factor)


class OutputVariables(Variables):
    '''Output variables consist of all variable data obtained as output from MMM (other than rho and rmin)'''

    def __init__(self):
        # Independent Variables
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
        '''Returns (list of str): all output variable names (other than rho and rmin)'''
        all_vars = self.get_variables()
        all_vars.remove('rho')
        all_vars.remove('rmin')
        return all_vars

    def get_etgm_vars(self):
        '''Returns (list of str): all ETGM model variables'''
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'ETGM' in var]

    def get_mtm_vars(self):
        '''Returns (list of str): all MTM model variables'''
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'MTM' in var]

    def get_dbm_vars(self):
        '''Returns (list of str): all DRIBM model variables'''
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'DBM' in var]

    def get_etg_vars(self):
        '''Returns (list of str): all Horton ETG model variables'''
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'ETG' in var and 'ETGM' not in var]

    def get_weiland_vars(self):
        '''Returns (list of str): all Weiland model variables'''
        output_vars = self.get_all_output_vars()
        return [var for var in output_vars if 'W20' in var]

    def save_all_vars(self, options, scan_factor=None):
        '''Saves output variables to a CSV (other than rho)'''

        # Put rmin at the front of the variable list
        var_list = self.get_variables()
        var_list.insert(0, var_list.pop(var_list.index('rmin')))
        var_list.remove('rho')

        data, header = self.get_data_as_array(var_list)
        self.save_to_csv(data, header, SaveType.OUTPUT, options, scan_factor)


class Variable:
    def __init__(self, name, cdfvar=None, smooth=None, label='', desc='', minvalue=None,
                 absminvalue=None, save_type=None, mmm_label='', units='', dimensions=None, values=None):
        # Public
        self.name = name
        self.cdfvar = cdfvar  # Name of variable as used in CDF's
        self.smooth = smooth  # None to disable smoothing, or n = positive integer
        self.label = label  # Plot label in LaTeX Format
        self.desc = desc  # Stores the long_name value from CDF's
        self.minvalue = minvalue  # minimum value variable is allowed to have
        self.absminvalue = absminvalue  # minimum value the absolute value of the variable is allowed to have
        self.save_type = save_type if save_type is not None else SaveType.NONE
        # Private
        self._units_label = ''
        self._units = ''
        self._dimensions = dimensions if dimensions is not None else ['', '']
        self._values = values

        # Call units setter to also set units_label
        self.units = units

    def __str__(self):
        return str(self.name)

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
        if units in UNITS_TO_UNITS_LABEL.keys():
            self._units_label = UNITS_TO_UNITS_LABEL[units]

    @property
    def dimensions(self):
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions):
        if type(dimensions) != list:
            raise ValueError(f'Variable dimensions must be {list} and not {type(dimensions)}')
        self._dimensions = dimensions

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        if type(values) != np.ndarray:
            raise ValueError(f'Variable values must be {np.ndarray} and not {type(values)}')
        self._values = values

    def set(self, **kwargs):
        '''Sets members using keyword arguments'''
        for key, value in kwargs.items():
            setattr(self, key, value)

    def apply_smoothing(self, input_points):
        '''
        Variable smoothing using a Gaussian filter

        The value of sigma needs to increase nearly linearly as the amount of
        input_points increases, to maintain the same level of smoothing.  To
        achieve this, the self.smooth value is multiplied by input_points.
        Additionally, the amount of smoothing applied also depends on how
        closely grouped values are along the x-axis.  Specifically, for a
        given value of sigma, loosely spaced values will be smoothed more
        than tightly clustered values.

        For this reason, using the uniform rho input option will result in
        uniform smoothing across a variable (since the values in rmin become
        more clustered as rmin approaches its maximum value).  If the uniform
        rho input option is not used, then variables will have stronger
        smoothing applied near rmin = 0, and weaker smoothing applied near
        rmin = 1.

        Parameters:
        * input_points (int): Number of radial points each variable has
        '''

        if self.smooth is not None:
            sigma = int(input_points * self.smooth / 100)
            self.values = scipy.ndimage.gaussian_filter(self.values, sigma=(sigma, 0))

    def set_minvalue(self, raise_exception=True):
        '''
        Sets the minimum or absolute minimum value for a variable

        Minimum values are used to handle variables that cannot take values
        below a minimum amount (such as negative Temperatures). Due to
        expected persisting errors following the interpolation process, we
        don't raise an exception if there is at most one nonphysical value
        along the radial dimension at any point in time.  Instead, these
        errors are silently fixed.  However, an exception is raised if
        multiple nonphysical values are detected along the radial dimension
        for any point in time.

        Absolute minimum values are used to handle variables that are allowed
        to be negative, but can't get too close to zero (due to divide by
        zero issues).  No exceptions are raised if variables go below their
        absolute minimum value, since these are not considered to be errors.
        '''

        if self.minvalue is not None:
            multiple_errors_per_timeval = (np.count_nonzero(self.values < self.minvalue, axis=0) > 1)
            if raise_exception and multiple_errors_per_timeval.any():
                idx_list = [i for i in np.where(multiple_errors_per_timeval)][0]
                raise ValueError(
                    f'Multiple Nonphysical values obtained for {self.name}\n'
                    f'    min value:     {self.values[:, idx_list].min()}\n'
                    f'    time indices:  {idx_list}\n'
                )
            # When an exception is not raised, fix the minimum value
            self.values[self.values < self.minvalue] = self.minvalue

        if self.absminvalue is not None:
            too_small = np.absolute(self.values) < self.absminvalue
            if too_small.any():
                value_signs = np.sign(self.values[too_small])
                value_signs[value_signs == 0] = 1  # np.sign(0) = 0, so set these to +1
                self.values[too_small] = self.absminvalue * value_signs

    def clamp_gradient(self, max):
        '''
        Clamps values between -max and max, and sets origin value to
        approximately 0

        Gradient values will later be clamped again within the MMM driver by
        the variable gmax = g_{max}, which has values lower than the clamping
        value used here.  This step of clamping first here just helps with
        making the presentation of input profiles a bit cleaner.
        '''

        self.values[0, :] = 1e-6
        self.values[self.values > max] = max
        self.values[self.values < -max] = -max

    def check_for_nan(self):
        '''Checks for nan values and raises an exception if any are found'''
        if np.isnan(self.values).any():
            raise ValueError(f'nan values found in variable {self.name}')


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
