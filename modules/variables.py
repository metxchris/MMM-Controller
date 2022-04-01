"""Handles all variables both needed for MMM input and produced as MMM output

The Variables class serves at the parent class to both InputVariables and
OutputVariables.  All variable data here will be defined in terms of rmin
(or rho).  The Variables class should never be instantiated directly; either
an InputVariables or OutputVariables class object should be instantiated.

InputVariables will have have two-dimensional value arrays when data is loaded
from a CDF (arrays radial and time values, in the order [position, time]).
However, value arrays will only be defined along the radial dimension when
data is loaded from a CSV.  InputVariables consist of all variables that are
used as direct input or calculation of direct input to MMM, as well as
additional variables that are calculated in MMM.  By calculating additional
variables that are not sent to MMM, we are able to plot and analyze their
values as needed.

OutputVariables values will always be one-dimensional arrays that are defined
on the radial dimension. OutputVariables consist of all variables that are
returned as output after running MMM.

In addition to storing values for each variable, other things such as variable
units, plot labels, minimum values, and the names of corresponding TRANSP
variables are stored here as well.  Each variable class also handles the
saving and loading of its data to CSVs, with the help of the utils class,
which provides the paths of directories for saving and loading of data.

The Variables class and its children are all coupled with the Options class,
which is used for storing options needed for plotting, checking values, and
both saving and loading variable data.  An Options object must be
instantiated and passed in when instantiating Variables objects.  When
loading variable data from CSVs, it is advised to load the options data
first, and then use the loaded options to load the variables.

Example Usage:
    # Instantiate and Load Options from CSV
    options = modules.options.Options()
    options.load(runid='120968A02', scan_num=1)

    # Instantiate InputVariables Object
    input_vars = variables.InputVariables(options)

    # Load Input Variables from CSV (both base and additional variables)
    input_vars.load_from_csv(SaveType.INPUT)
    input_vars.load_from_csv(SaveType.ADDITIONAL)

    # Save Input Variables to CSV
    input_vars.save()

    # Get List of variables with corresponding TRANSP values (in CDF)
    input_vars.get_cdf_variables()

    # Get all values of rmin at the first time index
    input_vars.rmin.values[:, 0]
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging

# 3rd Party Packages
import numpy as np
import scipy.ndimage

# Local Packages
import modules.constants as constants
import modules.utils as utils
from modules.enums import SaveType


_log = logging.getLogger(__name__)

# Used to create units labels to display on plots from units strings
_UNITS_TO_UNITS_LABEL = {
    'T*m': r'T$\,$m',
    'T*m^2': r'T$\,$m$^2$',
    'm^-3': r'm$^{-3}$',
    'm/s^2': r'm/s$^2$',
    'm^2/s': r'm$^2$/s',
    's^{-1}': r's$^{-1}$',
    'm^-1': r'm$^{-1}$',
    'keVm/s': r'keV$\,$m/s',
}


# Parent class for input and output variables
class Variables:
    def __init__(self, options):
        self.options = options

    def __str__(self):
        return str(self.get_nonzero_variables())

    def get_variables(self):
        '''Returns (list[str]): all variable names'''
        return [var for var in dir(self) if isinstance(getattr(self, var), Variable)]

    def get_nonzero_variables(self):
        '''Returns (list[str]): variable names with nonzero values'''
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).values is not None]

    def print_nonzero_variables(self):
        '''Prints various attributes of nonzero variables'''
        var_names = self.get_nonzero_variables()
        for v in var_names:
            print(f'{v}, '
                  f'{getattr(self, v).name}, '
                  f'{getattr(self, v).desc}, '
                  f'{getattr(self, v).units}, '
                  f'{getattr(self, v).values.shape}, '
                  f'{getattr(self, v).dimensions}')

    def set_radius_values(self):
        '''Sets rho from rmin'''
        if self.rmin.values.ndim == 2:
            self.rmina.values = self.rmin.values / self.rmin.values[-1, :]
            self.rho.values = np.tile(np.linspace(0, 1, self.rmin.values.shape[0]), (self.rmin.values.shape[1], 1)).T
        elif self.rmin.values.ndim == 1:
            # This is expected when loading data from rho files with rho = 0
            if self.rmin.values[-1] == 0:
                self.rmina.values = np.zeros_like(self.rmin.values)
                self.rho.values = np.zeros_like(self.rmin.values)
            else:
                self.rmina.values = self.rmin.values / self.rmin.values[-1]
                self.rho.values = np.linspace(0, 1, self.rmin.values.shape[0])

    def _get_data_as_array(self, var_list):
        '''
        Gets data from requested variables in array format

        Parameters:
        * var_list (list): List of names (str) of variables to get data for

        Returns:
        * data (np.ndarray): The requested data in a 2-dimensional array
        * header (str): The string of variable names corresponding to data in the array

        Raises:
        * ValueError: If the time index in options has not been initialized
        '''

        num_vars = len(var_list)
        data = np.zeros((self.options.input_points, num_vars), dtype=float)
        header = ','.join(var_list)

        if isinstance(self, InputVariables):
            if self.options.time_idx is None:
                raise ValueError('The time index has not been initialized')
            for i, var_name in enumerate(var_list):
                var_values = getattr(self, var_name).values
                if var_values is not None:
                    data[:, i] = var_values[:, self.options.time_idx]

        elif isinstance(self, OutputVariables):
            for i, var_name in enumerate(var_list):
                data[:, i] = getattr(self, var_name).values

        return data, header

    def _save_to_csv(self, data, header, save_type, scan_factor=None, rho_value=None):
        '''
        Saves data in np.ndarray format to a CSV

        Parameters:
        * data (np.ndarray): The data to save
        * header (str): The header to be saved to the CSV
        * save_type (SaveType): The SaveType of the data being saved
        * scan_factor (float): The scan_factor, if doing a parameter scan
        * rho_value (str | float): The rho value of the CSV to use (optional)
        '''

        dir_path, file_path = self._get_csv_save_path(save_type, scan_factor, rho_value)
        utils.create_directory(dir_path)
        np.savetxt(file_path, data, header=header, fmt='%.6e', delimiter=',')

        _log.info(f'\n\tSaved: {file_path}\n')

    def load_from_csv(self, save_type, scan_factor=None, rho_value=None):
        '''
        Loads data from a CSV into the current Variables subclass object

        Parameters:
        * save_type (SaveType): The SaveType of the data being saved
        * scan_factor (str | float): The scan_factor, if doing a parameter scan (optional)
        * rho_value (str | float): The rho value of the CSV to use (optional)
        '''

        __, file_path = self._get_csv_save_path(save_type, scan_factor, rho_value)
        self.load_from_file_path(file_path)

    def load_from_file_path(self, file_path):
        '''
        Loads data from a file into the current Variables subclass object

        Parameters:
        * file_path (str): The path of the file to load

        Raises:
        * ValueError: If no variable names are loaded
        '''

        # TODO: Add check if file exists
        data_array = np.genfromtxt(file_path, delimiter=',', dtype=float, names=True)
        var_names = data_array.dtype.names

        if not var_names:
            raise ValueError(f'No variable names were loaded from {file_path}')

        for var_name in var_names:
            if hasattr(self, var_name):
                getattr(self, var_name).values = data_array[var_name]

        if self.rmin.values is not None:
            self.set_radius_values()

    def _get_csv_save_path(self, save_type, scan_factor=None, rho_value=None):
        '''
        Gets the path where a CSV of variable data will be saved, based on the
        input parameter values

        Parameters:
        * save_type (SaveType): The SaveType of the data being saved
        * scan_factor (str | float): The scan_factor, if doing a parameter scan (optional)
        * rho_value (str | float): The rho value of the CSV to use (optional)

        Raises:
        * FileNotFoundError: If the file corresponding to the rho value cannot be found
        '''

        runid = self.options.runid
        scan_num = self.options.scan_num
        var_to_scan = self.options.var_to_scan

        if rho_value is not None:
            rho_str = rho_value if isinstance(rho_value, str) else f'{rho_value:{constants.RHO_VALUE_FMT}}'
            dir_path = utils.get_rho_path(runid, scan_num, var_to_scan)
            file_path = (f'{dir_path}\\{save_type.name.capitalize()} '
                         f'rho{constants.RHO_VALUE_SEPARATOR}{rho_str}.csv')
            if not utils.check_exists(file_path):
                raise FileNotFoundError(
                    f'Rho file not found for value {rho_str}\n'
                    f'Use utils.get_closest_rho function to find the correct rho value to load'
                )

            if scan_factor:
                _log.warning(f'\n\tThe scan_factor input parameter is not used when rho_value is specified')

        elif scan_factor is not None:
            scan_factor_str = (
                scan_factor if isinstance(scan_factor, str) else f'{scan_factor:{constants.SCAN_FACTOR_FMT}}'
            )
            dir_path = utils.get_var_to_scan_path(runid, scan_num, var_to_scan)
            file_path = (f'{dir_path}\\{save_type.name.capitalize()} {var_to_scan}'
                         f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}{scan_factor_str}.csv')

        else:
            dir_path = utils.get_scan_num_path(runid, scan_num)
            file_path = f'{dir_path}\\{runid} {save_type.name.capitalize()} Profiles.csv'

        return dir_path, file_path


class InputVariables(Variables):
    '''
    Input variables are defined as anything that isn't read as output from the
    MMM driver

    All members are defined using the Variable class.  See the Variable class
    definition for more info.  Please refer to the documentation provided
    with MMM for more information about the variables that are used as input
    to MMM.
    '''

    def __init__(self, options=None):

        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME', label=r'time', units='s')
        self.x = Variable('X', cdfvar='X', label=r'$x$')
        self.xb = Variable('XB', cdfvar='XB', label=r'$x_\mathrm{B}$')

        # CDF Variables needed for calculations
        self.aimp = Variable('Mean Mass of Impurities', cdfvar='AIMP', label=r'$\overline{M}_\mathrm{imp}$',
                             save_type=SaveType.INPUT, minvalue=1, smooth=1)
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT')
        self.bz = Variable('BZ', cdfvar='BZ')
        self.bzxr = Variable('BZXR', cdfvar='BZXR')
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
        self.omega = Variable('Toroidal Frequency', cdfvar='OMEGDATA', label=r'$\omega$',
                              minvalue=1e-6, save_type=SaveType.ADDITIONAL, units='1/s')
        self.q = Variable('Safety Factor', cdfvar='Q', label=r'$q$', minvalue=1e-6, smooth=1,
                          save_type=SaveType.INPUT)
        self.rho = Variable('Normalized Radius', label=r'$\rho$')
        self.rhochi = Variable('Radius', label=r'$\rho_\chi$')
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP', label=r'$R$',
                             save_type=SaveType.INPUT, units='m', minvalue=0)
        self.rmin = Variable('Minor Radius', cdfvar='RMNMP', label=r'$r$',
                             save_type=SaveType.INPUT, units=r'm', minvalue=0)
        self.rmina = Variable('Minor Radius (normalized)', label=r'$r/a$', units=r'', minvalue=0)
        self.te = Variable('Electron Temperature', cdfvar='TE', label=r'$T_\mathrm{e}$',
                           minvalue=1e-6, smooth=1, save_type=SaveType.INPUT, units='keV')
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI', label=r'$T_\mathrm{i}$',
                           minvalue=1e-6, smooth=1, save_type=SaveType.INPUT, units='keV')
        self.vpol = Variable('Poloidal Velocity', cdfvar='VPOLX_NC', label=r'$v_\theta$',
                             absminvalue=1e-6, smooth=3, save_type=SaveType.INPUT, units='m/s')
        self.vtor = Variable('Toroidal Velocity', cdfvar='VTOR_AVG', label=r'$v_\phi$',
                             absminvalue=1e-6, smooth=3, save_type=SaveType.INPUT, units='m/s')
        self.wexbs = Variable(r'ExB Shear Rate', label=r'$\omega_{E \times B}$', smooth=3,
                              save_type=SaveType.INPUT, units='s^{-1}', minvalue=1e-6)
        self.zimp = Variable('Mean Charge of Impurities', cdfvar='XZIMP', label=r'$\overline{Z}_\mathrm{imp}$',
                             smooth=1, save_type=SaveType.INPUT)

        # wexbs variables may show units of rad/s in TRANSP, but they are all actually in 1/s
        self.wexbsa = Variable(r'ExB Shear Rate', cdfvar='SREXBA', label=r'$\omega_{E \times B}$',
                               smooth=1, units='s^{-1}', minvalue=1e-6)
        self.wexbsmod = Variable(r'ExB Shear Rate', cdfvar='SREXBMOD', label=r'$\omega_{E \times B}$',
                                 smooth=1, units='s^{-1}', minvalue=1e-6)
        self.wexbsv2 = Variable(r'ExB Shear Rate', cdfvar='SREXBV2', label=r'$\omega_{E \times B}$',
                                smooth=1, units='s^{-1}', minvalue=1e-6)

        self.vtorx = Variable('Toroidal Velocity (Imp)', cdfvar='VTORX_NC', label=r'$v_\phi$',
                              absminvalue=1e-6, smooth=3, units='m/s')
        self.vtord = Variable('Toroidal Velocity (D+)', cdfvar='VTORD_NC', label=r'$v_\phi$',
                              absminvalue=1e-6, smooth=3, units='m/s')
        self.vtorh = Variable('Toroidal Velocity (H+)', cdfvar='VTORH_NC', label=r'$v_\phi$',
                              absminvalue=1e-6, smooth=3, units='m/s')
        self.vtoravg = Variable('Toroidal Velocity (avg)', cdfvar='VTOR_AVG', label=r'$v_\phi$',
                                absminvalue=1e-6, smooth=3, units='m/s')

        # Additional CDF variables
        self.betat = Variable('BETAT', cdfvar='BETAT')
        self.tepro = Variable('Electron Temperature', cdfvar='TEPRO')
        self.tipro = Variable('Thermal Ion Temperature', cdfvar='TIPRO')

        # Calculated Variables (some are also in the CDF)
        self.aimass = Variable('Mean Mass of Thermal Ions', label=r'$\overline{M}_\mathrm{i}$',
                               save_type=SaveType.INPUT, minvalue=1)
        self.ahyd = Variable('Mean Mass of Hydrogenic Ions', label=r'$\overline{M}_\mathrm{h}$',
                             save_type=SaveType.INPUT, minvalue=1)
        self.alphamhd = Variable('Alpha MHD', label=r'$\alpha_\mathrm{MHD}$',
                                 save_type=SaveType.ADDITIONAL)
        self.alphamhdunit = Variable('Alpha MHD', label=r'$\alpha_\mathrm{MHD,u}$',
                                     save_type=SaveType.ADDITIONAL)
        self.beta = Variable('Pressure Ratio', cdfvar='BTPL', label=r'$\beta$',
                             save_type=SaveType.ADDITIONAL, minvalue=0)
        self.betae = Variable('Electron Pressure Ratio', cdfvar='BTE', label=r'$\beta_\mathrm{\,e}$',
                              save_type=SaveType.ADDITIONAL, minvalue=0)
        self.betaeunit = Variable('Electron Pressure Ratio', label=r'$\beta_\mathrm{\,e,u}$',
                                  save_type=SaveType.ADDITIONAL, minvalue=0)
        self.betaepunit = Variable('Electron Beta Prime', label=r'$\beta^\prime_\mathrm{\,e,u}$',
                                   save_type=SaveType.ADDITIONAL)
        self.bftor = Variable('Toroidal Magnetic Flux', cdfvar='TRFLX', label=r'$\Psi_\mathrm{T}$',
                              save_type=SaveType.ADDITIONAL, minvalue=0)
        self.bpol = Variable('Poloidal Magnetic Field', cdfvar='BPOL', label=r'$B_\theta$',
                             save_type=SaveType.ADDITIONAL, units='T', minvalue=1e-6)
        self.btor = Variable('Toroidal Magnetic Field', cdfvar='', label=r'$B_\phi$',
                             save_type=SaveType.INPUT, units='T', absminvalue=1e-6)
        self.bunit = Variable('Magnetic Field', label=r'$B_\mathrm{u}$',
                              save_type=SaveType.INPUT,
                              units='T', absminvalue=1e-6)
        self.bunit_btor = Variable('Magnetic Field Ratio', label=r'$B_\mathrm{u} / B_\mathrm{\phi}$',
                                   save_type=SaveType.ADDITIONAL, units='', absminvalue=1e-6)
        self.csound = Variable('Sound Speed', label=r'$c_s$',
                               save_type=SaveType.ADDITIONAL, units='m/s')
        self.csound_a = Variable('Sound Frequency', label=r'$c_s / a$',
                                 save_type=SaveType.ADDITIONAL, units=r'$s^{-1}$')
        self.e_r_grp = Variable('Radial Electric Field (p)', cdfvar='ERPRESS', label=r'$E_\mathrm{r,p}$',
                                smooth=0)
        self.e_r_phi = Variable('Radial Electric Field (tor)', cdfvar='ERVTOR', label=r'$E_\mathrm{r,\phi}$',
                                smooth=0)
        self.e_r_tht = Variable('Radial Electric Field (pol)', cdfvar='ERVPOL', label=r'$E_\mathrm{r,\theta}$',
                                smooth=0)
        self.eps = Variable('Inverse Aspect Ratio', label=r'$\epsilon$',
                            save_type=SaveType.ADDITIONAL)
        self.epsilonne = Variable('2 gBu / gne', label=r'$\epsilon_\mathrm{ne,u}$',
                                  save_type=SaveType.ADDITIONAL)
        self.etae = Variable('Electron Gradient Ratio', cdfvar='ETAE', label=r'$\eta_\mathrm{\,e}$',
                             save_type=SaveType.ADDITIONAL)
        self.etai = Variable('Ion Gradient Ratio', label=r'$\eta_\mathrm{\,i}$')  # cdfvar='ETAI' is not gti/gni
        self.gmax = Variable('Max Gradient', label=r'$g_\mathrm{max}$',
                             save_type=SaveType.ADDITIONAL)
        self.gmaxunit = Variable('Max Gradient', label=r'$g_\mathrm{max,u}$',
                                 save_type=SaveType.ADDITIONAL)
        self.gne_threshold = Variable(r'Growth Rate Threshold', label=r'$g_\mathrm{ne}$')
        self.gte_threshold = Variable(r'Growth Rate Threshold', label=r'$g_\mathrm{Te}$')
        self.gxi = Variable('Flux Surface Vol. Avg.', cdfvar='GXI', units='m^-1', label=r'$\nabla \hat{\rho}$',
                            save_type=SaveType.INPUT)
        self.gyrfe = Variable('Electron Gyrofrequency', label=r'$\omega_\mathrm{ce}$',
                              save_type=SaveType.ADDITIONAL, units='s^{-1}')
        self.gyrfeunit = Variable('Electron Gyrofrequency', label=r'$\omega_\mathrm{ce,u}$',
                                  save_type=SaveType.ADDITIONAL, units='s^{-1}')
        self.gyrfi = Variable('Ion Gyrofrequency', label=r'$\omega_\mathrm{ci}$',
                              save_type=SaveType.ADDITIONAL, units='s^{-1}')
        self.gyrfiunit = Variable('Ion Gyrofrequency', label=r'$\omega_\mathrm{ci,u}$',
                                  save_type=SaveType.ADDITIONAL, units='s^{-1}')
        self.lare = Variable('Electron Larmor Radius', label=r'$\rho_\mathrm{e}$',
                             save_type=SaveType.ADDITIONAL, units='m')
        self.lareunit = Variable('Electron Larmor Radius', label=r'$\rho_\mathrm{e,u}$',
                                 save_type=SaveType.ADDITIONAL, units='m')
        self.loge = Variable('Electron Coulomb Logarithm', cdfvar='CLOGE', label=r'$\lambda_\mathrm{e}$',
                             save_type=SaveType.ADDITIONAL)
        self.logi = Variable('Ion Coulomb Logarithm', cdfvar='CLOGI', label=r'$\lambda_\mathrm{i}$')
        self.ni = Variable('Thermal Ion Density', cdfvar='NI', label=r'$n_\mathrm{i}$', units='m^-3',
                           save_type=SaveType.ADDITIONAL, minvalue=1e-6)
        self.ni2 = Variable('Thermal Ion Density', label=r'$n_\mathrm{i}$', units='m^-3', minvalue=1e-6)
        self.nh0 = Variable('Hydrogen Ion Density', cdfvar='NH', label=r'$n_\mathrm{h0}$', units='m^-3',
                            save_type=SaveType.ADDITIONAL)
        self.nh = Variable('Total Hydrogenic Ion Density', label=r'$n_\mathrm{h}$',
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
        self.rhosunit = Variable('Ion Larmor Radius', units='m', label=r'$\rho_\mathrm{s,u}$',
                                 save_type=SaveType.ADDITIONAL, minvalue=1e-6)
        self.shat = Variable('Effective Magnetic Shear', label=r'$\hat{s}_{\kappa}$',
                             save_type=SaveType.ADDITIONAL)
        self.shat_gxi = Variable('Effective Magnetic Shear', label=r'$\hat{s}_{\nabla \rho}}$',
                                 save_type=SaveType.ADDITIONAL)
        self.shear = Variable('Magnetic Shear', cdfvar='SHAT', label=r'$s$',
                              save_type=SaveType.ADDITIONAL)
        self.tau = Variable('Temperature Ratio', label=r'$\tau$',
                            save_type=SaveType.ADDITIONAL, minvalue=0)
        self.vpar = Variable('Parallel Velocity', label=r'$v_\parallel$', absminvalue=1e-6,
                             save_type=SaveType.INPUT, units='m/s')
        self.vthe = Variable('Electron Thermal Velocity', label=r'$v_{\mathrm{Te}}$',
                             save_type=SaveType.ADDITIONAL, units='m/s')
        self.vthi = Variable('Ion Thermal Velocity', label=r'$v_{\mathrm{Ti}}$',
                             save_type=SaveType.ADDITIONAL, units='m/s')
        self.zeff = Variable('Effective Charge', cdfvar='ZEFFP', label=r'$Z_\mathrm{eff}$',
                             save_type=SaveType.INPUT, minvalue=1)
        self.wbounce = Variable('Bounce Frequency', label=r'$\omega_\mathrm{be}$', units=r'$s^{-1}$',
                                save_type=SaveType.ADDITIONAL)
        self.wtransit = Variable('Transit Frequency', label=r'$\omega_\mathrm{te}$',
                                 save_type=SaveType.ADDITIONAL, units=r'$s^{-1}$')
        self.xetgm_const = Variable('Diffusivity Constant',
                                    label=r'${\rho_\mathrm{e,u}}^2 v_\mathrm{Te} / L_\mathrm{Te}$',
                                    save_type=SaveType.ADDITIONAL)

        # Calculated Gradients
        self.gbtor = Variable('Btor Gradient', label=r'$g_{\mathrm{B\phi}}$')
        self.gbunit = Variable('Bunit Gradient', label=r'$g_{\mathrm{Bu}}$',
                               save_type=SaveType.INPUT)
        self.gne = Variable('Electron Density Gradient', label=r'$g_{\mathrm{ne}}$',
                            save_type=SaveType.INPUT)
        self.gnh = Variable('Hydrogenic Ion Density Gradient', label=r'$g_{\mathrm{nh}}$',
                            save_type=SaveType.INPUT)
        self.gni = Variable('Thermal Ion Density Gradient', label=r'$g_{\mathrm{ni}}$',
                            save_type=SaveType.INPUT)
        self.gnz = Variable('Impurity Density Gradient', label=r'$g_{\mathrm{nz}}$',
                            save_type=SaveType.INPUT)
        self.gp_i = Variable('Ion Pressure Gradient', label=r'$g_{\mathrm{pi}}$')
        self.gq = Variable('Safety Factor Gradient', label=r'$g_{q}$',
                           save_type=SaveType.INPUT)
        self.gte = Variable('Electron Temperature Gradient', label=r'$g_{\mathrm{Te}}$',
                            save_type=SaveType.INPUT)
        self.gti = Variable('Thermal Ion Temperature Gradient', label=r'$g_{\mathrm{Ti}}$',
                            save_type=SaveType.INPUT)
        self.gvpar = Variable('Parallel Velocity Gradient', label=r'$g_{v\parallel}$',
                              save_type=SaveType.INPUT)
        self.gvpol = Variable('Poloidal Velocity Gradient', label=r'$g_{v\theta}$',
                              save_type=SaveType.INPUT)
        self.gvtor = Variable('Toroidal Velocity Gradient', label=r'$g_{v\phi}$',
                              save_type=SaveType.INPUT)

        super().__init__(options)  # Init parent class

    def set_x_values(self):
        '''
        Sets variable x from variable xb

        x is the grid between xb, and has one fewer point than xb.
        '''
        self.x.values = (self.xb.values[0:-1, :] + self.xb.values[1:, :]) / 2

    def get_vars_of_type(self, save_type):
        '''Returns (list of str): List of all variables with the specified save_type'''
        nonzero_variables = self.get_variables()
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

    def save_vars_of_type(self, save_type, scan_factor=None):
        '''
        Saves variable values of the specified save type to a CSV

        Parameters:
        * save_type (SaveType): The save type of the variables to save
        * time_idx (int): The index of the measurement time
        * runid (str): The runid of the CSV to use
        * scan_num (int): The scan number of the CSV to use
        * var_to_scan (str): The scanned variable of the CSV to use (optional)
        * scan_factor (scan_factor): The value of the scan factor (optional)
        '''

        # Put rmin at the front of the variable list
        var_list = self.get_vars_of_type(save_type)
        if 'rmin' not in var_list:
            var_list.insert(0, 'rmin')
        else:
            var_list.insert(0, var_list.pop(var_list.index('rmin')))

        data, header = self._get_data_as_array(var_list)
        self._save_to_csv(data, header, save_type, scan_factor)

    def save(self, scan_factor=None):
        '''
        Saves variable values of all relevant save types to a CSV

        Parameters:
        * time_idx (int): The index of the measurement time
        * runid (str): The runid of the CSV to use
        * scan_num (int): The scan number of the CSV to use
        * var_to_scan (str): The scanned variable of the CSV to use (optional)
        * scan_factor (scan_factor): The value of the scan factor (optional)
        '''

        self.save_vars_of_type(SaveType.INPUT, scan_factor)
        self.save_vars_of_type(SaveType.ADDITIONAL, scan_factor)


class OutputVariables(Variables):
    '''
    Output variables consist of all variable data obtained as output from MMM
    (as well as rho and rmin)

    Please refer to the documentation provided with MMM for more information
    about the variables that obtained as output from MMM.
    '''

    def __init__(self, options=None):
        # Independent Variables
        self.rho = Variable('rho', units='', label=r'$\rho$')
        self.rmin = Variable('Minor Radius', units='m', label=r'$r$', minvalue=0)
        self.rmina = Variable('rmina', label=r'$r/a$', units=r'', minvalue=0)
        # Total Fluxes
        self.fti = Variable('fti', units='keVm/s', label=r'$\Gamma_\mathrm{Ti}$')
        self.fdi = Variable('fdi', units='m^{-2}s^{-1}^', label=r'$\Gamma_\mathrm{Di}$')
        self.fte = Variable('fte', units='keVm/s', label=r'$\Gamma_\mathrm{Te}$')
        self.fdz = Variable('fdz', units='m^{-2}s^{-1}', label=r'$\Gamma_\mathrm{Dz}$')
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
        self.gmaW20ii = Variable('gmaW20ii', units='s^{-1}', label='gmaW20ii')
        self.omgW20ii = Variable('omgW20ii', units='s^{-1}', label='omgW20ii')
        self.gmaW20ie = Variable('gmaW20ie', units='s^{-1}', label='gmaW20ie')
        self.omgW20ie = Variable('omgW20ie', units='s^{-1}', label='omgW20ie')
        self.gmaW20ei = Variable('gmaW20ei', units='s^{-1}', label='gmaW20ei')
        self.omgW20ei = Variable('omgW20ei', units='s^{-1}', label='omgW20ei')
        self.gmaW20ee = Variable('gmaW20ee', units='s^{-1}', label='gmaW20ee')
        self.omgW20ee = Variable('omgW20ee', units='s^{-1}', label='omgW20ee')
        # DBM Components
        self.xtiDBM = Variable('xtiDBM', units='m^2/s', label='xtiDBM')
        self.xdiDBM = Variable('xdiDBM', units='m^2/s', label='xdiDBM')
        self.xteDBM = Variable('xteDBM', units='m^2/s', label='xteDBM')
        self.gmaDBM = Variable('gmaDBM', units='s^{-1}', label='gmaDBM')
        self.omgDBM = Variable('omgDBM', units='s^{-1}', label='omgDBM')
        # ETG Component
        self.xteETG = Variable('xteETG', units='m^2/s', label=r'$\chi_\mathrm{e, etg}$')
        self.gtecritETG = Variable(r'Critical $g_\mathrm{Te}$ (Jenko ETG)', units='',
                                   label=r'$g_\mathrm{Te, etg}^\mathrm{crit}$')
        # MTM Components
        self.xteMTM = Variable('xteMTM', units='m^2/s', label=r'$\chi_\mathrm{e, mtm}$')
        self.gmaMTM = Variable('gmaMTM', units='s^{-1}', label=r'$\gamma_\mathrm{mtm}$')
        self.omgMTM = Variable('omgMTM', units='s^{-1}', label=r'$\omega_\mathrm{mtm}$')
        self.kyrhosMTM = Variable('kyrhosMTM', units='', label=r'$k_y\rho_\mathrm{s}$')
        self.dbsqprf = Variable('dbsqprf', units='', label=r'$|\delta B/B|^2$')
        # ETGM Components
        self.xteETGM = Variable('Thermal Diffusivity', units='m^2/s', label=r'$\chi_\mathrm{e}$')
        self.xte2ETGM = Variable('Thermal Diffusivity', units='m^2/s', label=r'$\chi^{\ast}_\mathrm{e}$')
        self.gmaETGM = Variable('Growth Rate', units='s^{-1}', label=r'$\gamma_\mathrm{}$')
        self.omgETGM = Variable('Frequency', units='s^{-1}', label=r'$\omega_\mathrm{}$')
        self.kyrhoeETGM = Variable('Wave Number', units='', label=r'$k_y\rho_\mathrm{e}$')
        self.kyrhosETGM = Variable('Wave Number', units='', label=r'$k_y\rho_\mathrm{s}$')
        self.gaveETGM = Variable('Magnetic Field Curvature', units='', label=r'$\overline{G}$')
        self.alphaETGM = Variable('alphaMHD', units='', label=r'$\alpha_\mathrm{MHD}$')
        self.kpara2ETGM = Variable('kpara2ETGM', units='', label=r'$\langle k^{2}_\parallel \rangle$')
        self.fleETGM = Variable('fleETGM', units='', label=r'$\langle k^{2}_\perp \rho^{2}_\mathrm{e}\rangle$')
        self.phi2ETGM = Variable('Electrostatic Potential', units='', label=r'$|\hat{\phi}|^2$')
        self.Apara2ETGM = Variable('Electromagnetic Potential', units='', label=r'$|\hat{A}_{\!\parallel}\!|^2$')
        self.omegadETGM = Variable('omegadETGM', units='s^{-1}', label=r'$\omega_\mathrm{De}$')
        self.omegadiffETGM = Variable('omegadiffETGM', units='s^{-1}', label=r'$\omega - \omega_\mathrm{De}$')
        self.gammadiffETGM = Variable('gammadiffETGM', units='s^{-1}', label=r'$\gamma - \omega_\mathrm{De}$')
        self.omegasETGM = Variable('omegasETGM', units='s^{-1}', label=r'$\omega_{*\mathrm{e}}$')
        self.omegateETGM = Variable('omegateETGM', units='s^{-1}', label=r'$\omega_{\mathrm{Te}}$')
        self.omegasetaETGM = Variable('omegasetaETGM', units='s^{-1}',
                                      label=r'$\omega_{*\mathrm{e}} (1 + \eta_\mathrm{e})$')
        self.walfvenunit = Variable('Alfven Frequency', units='s^{-1}', label=r'$\omega_\mathrm{A}$')
        self.satETGM = Variable('Saturation Ratio', units='', label=r'$2\hat{\gamma}/|\hat{\phi}| R k_\mathrm{x}$')

        super().__init__(options)  # Init parent class

    def get_all_output_vars(self):
        '''Returns (list of str): all output variable names (other than rho and rmin)'''
        all_vars = self.get_variables()
        all_vars.remove('rho')
        all_vars.remove('rmin')
        all_vars.remove('rmina')
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

    def save(self, scan_factor=None):
        '''Saves output variables to a CSV (other than rho)'''

        # Put rmin at the front of the variable list
        var_list = self.get_variables()
        var_list.insert(0, var_list.pop(var_list.index('rmin')))
        var_list.remove('rho')

        data, header = self._get_data_as_array(var_list)
        self._save_to_csv(data, header, SaveType.OUTPUT, scan_factor)


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

        self.units = units  # Call units setter to also set units_label

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
        if units in _UNITS_TO_UNITS_LABEL.keys():
            self._units_label = _UNITS_TO_UNITS_LABEL[units]  # Set units_label in LaTeX format

    @property
    def dimensions(self):
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions):
        if not isinstance(dimensions, list):
            raise ValueError(f'Variable dimensions must be {list} and not {type(dimensions)}')
        self._dimensions = dimensions

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        if not isinstance(values, np.ndarray):
            raise ValueError(f'Variable values must be {np.ndarray} and not {type(values)}')
        self._values = values

    def set(self, **kwargs):
        '''Sets members using keyword arguments'''
        for key, value in kwargs.items():
            setattr(self, key, value)

    def apply_smoothing(self):
        '''
        Variable smoothing using a Gaussian filter

        The value of sigma needs to increase nearly linearly as the amount of
        input points increases, to maintain the same level of smoothing.  To
        achieve this, the self.smooth value is multiplied by the number of
        data points it has. Additionally, the amount of smoothing applied
        also depends on how closely grouped values are along the x-axis.
        Specifically, for a given value of sigma, loosely spaced values will
        be smoothed more than tightly clustered values.

        For this reason, using the uniform rho input option will result in
        uniform smoothing across a variable (since the values in rmin become
        more clustered as rmin approaches its maximum value).  If the uniform
        rho input option is not used, then variables will have stronger
        smoothing applied near rmin = 0, and weaker smoothing applied near
        rmin = 1.
        '''

        if self.smooth is not None:
            sigma = int(self.values.shape[0] * self.smooth / 100)
            if self.values.ndim == 2:
                self.values = scipy.ndimage.gaussian_filter(self.values, sigma=(sigma, 0))
            else:
                self.values = scipy.ndimage.gaussian_filter(self.values, sigma=sigma)

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

        Parameters:
        * raise_exception (bool): Possible exceptions will not be raised when false

        Raises:
        * ValueError: If multiple nonphysical values are found
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

    def clamp_values(self, clamp_value):
        '''Clamps values between -clamp_value and +clamp_value'''
        self.values[self.values > clamp_value] = clamp_value
        self.values[self.values < -clamp_value] = -clamp_value

    def set_origin_to_zero(self):
        '''
        Sets origin values to approximately zero

        Original values (rmin = 0) are multiplied by 1e-6 times the lowest
        absolute value (along the entire radius) at each time slice.  By
        doing this, we avoid potential division by zero errors when
        performing calculations.
        '''
        self.values[0, :] = 1e-6 * np.absolute(self.values).min(axis=0)

    def update_label(self, before, after):
        """
        Updates the label of the variable

        Parameters:
        * before (str): Inserted at the front the label
        * after (str): Appended to the end of the label
        """

        label_stripped = self.label.strip('$')
        self.label = f'${before}{label_stripped}{after}$'

    def check_for_nan(self):
        '''Checks for nan values and raises a ValueError if any are found'''
        if np.isnan(self.values).any():
            raise ValueError(f'nan values found in variable {self.name}')


# For testing purposes
if __name__ == '__main__':
    import modules.options
    options = modules.options.Options(runid='TEST', scan_num=25)

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
    # ivars.save(options)
    # ovars.save(options)

    ivars = InputVariables()
    # ivars.load_from_csv(SaveType.INPUT, scan_factor=1.5)
    # ivars.load_from_csv(SaveType.ADDITIONAL, scan_factor=1.5, rho_value=0.5)

    ivars.print_nonzero_variables()

    print(ivars.get_variables())
