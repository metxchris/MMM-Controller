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
import settings
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
    'm^2/s^2': r'm$^2$/s$^2$',
    's^{-1}': r's$^{-1}$',
    's^-1': r's$^{-1}$',
    'm^-1': r'm$^{-1}$',
    'm^-2': r'm$^{-2}$',
    'm^2': r'm$^2$',
    'MA/m^2': r'MA/m$^2$',
    'keVm/s': r'keV$\,$m/s',
    'm^-2s^-1': r'm$^{-2}$s$^{-1}$',
    'ohm*m': r'$\Omega\,$m',
    'A/m^2': r'A/m$^2$'
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
            if isinstance(getattr(self, v).values, np.ndarray):
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
        * TypeError: If the variable values is not an instance of np.ndarray
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
                    if not isinstance(var_values, np.ndarray):
                        raise TypeError(f'Variable {var_name} is not an instance of np.ndarray')
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
        if not utils.check_exists(file_path):
            if save_type != SaveType.ADDITIONAL:
                raise FileNotFoundError(
                    f'File not found: {file_path}\n'
                )

            file_path = None

        if file_path:
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
                if save_type != SaveType.ADDITIONAL:
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

        # EPM Variables
        self.gtf = Variable(
            'Fast Ion Temperature Gradient',
            label=r'$g_{\mathrm{Tf}}$',
        )

        self.gnf = Variable(
            'Fast Ion Density Gradient',
            label=r'$g_{\mathrm{nf}}$',
        )

        self.tf = Variable(
            'Fast Ion Temperature',
            label=r'$T_\mathrm{f}$',
            minvalue=1e-6,
            smooth=1,
            units='keV',
        )

        if settings.USE_EPM:
            self.gtf.save_type = SaveType.INPUT
            self.gnf.save_type = SaveType.INPUT
            self.tf.save_type = SaveType.INPUT

        # Fundamental CDF Variables
        self.time = Variable(
            'Time',
            cdfvar='TIME',
            label=r'time',
            units='s',
        )

        self.x = Variable(
            'X',
            cdfvar='X',
            label=r'$x$',
        )

        self.xb = Variable(
            'XB',
            cdfvar='XB',
            label=r'$x_\mathrm{B}$',
        )

        # CDF Variables needed for calculations
        self.az = Variable(
            'Mean Mass of Impurities',
            cdfvar='AIMP',
            label=r'$\overline{M}_\mathrm{imp}$',
            minvalue=1,
            smooth=1,
            save_type=SaveType.INPUT,
        )

        self.arat = Variable(
            'Aspect Ratio',
            cdfvar='ARAT',
        )

        self.bzxr = Variable(
            'BZXR',
            cdfvar='BZXR',
        )

        self.dbp = Variable(
            'dbp',
            label=r'$dB_\theta/d\rho$',
            save_type=SaveType.INPUT,
        )

        self.d2bp = Variable(
            'd2bp',
            label=r'$d^2B_\theta/d\rho^2$',
            save_type=SaveType.INPUT,
        )

        self.elong = Variable(
            'Elongation',
            cdfvar='ELONG',
            label=r'$\kappa$',
            smooth=0,
            save_type=SaveType.INPUT,
        )

        self.etanc = Variable(
            'NC Resistivity',
            cdfvar='ETA_NC',
            label=r'$\eta_{\mathrm{NC}}$',
            units='mohm',
            smooth=0,
            minvalue=0,
            save_type=SaveType.INPUT,
        )

        self.ne = Variable(
            'Electron Density',
            cdfvar='NE',
            label=r'$n_\mathrm{e}$',
            minvalue=1e-32,
            smooth=1,
            units='m^-3',
            default_values=1e-16,
            save_type=SaveType.INPUT,
        )

        self.nf = Variable(
            'Fast Ion Density',
            cdfvar='BDENS',
            label=r'$n_\mathrm{f}$',
            minvalue=1e-32,
            smooth=1,
            units='m^-3',
            default_values=1e-16,
            save_type=SaveType.INPUT,
        )

        self.nd = Variable(
            'Deuterium Ion Density',
            cdfvar='ND',
            label=r'$n_{\rm d}$',
            minvalue=1e-32,
            smooth=1,
            units='m^-3',
            default_values=1e-16,
            save_type=SaveType.ADDITIONAL,
        )

        self.nz = Variable(
            'Impurity Density',
            cdfvar='NIMP',
            label=r'$n_z$',
            minvalue=1e-32,
            smooth=1,
            units='m^-3',
            default_values=1e-16,
            save_type=SaveType.INPUT,
        )

        self.omega = Variable(
            'Toroidal Frequency',
            cdfvar='OMEGA',
            label=r'$\omega_\phi$',
            units='1/s',  # CDF may list units as rad/s, but assume 1/s (CDF is wrong)
        )

        self.q = Variable(
            'Safety Factor',
            cdfvar='Q',
            label=r'$q$',
            minvalue=1e-6,
            smooth=0,
            save_type=SaveType.INPUT,
        )

        self.rho = Variable(  # Rho is saved in every CSV
            'Normalized Radius',
            label=r'$\hat{\rho}$'
        )

        self.rhochi = Variable(
            'Radius',
            label=r'$\rho_\chi$'
        )

        self.rmaj = Variable(
            'Major Radius',
            cdfvar='RMJMP',
            label=r'$R$',
            units='m',
            minvalue=0,
            save_type=SaveType.INPUT,
        )

        self.rmin = Variable(
            'Minor Radius',
            cdfvar='RMNMP',
            label=r'$r$',
            units=r'm',
            minvalue=0,
            save_type=SaveType.INPUT,
        )

        self.rmina = Variable(
            'Minor Radius (normalized)',
            label=r'$r/a$',
            units=r'',
            minvalue=0,
        )

        self.te = Variable(
            'Electron Temperature',
            cdfvar='TE',
            label=r'$T_\mathrm{e}$',
            minvalue=1e-6,
            smooth=1,
            units='keV',
            save_type=SaveType.INPUT,
        )

        self.tepro = Variable(  # Data copied to te when experimental profiles are used
            'Electron Temperature',
            cdfvar='TEPRO',
            label=r'$T_\mathrm{e,exp}$',
            minvalue=1e-6,
            smooth=0,
            units='keV',
        )

        self.ti = Variable(
            'Thermal Ion Temperature',
            cdfvar='TI',
            label=r'$T_\mathrm{i}$',
            minvalue=1e-6,
            smooth=1,
            units='keV',
            save_type=SaveType.INPUT,
        )

        self.tipro = Variable(  # Data copied to ti when experimental profiles are used
            'Thermal Ion Temperature',
            cdfvar='TIPRO',
            label=r'$T_\mathrm{i,exp}$',
            minvalue=1e-6,
            smooth=0,
            units='keV',
        )

        self.vpol = Variable(
            'Poloidal Velocity',
            cdfvar='VPOLX_NC',
            label=r'$v_\theta$',
            absminvalue=1e-6,
            smooth=3,
            units='m/s',
            save_type=SaveType.INPUT,
        )

        self.vtor = Variable(
            'Toroidal Velocity',
            cdfvar='VTOR_AVG',  # CDF variable only for comparison with MMM definition
            label=r'$v_\phi$',
            absminvalue=1e-6,
            smooth=3,
            units='m/s',
            save_type=SaveType.INPUT,
        )

        self.wexb = Variable(
            'ExB Shear Rate',
            label=r'$\omega_{E \times B}$',
            smooth=3,
            units='s^{-1}',
            minvalue=1e-6,
            save_type=SaveType.INPUT,
        )

        self.zz = Variable(
            'Mean Charge of Impurities',
            cdfvar='XZIMP',
            label=r'$\overline{Z}_\mathrm{imp}$',
            smooth=1,
            save_type=SaveType.INPUT,
        )

        self.ebeam = Variable(
            'Fast Ion Beam Energy',
            cdfvar='EBEAM_D',
            label=r'EBEAM$\_$D',
            minvalue=1e-4,
            units='keV',
        )

        self.pmhdf = Variable(
            'pmhdf',
            cdfvar='PMHDF_IN',
            label=r'PMHDF$\_$IN',
            units='Pa',
        )

        self.tmhdf = Variable(
            'tmhdf',
            label=r'$T_{\rm f}$',
            units='keV',
        )

        self.btbe = Variable(
            'btbe',
            cdfvar='BTBE',
            label=r'BTBE',
            units='',
        )

        self.tbtbe = Variable(
            'tbtbe',
            label=r'$T_{\rm f, BTBE}$',
            units='keV',
        )

        self.wexbsa = Variable(
            'ExB Shear Rate',
            cdfvar='SREXBA',
            label=r'$\omega_{E \times B}$',
            smooth=1,
            units='s^{-1}',  # CDF may list units as rad/s, but assume 1/s (CDF is wrong)
            minvalue=1e-6,
        )

        self.wexbsmod = Variable(
            'ExB Shear Rate',
            cdfvar='SREXBMOD',
            label=r'$\omega_{E \times B}$',
            smooth=1,
            units='s^{-1}',  # CDF may list units as rad/s, but assume 1/s (CDF is wrong)
            minvalue=1e-6,
        )

        self.wexbsv2 = Variable(
            'ExB Shear Rate',
            cdfvar='SREXBV2',
            label=r'$\omega_{E \times B}$',
            smooth=1,
            units='s^{-1}',  # CDF may list units as rad/s, but assume 1/s (CDF is wrong)
            minvalue=1e-6,
        )

        self.ai = Variable(
            'Mean Mass of Thermal Ions',
            label=r'$\overline{M}_\mathrm{i}$',
            units='u',
            minvalue=1,
            save_type=SaveType.INPUT,
        )

        self.ah = Variable(  # Hydrogen + Deuterium + Tritium
            'Mean Mass of Hydrogenic Ions',
            label=r'$\overline{M}_\mathrm{h}$',
            units='u',
            minvalue=1,
            save_type=SaveType.INPUT,
        )

        self.alphamhd = Variable(
            'Alpha MHD',
            label=r'$\alpha_\mathrm{MHD}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.alphamhdu = Variable(
            'Alpha MHD',
            label=r'$\alpha_\mathrm{MHD,u}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.beta = Variable(
            'Pressure Ratio',
            cdfvar='BTPL',
            label=r'$\beta$',
            minvalue=0,
        )

        self.betae = Variable(
            'Electron Pressure Ratio',
            cdfvar='BTE',
            label=r'$\beta_\mathrm{\,e}$',
            minvalue=0,
            save_type=SaveType.ADDITIONAL,
        )

        self.betaeu = Variable(
            'Electron Pressure Ratio',
            label=r'$\beta_\mathrm{\,e,u}$',
            minvalue=0,
            save_type=SaveType.ADDITIONAL,
        )

        self.betapu = Variable(
            'Beta Prime',
            label=r'$\beta^\prime$',
            save_type=SaveType.ADDITIONAL,
        )

        self.bftor = Variable(
            'Toroidal Magnetic Flux',
            cdfvar='TRFLX',
            label=r'$\Psi_\mathrm{T}$',
            minvalue=0,
        )

        self.bpol = Variable(
            'Poloidal Magnetic Field',
            cdfvar='BPOL',
            label=r'$B_\theta$',
            units='T',
            minvalue=1e-32,
            save_type=SaveType.ADDITIONAL,
        )

        self.btor = Variable(
            'Toroidal Magnetic Field',
            label=r'$B_\phi$',
            units='T',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.bu = Variable(
            'Magnetic Field',
            label=r'$B_\mathrm{u}$',
            units='T',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.csound = Variable(
            'Sound Speed',
            label=r'$c_s$',
            units='m/s',
            save_type=SaveType.ADDITIONAL,
        )

        self.csound_a = Variable(
            'Sound Frequency',
            label=r'$c_s / a$',
            units=r'$s^{-1}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.curoh = Variable(
            'OH Current',
            label=r'$I_\mathrm{OH}$',
            units=r'$MA$',
        )

        self.curohrho = Variable(
            'OH Current',
            label=r'$I_\mathrm{OH}$',
            units=r'$MA$',
        )

        self.curdoh = Variable(
            'OH Current Density',
            cdfvar='CUROH',
            label=r'$J_\mathrm{OH}$',
            units=r'MA/m^2',
            default_values=0,
        )

        self.curlh = Variable(
            'LH Current',
            label=r'$I_\mathrm{LH}$',
            units=r'$MA$',
        )

        self.curdlh = Variable(
            'LH Current Density',
            cdfvar='LHCUR',
            label=r'$J_\mathrm{LH}$',
            units=r'MA/m^2',
            default_values=0,
        )

        self.e_r_grp = Variable(
            'Radial Electric Field (p)',
            cdfvar='ERPRESS',  # MMM Calculations don't match CDF variables
            label=r'$E_\mathrm{r,p}$',
            smooth=0,
        )

        self.e_r_phi = Variable(
            'Radial Electric Field (tor)',
            cdfvar='ERVTOR',  # MMM Calculations don't match CDF variables
            label=r'$E_\mathrm{r,\phi}$',
            smooth=0,
        )

        self.e_r_tht = Variable(
            'Radial Electric Field (pol)',
            cdfvar='ERVPOL',  # MMM Calculations don't match CDF variables
            label=r'$E_\mathrm{r,\theta}$',
            smooth=0,
        )

        self.eps = Variable(
            'Inverse Aspect Ratio',
            label=r'$\epsilon$',
            save_type=SaveType.ADDITIONAL,
        )

        self.epsne = Variable(
            '2 gBu / gne',
            label=r'$\epsilon_\mathrm{ne,u}$',
        )

        self.etae = Variable(
            'Electron Gradient Ratio',
            cdfvar='ETAE',
            label=r'$\eta_\mathrm{\,e}$',
            contour_max=20,
            contour_min=-20,
            save_type=SaveType.ADDITIONAL,
        )

        self.etai = Variable(
            'Ion Gradient Ratio',
            cdfvar='ETAI',  # This is not the same definition as our gti/gni
            label=r'$\eta_\mathrm{\,i}$',
        )

        self.gne_threshold = Variable(
            'Growth Rate Threshold',
            label=r'$g_\mathrm{ne}$',
        )

        self.gte_threshold = Variable(
            'Growth Rate Threshold',
            label=r'$g_\mathrm{Te}$',
        )

        self.gxi = Variable(
            'Flux Surface Vol. Avg.',
            cdfvar='GXI',
            units='m^-1',
            label=r'$\nabla \hat{\rho}$',
            # save_type=SaveType.INPUT,
        )

        self.gyrfe = Variable(
            'Electron Gyrofrequency',
            label=r'$\omega_\mathrm{ce}$',
            units='s^{-1}',
        )

        self.gyrfeu = Variable(
            'Electron Gyrofrequency',
            label=r'$\omega_\mathrm{ce,u}$',
            units='s^{-1}',
        )

        self.gyrfi = Variable(
            'Ion Gyrofrequency',
            label=r'$\omega_\mathrm{ci}$',
            units='s^{-1}',
            save_type=SaveType.ADDITIONAL,
        )

        self.gyrfiu = Variable(
            'Ion Gyrofrequency',
            label=r'$\omega_\mathrm{ci,u}$',
            units='s^{-1}',
            save_type=SaveType.ADDITIONAL,
        )

        self.lare = Variable(
            'Electron Larmor Radius',
            label=r'$\rho_\mathrm{e}$',
            units='m',
        )

        self.lareu = Variable(
            'Electron Larmor Radius',
            label=r'$\rho_\mathrm{e,u}$',
            units='m',
        )

        self.loge = Variable(
            'Electron Coulomb Logarithm',
            cdfvar='CLOGE',
            label=r'$\lambda_\mathrm{e}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.logi = Variable(
            'Ion Coulomb Logarithm',
            cdfvar='CLOGI',
            label=r'$\lambda_\mathrm{i}$',
        )

        self.ni = Variable(
            'Thermal Ion Density',
            cdfvar='NI',
            label=r'$n_\mathrm{i}$',
            units='m^-3',
            minvalue=1e-32,
            default_values=0,
            save_type=SaveType.ADDITIONAL,
        )

        self.ni2 = Variable(
            'Thermal Ion Density',
            label=r'$n_\mathrm{i}$',
            units='m^-3',
            minvalue=1e-32,
        )

        self.nh0 = Variable(
            'Hydrogen Ion Density',
            cdfvar='NH',
            label=r'$n_\mathrm{h0}$',
            units='m^-3',
            default_values=0,
            save_type=SaveType.ADDITIONAL,
        )

        self.nh = Variable(
            'Total Hydrogenic Ion Density',
            label=r'$n_\mathrm{h}$',
            units='m^-3',
            minvalue=1e-32,
            default_values=0,
            save_type=SaveType.INPUT,
        )

        self.nuei = Variable(
            'Electron Collision Frequency',
            label=r'$\nu_\mathrm{ei}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.nuste = Variable(
            'Electron Collisionality',
            cdfvar='NUSTE',
            label=r'$\nu^{*}_\mathrm{e}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.p = Variable(
            'Plasma Pressure',
            cdfvar='PPLAS',
            label=r'$p$',
            save_type=SaveType.ADDITIONAL,
            minvalue=1e-100
        )

        # self.pav = Variable('Plasma Vol Avg', cdfvar='PTOWB', label=r'$pav$')
        self.pcur = Variable(
            'Measured Plasma Current',
            cdfvar='PCUR',
            label=r'$I_\mathrm{measured}$',
        )

        self.rhos = Variable(
            'Ion Larmor Radius',
            units='m',
            label=r'$\rho_\mathrm{s}$',
            minvalue=1e-32,
            save_type=SaveType.ADDITIONAL,
        )

        self.rhosu = Variable(
            'Ion Larmor Radius', 
            units='m',
            label=r'$\rho_\mathrm{s,u}$',
            save_type=SaveType.ADDITIONAL,
            minvalue=1e-32,
        )

        self.shat = Variable(
            'Effective Magnetic Shear',
            label=r'$\hat{s}_{\kappa}$',
            contour_max=20,
            contour_min=-20,
            save_type=SaveType.ADDITIONAL,
        )

        self.shat_gxi = Variable(
            'Effective Magnetic Shear',
            label=r'$\hat{s}$',
            contour_max=20,
            contour_min=-20,
            save_type=SaveType.ADDITIONAL,
        )

        self.shear = Variable(
            'Magnetic Shear',
            cdfvar='SHAT',
            label=r'$s$',
            save_type=SaveType.ADDITIONAL,
        )

        self.te_ti = Variable(
            'Temperature Ratio',
            label=r'$T_\mathrm{e}/T_\mathrm{i}$',
            minvalue=0,
        )

        self.ti_te = Variable(
            'Temperature Ratio',
            label=r'$T_\mathrm{i}/T_\mathrm{e}$',
            minvalue=0,
            save_type=SaveType.ADDITIONAL,
        )

        self.tf_te = Variable(
            'Temperature Ratio',
            label=r'$T_\mathrm{f}/T_\mathrm{e}$',
            minvalue=0,
        )

        self.vpar = Variable(
            'Parallel Velocity',
            label=r'$v_\parallel$',
            absminvalue=1e-32,
            units='m/s',
            save_type=SaveType.INPUT,
        )

        self.vthe = Variable(
            'Electron Thermal Velocity',
            label=r'$v_{\mathrm{Te}}$',
            units='m/s',
            save_type=SaveType.ADDITIONAL,
        )

        self.vthi = Variable(
            'Ion Thermal Velocity',
            label=r'$v_{\mathrm{Ti}}$',
            units='m/s',
            save_type=SaveType.ADDITIONAL,
        )

        self.zeff = Variable(
            'Effective Charge',
            cdfvar='ZEFFP',
            label=r'$Z_\mathrm{eff}$',
            save_type=SaveType.INPUT,
            minvalue=0.999,
        )

        self.wbe = Variable(
            'Bounce Frequency',
            label=r'$\omega_\mathrm{be}$',
            units=r'$s^{-1}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.wte = Variable(
            'Transit Frequency',
            label=r'$\omega_\mathrm{te}$',
            units=r'$s^{-1}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.xetgm_const = Variable(
            'Diffusivity Constant',
            label=r'${\rho^2_\mathrm{e,u}} v_\mathrm{Te} / L_\mathrm{Te}$',
            save_type=SaveType.ADDITIONAL,
        )

        self.vA = Variable(
            'Alfven Velocity',
            label=r'$v_{\rm A}$',
            units='m/s',
            save_type=SaveType.ADDITIONAL,
        )

        self.vei_nc = Variable(
            'Collision Frequency NC',
            label=r'$\nu_\mathrm{ei,NC}$',
            minvalue=0,
        )

        # Calculated Gradients
        self.gbtor = Variable(
            'Btor Gradient',
            label=r'$g_{\mathrm{B\phi}}$',
            absminvalue=1e-32,
        )

        self.gbu = Variable(
            'bu Gradient',
            label=r'$g_{\mathrm{Bu}}$',
            absminvalue=1e-32,
            # save_type=SaveType.INPUT,
        )

        self.gne = Variable(
            'Electron Density Gradient',
            label=r'$g_{\mathrm{ne}}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gq2 = Variable(
            'gq2',
            label=r'$g_{\mathrm{q}}$',
            absminvalue=1e-32,
        )

        self.gnh = Variable(
            'Hydrogenic Ion Density Gradient',
            label=r'$g_{\mathrm{nh}}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gni = Variable(
            'Thermal Ion Density Gradient',
            label=r'$g_{\mathrm{ni}}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gnz = Variable(
            'Impurity Density Gradient',
            label=r'$g_{\mathrm{nz}}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gq = Variable(
            'Safety Factor Gradient',
            label=r'$g_{q}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gte = Variable(
            'Electron Temperature Gradient',
            label=r'$g_{\mathrm{Te}}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gti = Variable(
            'Thermal Ion Temperature Gradient',
            label=r'$g_{\mathrm{Ti}}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gvpar = Variable(
            'Parallel Velocity Gradient',
            label=r'$g_{v\parallel}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gvpol = Variable(
            'Poloidal Velocity Gradient',
            label=r'$g_{v\theta}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gvtor = Variable(
            'Toroidal Velocity Gradient',
            label=r'$g_{v\phi}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        self.gelong = Variable(
            'Elongation Gradient',
            label=r'$g_{\kappa}$',
            absminvalue=1e-32,
            save_type=SaveType.INPUT,
        )

        # CDF MMM Variables
        self.conde = Variable(
            'CONDE',
            cdfvar='CONDE',
            units='m^2/s',
            label=r'$\chi_{\mathrm{e}}$',
            default_values=0,
            contour_max=100,
        )

        self.condi = Variable(
            'CONDI',
            cdfvar='CONDI',
            units='m^2/s',
            label=r'$\chi_{\mathrm{i}}$',
            default_values=0,
        )

        self.condepr = Variable(
            'CONDEPR',
            cdfvar='CONDEPR',
            units='m^2/s',
            label=r'$\chi_{\mathrm{e}}$',
            default_values=0,
        )

        self.condipr = Variable(
            'CONDIPR',
            cdfvar='CONDIPR',
            units='m^2/s',
            label=r'$\chi_{\mathrm{i}}$',
            default_values=0,
        )

        self.condewnc = Variable(
            'CONDEWNC',
            cdfvar='CONDEWNC',
            units='m^2/s',
            label=r'$\chi_{\mathrm{e}}$',
            default_values=0,
        )

        self.condiwnc = Variable(
            'CONDIWNC',
            cdfvar='CONDIWNC',
            units='m^2/s',
            label=r'$\chi_{\mathrm{i}}$',
            default_values=0,
        )

        self.xkemmm07 = Variable(
            'XKEMMM07',
            cdfvar='XKEMMM07',
            units='m^2/s',
            label=r'$\chi_{\mathrm{e}}$',
            default_values=0,
        )

        self.xkimmm07 = Variable(
            'XKIMMM07',
            cdfvar='XKIMMM07',
            units='m^2/s',
            label=r'$\chi_{\mathrm{i}}$',
            default_values=0,
        )

        self.xkepaleo = Variable(
            'Electron Thermal Diffusivity',
            cdfvar='XKEPALEO',
            units='m^2/s',
            label=r'$\chi_{\mathrm{e}}$',
            default_values=0,
        )

        self.xke = Variable(
            'Electron Thermal Diffusivity',
            units='m^2/s',
            label=r'$\chi_{\mathrm{e}}$',
            default_values=0,
            contour_max=200,
        )

        self.xki = Variable(
            'Ion Thermal Diffusivity',
            units='m^2/s',
            label=r'$\chi_{\mathrm{i}}$',
            default_values=0,
            contour_max=200,
        )

        self.drmin = Variable(
            'dr',
            label=r'$dr$',
            units='m',
        )

        self.darea = Variable(
            'dS',
            cdfvar='DAREA',
            label=r'$DAREA$',
            units='m^2',
        )

        self.icur = Variable(
            'Total plasma current',
            label=r'I$_{\rm p}$',
            units='A',
        )

        self.jcur = Variable(
            'Total plasma current density',
            cdfvar='CUR', label=r'j$_{\rm p}$',
            units='A/m^2',
        )

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

    def use_experimental_profiles(self):
        '''Attempts to use experimental temperature and safety factor profiles in place of calculated profiles'''
        if self.tepro.values is not None and (self.tepro.values > 1.001 * self.tepro.default_values).any():
            self.te.values = self.tepro.values

        if self.tipro.values is not None and (self.tipro.values > 1.001 * self.tipro.default_values).any():
            self.ti.values = self.tipro.values

        if self.qpro.values is not None and (self.qpro.values > 1.001 * self.qpro.default_values).any():
            self.q.values = self.qpro.values

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

        if settings.SAVE_ADDITIONAL_VARIABLES:
            self.save_vars_of_type(SaveType.ADDITIONAL, scan_factor)

    def reduce_memory(self):
        '''
        Remove variables listed here to reduce memory consumption

        Some variables are needed for calculations, but are no longer
        essential after calculations have been completed. This method should
        typically be called when running scans and not when plotting data.
        '''

        self.btbe = None
        self.conde = None
        self.condepr = None
        self.condewnc = None
        self.condi = None
        self.condipr = None
        self.condiwnc = None
        self.curdlh = None
        self.curdoh = None
        self.curlh = None
        self.curoh = None
        self.curohrho = None
        self.e_r_grp = None
        self.e_r_phi = None
        self.e_r_tht = None
        self.gne_threshold = None
        self.gte_threshold = None
        self.logi = None
        self.pcur = None
        self.pmhdf = None
        self.tbtbe = None
        self.tepro = None
        self.tipro = None
        self.tmhdf = None
        self.wexbsa = None
        self.wexbsmod = None
        self.wexbsv2 = None
        self.xke = None
        self.xkemmm07 = None
        self.xkepaleo = None
        self.xki = None
        self.xkimmm07 = None


class OutputVariables(Variables):
    '''
    Output variables consist of all variable data obtained as output from MMM
    (as well as rho and rmin)

    Please refer to the documentation provided with MMM for more information
    about the variables that obtained as output from MMM.
    '''

    def __init__(self, options=None):
        # Independent Variables
        self.rho = Variable('rho', units='', label=r'$\hat{\rho}$')
        self.rmin = Variable('Minor Radius', units='m', label=r'$r$', minvalue=0)
        self.rmina = Variable('rmina', label=r'$r/a$', units=r'', minvalue=0)
        # Total Fluxes
        self.fti = Variable('fti', units='keVm/s', label=r'$\Gamma_\mathrm{i}$', contour_min=0)
        self.fde = Variable('fde', units='m^-2s^-1', label=r'$\Gamma_\mathrm{n}$', contour_min=0)
        self.fte = Variable('fte', units='keVm/s', label=r'$\Gamma_\mathrm{e}$', contour_min=0)
        self.fdz = Variable('fdz', units='m^-2s^-1', label=r'$\Gamma_\mathrm{z}$', contour_min=0)
        self.fvt = Variable('fvt', units='m^2/s^2', label=r'$\Gamma_\phi$', contour_min=0)
        self.fvp = Variable('fvp', units='m^2/s^2', label=r'$\Gamma_\theta$', contour_min=0)
        # Convective Velocities
        self.vci = Variable('vci', units='m/s', label=r'$V_{\mathrm{i}}$')
        self.vch = Variable('vch', units='m/s', label=r'$V_{\mathrm{n}}$')
        self.vce = Variable('vce', units='m/s', label=r'$V_{\mathrm{e}}$')
        self.vcz = Variable('vcz', units='m/s', label=r'$V_{\mathrm{z}}$', contour_max=1e5, contour_min=-1e5)
        self.vct = Variable('vct', units='m/s', label=r'$V_{\phi}$', contour_max=1e5, contour_min=-1e5)
        self.vcp = Variable('vcp', units='m/s', label=r'$V_{\theta}$', contour_max=1e5, contour_min=-1e5)
        # Total Diffusivities
        self.xti = Variable('xti', units='m^2/s', label=r'$\chi_\mathrm{i}$')
        self.xde = Variable('xde', units='m^2/s', label=r'$\chi_\mathrm{n}$')
        self.xdi = Variable('xdi', units='m^2/s', label=r'$\chi_\mathrm{n}$')  # TODO: Remove eventually
        self.xte = Variable('xte', units='m^2/s', label=r'$\chi_\mathrm{e}$')
        self.xdz = Variable('xdz', units='m^2/s', label=r'$\chi_\mathrm{z}$')
        self.xvt = Variable('xvt', units='m^2/s', label=r'$\chi_\phi$')
        self.xvp = Variable('xvp', units='m^2/s', label=r'$\chi_\theta$')
        # Weiland Components
        self.xtiW20 = Variable('xtiW20', units='m^2/s', label=r'$\chi_\mathrm{i, w20}$')
        self.xdeW20 = Variable('xdeW20', units='m^2/s', label=r'$\chi_\mathrm{n, w20}$')
        self.xteW20 = Variable('xteW20', units='m^2/s', label=r'$\chi_\mathrm{e, w20}$')
        self.gmaW20ii = Variable('gmaW20ii', units='s^{-1}', label=r'$\gamma_\mathrm{ii, w20}$',
                                 contour_max=5e7, contour_min=-5e7)
        self.omgW20ii = Variable('omgW20ii', units='s^{-1}', label=r'$\omega_\mathrm{ii, w20}$',
                                 contour_max=5e6, contour_min=-5e6)
        self.gmaW20ie = Variable('gmaW20ie', units='s^{-1}', label=r'$\gamma_\mathrm{ie, w20}$',
                                 contour_max=5e7, contour_min=-5e7)
        self.omgW20ie = Variable('omgW20ie', units='s^{-1}', label=r'$\omega_\mathrm{ie, w20}$',
                                 contour_max=5e6, contour_min=-5e6)
        self.gmaW20ei = Variable('gmaW20ei', units='s^{-1}', label=r'$\gamma_\mathrm{ei, w20}$',
                                 contour_max=5e7, contour_min=-5e7)
        self.omgW20ei = Variable('omgW20ei', units='s^{-1}', label=r'$\omega_\mathrm{ei, w20}$',
                                 contour_max=5e6, contour_min=-5e6)
        self.gmaW20ee = Variable('gmaW20ee', units='s^{-1}', label=r'$\gamma_\mathrm{ee, w20}$',
                                 contour_max=5e7, contour_min=-5e7)
        self.omgW20ee = Variable('omgW20ee', units='s^{-1}', label=r'$\omega_\mathrm{ee, w20}$',
                                 contour_max=5e6, contour_min=-5e6)
        self.gmaW20e = Variable('gmaW20e', units='s^{-1}', label=r'$\gamma_\mathrm{e, w20}$',
                                 contour_max=5e5, contour_min=-5e5)
        self.omgW20e = Variable('omgW20e', units='s^{-1}', label=r'$\omega_\mathrm{e, w20}$',
                                 contour_max=5e5, contour_min=-5e5)
        self.gmaW20i = Variable('gmaW20i', units='s^{-1}', label=r'$\gamma_\mathrm{i, w20}$',
                                 contour_max=5e5, contour_min=-5e5)
        self.omgW20i = Variable('omgW20i', units='s^{-1}', label=r'$\omega_\mathrm{i, w20}$',
                                 contour_max=5e5, contour_min=-5e5)
        self.gaveW20i = Variable('gaveW20i', units='', label=r'$\overline{G}_{\rm i}$',
                                 contour_max=2, contour_min=0)
        self.gaveW20e = Variable('gaveW20e', units='', label=r'$\overline{G}_{\rm e}$',
                                 contour_max=2, contour_min=0)
        self.kyrhosW20i = Variable('kyrhosW20i', units='', label=r'$k_y\rho_\mathrm{s, i}$',
                                  )
        self.kyrhosW20e = Variable('kyrhosW20e', units='', label=r'$k_y\rho_\mathrm{s, e}$',
                                  )
        self.kparaW20i = Variable('kparaW20i', units='1/m', label=r'$\langle k_{\parallel}\! \rangle_{\rm i}$',
                                  )
        self.kparaW20e = Variable('kparaW20e', units='1/m', label=r'$\langle k_{\parallel}\! \rangle_{\rm e}$',
                                  )
        # DBM Components
        self.Apara2DBM = Variable('Apara2DBM', units='', label=r'$|\hat{A}_{\!\parallel}\!|^2$')
        self.gaveDBM = Variable('gaveDBM', units='', label=r'$\overline{G}$',contour_max=5, contour_min=-5)
        self.gmaDBM = Variable('gmaDBM', units='s^{-1}', label=r'$\gamma_\mathrm{dbm}$', contour_max=5e6, contour_min=-5e6)
        self.kyrhosDBM = Variable('kyrhosDBM', units='', label=r'$k_y\rho_\mathrm{s}$')
        self.kparaDBM = Variable('kparaDBM', units='m^-1', label=r'$\langle k_\parallel \rangle$')
        self.omgDBM = Variable('omgDBM', units='s^{-1}', label=r'$\omega_\mathrm{dbm}$', contour_max=5e6, contour_min=-5e6)
        self.phi2DBM = Variable('phi2DBM', units='', label=r'$|\hat{\phi}|^2$')
        self.satDBM = Variable('satDBM', units='', label=r'$2\hat{\gamma}/|\hat{\phi}| R k_\mathrm{x}$')
        self.xde2DBM = Variable('xde2DBM', units='m^2/s', label=r'$\chi^*_\mathrm{n, dbm}$')
        self.xdeDBM = Variable('xdeDBM', units='m^2/s', label=r'$\chi_\mathrm{n, dbm}$')
        self.xte2DBM = Variable('xte2DBM', units='m^2/s', label=r'$\chi^*_\mathrm{e, dbm}$')
        self.xteDBM = Variable('xteDBM', units='m^2/s', label=r'$\chi_\mathrm{e, dbm}$')
        self.xti2DBM = Variable('xti2DBM', units='m^2/s', label=r'$\chi^*_\mathrm{i, dbm}$')
        self.xtiDBM = Variable('xtiDBM', units='m^2/s', label=r'$\chi_\mathrm{i, dbm}$')
        # EPM Component
        self.gmaEPM = Variable('gmaEPM', units='s^{-1}', label=r'$\gamma_\mathrm{epm}$', contour_max=5e7, contour_min=-5e7)
        self.kyrhosEPM = Variable('kyrhosEPM', units='', label=r'$k_y\rho_\mathrm{s}$')
        self.omgEPM = Variable('omgEPM', units='s^{-1}', label=r'$\omega_\mathrm{epm}$', contour_max=5e7, contour_min=-5e7)
        self.xdeEPM = Variable('xdeEPM', units='m^2/s', label=r'$\chi_\mathrm{n, epm}$')
        self.xteEPM = Variable('xteEPM', units='m^2/s', label=r'$\chi_\mathrm{e, epm}$')
        self.xtiEPM = Variable('xtiEPM', units='m^2/s', label=r'$\chi_\mathrm{i, epm}$')
        self.nEPM = Variable('nEPM', units='', label=r'$n$')
        self.gaveEPM = Variable('gaveEPM', units='', label=r'$\overline{G}$', contour_max=5, contour_min=-5)
        self.kparaEPM = Variable('kparaEPM', units='m^-1', label=r'$\langle k_\parallel \rangle$')
        self.errorEPM = Variable('errorEPM', units='', label=r'Error$_{\rm epm}$')
        self.wdfEPM = Variable('wdfEPM', units='s^-1', label=r'$\omega_{\rm Df}$')
        self.wdeEPM = Variable('wdfEPM', units='s^-1', label=r'$\omega_{\rm De}$')
        self.wseEPM = Variable('wdfEPM', units='s^-1', label=r'$\omega_{\rm *e}$')
        # ETG Component
        self.xteETG = Variable('xteETG', units='m^2/s', label=r'$\chi_\mathrm{e, etg}$')
        self.gtecETG = Variable(r'Critical $g_\mathrm{Te}$ (Jenko ETG)', units='',
                                   label=r'$g_\mathrm{Te, crit}$')
        # MTM Components
        self.xteMTM = Variable('xteMTM', units='m^2/s', label=r'$\chi_\mathrm{e, mtm}$',
                               contour_max=25, contour_min=0)
        self.gmaMTM = Variable('gmaMTM', units='s^{-1}', label=r'$\gamma_\mathrm{mtm}$',
                               contour_max=5e7, contour_min=-5e7)
        self.gmanMTM = Variable('gmanMTM', units='', label=r'$\gamma\, a / c_\mathrm{s}$')
        self.omgMTM = Variable('omgMTM', units='s^{-1}', label=r'$\omega_\mathrm{mtm}$',
                               contour_max=3e6, contour_min=-3e6)
        self.omgnMTM = Variable('omgnMTM', units='', label=r'$\omega\, a / c_\mathrm{s}$')
        self.kyrhosMTM = Variable('kyrhosMTM', units='', label=r'$k_y\rho_\mathrm{s}$')
        self.dbsqMTM = Variable('dbsqMTM', units='', label=r'$|\delta B/B|^2$')
        # ETGM Components
        self.xteETGM = Variable('Thermal Diffusivity', units='m^2/s', label=r'$\chi_\mathrm{e, etgm}$',
                                contour_max=12.6, contour_min=-2e1)
        self.xte2ETGM = Variable('Thermal Diffusivity', units='m^2/s', label=r'$\chi^{\ast}_\mathrm{e, etgm}$',
                                 contour_max=1e1, contour_min=-1e1)
        self.gmaETGM = Variable('Growth Rate', units='s^{-1}', label=r'$\gamma_\mathrm{etgm}$',
                                contour_max=5e7, contour_min=-5e7)
        self.gmanETGM = Variable('Growth Rate', units='', label=r'$\gamma\ (c_{\rm s} / a)$',
                                contour_max=5e7, contour_min=-5e7)
        self.omgETGM = Variable('Frequency', units='', label=r'$\omega_\mathrm{etgm}$',
                                contour_max=3e6, contour_min=-3e6)
        self.omgnETGM = Variable('Frequency', units='', label=r'$\omega_{\rm r}\ (c_{\rm s} / a)$',
                                contour_max=3e6, contour_min=-3e6)
        self.kyrhoeETGM = Variable('Wave Number', units='', label=r'$k_y\rho_\mathrm{e}$')
        self.kyrhosETGM = Variable('Wave Number', units='', label=r'$k_y\rho_\mathrm{s}$')
        self.gaveETGM = Variable('Magnetic Field Curvature', units='', label=r'$\overline{G}$',
                                 contour_max=5, contour_min=-5)
        self.alphaETGM = Variable('alphaMHD', units='', label=r'$\alpha_\mathrm{MHD}$')
        self.kpara2ETGM = Variable('kpara2ETGM', units='m^-2', label=r'$\langle k^{2}_\parallel \rangle$')
        self.kparaETGM = Variable('kparaETGM', units='m^-1', label=r'$\langle k_\parallel \rangle$')
        self.fleETGM = Variable('fleETGM', units='', label=r'$\langle k^{2}_\perp \rho^{2}_\mathrm{e}\rangle$')
        self.phi2ETGM = Variable('Electrostatic Potential', units='', label=r'$|\hat{\phi}|^2$')
        self.Apara2ETGM = Variable('Electromagnetic Potential', units='', label=r'$|\hat{A}_{\!\parallel}\!|^2$')
        self.wdeETGM = Variable('wdeETGM', units='s^{-1}', label=r'$\omega_\mathrm{De}$',
                                   contour_max=3e6, contour_min=-3e6)
        self.wde_gaveETGM = Variable('wde_gaveETGM', units='s^{-1}', label=r'$\omega_\mathrm{De} / \overline{G}$',
                                        contour_max=3e6, contour_min=-3e6)
        self.omgdiffETGM = Variable('omgdiffETGM', units='s^{-1}', label=r'$\omega - \omega_\mathrm{De}$',
                                      contour_max=1e6, contour_min=-1e6)
        self.gmadiffETGM = Variable('gmadiffETGM', units='s^{-1}', label=r'$\gamma - \omega_\mathrm{De}$',
                                      contour_max=1e6, contour_min=-1e6)
        self.wseETGM = Variable('wseETGM', units='s^{-1}', label=r'$\omega_{*\mathrm{e}}$',
                                   contour_max=1e7, contour_min=-1e7)
        self.wteETGM = Variable('wteETGM', units='s^{-1}', label=r'$\omega_{\mathrm{Te}}$',
                                    contour_max=1e7, contour_min=-1e7)
        self.wsetaETGM = Variable('wsetaETGM', units='s^{-1}',
                                      label=r'$\omega_{*\mathrm{e}} (1 + \eta_\mathrm{e})$',
                                      contour_max=1e7, contour_min=-1e7)
        self.waETGM = Variable('Alfven Frequency', units='s^{-1}', label=r'$\omega_\mathrm{A}$')
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
    def __init__(self, name, cdfvar=None, smooth=None, label='', desc='', minvalue=None, absminvalue=None,
                 save_type=None, default_values=1.e-16, contour_max=np.inf, contour_min=-np.inf,
                 reduce_memory=False, units='', dimensions=None, values=None):
        # Public
        self.name = name
        self.cdfvar = cdfvar  # Name of variable as used in CDF's
        self.smooth = smooth  # None to disable smoothing, or n = positive integer
        self.label = label  # Plot label in LaTeX Format
        self.desc = desc  # Stores the long_name value from CDF's
        self.minvalue = minvalue  # minimum value variable is allowed to have
        self.absminvalue = absminvalue  # minimum value the absolute value of the variable is allowed to have
        self.save_type = save_type if save_type is not None else SaveType.NONE
        self.default_values = default_values  # values to use if variable not in CDF
        self.contour_max = contour_max  # maximum value to display on contour plots
        self.contour_min = contour_min  # minimum value to display on contour plots
        self.reduce_memory = reduce_memory  # delete this variable after calculations to reduce memory usage
        # Private
        self._units_label = ''
        self._units = ''
        # self._dimensions = dimensions if dimensions is not None else ['', '']
        self._dimensions = dimensions
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
        return self._values if self._values is not None else self.default_values

    @values.setter
    def values(self, values):
        if not isinstance(values, np.ndarray):
            raise ValueError(f'Variable values must be type {np.ndarray} and not {type(values)}')
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

        if self.smooth is not None and isinstance(self.values, np.ndarray):
            sigma = int(self.values.shape[0] * self.smooth / 100)
            if self.values.ndim == 2:
                self.values = scipy.ndimage.gaussian_filter(self.values, sigma=(sigma, 0))
            else:
                self.values = scipy.ndimage.gaussian_filter(self.values, sigma=sigma)

    def set_minvalue(self, ignore_exceptions=False):
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
        * ignore_exceptions (bool): Possible exceptions will be ignored when True

        Raises:
        * ValueError: If multiple nonphysical values are found
        '''

        if self.minvalue is not None and isinstance(self.values, np.ndarray):
            multiple_errors_per_timeval = (np.count_nonzero(self.values < self.minvalue, axis=0) > 1)
            if not ignore_exceptions and multiple_errors_per_timeval.any():
                idx_list = [i for i in np.where(multiple_errors_per_timeval)][0]
                raise ValueError(
                    f'Multiple Nonphysical values obtained for {self.name}\n'
                    f'    min value:     {self.values[:, idx_list].min()}\n'
                    f'    threshold:     {self.minvalue}\n'
                    f'    time indices:  {idx_list}\n'
                )
            # When an exception is not raised, fix the minimum value
            self.values[self.values < self.minvalue] = self.minvalue

        if self.absminvalue is not None and isinstance(self.values, np.ndarray):
            too_small = np.absolute(self.values) < self.absminvalue
            if too_small.any():
                value_signs = np.sign(self.values[too_small])
                value_signs[value_signs == 0] = 1  # np.sign(0) = 0, so set these to +1
                self.values[too_small] = self.absminvalue * value_signs

    def clamp_values(self, clamp_value):
        '''Clamps values between -clamp_value and +clamp_value'''
        self.values[self.values > clamp_value] = clamp_value
        self.values[self.values < -clamp_value] = -clamp_value

    def get_sign(self):
        """Get sign of variable values"""
        sign = np.sign(self.values)
        sign[sign==0] = 1
        return sign

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

    def check_for_nan(self, ignore_exceptions=False):
        '''Checks for nan values and raises a ValueError if any are found'''
        if np.isnan(self.values).any() and not ignore_exceptions:
            raise ValueError(f'nan values found in variable {self.name}')

    def reduce_memory(self):
        if self.delete:
            del self


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
