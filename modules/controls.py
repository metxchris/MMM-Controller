"""Handles values of controls used in the header of the MMM input file

A control is any value in the MMM input file that is not dependent on rmin,
the minor radius of the plasma.  Controls are generally used for switches or
scalar coefficients by MMM.  Some controls may also be scanned over, such as
kyrhoe in any of the component models. The InputControls class stores all
control values needed by MMM, and also handles generating the header of the
MMM input file, as used by the MMM wrapper.

The InputControls class is coupled with the Options class, which is used for
storing options needed for plotting, checking values, and both saving and
loading control data.  An Options object must be instantiated and passed in
when instantiating InputControls objects.  When loading control data from
CSVs, it is advised to load the options data first, and then use the loaded
options to load the controls.

Example Usage:
    # Instantiate Options for New InputControls
    options = modules.options.Options(
        runid='132017T01',
        shot_type=ShotType.NSTX,
        input_points=101,
    )

    # Instantiate InputControls and Set Values
    controls = InputControls(
        options,
        input_points=101,
        cmodel_weiland=1,
        cmodel_dribm=0,
        cmodel_etg=0,
    )

    # Get MMM Header
    controls.get_mmm_header()

    # Load Options for Loading InputControls
    options = modules.options.Options()
    options.load(runid='132017T01', scan_num=1)

    # Load Controls
    controls = InputControls(options)
    controls.load_from_csv()

    # Save Controls
    controls.save()
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.utils as utils
import modules.constants as constants
import modules.headers.controls123 as controls123
import modules.headers.controls114 as controls114
import modules.headers.controls113 as controls113
import modules.headers.controls107 as controls107
import modules.headers.controls105 as controls105
import modules.headers.controls90 as controls90
import modules.headers.controlsold as controlsold
from modules.enums import SaveType
import settings


_log = logging.getLogger(__name__)


class InputControls:
    '''
    Input Controls for the MMM input file

    Please refer to the documentation for MMM for more information about the
    various controls used here.

    Notes:
    * Controls defined here are will be placed into the header of the MMM input file
    * Controls with vtype=int are expected as Fortran Integer types in the input file
    * Controls with vtype=float are expected as Fortran Real types in the input file
    * The value defined in each Control object is the default value of that control

    Raises:
    * ValueError: If keyword arguments (kwargs) are provided on init and options are not
    * ValueError: If input_points.values is None
    * TypeError: If input_points.values is of type np.ndarray
    '''

    def __init__(self, options=None, **kwargs):
        if kwargs and not options:
            # options is only allowed to be None when instantiating the class solely for membership checks
            raise ValueError('options must be provided when providing keyword arguments')

        self.input_points = Control('Number of radial points', int, None)
        
        # Switches for component models
        self.cmodel_weiland = Control('Weiland', float, 1)
        self.cmodel_dribm = Control('DRIBM', float, 1)
        self.cmodel_epm = Control('EPM', float, 1)
        self.cmodel_etg = Control('ETG', float, 1)
        self.cmodel_etgm = Control('ETGM', float, 1)
        self.cmodel_mtm = Control('MTM', float, 1)

        # MMM Options
        self.mmm_separate_conv_vel = Control('[Separate Diffusivity and Convective Velocity] 0: False, 1: True', int, 1)
        self.mmm_convert_negative_chi = Control('[Convert Negative Diffusivity to Convective Velocity] 0: False, 1: True', int, 1)
        self.mmm_allow_negative_chi = Control('[Allow Negative Diffusivity Output] 0: False, 1: True', int, 1)
        self.mmm_use_solver_grads = Control('[Convert Diffusivity Using Solver Gradients] 0: False, 1: True', int, 1)
        self.mmm_limit_small_grads = Control('[Limit Diffusivity From Small Gradients] 0: False, 1: True', int, 1)
        self.mmm_omp_threads = Control('[OpenMP Threads] 0: auto, 1+: thread count', int, 0)

        self.mmm_xti_max = Control('max |xti|', float, 200)
        self.mmm_xde_max = Control('max |xde|', float, 200)
        self.mmm_xte_max = Control('max |xte|', float, 200)
        self.mmm_xdz_max = Control('max |xdz|', float, 200)
        self.mmm_xvt_max = Control('max |xvt|', float, 200)
        self.mmm_xvp_max = Control('max |xvp|', float, 200)
        self.mmm_vti_max = Control('max |vti|', float, 5e6)
        self.mmm_vde_max = Control('max |vde|', float, 5e6)
        self.mmm_vte_max = Control('max |vte|', float, 5e6)
        self.mmm_vdz_max = Control('max |vdz|', float, 5e6)
        self.mmm_vvt_max = Control('max |vvt|', float, 5e6)
        self.mmm_vvp_max = Control('max |vvp|', float, 5e6)

        # Weiland options
        self.w20_shear_def = Control('[Effective Shear Definition] 0: use elong; 1: use grad rho', int, 0)
        self.w20_negative_chi = Control('[Minimum Diffusivity] 0: Zero, 1: Negative of max', int, 1)

        self.w20_exbs = Control('ExB shear coefficient', float, 1)
        self.w20_vconv = Control('Convective velocity coefficient', float, 1)
        self.w20_vconv_mult = Control('Max convective velocity coefficient', float, 1)
        self.w20_diff_mult = Control('Max diffusivity coefficient', float, 0.5)
        self.w20_kyrhos = Control('kyrhos', float, 0.316, label=r'$k_\mathrm{y}\rho_\mathrm{s}$')
        self.w20_xti_max = Control('max |xti|', float, 1e6)
        self.w20_xde_max = Control('max |xde|', float, 1e6)
        self.w20_xte_max = Control('max |xte|', float, 1e6)
        self.w20_xdz_max = Control('max |xdz|', float, 1e6)
        self.w20_xvt_max = Control('max |xvt|', float, 1e6)
        self.w20_xvp_max = Control('max |xvp|', float, 1e6)
        self.w20_xte_min = Control('Lower Bound', float, -100)
        self.w20_xti_min = Control('Lower Bound', float, -100)
        self.w20_xde_min = Control('Lower Bound', float, -100)
        self.w20_xdz_min = Control('Lower Bound', float, -100)
        self.w20_xvp_min = Control('Lower Bound', float, -100)
        self.w20_xvt_min = Control('Lower Bound', float, -100)
        self.w20_gmult = Control('gmult', float, 1, label='guess factor', )
        self.w20_mpsf = Control('Convective velocity coefficient', float, 1)
        
        # EPM options
        self.epm_direction = Control('[Search Direction] 0: Any, 1: Electron, -1: Ion', int, 0)
        self.epm_n_start = Control('Toroidal wave number start value (min 1)', int, 1, label=r'$n$')
        self.epm_n_end = Control('Toroidal wave number end value ', int, 10)
        self.epm_n_step = Control('Toroidal wave number step value (min 1)', int, 1)
        self.epm_sum_modes = Control('[Mode Handling] 1: Sum modes, 0: Most unstable mode', int, 0)

        self.epm_exbs = Control('ExB shear coefficient', float, 1)
        self.epm_xti_min = Control('Lower bound of ion thermal diffusivity', float, -100)
        self.epm_xti_max = Control('Upper bound of ion thermal diffusivity', float, 100)
        self.epm_xde_min = Control('Lower bound of electron particle diffusivity', float, -100)
        self.epm_xde_max = Control('Upper bound of electron particle diffusivity', float, 100)
        self.epm_xte_min = Control('Lower bound of electron thermal diffusivity', float, -100)
        self.epm_xte_max = Control('Upper bound of electron thermal diffusivity', float, 100)
        self.epm_chi_cal = Control('Diffusivity calibration', float, 1)
 
        # DRIBM options
        self.dribm_direction = Control('[Search Direction] 0: Any, 1: Electron, -1: Ion', int, 0)
        self.dribm_kyrhos_scan = Control('Number of scanned ky modes (min 50), 0: Disable ky scan', int, 20)
        self.dribm_kyrhos_layers = Control('kyrhos layers', int, 10)
        self.dribm_diffusivity_type = Control('[DBM Model] 0: Default, 1: Alternate', int, 0)
        self.dribm_sat_expo = Control('Saturation ratio exponent (min 0)', int, 2)
        self.dribm_sum_modes = Control('[Mode Handling] 1: Sum modes, 0: Most unstable mode', int, 0)
        self.dribm_negative_chi = Control('[Minimum Diffusivity] 0: Zero, 1: Negative of max', int, 0)
        self.dribm_kyrhos_type = Control('[kyrhos Increments] 1: Exponential, 0: Linear', int, 1)

        self.dribm_exbs = Control('ExB shear coefficient', float, 1)
        self.dribm_kyrhos_min = Control('lower limit of kyrhos scan', float, 1, label=r'$k_\mathrm{y}\rho_\mathrm{s}$')
        self.dribm_kyrhos_max = Control('upper limit of kyrhos scan', float, 50, label=r'$(k_\mathrm{y}\rho_\mathrm{s})_\mathrm{max}$')
        self.dribm_kxoky = Control('kx / ky', float, 0.1)
        self.dribm_diff_mult = Control('Max diffusivity coefficient', float, 0.5)
        self.dribm_xti_max = Control('max |xti|', float, 100)
        self.dribm_xde_max = Control('max |xde|', float, 100)
        self.dribm_xte_max = Control('max |xte|', float, 100)
        self.dribm_xti_min = Control('Lower bound of ion thermal diffusivity', float, -100)
        self.dribm_xde_min = Control('Lower bound of electron particle diffusivity', float, -100)
        self.dribm_xte_min = Control('Lower bound of electron thermal diffusivity', float, -100)
        self.dribm_chi_max_cal = Control('Diffusivity calibration (max)', float, 1)
        self.dribm_chi_sum_cal = Control('Diffusivity calibration (sum)', float, 1)
        self.dribm_chi2_max_cal = Control('Alternate Diffusivity calibration (max)', float, 1)
        self.dribm_chi2_sum_cal = Control('Alternate Diffusivity calibration (sum)', float, 1)
        self.dribm_kpsh_mult = Control('DRIBM extra real', float, 1, label='zkpsh')
        self.dribm_ti_te_mult = Control('DRIBM extra real', float, 1, label=r'$T_\mathrm{i} / T_\mathrm{e}$')
        self.dribm_vei_mult = Control('DRIBM extra real', float, 1, label=r'$\nu_\mathrm{ei}$')
        
        # MTM options
        self.mtm_kyrhos_loops = Control('kyrhos scans for sum', int, 200)
        self.mtm_kyrhos_layer_loops = Control('kyrhos scans per layer', int, 7)
        self.mtm_negative_chi = Control('[Minimum Diffusivity] 0: Zero, 1: Negative of max', int, 0)
        self.mtm_kyrhos_type = Control('[kyrhos Increments] 1: Exponential, 0: Linear', int, 1)
        self.mtm_kyrhos_layers = Control('kyrhos layers', int, 10)

        self.mtm_ky_kx = Control('ky / kx', float, 0.2, label=r'$k_y/k_x$')
        self.mtm_cf = Control('calibration factor', float, 1)
        self.mtm_kyrhos_min = Control('lower limit of kyrhos scan', float, 0.005, label=r'$k_\mathrm{y}\rho_\mathrm{s}$')
        self.mtm_kyrhos_max = Control('upper limit of kyrhos scan', float, 10, label=r'$(k_\mathrm{y}\rho_\mathrm{s})_\mathrm{max}$')
        self.mtm_diff_mult = Control('Max diffusivity coefficient', float, 0.25)
        self.mtm_lbetd = Control('Lower bound of electron thermal diffusivity', float, 0)
        self.mtm_ubetd = Control('Upper bound of electron thermal diffusivity', float, 100)
        self.mtm_xte_max = Control('max |xte|', float, 100)
        self.mtm_xte_min = Control('min xte', float, -100)
        
        # ETG options
        self.etg_jenko_threshold = Control('[Jenko threshold] 1: Electrostatic, 2: Electromagnetic', int, 1)
        self.etg_restrict_chi = Control('[Diffusivity Values] 0: No restrictions, 1: Strictly non-negative', int, 0)

        self.etg_cees_scale = Control('CEES scale', float, 0.06)
        self.etg_ceem_scale = Control('CEEM scale', float, 0.06)
        self.etg_xte_min = Control('CEEM scale', float, 0)
        self.etg_xte_max = Control('max |xte|', float, 100)

        # ETGM options
        self.etgm_direction = Control('[Search Direction] 0: All, 1: Electron, -1: Ion', int, 0)
        self.etgm_sat_expo = Control('Saturation ratio exponent (min 0)', int, 2)
        self.etgm_sum_modes = Control('Sum modes, 0: Use most unstable mode', int, 0)
        self.etgm_kyrhos_type = Control('1: Exponential kyrhos increments, 0: Linear kyrhos increments', int, 1)
        self.etgm_kyrhos_scan = Control('Number of kyrhos scan loops (min 50), 0: disable kyrhos scan', int, 50)
        self.etgm_diffusivity_type = Control('[ETGM Model] 0: Default, 1: Alternate, 2: Horton', int, 0)
        self.etgm_negative_chi = Control('[Minimum Diffusivity] 0: Zero, 1: Negative of max', int, 0)
        self.etgm_jenko_threshold = Control('[Jenko threshold] 1: Electrostatic, 2: Electromagnetic', int, 1)
        self.etgm_kyrhos_layers = Control('kyrhos layers', int, 10)
        self.etgm_empty_int = Control('empty', int, 1)

        self.etgm_exbs = Control('ExB shear coefficient', float, 0)
        self.etgm_kyrhos_min = Control('lower limit of kyrhos scan', float, 1, label=r'$k_\mathrm{y}\rho_\mathrm{s}$')
        self.etgm_kyrhos_max = Control('upper limit of kyrhos scan', float, 50)
        self.etgm_kxoky = Control('kx / ky', float, 0.1)
        self.etgm_lbetd = Control('Lower bound of electron thermal diffusivity', float, -100)
        self.etgm_ubetd = Control('Upper bound of electron thermal diffusivity', float, 100)
        self.etgm_diff_mult = Control('Max diffusivity coefficient', float, 0.25)
        self.etgm_xte_max = Control('max |xte|', float, 100)
        self.etgm_xte_min = Control('min xte', float, 100)
        self.etgm_xte_max_cal = Control('xte calibration (max)', float, 1)
        self.etgm_xte_sum_cal = Control('xte calibration (sum)', float, 1)
        self.etgm_xte2_max_cal = Control('xte2 calibration (max)', float, 1)
        self.etgm_xte2_sum_cal = Control('xte2 calibration (max)', float, 1)

        self.etgm_disable_geometry = Control('1: Theta, alpha = 0, 0: Default G_ave', int, 0)
        self.etgm_electrostatic = Control('1: Electrostatic, 0: Electromagnetic', int, 0)
        self.etgm_alpha_mult = Control('alphaMHD mult', float, 1, label=r'$\alpha_\mathrm{MHD,u}$')
        self.etgm_betae_mult = Control('betae mult', float, 1, label=r'$\beta_\mathrm{e,u}$')
        self.etgm_nuei_mult = Control('nuei mult', float, 1, label=r'$\nu_\mathrm{ei}$')
        self.etgm_vthe_mult = Control('vthe mult', float, 1, label=r'$v_\mathrm{Te}$')
        self.etgm_betaep_mult = Control('beta_prime mult', float, 1, label=r'$\beta^\prime$')
        self.etgm_extra_mult = Control('extra', float, 1)

        # Verbose level
        self.lprint = Control('Verbose Level', values=0, vtype=int)

        self.input_points.values = options.input_points if options else None
        self.options = options

        self.set(**kwargs)
        self.set_wexb_factor()
        self.set_allow_negative_chi()
        self.set_cmodel_values()

        if not self.input_points.values and self.options.input_points:
            self.input_points.values = self.options.input_points

        # Error checking
        if not self.input_points.values:
            raise ValueError('input_points must be set to generate the MMM header')

        if isinstance(self.input_points.values, np.ndarray):
            raise TypeError('Unable to create MMM header for controls loaded with array values')

    def set(self, **kwargs):
        '''Sets specified control values'''
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise ValueError(f'Invalid control specified: {key}')
            getattr(self, key).values = value

        self.verify_values()

    def set_wexb_factor(self):
        self.w20_exbs.values = self.options.wexb_factor
        self.dribm_exbs.values = self.options.wexb_factor
        self.etgm_exbs.values = self.options.wexb_factor
        self.epm_exbs.values = self.options.wexb_factor

    def set_allow_negative_chi(self):
        self.mmm_allow_negative_chi.values = self.options.allow_negative_chi
        self.w20_negative_chi.values = self.options.allow_negative_chi
        self.dribm_negative_chi.values = self.options.allow_negative_chi
        self.etgm_negative_chi.values = self.options.allow_negative_chi
        self.mtm_negative_chi.values = self.options.allow_negative_chi

    def set_cmodel_values(self):
        self.cmodel_weiland.values = self.options.cmodel_w20
        self.cmodel_dribm.values = self.options.cmodel_dbm
        self.cmodel_epm.values = self.options.cmodel_epm
        self.cmodel_etg.values = self.options.cmodel_etg
        self.cmodel_etgm.values = self.options.cmodel_etgm
        self.cmodel_mtm.values = self.options.cmodel_mtm

    def verify_values(self):
        '''Verifies that certain control values are correct and fixes them if needed'''
        ...

    def get_mmm_header(self):
        '''
        Gets the header for the MMM input file

        Raises:
        * ValueError: If no header is defined for settings.MMM_HEADER_VERSION
        '''

        if settings.MMM_HEADER_VERSION in ['#123', '#129']:
            return controls123.get_mmm_header(self)
        if settings.MMM_HEADER_VERSION in ['#114', '#117']:
            return controls114.get_mmm_header(self)
        if settings.MMM_HEADER_VERSION in ['#113']:
            return controls113.get_mmm_header(self)
        if settings.MMM_HEADER_VERSION in ['#107', '#111']:
            return controls107.get_mmm_header(self)
        if settings.MMM_HEADER_VERSION == '#105':
            return controls105.get_mmm_header(self)
        if settings.MMM_HEADER_VERSION == '#90':
            return controls90.get_mmm_header(self)
        if settings.MMM_HEADER_VERSION == 'old':
            return controlsold.get_mmm_header(self)

        raise ValueError(f'No header rule defined for {settings.MMM_HEADER_VERSION}')

    def get_scanned_control(self):
        return getattr(self, self.options.var_to_scan)

    def get_keys(self):
        '''Returns (list): All keys of input controls'''
        return [o for o in dir(self) if isinstance(getattr(self, o), Control)]

    def get_key_values_pairs(self):
        '''Returns (list): All key-values pairs of input controls'''
        keys = self.get_keys()
        return [f'{o}, {getattr(self, o).values}' for o in keys]

    def print_key_values_pairs(self):
        '''Prints: All key-values pairs of input controls'''
        kvps = self.get_key_values_pairs()
        for kvp in kvps:
            print(kvp)

    def save(self, scan_factor=None):
        '''
        Saves InputControls data to CSV

        If scan_factor is specified, then var_to_scan must also be specified

        Parameters:
        * scan_factor (float): The value of the scan factor (Optional)
        '''

        runid = self.options.runid
        scan_num = self.options.scan_num
        var_to_scan = self.options.var_to_scan

        if scan_factor:
            scan_factor_str = f'{scan_factor:{constants.SCAN_FACTOR_FMT}}'
            save_dir = utils.get_var_to_scan_path(runid, scan_num, var_to_scan)
            file_name = (f'{save_dir}\\{SaveType.CONTROLS.name.capitalize()} {var_to_scan}'
                         f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}{scan_factor_str}.csv')
        else:
            save_dir = utils.get_scan_num_path(runid, scan_num)
            file_name = f'{save_dir}\\{SaveType.CONTROLS.name.capitalize()}.csv'

        control_data = self.get_key_values_pairs()

        with open(file_name, 'w') as f:
            for data in control_data:
                f.write(f'{data}\n')

        _log.info(f'\n\tSaved: {file_name}\n')

    def load_from_csv(self, scan_factor=None, use_rho=False):
        '''
        Loads Controls data from a CSV into the current Controls object

        If either parameter use_rho or scan_factor are specified, then
        var_to_scan must also be specified.

        Parameters:
        * scan_factor (float): The scan_factor, if doing a variable scan (optional)
        * use_rho (bool): True if the CSV to load is in the rho folder (optional)
        '''

        runid = self.options.runid
        scan_num = self.options.scan_num
        var_to_scan = self.options.var_to_scan
        controls_name = SaveType.CONTROLS.name.capitalize()

        if use_rho:
            dir_path = utils.get_rho_path(runid, scan_num, var_to_scan)
            control_files = utils.get_files_in_dir(dir_path, f'{controls_name}*', show_warning=False)

            if scan_factor:
                _log.warning(f'\n\tThe scan_factor input parameter is not used when use_rho is True')

        elif scan_factor is not None:
            dir_path = utils.get_var_to_scan_path(runid, scan_num, var_to_scan)
            control_files = utils.get_files_in_dir(dir_path, f'{controls_name}*', show_warning=False)
            scan_factor_str = f'{scan_factor:{constants.SCAN_FACTOR_FMT}}'
            control_files = [file for file in control_files if scan_factor_str in file]

        else:
            dir_path = utils.get_scan_num_path(runid, scan_num)
            control_files = utils.get_files_in_dir(dir_path, f'{controls_name}*')

        # There may not be a Controls.CSV file to load in some cases, which is expected for
        # plotting scanned variables when var_to_scan is not a Control
        if control_files:
            if use_rho:
                self._load_from_np_csv(control_files[0])
            else:
                self._load_from_simple_csv(control_files[0])

    def _load_from_simple_csv(self, file_name):
        '''
        Loads a simple CSV where each line is a single (key, value) pair

        Parameters:
        * file_name (str): The name and path of the file to open
        '''

        with open(file_name, 'r') as file:
            for line in file:
                key, value = line.replace('\n', '').split(',')
                if hasattr(self, key):
                    getattr(self, key).values = float(value)

    def _load_from_np_csv(self, file_name):
        '''
        Loads a traditional CSV saved by Numpy where the first row contains all the keys, and values are in columns

        Parameters:
        * file_name (str): The name and path of the file to open
        '''

        data_array = np.genfromtxt(file_name, delimiter=',', dtype=float, names=True)
        control_names = data_array.dtype.names
        for name in control_names:
            if hasattr(self, name):
                getattr(self, name).values = data_array[name]


class Control:
    def __init__(self, name, vtype, values, label='', units_label=''):
        self.name = name
        self.vtype = vtype
        self.values = values
        self.label = label
        self.units_label = units_label

    def get_value_str(self):
        return int(self.values) if self.vtype is int else self.values

    def get_input_line(self):
        return f'{self.get_value_str()}  ! {self.name}\n'


'''
For testing purposes:
* There need to be existing folders corresponding to the runid and scan_num when saving controls
'''
if __name__ == '__main__':
    settings.MMM_HEADER_VERSION = '#113'

    import modules.options
    options = modules.options.Options(runid='TEST', scan_num=373, input_points=51)

    '''Print sample MMM Header from user-specified Options'''
    controls = InputControls(options)
    print(f'Example MMM Header:\n{"-"*50}')
    print(controls.get_mmm_header())

    '''Print Controls values loaded from CSV'''
    # controls = InputControls(options)
    # controls.load_from_csv()
    # print(f'Controls Values Loaded From CSV:\n{"-"*50}')
    # controls.print_key_values_pairs()

    '''Save Controls to CSV'''
    # controls.save()
