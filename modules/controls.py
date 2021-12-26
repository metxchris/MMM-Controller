# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.utils as utils
import modules.constants as constants
from modules.enums import ShotType, SaveType


class InputControls:
    '''
    Input Controls for the MMM input file

    Notes:
    * Controls defined here are will be placed into the header of the MMM input file
    * Controls with vtype=int are expected as Fortran Integer types in the input file
    * Controls with vtype=float are expected as Fortran Real types in the input file
    * The value defined in each Control object is the default value of that control
    '''

    def __init__(self, options=None):
        self.npoints = Control('npoints', 'Number of radial points', values=51, vtype=int)
        # Switches for component models
        self.cmodel_weiland = Control('cmodel_weiland', 'Weiland', values=1, vtype=float)
        self.cmodel_dribm = Control('cmodel_dribm', 'DRIBM', values=1, vtype=float)
        self.cmodel_etg = Control('cmodel_etg', 'ETG', values=1, vtype=float)
        self.cmodel_etgm = Control('cmodel_etgm', 'ETGM', values=1, vtype=float)
        self.cmodel_mtm = Control('cmodel_mtm', 'MTM', values=1, vtype=float)
        # Weiland options
        self.weiland_exbs = Control('weiland_exbs', 'ExB shear coefficient', values=1, vtype=float)
        self.weiland_mpsf = Control('weiland_mpsf', 'Momentum pinch scaling factor', values=1, vtype=float)
        self.weiland_lbetd = Control('weiland_lbetd', 'Lower bound of electron thermal diffusivity', values=0, vtype=float)
        self.weiland_ubetd = Control('weiland_ubetd', 'Upper bound of electron thermal diffusivity', values=100, vtype=float)
        self.weiland_lbitd = Control('weiland_lbitd', 'Lower bound of ion thermal diffusivity', values=0, vtype=float)
        self.weiland_ubitd = Control('weiland_ubitd', 'Upper bound of ion thermal diffusivity', values=100, vtype=float)
        # DRIBM options
        self.dribm_exbs = Control('dribm_exbs', 'ExB shear coefficient', values=1, vtype=float)
        self.dribm_kyrhos = Control('dribm_kyrhos', 'kyrhos', values=0.1, vtype=float, label=r'$k_y \rho_s$')
        # MTM options
        self.mtm_ky_kx = Control('mtm_ky_kx', 'ky/kx for MTM', values=0.2, vtype=float, label=r'$k_y/k_x$')
        self.mtm_cf = Control('mtm_cf', 'calibration factor', values=1.0, vtype=float)
        # ETG options
        self.etg_jenko_threshold = Control('etg_jenko_threshold', 'Jenko threshold', values=2, vtype=int)
        self.etg_cees_scale = Control('etg_cees_scale', 'CEES scale', values=0.06, vtype=float)
        self.etg_ceem_scale = Control('etg_ceem_scale', 'CEEM scale', values=0.06, vtype=float)
        # ETGM options
        self.etgm_cl = Control('etgm_cl', 'Collisionless limit', values=1, vtype=int)
        self.etgm_exbs = Control('etgm_exbs', 'ExB shear coefficient', values=0.0, vtype=float)
        self.etgm_kyrhoe = Control('etgm_kyrhoe', 'kyrhoe', values=0.25, vtype=float, label=r'$k_y \rho_e$')
        self.etgm_kyrhos = Control('etgm_kyrhos', 'kyrhos', values=0.33, vtype=float, label=r'$k_y \rho_s$')
        # Verbose level
        self.lprint = Control('lprint', 'Verbose Level', values=0, vtype=int)

        # Private members are not input controls, but are used to help verify control values
        self._shot_type = ShotType.NONE

        # Update values from options, if available
        if options is not None:
            self.update_from_options(options)

    def update_from_options(self, options):
        '''Updates any controls dependent on options values, then verifies values'''
        self.npoints.values = options.input_points
        self._shot_type = options.shot_type
        self.verify_values()

    def set(self, **kwargs):
        '''Sets all control values, then verifies values'''
        for key, value in kwargs.items():
            getattr(self, key).values = value

        self.verify_values()

    def verify_values(self):
        '''Verifies that certain values are correct and fixes them if needed'''

        # Note: The DRIBM model is currently disabled for NSTX discharges
        if self._shot_type == ShotType.NSTX:
            self.cmodel_dribm.values = 0

        self.etgm_kyrhos.values = max(1e-6, self.etgm_kyrhos.values)
        self.etgm_kyrhoe.values = max(0, self.etgm_kyrhoe.values)

    def get_mmm_header(self):
        if type(self.npoints.values) is np.ndarray:
            raise ValueError('Unable to create MMM header for controls loaded with array values')

        return MMM_HEADER.format(
            npoints=self.npoints.get_value_str(),
            cmodel_weiland=self.cmodel_weiland.get_value_str(),
            cmodel_dribm=self.cmodel_dribm.get_value_str(),
            cmodel_etg=self.cmodel_etg.get_value_str(),
            cmodel_etgm=self.cmodel_etgm.get_value_str(),
            cmodel_mtm=self.cmodel_mtm.get_value_str(),
            weiland_exbs=self.weiland_exbs.get_value_str(),
            weiland_mpsf=self.weiland_mpsf.get_value_str(),
            weiland_lbetd=self.weiland_lbetd.get_value_str(),
            weiland_ubetd=self.weiland_ubetd.get_value_str(),
            weiland_lbitd=self.weiland_lbitd.get_value_str(),
            weiland_ubitd=self.weiland_ubitd.get_value_str(),
            dribm_exbs=self.dribm_exbs.get_value_str(),
            dribm_kyrhos=self.dribm_kyrhos.get_value_str(),
            mtm_ky_kx=self.mtm_ky_kx.get_value_str(),
            mtm_cf=self.mtm_cf.get_value_str(),
            etg_jenko_threshold=self.etg_jenko_threshold.get_value_str(),
            etg_cees_scale=self.etg_cees_scale.get_value_str(),
            etg_ceem_scale=self.etg_ceem_scale.get_value_str(),
            etgm_cl=self.etgm_cl.get_value_str(),
            etgm_exbs=self.etgm_exbs.get_value_str(),
            etgm_kyrhoe=self.etgm_kyrhoe.get_value_str(),
            etgm_kyrhos=self.etgm_kyrhos.get_value_str(),
            lprint=self.lprint.get_value_str(),
        )

    def get_keys(self):
        '''Returns (list): All keys of input controls'''
        return [o for o in dir(self) if not callable(getattr(self, o)) and not o.startswith("_")]

    def get_key_values_pairs(self):
        '''Returns (list): All key-values pairs of input controls'''
        keys = self.get_keys()
        return [str(o) + ', ' + str(getattr(self, o).values).replace('\n', '') for o in keys]

    def print_key_values_pairs(self):
        '''Prints: All key-values pairs of input controls'''
        kvps = self.get_key_values_pairs()
        for kvp in kvps:
            print(kvp)

    def save_to_csv(self, options, scan_factor=None):
        '''
        Saves InputControls data to CSV

        Parameters:
        * options (OptionsData): Options.instance
        * scan_factor (float): The value of the scan factor (Optional)
        '''

        if scan_factor is not None:
            scan_factor_str = constants.SCAN_FACTOR_FMT_STR.format(scan_factor)
            save_dir = utils.get_var_to_scan_path(options.runid, options.scan_num, options.var_to_scan)
            file_name = (f'{save_dir}\\{SaveType.CONTROLS.name.capitalize()} {options.var_to_scan}'
                         f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}{scan_factor_str}.csv')
        else:
            save_dir = utils.get_scan_num_path(options.runid, options.scan_num)
            file_name = f'{save_dir}\\{SaveType.CONTROLS.name.capitalize()}.csv'

        control_data = self.get_key_values_pairs()

        f = open(file_name, 'w')
        for data in control_data:
            f.write(f'{data}\n')

        print(f'Input controls saved to \n    {file_name}\n')

    def load_from_csv(self, runid, scan_num, var_to_scan=None, scan_factor=None, use_rho=False):
        '''
        Loads Controls data from a CSV into the current Controls object

        Parameters:
        * runid (str): The runid of the CSV to use
        * scan_num (int): The scan number of the CSV to use
        * var_to_scan (str): The scanned variable of the CSV to use (optional)
        * scan_factor (float): The scan_factor, if doing a parameter scan (optional)
        * use_rho (bool): True if the CSV to load is in the rho folder (optional)
        '''

        controls_name = SaveType.CONTROLS.name.capitalize()

        if use_rho:
            dir_path = utils.get_rho_path(runid, scan_num, var_to_scan)
            control_files = utils.get_files_in_dir(dir_path, f'{controls_name}*', show_warning=False)

        elif scan_factor is not None:
            dir_path = utils.get_var_to_scan_path(runid, scan_num, var_to_scan)
            control_files = utils.get_files_in_dir(dir_path, f'{controls_name}*', show_warning=False)
            scan_factor_str = constants.SCAN_FACTOR_FMT_STR.format(scan_factor)
            control_files = [file for file in control_files if scan_factor_str in file]

        else:
            dir_path = utils.get_scan_num_path(runid, scan_num)
            control_files = utils.get_files_in_dir(dir_path, f'{controls_name}*', show_warning=False)

        # There may not be a Controls.CSV file to load in some cases, which is expected when
        # plotting scanned variables, if the scanned variable was not an InputControl
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

        file = open(file_name, 'r')
        for line in file:
            key, value = line.replace('\n', '').split(',')
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
            getattr(self, name).values = data_array[name]


class Control:
    def __init__(self, name, desc, values, vtype, label='', units_label=''):
        self.name = name
        self.desc = desc
        self.values = values
        self.vtype = vtype
        self.label = label
        self.units_label = units_label

    def get_value_str(self):
        return int(self.values) if self.vtype is int else f'{constants.INPUT_CONTROL_VALUE_FMT_STR.format(self.values)}D0'


# Header for MMM input file
MMM_HEADER = (
    '''&testmmm_input_control
    npoints = {npoints}    ! Number of radial points
    input_kind = 1
    /
    &testmmm_input_1stkind
    ! This is a sample input file of the first kind
    ! of an NSTX discharge

    !.. Switches for component models
    !   1D0 - ON, 0D0 - OFF
    cmodel  =
       {cmodel_weiland}    ! Weiland
       {cmodel_dribm}    ! DRIBM
       {cmodel_etg}    ! ETG
       {cmodel_mtm}    ! MTM
       {cmodel_etgm}    ! ETGM

    !.. Weiland real options
    cW20 =
       {weiland_exbs}     ! ExB shear coefficient
       {weiland_mpsf}     ! Momentum pinch scaling factor
       {weiland_lbetd}     ! Lower bound of electron thermal diffusivity
       {weiland_ubetd}     ! Upper bound of electron thermal diffusivity
       {weiland_lbitd}     ! Lower bound of ion thermal diffusivity
       {weiland_ubitd}     ! Upper bound of ion thermal diffusivity

    !.. DRIBM real options
    cDBM =
       {dribm_exbs}     ! ExB shear coefficient
       {dribm_kyrhos}   ! kyrhos

    !.. MTM real options
    cMTM =
       {mtm_ky_kx}   ! ky/kx for MTM
       {mtm_cf}   ! calibration factor

    !.. ETG integer options
    lETG =
       {etg_jenko_threshold}       ! Jenko threshold
               ! applied to both electrostatic and electromagnetic regimes

    !.. ETG real options
    cETG =
       {etg_cees_scale}    ! CEES scale
       {etg_ceem_scale}    ! CEEM scale

    !.. ETGM integer options
    lETGM =
       {etgm_cl}      ! Collisionless limit

    !.. ETGM real options
    cETGM =
       {etgm_exbs}     ! ExB shear coefficient
       {etgm_kyrhos}   ! kyrhos
       {etgm_kyrhoe}   ! kyrhoe

    lprint   = {lprint}      ! Verbose level\n\n''')


'''
For testing purposes:
* There need to be existing folders corresponding to the runid and scan_num when saving controls
'''
if __name__ == '__main__':
    from modules.options import Options

    '''Print sample MMM Header from user-specified Options'''
    Options.instance.set(
        runid='TEST',
        shot_type=ShotType.NSTX,
        input_points=51,
        scan_num=1)
    controls = InputControls(Options.instance)
    print(f'Example MMM Header:\n{"-"*50}')
    print(controls.get_mmm_header())

    '''Print Controls values from saved Options CSV'''
    Options.instance.load_options(
        runid='TEST',
        scan_num=5,
    )
    controls = InputControls(Options.instance)
    controls.load_from_csv(
        Options.instance.runid,
        Options.instance.scan_num,
        Options.instance.var_to_scan,
        scan_factor=3,
        use_rho=True,
    )
    print(f'Controls Values Loaded From CSV:\n{"-"*50}')
    controls.print_key_values_pairs()

    # controls.save_controls(Options.instance)
