"""Stores options needed to run modules in MMM Explorer

Some options are set at runtime by the user, and others may be set while
various code executes.  Since options values are needed in the majority of
modules used by MMM Explorer, the Options class is meant to be accessed by
the public instance of it that's created at the bottom of the module.  This
lets us reference options values without having to pass the object around.

Options are saved using the pickle module, which serializes the entire object
in binary, and is unreadable by the user.  Consequently, we also save a CSV
of all options values stored in the pickle file, which allows the user to see
what options values were stored for each scan; this CSV is otherwise not used
by MMM Explorer.

Example Usage:
    import modules.options as options

    # Set options values
    options.instance.set(
        runid='120968A02',
        shot_type=ShotType.NSTX,
        input_time=0.5,
        input_points=51,
    )

    # Save options (to pickle file)
    options.instance.save()

    # Load options (from pickle file)
    options.instance.load('120968A02', 1)
"""

# Standard Packages
import pickle
import inspect

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.utils as utils
import modules.datahelper as datahelper
import modules.constants as constants
from modules.enums import ShotType, ScanType


class Options:
    '''
    Stores options for MMM Controller

    Properties and Members:
    * apply_smoothing (bool): kill-switch to disable smoothing of all variables
    * input_points (int): the amount of radial points each variable is interpolated to when sent to MMM
    * input_time (float): the time to check the CDF for values
    * runid (str): the Runid in the CDF, usually also the name of the CDF
    * scan_factor_str (str): the string of the scan factor, rounded for better visual presentation
    * scan_num (int): the number identifying where data is stored within the ./output/runid/ directory
    * scan_range (np.ndarray): the range of factors to multiply the var_to_scan by
    * scan_type (ScanType): the type of the scan
    * shot_type (ShotType): the shot type of the CDF
    * temperature_profiles (bool): replace temperature variables with experimental profiles
    * time_str (str): the string of the measurement time, rounded for better visual presentation
    * time_idx (int): the index of the CDF time value that is closest to input_time
    * uniform_rho (bool): interpolate input variables onto a rho of uniform spacing
    * var_to_scan (str): the variable to scan; syntax must match a member of InputVariables
    '''

    # Private members (each has a property)
    _input_points = None
    _runid = None
    _scan_range = None
    _var_to_scan = None

    # Public members
    apply_smoothing = False
    input_time = None
    scan_num = None
    scan_type = ScanType.NONE
    shot_type = ShotType.NONE
    temperature_profiles = False
    time_str = None
    time_idx = None
    uniform_rho = False

    # Properties
    @property
    def input_points(self):
        return self._input_points

    @input_points.setter
    def input_points(self, points):
        if points:
            self._input_points = max(points, 5)

    @property
    def runid(self):
        return self._runid

    @runid.setter
    def runid(self, runid):
        if not isinstance(runid, str):
            runid = str(runid)
        self._runid = runid.strip()

    @property
    def scan_range(self):
        return self._scan_range

    @scan_range.setter
    def scan_range(self, scan_range):
        if scan_range is not None:
            if not isinstance(scan_range, np.ndarray):
                raise TypeError(f'scan_range must be {np.ndarray} or {None} and not {type(scan_range)}')
            too_small = np.absolute(scan_range) < constants.ABSMIN_SCAN_FACTOR_VALUE
            if too_small.any():
                value_signs = np.sign(scan_range[too_small])
                value_signs[value_signs == 0] = 1  # np.sign(0) = 0, so set these to +1
                scan_range[too_small] = constants.ABSMIN_SCAN_FACTOR_VALUE * value_signs
        self._scan_range = scan_range

    @property
    def var_to_scan(self):
        return self._var_to_scan

    @var_to_scan.setter
    def var_to_scan(self, var_to_scan):
        self._var_to_scan = var_to_scan
        self.scan_type = datahelper.get_scan_type(var_to_scan)

    # Methods
    def get_keys(self):
        '''Returns (list): Names of properties and public members'''
        kvps = self.get_key_value_pairs()
        return [kvp[0] for kvp in kvps]

    def get_key_value_pairs(self):
        '''Returns (list): Names and values of properties and public members'''
        attributes = inspect.getmembers(self, lambda m: not inspect.isroutine(m))
        return [a for a in attributes if not a[0].startswith('_')]

    def set(self, **kwargs):
        '''Sets specified options values'''
        for key, value in kwargs.items():
            setattr(self, key, value)

    def load(self, runid, scan_num):
        '''
        Loads Options object from a pickle file

        Parameters:
        * runid (str): The run id to load options from
        * scan_num (int): The scan number to load options from
        '''

        pickle_path = utils.get_options_path(runid, scan_num)

        with open(pickle_path, 'rb') as handle:
            loaded_options = pickle.load(handle)

        # Setting options values one-by-one will not break any existing references to Options
        options_to_set = loaded_options.get_keys()
        for option in options_to_set:
            setattr(self, option, getattr(loaded_options, option))

    def save(self):
        '''
        Saves Options object to a pickle file

        A CSV is also saved to make saved options viewable by the user, but is otherwise not used
        '''

        pickle_path = utils.get_options_path(self.runid, self.scan_num)

        with open(pickle_path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Options are also written to a CSV to make viewing their values easier
        csv_path = pickle_path.replace('.pickle', '.csv')
        options_values = self.get_key_value_pairs()
        with open(csv_path, 'w') as handle:
            for o in options_values:
                handle.write(f'{o[0]}, {o[1]}\n')

        print(f'Options saved to {pickle_path}\n')

    def set_measurement_time(self, tvar):
        '''Find the index of the measurement time closest to the input_time, then store that value and its index'''
        self.time_idx = np.argmin(np.abs(tvar.values - self.input_time))
        self.time_str = f'{tvar.values[self.time_idx]:{constants.TIME_VALUE_FMT}}'

    def find_scan_factor(self, scan_factor):
        '''Returns (float or None): Value in scan_range closest to the specified scan_factor'''
        if scan_factor is None:
            return scan_factor
        elif self.scan_range is None:
            raise ValueError('Cannot find scan_factor value when scan_range is None')
        return self.scan_range[np.argmin(np.abs(self.scan_range - scan_factor))]


# Store a public instance of the Options class
instance = Options()
