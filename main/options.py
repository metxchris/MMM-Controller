# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import pickle

# Local Packages
from main import utils
from main.enums import ShotType
from main.variables import InputVariables


class _Options:
    '''
    Options for MMM Controller (private class)

    Properties:
    * apply_smoothing (bool): killswitch to disable smoothing of all variables
    * input_points (int): the amount of radial points each variable is interpolated to when sent to MMM
    * input_time (float): the time to check the CDF for values
    * reject_outliers (bool): replaces outliers with values of 0
    * runid (str): the Runid in the CDF, usually also the name of the CDF
    * scan_factor_str (str): the string of the scan factor, rounded for better visual presentation
    * scan_range (np.ndarray): the range of factors to multiply the var_to_scan by
    * shot_type (ShotType): the shot type of the CDF
    * temperature_profiles (bool): replace temperature variables with experimental profiles
    * time_str (str): the string of the measurement time, rounded for better visual presentation
    * time_idx (int): the index of the CDF time value that is closest to input_time
    * uniform_rho (bool): interpolate input variables onto a rho of uniform spacing
    * var_to_scan (str): the variable to scan; syntax must match a member of InputVariables
    '''

    # Private members (each has a property)
    _apply_smoothing = False
    _input_points = None
    _input_time = None
    _reject_outliers = False
    _runid = None
    _scan_factor_str = None
    _scan_range = None
    _shot_type = ShotType.NONE
    _temperature_profiles = False
    _time_str = None
    _time_idx = None
    _uniform_rho = False
    _var_to_scan = None

    # Properties
    @property
    def apply_smoothing(self):
        return self._apply_smoothing
    @apply_smoothing.setter
    def apply_smoothing(self, apply_smoothing):
        if type(apply_smoothing) is bool:
            self._apply_smoothing = apply_smoothing
        else:
            raise TypeError(f'apply_smoothing must be a bool and not {type(apply_smoothing)}')

    @property
    def input_points(self):
        return self._input_points
    @input_points.setter
    def input_points(self, points):
        if points is not None:
            self._input_points = max(points, 5)

    @property
    def input_time(self):
        return self._input_time
    @input_time.setter
    def input_time(self, input_time):
        self._input_time = input_time

    @property
    def reject_outliers(self):
        return self._reject_outliers
    @reject_outliers.setter
    def reject_outlilers(self, reject_outlilers):
        self._reject_outliers = reject_outlilers

    @property
    def runid(self):
        return self._runid
    @runid.setter
    def runid(self, runid):
        if type(runid) is str:
            self._runid = runid.strip()

    @property
    def scan_factor_str(self):
        return self._scan_factor_str
    @scan_factor_str.setter
    def scan_factor_str(self, scan_factor_str):
        if scan_factor_str is not None:
            self._scan_factor_str = '{:.3f}'.format(scan_factor_str)

    @property
    def scan_range(self):
        return self._scan_range
    @scan_range.setter
    def scan_range(self, scan_range):
        if type(scan_range) is not np.ndarray and scan_range is not None:
            raise TypeError(f'scan_range must be type np.ndarray and not {type(scan_range)}')
        self._scan_range = scan_range

    @property
    def shot_type(self):
        return self._shot_type
    @shot_type.setter
    def shot_type(self, shot_type):
        self._shot_type = shot_type

    @property
    def temperature_profiles(self):
        return self._temperature_profiles
    @temperature_profiles.setter
    def temperature_profiles(self, temperature_profiles):
        self._temperature_profiles = temperature_profiles

    @property
    def time_idx(self):
        return self._time_idx
    @time_idx.setter
    def time_idx(self, idx):
        self._time_idx = idx

    @property
    def time_str(self):
        return self._time_str
    @time_str.setter
    def time_str(self, time_str):
        self._time_str = time_str

    @property
    def uniform_rho(self):
        return self._uniform_rho
    @uniform_rho.setter
    def uniform_rho(self, uniform_rho):
        self._uniform_rho = uniform_rho

    @property
    def var_to_scan(self):
        return self._var_to_scan
    @var_to_scan.setter
    def var_to_scan(self, var_to_scan):
        if var_to_scan is None or hasattr(InputVariables(), var_to_scan):
            self._var_to_scan = var_to_scan
        else:
            raise ValueError(f'Variable {var_to_scan} is not defined under InputVariables')

    # Methods
    def get_options(self):
        return [o[1:] for o in dir(self) if not callable(getattr(self, o)) and o.startswith("_") and not o.startswith("__")]

    def get_options_values(self):
        options = self.get_options()
        return [str(o) + ',' + str(getattr(self, o)).replace('\n', '') for o in options]

    def set_options(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f'Error: Options does not have attribute {key}')

    def load_options(self, file_path):
        '''Loads _Options object from a pickle file'''

        pickle_path = file_path.replace('.csv', '.pickle')
        with open(pickle_path, 'rb') as handle:
            loaded_options = pickle.load(handle)

        # Setting options values one-by-one will not break any existing references to Options
        options_to_set = loaded_options.get_options()
        for o in options_to_set:
            option_value = getattr(loaded_options, o)
            setattr(self, o, option_value)

    def save_options(self, file_path, csv_path=''):
        '''Saves _Options object to a pickle file, as well as a CSV'''

        pickle_path = file_path.replace('.csv', '.pickle')
        with open(file_path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Options can also be written to a CSV to make viewing their values easier
        if csv_path != '':
            options_values = self.get_options_values()
            f = open(csv_path, 'w')
            for option in options_values:
                f.write(f'{option}\n')

    def set_measurement_time(self, tvar):
        '''Find the index of the measurement time closest to the input_time, then store that value and its index'''
        self._time_idx = np.argmin(np.abs(tvar.values - self.input_time))
        self._time_str = "{:.3f}".format(tvar.values[self.time_idx])


class Options:
    '''Stores a public instance of the _Options class (Singleton pattern)'''

    instance = _Options()


if __name__ == '__main__':
    ...
