# Standard Packages
import sys; sys.path.insert(0, '../')

# Local Packages
import plotting.modules.profiles as profiles
from main.variables import InputVariables, OutputVariables
from main.enums import SaveType, ProfileType
from main.options import Options


def main(runid, scan_num, input_scan_factor):
    '''
    Creates profile plots from data saved in CSVs.  The value closest to input_scan_factor in the saved
    scan_range will be used as the value for scan_factor

    Parameters:
    runid (str): The name of the CDF or run to reference
    scan_number (int): The number of the scan to reference
    scan_factor (float): The value of the scan factor file to reference
    '''

    Options.instance.load_options(runid, scan_num)
    scan_factor = Options.instance.find_scan_factor(input_scan_factor)
    args = (runid, scan_num, Options.instance.var_to_scan, scan_factor)

    input_vars = InputVariables()
    input_vars.load_from_csv(SaveType.INPUT, *args)
    input_vars.load_from_csv(SaveType.ADDITIONAL, *args)

    output_vars = OutputVariables()
    output_vars.load_from_csv(SaveType.OUTPUT, *args)

    profiles.plot_profiles(ProfileType.INPUT, input_vars, scan_factor=scan_factor)
    profiles.plot_profiles(ProfileType.ADDITIONAL, input_vars, scan_factor=scan_factor)
    profiles.plot_profiles(ProfileType.OUTPUT, output_vars, scan_factor=scan_factor)


if __name__ == '__main__':
    # Runid and Scan Number (uncomment the line you wish to use)
    # runid, scan_num = '120982A09', 1
    runid, scan_num = 'TEST', 45

    # Scan Factor (var_to_scan will be read from the saved options file)
    input_scan_factor = 2.55

    main(runid, scan_num, input_scan_factor)
