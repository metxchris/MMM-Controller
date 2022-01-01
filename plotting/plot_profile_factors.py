# Standard Packages
import sys; sys.path.insert(0, '../')

# Local Packages
import modules.options
import plotting.modules.profiles as profiles
from modules.variables import InputVariables, OutputVariables
from modules.enums import SaveType, ProfileType


def main(runid, scan_num, input_scan_factor, save_types):
    '''
    Creates profile plots from data saved in CSVs.  The value closest to input_scan_factor in the saved
    scan_range will be used as the value for scan_factor

    Parameters:
    runid (str): The name of the CDF or run to reference
    scan_number (int): The number of the scan to reference
    scan_factor (float): The value of the scan factor file to reference
    save_types (list[Savetype]): The save types to plot profiles of
    '''

    options = modules.options.Options()
    options.load(runid, scan_num)
    scan_factor = options.find_scan_factor(input_scan_factor)

    input_vars = InputVariables(options)
    output_vars = OutputVariables(options)

    for save_type in save_types:
        if save_type == SaveType.INPUT:
            input_vars.load_from_csv(SaveType.INPUT, scan_factor)
            profiles.plot_profiles(ProfileType.INPUT, input_vars, scan_factor=scan_factor)
        elif save_type == SaveType.ADDITIONAL:
            input_vars.load_from_csv(SaveType.ADDITIONAL, scan_factor)
            profiles.plot_profiles(ProfileType.ADDITIONAL, input_vars, scan_factor=scan_factor)
        elif save_type == SaveType.OUTPUT:
            output_vars.load_from_csv(SaveType.OUTPUT, scan_factor)
            profiles.plot_profiles(ProfileType.OUTPUT, output_vars, scan_factor=scan_factor)


if __name__ == '__main__':
    # Runid and Scan Number (uncomment the line you wish to use)
    # runid, scan_num = '120982A09', 1
    runid, scan_num = '138536A01', 2

    # Scan Factor (var_to_scan will be read from the saved options file)
    input_scan_factor = 2.5

    save_types = [SaveType.INPUT]
    # save_types = [SaveType.INPUT, SaveType.ADDITIONAL, SaveType.OUTPUT]

    main(runid, scan_num, input_scan_factor, save_types)
