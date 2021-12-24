# Standard Packages
import sys; sys.path.insert(0, '../')

# Local Packages
import plotting.modules.profiles as profiles
import modules.options as options
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
    save_types (list of Savetype): The save types to plot profiles of
    '''

    options.instance.load_options(runid, scan_num)
    scan_factor = options.instance.find_scan_factor(input_scan_factor)
    args = (runid, scan_num, options.instance.var_to_scan, scan_factor)

    input_vars = InputVariables()
    output_vars = OutputVariables()

    for save_type in save_types:
        if save_type == SaveType.INPUT:
            input_vars.load_from_csv(SaveType.INPUT, *args)
            profiles.plot_profiles(ProfileType.INPUT, input_vars, scan_factor=scan_factor)
        elif save_type == SaveType.ADDITIONAL:
            input_vars.load_from_csv(SaveType.ADDITIONAL, *args)
            profiles.plot_profiles(ProfileType.ADDITIONAL, input_vars, scan_factor=scan_factor)
        elif save_type == SaveType.OUTPUT:
            output_vars.load_from_csv(SaveType.OUTPUT, *args)
            profiles.plot_profiles(ProfileType.OUTPUT, output_vars, scan_factor=scan_factor)


if __name__ == '__main__':
    # Runid and Scan Number (uncomment the line you wish to use)
    # runid, scan_num = '120982A09', 1
    runid, scan_num = 'TEST', 180

    # Scan Factor (var_to_scan will be read from the saved options file)
    input_scan_factor = 2.5

    save_types = [SaveType.INPUT]
    # save_types = [SaveType.INPUT, SaveType.ADDITIONAL, SaveType.OUTPUT]

    main(runid, scan_num, input_scan_factor, save_types)
