# Standard Packages
import os
import glob
from os.path import exists, dirname
import sys
sys.path.insert(0, '../')
# Local Packages
import pdftk, output, temp, cdfs

# Returns the path to the CDF folder
def get_cdf_path(file_name):
    return "{0}\\{1}.CDF".format(dirname(cdfs.__file__), file_name)

# Returns the path to the temp folder
def get_temp_path(file_name=''):
    return '{0}\\{1}'.format(dirname(temp.__file__), file_name)

# Returns the path to the output folder
def get_output_path(file_name=''):
    return '{0}\\{1}'.format(dirname(output.__file__), file_name)

# Returns the path to the output folder
def get_mmm_path(file_name=''):
    return '{0}\\{1}'.format(dirname(mmm.__file__), file_name)

# Returns the path to the pdftk executable
def get_pdftk_path():
    return '{0}\\pdftk.exe'.format(dirname(pdftk.__file__))

# checks if output dir exists and creates it if needed
def create_directory(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

# Returns original file_path if no duplicate files exist
# Otherwise appends (#) to the end of the file name to avoid overwritting a file
# TODO: this fails when the path starts with any '.'
def check_filename(file_path):
    if not os.path.exists(file_path):
        return file_path

    for i in range(2, 1000):
        path_split = file_path.split('.')
        new_file_path = '{0} ({1}).{2}'.format(path_split[0], i, path_split[1])
        if not os.path.exists(new_file_path):
            return new_file_path

    # Throw an exception if this many duplicate files exist
    raise NameError('Too many duplicate files exist to save {0}'.format(file_path))

# Opens the output pdf (likely only works on Windows)
def open_file(file_path):
    os.startfile(file_path)

# Clears temporary files from the temp folder
def clear_temp_folder():
    # Clear individual pdf sheets
    temp_files = get_temp_path('*.pdf')
    for file in glob.glob(temp_files):
        os.remove(file)

    # Clear input and output files for MMM Driver
    if os.path.exists(get_temp_path('input')):
        os.remove(get_temp_path('input'))
    if os.path.exists(get_temp_path('output')):
        os.remove(get_temp_path('output'))

# Merge pdf sheets using pdftk in the temp folder into a single pdf and place in the output folder
def merge_input_profile_sheets(input_options):
    create_directory(get_output_path(input_options.runid))

    merged_name = '{0}\\{1} Input Profiles.pdf'.format(input_options.runid, input_options.runid)
    output_file = check_filename(get_output_path(merged_name))
    temp_path = get_temp_path()
    pdftk_path = get_pdftk_path()
    
    # Shell command to use pdftk.exe
    os.system('cd {0} & {1} *.pdf cat output \"{2}\"'.format(temp_path, pdftk_path, output_file))

    return output_file

if __name__ == '__main__':
    clear_temp_folder()
