# Standard Packages
import os
import glob
from os.path import exists, dirname
import sys
sys.path.insert(0, '../')
# Local Packages
import pdftk, output, temp 

# Returns the path to the temp folder
def get_temp_path(file_name=''):
    return '{0}\\{1}'.format(dirname(temp.__file__), file_name)

# Returns the path to the pdftk executable
def get_pdftk_path():
    return '{0}\\pdftk.exe'.format(dirname(pdftk.__file__))

# Returns the path to the output folder
def get_output_path(file_name):
    return '{0}\\{1}'.format(dirname(output.__file__), file_name)

# checks if output dir exists and creates it if needed
def create_output_dir(input_options):
    output_dir = get_output_path(input_options.runid)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

# Returns original file_path if no duplicate files exist
# Otherwise appends (#) to the end of the file name to avoid overwritting a file
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
def open_output_pdf(file_path):
    os.startfile(file_path)

# Deletes any pdf from the temp folder
def clear_temp_folder():
    temp_files = get_temp_path('*.pdf')
    for file in glob.glob(temp_files):
        os.remove(file)

# Merge pdf sheets using pdftk in the temp folder into a single pdf and place in the output folder
def merge_input_profile_sheets(input_options):
    create_output_dir(input_options)

    merged_name = '{0}\\{1} Input Profiles.pdf'.format(input_options.runid, input_options.runid)
    output_file = check_filename(get_output_path(merged_name))
    temp_path = get_temp_path()
    pdftk_path = get_pdftk_path()
    
    # Shell command to use pdftk.exe
    os.system('cd {0} & {1} *.pdf cat output \"{2}\"'.format(temp_path, pdftk_path, output_file))

    open_output_pdf(output_file)

if __name__ == '__main__':
    clear_temp_folder()
