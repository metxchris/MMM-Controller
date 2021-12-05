# Standard Packages
import os
import glob
import sys
sys.path.insert(0, '../')
# Local Packages
import pdftk, output, temp, cdfs

def get_cdf_path(file_name):
    '''Returns the path to the CDF folder (str)'''
    return f'{os.path.dirname(cdfs.__file__)}\\{file_name}.CDF'

def get_temp_path(file_name=''):
    '''Returns the path to the temp folder (str)'''
    return f'{os.path.dirname(temp.__file__)}\\{file_name}'

def get_output_path(file_name=''):
    '''Returns the path to the output folder (str)'''
    return f'{os.path.dirname(output.__file__)}\\{file_name}'

def get_pdftk_path():
    '''Returns the path to the pdftk executable (str)'''
    return f'{os.path.dirname(pdftk.__file__)}\\pdftk.exe'

def create_directory(dir_name):
    '''
    Checks if output dir exists and creates it if needed.

    Parameters:
    * dir_name (str): Path of directory
    '''

    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def check_filename(file_path):
    '''
    Checks if file exists and returns a new file path if the checked file exists.

    A number in the form of (#) is appended to the end of the file path in the event that
    the file already exists.  This is done so that files are not overwritten in an existing directory.
    An exception is raised if too many duplicate files already exist.

    TODO: this fails when the path starts with any '.'

    Parameters:
    * file_path (str): Path to file

    Returns:
    * file_path (str): Path to file that does not exist
    '''

    if os.path.exists(file_path):
        num_range = range(2, 1000)
        for i in num_range:
            path_split = file_path.split('.')
            new_file_path = f'{path_split[0]} ({i}).{path_split[1]}'
            if not os.path.exists(new_file_path):
                file_path = new_file_path
                break

        if i == max(num_range):
            raise NameError(f'Too many duplicate files exist to save {file_path}')

    return file_path

def check_dirname(dir_path):
    '''
    Checks if directory exists and returns a new directory name if the checked directory exists.

    A number in the form of (#) is appended to the end of the directory name in the event that
    the directory already exists.  This is done so that files are not overwritten in an existing directory.
    An exception is raised if too many duplicate directories already exist.

    Parameters:
    * dir_path (str): Path to directory

    Returns:
    * dir_path (str): Path to directory that does not exist yet
    '''

    if os.path.isdir(dir_path):
        num_range = range(2, 1000)
        for i in num_range:
            new_dir_path = f'{dir_path} ({i})'
            if not os.path.isdir(new_dir_path):
                dir_path = new_dir_path
                break

        if i == max(num_range):
            raise NameError(f'Too many duplicate directories exist to save directory {dir_path}')

    return dir_path    

def open_file(file_path):
    '''
    Opens the specified file (likely only works on Windows)

    Parameters:
    * file_path (str): Path of file to open
    '''

    os.startfile(file_path)

def clear_folder(dir_path, file_type):
    '''
    Clears all files in a specified folder of the specified file_type.

    Parameters:
    * dir_path (str): Path of directory
    * file_type (str): Type of file to clear (include * in the string)
    '''

    folder = f'{dir_path}\\{file_type}'
    for file in glob.glob(folder):
        os.remove(file)
    print(f'Cleared all files of type {file_type} from {dir_path}\n')

def clear_temp_folder():
    '''Clears temporary files from the temp folder.'''

    # Clear individual pdf sheets
    temp_files = get_temp_path('*.pdf')
    for file in glob.glob(temp_files):
        os.remove(file)

    # Clear input and output files for MMM Driver
    if os.path.exists(get_temp_path('input')):
        os.remove(get_temp_path('input'))
    if os.path.exists(get_temp_path('output')):
        os.remove(get_temp_path('output'))

def get_files_in_dir(dir_path, file_type=''):
    '''
    Lists all files in dir_path of file_type.

    Parameters:
    * dir_path (str): Path of directory
    * file_type (str): Type of file to search for (include * in the string)

    Returns:
    * file_names (list): List of file names
    '''

    files = f'{dir_path}\\{file_type}'
    file_names = [file for file in glob.glob(files)]

    if len(file_names) == 0:
        print(f'*** Warning: No files found for {files}')

    return file_names

def merge_profile_sheets(runid, profile_type):
    '''
    Merge PDF sheets using Pdftk in the temp folder into a single PDF, then place the merged PDF in the output folder.

    Pdftk is a 3rd party executable that is used to merge individual PDF sheets into one PDF, 
    and is called using a shell command.

    Parameters:
    * profile_type (str): The type of profile to merge (Input, Output, etc.)

    Returns:
    * output_file (str): Path to merged PDF
    '''

    create_directory(get_output_path(runid))

    merged_name = f'{runid}\\{runid} {profile_type} Profiles.pdf'
    output_file = check_filename(get_output_path(merged_name))
    temp_path = get_temp_path()
    pdftk_path = get_pdftk_path()
    
    # Shell command to use pdftk.exe 
    # TODO: Replace os.system with subprocess.run()
    os.system(f'cd {temp_path} & {pdftk_path} *{profile_type}*.pdf cat output \"{output_file}\"')
    print(f'Profiles saved to {output_file}\n')

    return output_file

# For testing purposes
if __name__ == '__main__':

    clear_temp_folder()
