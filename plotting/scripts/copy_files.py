#!/usr/bin/python3

"""Copies figures or figure data from the plotting folder to another folder

"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')

# Local Packages
import modules.utils as utils


def copy_files(target_folder, source_folder, target_type, source_type):
    target_files = utils.get_files_in_dir(target_folder, f'*{target_type}')
    target_filenames = [name.split("\\")[-1].replace(target_type, '') for name in target_files]
    source_files = [f'{source_folder}\\{name}{source_type}' for name in target_filenames]

    for target, source in zip(target_files, source_files):
        if utils.check_exists(source):
            utils.copy_file(source, target)
        else:
            print(f'\nfile not found: {source}')


if __name__ == '__main__':
    runid = '138536A01'
    target_type = '.pdf'
    source_type = '.pdf'

    target_base = 'C:\\Users\\metxc\\Documents\\MMM Papers\\My Presentations\\ETGM Paper\\Figures'
    source_base = 'C:\\Users\\metxc\\Documents\\MMM-Package\\plotting\\output'

    copy_files(f'{target_base}\\Singles', f'{source_base}\\singles\\{runid}\\figures', target_type, source_type)
    copy_files(f'{target_base}\\Outputs', f'{source_base}\\singles\\{runid}\\figures', target_type, source_type)
    copy_files(f'{target_base}\\Contours', f'{source_base}\\contours\\{runid}\\figures', target_type, source_type)
