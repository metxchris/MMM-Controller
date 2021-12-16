# Standard Packages
import sys
sys.path.insert(0, '../')

# Local Packages
from main import utils
from main.enums import ShotType
from main.options import Options


MMM_LABELS = {
    'rmin': '! Half-width of the magnetic surface, r [m]',
    'rmaj': '! Major radius to geometric center of the magnetic surface, R [m]',
    'elong': '! Elongation profile',
    'ne': '! Electron density profile [m^-3]',
    'nh': '! Hydrogenic thermal particle density profile [m^-3]',
    'nz': '! Mean impurity density profile [m^-3]',
    'nf': '! Density from fast (non-thermal) ions profile [m^-3]',
    'zeff': '! Z_effective profile',
    'te': '! Electron temperature profile [keV]',
    'ti': '! Temperature of thermal ions profile [keV]',
    'q': '! Magnetic q profile',
    'btor': '! Toroidal magnetic field profile [Tesla]',
    'zimp': '! The profile of mean charge of impurities',
    'aimp': '! The profile of mean atomic mass of impurities',
    'ahyd': '! The profile of mean atomic mass of hydrogenic ions',
    'aimass': '! The profile of mean atomic mass of thermal ions',
    'wexbs': '! ExB shearing rate profile [rad/s]',
    'gne': '! Normalized electron density gradient profile',
    'gni': '! Normalized ion density gradient profile',
    'gnh': '! Normalized hydrogenic particle density gradient profile',
    'gnz': '! Normalized impurity density gradient profile',
    'gte': '! Normalized electron temperature gradient profile',
    'gti': '! Normalized ion temperature gradient profile',
    'gq': '! Normalized q gradient profile',
    'gvtor': '! Normalized toroidal velocity gradient profile',
    'vtor': '! Toroidal velocity profile [m/s]',
    'gvpol': '! Normalized poloidal velocity gradient profile',
    'vpol': '! Poloidal velocity profile [m/s]',
    'gvpar': '! Normalized parallel velocity gradient profile',
    'vpar': '! Parallel velocity profile [m/s]'
}


# Writes the input file used by the MMM driver
def write_input_file(input_vars, controls):
    input_options = Options.instance
    file_name = utils.get_temp_path('input')
    f = open(file_name, 'w')

    # Write MMM input file header
    f.write(controls.get_mmm_header())

    # Loop through MMM variables and write input file values
    for var_name in MMM_LABELS.keys():
        var = getattr(input_vars, var_name)

        # Write label and variable
        f.write(f'{MMM_LABELS[var_name]}\n')
        f.write(f'{var_name} = \n')

        # Write all values for variable
        values = var.values[:, input_options.time_idx]
        for value in values:
            # Writes values with 8 decimal places in exponential notation
            f.write('   {:.8e}\n'.format(value))
        f.write('\n')

    f.write('/\n')  # Needed to end the input file


if __name__ == '__main__':
    ...
