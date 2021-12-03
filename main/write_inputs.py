# Standard Packages
import sys
sys.path.insert(0, '../')
# Local Packages
from main import utils

MMM_LABELS = {
    'rmin'     :'! Half-width of the magnetic surface, r [m]',
    'rmaj'     :'! Major radius to geometric center of the magnetic surface, R [m]',
    'elong'    :'! Elongation profile',
    'ne'       :'! Electron density profile [m^-3]',
    'nh'       :'! Hydrogenic thermal particle density profile [m^-3]',
    'nz'       :'! Mean impurity density profile [m^-3]',
    'nf'       :'! Density from fast (non-thermal) ions profile [m^-3]',
    'zeff'     :'! Z_effective profile',
    'te'       :'! Electron temperature profile [keV]',
    'ti'       :'! Temperature of thermal ions profile [keV]',
    'q'        :'! Magnetic q profile',
    'btor'     :'! Toroidal magnetic field profile [Tesla]',
    'zimp'     :'! The profile of mean charge of impurities',
    'aimp'     :'! The profile of mean atomic mass of impurities',
    'ahyd'     :'! The profile of mean atomic mass of hydrogenic ions',
    'aimass'   :'! The profile of mean atomic mass of thermal ions',
    'wexbs'    :'! ExB shearing rate profile [rad/s]',
    'gne'      :'! Normalized electron density gradient profile',
    'gni'      :'! Normalized ion density gradient profile',
    'gnh'      :'! Normalized hydrogenic particle density gradient profile',
    'gnz'      :'! Normalized impurity density gradient profile',
    'gte'      :'! Normalized electron temperature gradient profile',
    'gti'      :'! Normalized ion temperature gradient profile',
    'gq'       :'! Normalized q gradient profile',
    'gvtor'    :'! Normalized toroidal velocity gradient profile',
    'vtor'     :'! Toroidal velocity profile [m/s]',
    'gvpol'    :'! Normalized poloidal velocity gradient profile',
    'vpol'     :'! Poloidal velocity profile [m/s]',
    'gvpar'    :'! Normalized parallel velocity gradient profile',
    'vpar'     :'! Parallel velocity profile [m/s]'
    }

# Header for MMM input file.  Add additional formatting options as needed
MMM_HEADER = '''&testmmm_input_control
 npoints = {npoints}    ! Number of radial points
 input_kind = 1
/
&testmmm_input_1stkind
! This is a sample input file of the first kind
! of an NSTX discharge

!.. Switches for component models
!   1D0 - ON, 0D0 - OFF
cmodel  =
   1D0     ! Weiland
   1D0     ! DRIBM
   1D0     ! ETG
   1D0     ! ETGM
   1D0     ! MTM  
    
!.. Weiland real options
cW20 =
   1D0     ! ExB shear coefficient
   1D0     ! Momentum pinch scaling factor
   0D0     ! Lower bound of electron thermal diffusivity
   1D2     ! Upper bound of electron thermal diffusivity
   0D0     ! Lower bound of ion thermal diffusivity
   1D2     ! Upper bound of ion thermal diffusivity

!.. DRIBM real options
cDBM =
   1D0     ! ExB shear coefficient
   0.1D0   ! kyrhos
   
!.. MTM real options
cMTM =
   0.2D0   ! ky/kx for MTM
   1.0D0   ! calibration factor
   

!.. ETG integer options
lETG =
   2       ! Jenko threshold
           ! applied to both electrostatic and electromagnetic regimes

!.. ETG real options
cETG =
   6D-2    ! CEES scale
   6D-2    ! CEEM scale
   
!.. ETGM integer options
lETGM =
   1      ! Collisionless limit

!.. ETGM real options
cETGM =
   0.0D0     ! ExB shear coefficient
   0.330D0   ! kyrhos
   0.250D0   ! kyrhoe
   
lprint   = 0      ! Verbose level\n\n'''

# Writes the input file used by the MMM driver
def write_input_file(input_vars, input_options):
    file_name = utils.get_temp_path('input')
    f = open(file_name, 'w')

    # Write mmm input file header
    f.write(MMM_HEADER.format(npoints=input_options.input_points))

    # Loop through mmm variables and write input file values
    for var_name in MMM_LABELS.keys():
        var = getattr(input_vars, var_name)

        # Write label and variable
        f.write('{0}\n'.format(MMM_LABELS[var_name]))
        f.write('{0} = \n'.format(var_name))

        # Write all values for variable
        values = var.values[:, input_options.time_idx]
        for value in values:
            # Writes values with 12 decimal places in exponential notation
            f.write('   {:.12e}\n'.format(value))
        f.write('\n')

    f.write('/\n') # Needed to end the input file

if __name__ == '__main__':
    pass
