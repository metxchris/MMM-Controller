# Values for Formatting
TIME_VALUE_SIGFIGS = 3

# Number Formatters
INPUT_VARIABLE_VALUE_FMT = '.8e'
INPUT_CONTROL_VALUE_FMT = '.8f'
SCAN_FACTOR_DISPLAY_FMT = '.3g'
TIME_DISPLAY_FMT = f'.{TIME_VALUE_SIGFIGS}g'
SCAN_FACTOR_FMT = '0>7.3f'
TIME_VALUE_FMT = f'.{TIME_VALUE_SIGFIGS}f'
RHO_VALUE_FMT = '.3f'
SHEET_NUM_FMT = '03d'

# Strings
SCAN_FACTOR_VALUE_SEPARATOR = ' x= '
RHO_VALUE_SEPARATOR = ' = '

# Values
ABSMIN_SCAN_FACTOR_VALUE = 1e-12
MAX_GRADIENT = 1e3

# Misc
INTERP_TYPE = 'quadratic'  # scipy.interpolate.interp1d(kind=INTERP_TYPE)
                           # linear, nearest, nearest-up, zero, slinear, quadratic, cubic, previous, next.

# Physical Constants
PI = 3.1415926535
ZCC = 299792458                                               # Speed of light                 [m/s]
ZCE = 1.602176565 * 10**(-19)                                 # Electron charge                [Coulomb]
ZCEPS0 = 8.854187817 * 10**(-12)                              # Vacuum electrical permittivity [farads/m]
ZCMU0 = 4 * PI * 10**(-7)                                     # Vacuum magnetic permeability   [henrys/m]
ZCME = 9.10938215 * 10**(-31)                                 # Electron mass                  [kg]
ZCMP = 1.672621777 * 10**(-27)                                # Proton mass                    [kg]
ZCKB = ZCE * 10**3                                            # Energy conversion factor       [J/keV]
ZCF = ((4 * PI**(1 / 2) / 3) * (ZCE / (4 * PI * ZCEPS0))**2   # Collision frequency factor     TODO: units?
       * (ZCE / ZCKB) * (ZCE / ZCME * ZCE / ZCKB)**(1 / 2))
