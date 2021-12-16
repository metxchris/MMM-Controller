# Format strings
SCAN_FACTOR_FMT_STR = '{:0>7.3f}'
RHO_VALUE_FMT_STR = '{:.3f}'
SHEET_NUM_FMT_STR = '{:03d}'
TIME_FMT_STR = '{:.3f}'

# Unit string pairs for creating labels.  First item is the search string, second item is the replacement string
UNIT_STRINGS = [
    ['N/M**3', r'$\mathrm{m}^{-3}$'],
    ['M**2/SEC', r'$\mathrm{m}^{2}/s$'],
    ['M/SEC', r'm/s'],
    ['M', r'm'],
    ['SEC**-1', r's$^{-1}$'],
    ['MAMPS', r'MA'],
    ['RAD/SEC', r'rad/s'],
    ['PASCALS', r'Pa'],
    ['SECONDS', r's'],
    ['TESLA', r'T'],
    ['TESLA*M', r'T\,m'],
    ['EV', r'eV'],
    ['kEV', r'keV'],
    ['m/s^2', r'm/s$^2$'],
    ['m^2/s', r'm$^2$/s'],
    ['s^-1', r's$^{-1}$'],
]

# Physical Constants
PI = 3.1415926535
ZCE = 1.602176565 * 10**(-19)                                 # Electron charge                [Coulomb]
ZCEPS0 = 8.854187817 * 10**(-12)                              # Vacuum electrical permittivity [farads/m]
ZCMU0 = 4 * PI * 10**(-7)                                     # Vacuum magnetic permeability   [henrys/m]
ZCME = 9.10938215 * 10**(-31)                                 # Electron mass                  [kg]
ZCMP = 1.672621777 * 10**(-27)                                # Proton mass                    [kg]
ZCKB = ZCE * 10**3                                            # Energy conversion factor       [J/keV]
ZCF = ((4 * PI**(1 / 2) / 3) * (ZCE / (4 * PI * ZCEPS0))**2   # Collision frequency factor     TODO: units?
       * (ZCE / ZCKB) * (ZCE / ZCME * ZCE / ZCKB)**(1 / 2))
