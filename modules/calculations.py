"""Calculates variable values

Variable Classifications:
* Base variable: A non-gradient variable that is either used as input to MMM,
  or is needed to calculate an MMM input variable.  Base variables do not
  depend on values of gradient or additional variables.
* Gradient variable: A variable that is a normalized gradient, and may be used
  as an input to MMM.  These variables depend on values of base variables,
  but not on values of additional variables.
* Additional variable: A variable that is not needed for MMM input nor input
  calculations. Additional variables can depend on both base and gradient
  values.

Calculations for base variables try to match the definition used by TRANSP,
whenever possible.  In some cases, perfect matches cannot be obtained, since
TRANSP doesn't save every variable we need for calculations.  For example,
our vpar calculation is just an approximation for the TRANSP value of vpar,
although the difference should be fairly negligible.

Calculations for additional variables come directly from the MMM source files,
and these sources will be noted as necessary in the calculation functions
below.  The names of variables used here will either directly match or
closely resemble the names of variables in MMM.
"""

# Standard Packages
import sys
import copy
import inspect

# 3rd Party Packages
import numpy as np
from scipy.interpolate import interp1d  # TODO: use Akima1DInterpolator?

# Local Packages
import modules.options as options
import modules.constants as constants


def ahyd(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values

    ahyd = (nh0 + 2 * nd) / (nh0 + nd)

    calc_vars.ahyd.set(values=ahyd, units='')


def aimass(calc_vars):
    '''Mean Atomic Mass of Thermal Ions'''
    ahyd = calc_vars.ahyd.values
    aimp = calc_vars.aimp.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values

    aimass = (ahyd * nh + aimp * nz) / (nh + nz)

    calc_vars.aimass.set(values=aimass, units='')


def alphamhd(calc_vars):
    '''Alpha MHD (Weiland Definition)'''
    betae = calc_vars.betae.values
    gne = calc_vars.gne.values
    gni = calc_vars.gni.values
    gte = calc_vars.gte.values
    gti = calc_vars.gti.values
    q = calc_vars.q.values
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    alphamhd = q**2 * betae * (gne + gte + ti / te * (gni + gti))

    calc_vars.alphamhd.set(values=alphamhd, units='')


def beta(calc_vars):
    '''Beta'''
    zcmu0 = constants.ZCMU0
    btor = calc_vars.btor.values
    p = calc_vars.p.values

    beta = 2 * zcmu0 * p / btor**2

    calc_vars.beta.set(values=beta, units='')


def betae(calc_vars):
    '''Electron Beta'''
    zckb = constants.ZCKB
    zcmu0 = constants.ZCMU0
    btor = calc_vars.btor.values
    ne = calc_vars.ne.values
    te = calc_vars.te.values

    betae = 2 * zcmu0 * ne * te * zckb / btor**2

    calc_vars.betae.set(values=betae, units='')


def bpol(calc_vars):
    '''Poloidal Magnetic Field'''
    btor = calc_vars.btor.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values

    bpol = rmin / rmaj * btor / q

    calc_vars.bpol.set(values=bpol, units=calc_vars.btor.units)


def btor(calc_vars):
    '''Toroidal Magnetic Field'''
    bz = calc_vars.bz.values
    raxis = calc_vars.rmaj.values[0, :]
    rmaj = calc_vars.rmaj.values

    btor = raxis / rmaj * bz

    calc_vars.btor.set(values=btor, units=calc_vars.bz.units)


def eps(calc_vars):
    '''Inverse Aspect Ratio'''
    arat = calc_vars.arat.values

    eps = 1 / arat

    calc_vars.eps.set(values=eps, units='')


def etae(calc_vars):
    '''etae = gte / gne'''
    gte = calc_vars.gte.values
    gne = calc_vars.gne.values

    etae = gte / gne

    calc_vars.etae.set(values=etae, units='')


def etai(calc_vars):
    '''
    etae = gti / gni

    Note that TRANSP appears to use an entirely different definition of ni
    when calculating gni, than it uses for the values of ni itself.  As such,
    our calculation of etai will generally not match with TRANSP values.
    '''

    gti = calc_vars.gti.values
    gni = calc_vars.gni.values

    etai = gti / gni

    calc_vars.etai.set(values=etai, units='')


def gmax(calc_vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    eps = calc_vars.eps.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    gyrfi = calc_vars.gyrfi.values
    vthi = calc_vars.vthi.values

    gmax = rmaj / (vthi / gyrfi * q / eps)

    calc_vars.gmax.set(values=gmax, units='')


def gyrfi(calc_vars):
    '''Ion Gyrofrequency'''
    zce = constants.ZCE
    zcmp = constants.ZCMP
    aimass = calc_vars.aimass.values
    btor = calc_vars.btor.values

    gyrfi = zce * btor / (zcmp * aimass)

    calc_vars.gyrfi.set(values=gyrfi, units='s^-1')


def gave(calc_vars):
    '''Average Magnetic Surface Curvature'''
    shear = calc_vars.shear.values
    alphamhd = calc_vars.alphamhd.values

    gave = 2 / 3 + 5 / 9 * shear - 5 / 12 * alphamhd

    calc_vars.gave.set(values=gave, units='')


def loge(calc_vars):
    '''Electron Coulomb Logarithm'''

    # TODO: Need to add equations for different TE ranges
    ne = calc_vars.ne.values
    te = calc_vars.te.values

    # NRL Plasma Formulary Definition
    loge = 37.8 - np.log(ne**(1 / 2) / te)

    # TRANSP definition (equivalent)
    # zeff = calc_vars.zeff.values
    # loge = 39.23 - np.log(zeff*ne**(1 / 2) / te)

    calc_vars.loge.set(values=loge, units='')


def p(calc_vars):
    '''Plasma Pressure'''
    zckb = constants.ZCKB
    ne = calc_vars.ne.values
    ni = calc_vars.ni.values
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    p = (ne * te + ni * ti) * zckb

    calc_vars.p.set(values=p, units='Pa')


def nh0(calc_vars):
    '''Hydrogen Ion Density'''
    nd = calc_vars.nd.values
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nz = calc_vars.nz.values
    zimp = calc_vars.zimp.values

    nh0 = ne - zimp * nz - nf - nd

    calc_vars.nh0.set(values=nh0, units=calc_vars.ne.units)


def nh(calc_vars):
    '''Total Hydrogenic Ion Density'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values

    nh = nh0 + nd

    calc_vars.nh.set(values=nh, units=calc_vars.nd.units)


def ni(calc_vars):
    '''Thermal Ion Density'''
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values

    # TRANSP Definition
    ni = nh + nz

    calc_vars.ni.set(values=ni, units=calc_vars.ne.units)


def nuei(calc_vars):
    '''Collision Frequency'''
    zcf = constants.ZCF
    ne = calc_vars.ne.values
    te = calc_vars.te.values
    zeff = calc_vars.zeff.values
    loge = calc_vars.loge.values

    nuei = zcf * 2**(1 / 2) * ne * loge * zeff / te**(3 / 2)

    calc_vars.nuei.set(values=nuei, units='s^-1')


def nuei2(calc_vars):
    '''OLD NOTE: Not sure what to call this, but it leads to the approx the correct NUSTI'''
    zcf = constants.ZCF
    ni = calc_vars.ni.values
    ti = calc_vars.ti.values
    zeff = calc_vars.zeff.values
    loge = calc_vars.loge.values

    nuei2 = zcf * 2**(1 / 2) * ni * loge * zeff / ti**(3 / 2)

    calc_vars.nuei2.set(values=nuei2, units='s^-1')


def nuste(calc_vars):
    '''
    Electron Collisionality

    This is in approximate agreement with NUSTE in TRANSP.  One source of the
    disagreement is likely because the modmmm7_1.f90 Coulomb logarithm
    (loge) does not match perfectly with the TRANSP version (CLOGE).  However,
    our Coulomb logarithm definition follows from the NRL plasma formulary,
    and we feel that it's use is correct here.
    '''

    eps = calc_vars.eps.values
    nuei = calc_vars.nuei.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    vthe = calc_vars.vthe.values

    nuste = nuei * eps**(-3 / 2) * q * rmaj / vthe

    calc_vars.nuste.set(values=nuste, units='')


def nusti(calc_vars):
    '''
    Ion Collisionality

    OLD NOTE: This is approximately correct, but agreement is also somewhat
    time-dependent.  We likely need to use the Coulomb logarithm for ions in
    nuei instead of the Coulomb logarithm for electrons.
    '''

    zcme = constants.ZCME
    zcmp = constants.ZCMP
    eps = calc_vars.eps.values
    q = calc_vars.q.values
    nuei2 = calc_vars.nuei2.values
    rmaj = calc_vars.rmaj.values
    vthi = calc_vars.vthi.values

    nusti = nuei2 * eps**(-3 / 2) * q * rmaj / (2 * vthi) * (zcme / zcmp)**(1 / 2)

    calc_vars.nusti.set(values=nusti, units='')


def shat(calc_vars):
    '''
    Effective Magnetic Shear

    TRANSP uses a different definition for shat than we use here, so we expect
    our values to be similar, but not a direct match.
    '''

    elong = calc_vars.elong.values
    shear = calc_vars.shear.values

    shat = (2 * shear - 1 + (elong * (shear - 1))**2)**(1 / 2)
    shat[shat < 0] = 0

    calc_vars.shat.set(values=shat, units='')


def shear(calc_vars):
    '''Magnetic Shear'''
    gq = calc_vars.gq.values
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values

    shear = gq * rmin / rmaj

    calc_vars.shear.set(values=shear, units='')


def tau(calc_vars):
    '''Temperature Ratio te / ti'''
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    tau = te / ti

    calc_vars.tau.set(values=tau, units='')


def vthe(calc_vars):
    '''Thermal Velocity of Electrons'''
    zckb = constants.ZCKB
    zcme = constants.ZCME
    te = calc_vars.te.values

    vthe = (2 * zckb * te / zcme)**(1 / 2)

    calc_vars.vthe.set(values=vthe, units='m/s')


def vthi(calc_vars):
    '''Thermal Velocity of Ions'''
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    aimass = calc_vars.aimass.values
    ti = calc_vars.ti.values

    vthi = (zckb * ti / (zcmp * aimass))**(1 / 2)

    calc_vars.vthi.set(values=vthi, units='m/s')


def vpar(calc_vars):
    '''Parallel Velocity'''
    bpol = calc_vars.bpol.values
    btor = calc_vars.btor.values
    vpol = calc_vars.vpol.values
    vtor = calc_vars.vtor.values

    vpar = vtor + vpol * bpol / btor

    calc_vars.vpar.set(values=vpar, units=calc_vars.vtor.units)


def zeff(calc_vars):
    '''Effective Charge'''
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    zimp = calc_vars.zimp.values

    zeff = (nh + nf + zimp**2 * nz) / ne

    calc_vars.zeff.set(values=zeff, units='')


def calculate_gradient(gvar_name, var_name, drmin, calc_vars):
    '''
    Calculates the normalized gradient

    Parameters:
    * gvar_name (str): The name of the variable to store the gradient result
    * var_name (str): The name of the variable to take the gradient of
    * drmin (np.ndarray): Differential rmin
    * calc_vars (InputVariables): Object containing variable data
    '''

    rmaj = calc_vars.rmaj.values
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]  # includes origin

    # get variables related to the gradient from variable names
    gvar = getattr(calc_vars, gvar_name)
    var = getattr(calc_vars, var_name)

    # partial derivative along x-axis
    dxvar = np.diff(var.values, axis=0) / drmin

    # interpolate from x to xb
    set_interp = interp1d(x, dxvar, kind='cubic', fill_value="extrapolate", axis=0)
    dxvar = set_interp(xb)

    # take gradient
    gradient_values = rmaj * dxvar / var.values
    gvar.set(values=gradient_values, units='')

    opts = options.instance
    if opts.apply_smoothing:
        gvar.apply_smoothing(opts.input_points)

    gvar.clamp_gradient(100)  # TODO: should we be doing this?
    gvar.set_minvalue()

    if opts.reject_outliers:
        gvar.reject_outliers()

    gvar.check_for_nan()


def calculate_variable(var_function, calc_vars):
    '''
    Calculates the specified variable

    Values are saved to the input variables object by reference, so no return
    statements are needed.  After these values are calculated, we also apply
    optional smoothing, check variable minimum values, and check for nan
    values in the event that any errors occurred.  Note that the name of
    var_function must match the name of the variable in InputVariables in
    order for the calculation to be correctly called.

    This function was implemented in such a way that each variable can be
    checked as needed after completing it's calculation as if a callback
    function was passed to the calculation directly.

    Parameters:
    * var_function (function): Local function of the variable to calculate
    * calc_vars (InputVariables): Object containing variable data
    '''

    var_function(calc_vars)

    # Get the variable name corresponding to var_function by reflection
    var_name = var_function.__name__

    if options.instance.apply_smoothing:
        getattr(calc_vars, var_name).apply_smoothing(options.instance.input_points)

    getattr(calc_vars, var_name).set_minvalue()

    if options.instance.reject_outliers:
        getattr(calc_vars, var_name).reject_outliers()

    getattr(calc_vars, var_name).check_for_nan()


def calculate_base_variables(calc_vars):
    '''
    Calculates base variables

    Since the input parameter calc_vars is a reference, no return value for
    calculations are needed. Base calculations do not depend on any gradient
    variables or additional variables.

    Parameters:
    * calc_vars (InputVariables): Object containing variable data
    '''

    # Each use of the following calculate_variable functions are passing in
    # the local function of the variable to be calculated, which shares the
    # same name as the variable it calculates. Note that calculation order
    # matters here.

    calculate_variable(nh0, calc_vars)
    calculate_variable(nh, calc_vars)
    calculate_variable(ni, calc_vars)
    calculate_variable(ahyd, calc_vars)
    calculate_variable(aimass, calc_vars)
    calculate_variable(zeff, calc_vars)
    calculate_variable(btor, calc_vars)
    calculate_variable(bpol, calc_vars)
    calculate_variable(vpar, calc_vars)


def calculate_gradient_variables(calc_vars):
    '''
    Calculates gradient variables

    Since the input parameter calc_vars is a reference, no return value for
    calculations are needed.  Gradient calculations do not depend on any
    additional variables.

    Parameters:
    * calc_vars (InputVariables): Object containing variable data
    '''

    # Each use of the following calculate_variable functions are passing in
    # the local function of the variable to be calculated, which shares the
    # same name as the variable it calculates. Additionally, the sign on
    # drmin (differential rmin) sets the sign of the gradient equation

    drmin = np.diff(calc_vars.rmin.values, axis=0)
    calculate_gradient('gq', 'q', drmin, calc_vars)  # Only gq has positive drmin
    calculate_gradient('gne', 'ne', -drmin, calc_vars)
    calculate_gradient('gnh', 'nh', -drmin, calc_vars)
    calculate_gradient('gni', 'ni', -drmin, calc_vars)
    calculate_gradient('gnz', 'nz', -drmin, calc_vars)
    calculate_gradient('gte', 'te', -drmin, calc_vars)
    calculate_gradient('gti', 'ti', -drmin, calc_vars)
    calculate_gradient('gvpar', 'vpar', -drmin, calc_vars)
    calculate_gradient('gvpol', 'vpol', -drmin, calc_vars)
    calculate_gradient('gvtor', 'vtor', -drmin, calc_vars)


def calculate_additional_variables(calc_vars):
    '''
    Calculates additional variables

    Since the input parameter calc_vars is a reference, no return value for
    calculations are needed.  Additional variables are not used for
    constructing the MMM input file, and often mirror calculations made
    within MMM so that their values may be plotted.

    Parameters:
    * calc_vars (InputVariables): Object containing variable data
    '''

    # Each use of the following calculate_variable functions are passing in
    # the local function of the variable to be calculated, which shares the
    # same name as the variable it calculates. Note that calculation order
    # matters here.

    calculate_variable(tau, calc_vars)
    calculate_variable(eps, calc_vars)
    calculate_variable(p, calc_vars)
    calculate_variable(beta, calc_vars)
    calculate_variable(betae, calc_vars)
    calculate_variable(loge, calc_vars)
    calculate_variable(nuei, calc_vars)
    calculate_variable(nuei2, calc_vars)
    calculate_variable(vthe, calc_vars)
    calculate_variable(vthi, calc_vars)
    calculate_variable(nuste, calc_vars)
    calculate_variable(nusti, calc_vars)
    calculate_variable(gyrfi, calc_vars)
    calculate_variable(gmax, calc_vars)
    calculate_variable(shear, calc_vars)
    calculate_variable(shat, calc_vars)
    calculate_variable(alphamhd, calc_vars)
    calculate_variable(gave, calc_vars)
    calculate_variable(etae, calc_vars)
    calculate_variable(etai, calc_vars)


def calculate_new_variables(cdf_vars):
    '''
    Calculates new variables needed for MMM and data display

    A deepcopy of the variables object containing TRANSP data is made, since
    some calculations may overwrite variable values obtained from the CDF;
    this allows us to later compare calculated variables with CDF variables.
    This function should only be called when calculations need to be stored
    in a new variables object.

    Parameters:
    * cdf_vars (InputVariables): Variables object containing data from a CDF

    Returns:
    * calc_vars (InputVariables): Variables object containing calculation results
    '''

    calc_vars = copy.deepcopy(cdf_vars)
    calculate_base_variables(calc_vars)
    calculate_gradient_variables(calc_vars)
    calculate_additional_variables(calc_vars)

    return calc_vars


def get_calculated_vars():
    '''Returns (list of str): function names of calculated non-gradient variables'''
    return ([
        o[0] for o in inspect.getmembers(sys.modules[__name__])
        if inspect.isfunction(o[1]) and 'calculate' not in o[0]
    ])
