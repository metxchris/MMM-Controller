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

TODO:
* Consider replacing interp1d with Akima1DInterpolator, since TRANSP
  apparently uses this method of interpolation.
"""

# Standard Packages
import sys
import inspect
import functools

# 3rd Party Packages
import numpy as np
from scipy.interpolate import interp1d

# Local Packages
import modules.constants as constants
import modules.datahelper as datahelper


_gradients = set()  # Stores the names of calculated gradient variables


def gradient(gvar_name, var_name, drmin, calc_vars):
    '''
    Calculates the normalized gradient

    After the gradient value is calculated, optional smoothing is applied, and
    then the gradient is checked for min and nan values.  The overall sign of
    the gradient equation is determined by the sign given for drmin.

    Parameters:
    * gvar_name (str): The name of the variable to store the gradient result in
    * var_name (str): The name of the variable to take the gradient of
    * drmin (np.ndarray): Differential rmin
    * calc_vars (InputVariables): Object containing variable data
    '''

    _gradients.add(gvar_name)
    rmaj = calc_vars.rmaj.values
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]  # includes origin

    # get variables related to the gradient from variable names
    gvar = getattr(calc_vars, gvar_name)
    var = getattr(calc_vars, var_name)

    # partial derivative along radial dimension
    dxvar = np.diff(var.values, axis=0) / drmin

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, dxvar, kind='cubic', fill_value="extrapolate", axis=0)
    dxvar = set_interp(xb)

    # take gradient
    gradient_values = rmaj * dxvar / var.values
    gvar.set(values=gradient_values, units='')

    if calc_vars.options.apply_smoothing:
        gvar.apply_smoothing()

    gvar.set_origin_to_zero()
    gvar.clamp_values(constants.MAX_GRADIENT)
    gvar.set_minvalue()
    gvar.check_for_nan()


def calculation(func):
    '''
    Decorator function that wraps each non-gradient variable calculation

    In addition to storing the result of each calculation function to the
    corresponding variable object, the calculation decorator adds additional
    functionality at the end of each variable calculation as well.  In
    particular, optional smoothing is applied to the variable, and then the
    variable is checked for min and nan values. The units of each calculation
    are as specified for each corresponding variable in the InputVariables
    class.

    Note: The name of the variable functions must match the name of the
    variable in the InputVariables class in order for the calculations to
    work.

    Parameters:
    * func (function): Function of the variable to calculate
    '''

    @functools.wraps(func)  # Preserves the name of functions decorated with @calculation
    def wrapper(calc_vars):
        var = getattr(calc_vars, func.__name__)  # Get the variable corresponding to func
        var.values = func(calc_vars)  # Do the calculation

        if calc_vars.options.apply_smoothing:
            var.apply_smoothing()

        var.set_minvalue()
        var.check_for_nan()

        return func

    wrapper.calculation = True
    return wrapper


@calculation
def ahyd(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values

    return (nh0 + 2 * nd) / (nh0 + nd)


@calculation
def aimass(calc_vars):
    '''Mean Atomic Mass of Thermal Ions'''
    ahyd = calc_vars.ahyd.values
    aimp = calc_vars.aimp.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values

    return (ahyd * nh + aimp * nz) / (nh + nz)


@calculation
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

    return q**2 * betae * (gne + gte + ti / te * (gni + gti))


@calculation
def alphamhdunit(calc_vars):
    '''Alpha MHD (Weiland Definition)'''
    betaeunit = calc_vars.betaeunit.values
    gmaxunit = calc_vars.gmaxunit.values
    gne = np.maximum(np.minimum(calc_vars.gne.values, gmaxunit), -gmaxunit)
    gni = np.maximum(np.minimum(calc_vars.gni.values, gmaxunit), -gmaxunit)
    gte = np.maximum(np.minimum(calc_vars.gte.values, gmaxunit), -gmaxunit)
    gti = np.maximum(np.minimum(calc_vars.gti.values, gmaxunit), -gmaxunit)
    q = calc_vars.q.values
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    return q**2 * betaeunit * (gne + gte + ti / te * (gni + gti))


@calculation
def beta(calc_vars):
    '''Beta'''
    zcmu0 = constants.ZCMU0
    btor = calc_vars.btor.values
    p = calc_vars.p.values

    return 2 * zcmu0 * p / btor**2


@calculation
def betae(calc_vars):
    '''Electron Beta'''
    zckb = constants.ZCKB
    zcmu0 = constants.ZCMU0
    btor = calc_vars.btor.values
    ne = calc_vars.ne.values
    te = calc_vars.te.values

    return 2 * zcmu0 * ne * te * zckb / btor**2


@calculation
def betaeunit(calc_vars):
    '''Electron Beta'''
    zckb = constants.ZCKB
    zcmu0 = constants.ZCMU0
    bunit = calc_vars.bunit.values
    ne = calc_vars.ne.values
    te = calc_vars.te.values

    return 2 * zcmu0 * ne * te * zckb / bunit**2


@calculation
def bpol(calc_vars):
    '''Poloidal Magnetic Field'''
    btor = calc_vars.btor.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values

    return rmin / rmaj * btor / q


@calculation
def btor(calc_vars):
    '''Toroidal Magnetic Field'''
    bzxr = calc_vars.bzxr.values
    raxis = calc_vars.rmaj.values[0, :]
    rbound = calc_vars.rmaj.values[-1, :]
    rmaj = calc_vars.rmaj.values

    bz = bzxr / rbound
    # bz = calc_vars.bz.values

    return raxis / rmaj * bz


@calculation
def bunit(calc_vars):
    '''Magnetic Field (unit)'''
    btor0 = calc_vars.btor.values[0, :]
    rhochi = calc_vars.rhochi.values
    rmin = calc_vars.rmin.values
    x = calc_vars.x.values[:, 0]  # same for all time values
    xb = calc_vars.xb.values[:, 0]  # same for all time values

    drho_drmin = np.diff(rhochi, axis=0) / np.diff(rmin, axis=0)

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, drho_drmin, kind='cubic', fill_value="extrapolate", axis=0)
    dxrho = set_interp(xb)

    bunit = np.empty_like(dxrho)
    bunit[1:, :] = btor0 * rhochi[1:, :] / rmin[1:, :] * dxrho[1:, :]
    bunit[0, :] = bunit[1, :]

    return bunit


@calculation
def bunit_btor(calc_vars):
    '''Toroidal Magnetic Field'''
    bunit = calc_vars.bunit.values
    btor = calc_vars.btor.values

    return bunit / btor


@calculation
def csound(calc_vars):
    '''Sound Speed'''
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    aimass = calc_vars.aimass.values
    te = calc_vars.te.values

    return (zckb * te / (zcmp * aimass))**(1 / 2)


@calculation
def csound_a(calc_vars):
    '''Sound Frequency'''
    csound = calc_vars.csound.values
    amin = calc_vars.rmin.values[-1, :]

    return csound / amin


@calculation
def eps(calc_vars):
    '''Inverse Aspect Ratio'''
    arat = calc_vars.arat.values

    return 1 / arat


@calculation
def etae(calc_vars):
    '''etae = gte / gne'''
    gte = calc_vars.gte.values
    gne = calc_vars.gne.values

    return gte / gne


@calculation
def etai(calc_vars):
    '''
    etae = gti / gni

    Note that TRANSP appears to use an entirely different definition of ni
    when calculating gni, than it uses for the values of ni itself.  As such,
    our calculation of etai will generally not match with TRANSP values.
    '''

    gti = calc_vars.gti.values
    gni = calc_vars.gni.values

    return gti / gni


@calculation
def fle_etgm(calc_vars):
    '''Effective Charge'''
    alpha = calc_vars.alphamhd.values
    shear = calc_vars.shat_gxi.values

    kxoky = 0
    kyrhoe = 0.25

    return (
        (1 + kxoky**2) * kyrhoe**2 * (1 + 0.789868 * shear**2 - (10 / 9) * shear * alpha + (5 / 12) * alpha**2)
    )


@calculation
def gave(calc_vars):
    '''Average Magnetic Surface Curvature'''
    shear = calc_vars.shear.values
    # shear = calc_vars.shat_gxi.values
    alphamhdunit = calc_vars.alphamhdunit.values
    q = calc_vars.q.values

    return 2 / 3 + 5 / 9 * shear - 5 / 12 * alphamhdunit - alphamhdunit / (2 * q**2)


@calculation
def gave_noss(calc_vars):
    '''Average Magnetic Surface Curvature'''
    shear = calc_vars.shat_gxi.values
    alphamhdunit = calc_vars.alphamhdunit.values

    return 2 / 3 + 5 / 9 * shear - 5 / 12 * alphamhdunit


@calculation
def gmax(calc_vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    eps = calc_vars.eps.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    gyrfi = calc_vars.gyrfi.values
    vthi = calc_vars.vthi.values

    return rmaj / (vthi / gyrfi * q / eps)


@calculation
def gmaxunit(calc_vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    eps = calc_vars.eps.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    gyrfiunit = calc_vars.gyrfiunit.values
    vthi = calc_vars.vthi.values

    return 2 * rmaj / (vthi / gyrfiunit * q / eps)


@calculation
def gmaxunit_gte(calc_vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    gmaxunit = calc_vars.gmaxunit.values
    gte = calc_vars.gte.values

    gmaxunit_gte = np.minimum(gmaxunit / np.absolute(gte), 10)
    gmaxunit_gte[0, :] = gmaxunit_gte[1, :]

    return gmaxunit_gte


@calculation
def gmaxunit_gne(calc_vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    gmaxunit = calc_vars.gmaxunit.values
    gne = calc_vars.gne.values

    gmaxunit_gne = np.minimum(gmaxunit / np.absolute(gne), 10)
    gmaxunit_gne[0, :] = gmaxunit_gne[1, :]

    return gmaxunit_gne



@calculation
def gmaxunit_gnh(calc_vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    gmaxunit = calc_vars.gmaxunit.values
    gnh = calc_vars.gnh.values

    gmaxunit_gnh = np.minimum(gmaxunit / np.absolute(gnh), 10)
    gmaxunit_gnh[0, :] = gmaxunit_gnh[1, :]

    return gmaxunit_gnh



@calculation
def gyrfi(calc_vars):
    '''Ion Gyrofrequency'''
    zce = constants.ZCE
    zcmp = constants.ZCMP
    aimass = calc_vars.aimass.values
    btor = calc_vars.btor.values

    return zce * btor / (zcmp * aimass)


@calculation
def gyrfiunit(calc_vars):
    '''Ion Gyrofrequency'''
    zce = constants.ZCE
    zcmp = constants.ZCMP
    aimass = calc_vars.aimass.values
    bunit = calc_vars.bunit.values

    return zce * bunit / (zcmp * aimass)


@calculation
def gyrfe(calc_vars):
    '''Electron Gyrofrequency'''
    zce = constants.ZCE
    zcme = constants.ZCME
    btor = calc_vars.btor.values

    return zce * btor / zcme


@calculation
def gyrfeunit(calc_vars):
    '''Electron Gyrofrequency'''
    zce = constants.ZCE
    zcme = constants.ZCME
    bunit = calc_vars.bunit.values

    return zce * bunit / zcme


@calculation
def lare(calc_vars):
    '''Electron Gyroradius'''
    vthe = calc_vars.vthe.values
    gyrfe = calc_vars.gyrfe.values

    return vthe / gyrfe


@calculation
def lareunit(calc_vars):
    '''Electron Gyroradius'''
    vthe = calc_vars.vthe.values
    gyrfeunit = calc_vars.gyrfeunit.values

    return vthe / gyrfeunit


@calculation
def loge(calc_vars):
    '''Electron Coulomb Logarithm'''

    # TODO: Need to add equations for different TE ranges
    ne = calc_vars.ne.values
    te = calc_vars.te.values

    # NRL Plasma Formulary Definition
    # return 37.8 - np.log(ne**(1 / 2) / te)

    # TRANSP definition (equivalent)
    zeff = calc_vars.zeff.values
    return 39.23 - np.log(zeff * ne**(1 / 2) / te)


@calculation
def p(calc_vars):
    '''Plasma Pressure'''
    zckb = constants.ZCKB
    ne = calc_vars.ne.values
    ni = calc_vars.ni.values
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    return (ne * te + ni * ti) * zckb


@calculation
def nh0(calc_vars):
    '''Hydrogen Ion Density'''
    nd = calc_vars.nd.values
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nz = calc_vars.nz.values
    zimp = calc_vars.zimp.values

    return ne - zimp * nz - nf - nd


@calculation
def nh(calc_vars):
    '''Total Hydrogenic Ion Density'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values

    return nh0 + nd


@calculation
def ni(calc_vars):
    '''Thermal Ion Density'''
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values

    # TRANSP Definition
    return nh + nz


@calculation
def nuei(calc_vars):
    '''Collision Frequency'''
    zcf = constants.ZCF
    ne = calc_vars.ne.values
    te = calc_vars.te.values
    zeff = calc_vars.zeff.values
    loge = calc_vars.loge.values

    return zcf * 2**(1 / 2) * ne * loge * zeff / te**(3 / 2)


@calculation
def nuei2(calc_vars):
    '''OLD NOTE: Not sure what to call this, but it leads to the approx the correct NUSTI'''
    zcf = constants.ZCF
    ni = calc_vars.ni.values
    ti = calc_vars.ti.values
    zeff = calc_vars.zeff.values
    loge = calc_vars.loge.values

    return zcf * 2**(1 / 2) * ni * loge * zeff / ti**(3 / 2)


@calculation
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

    return nuei * eps**(-3 / 2) * q * rmaj / vthe


@calculation
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

    return nuei2 * eps**(-3 / 2) * q * rmaj / (2 * vthi) * (zcme / zcmp)**(1 / 2)


@calculation
def rhochi(calc_vars):
    '''
    Rho, derived from magnetic flux
    '''

    pi = constants.PI
    btor0 = calc_vars.btor.values[0, :]
    bftor = calc_vars.bftor.values

    return (bftor / (pi * btor0))**(1 / 2)


@calculation
def shat(calc_vars):
    '''
    Effective Magnetic Shear

    TRANSP uses a different definition for shat than we use here, so we expect
    our values to be similar, but not a direct match.
    '''

    elong = calc_vars.elong.values
    shear = calc_vars.shear.values

    return np.maximum(2 * shear - 1 + (elong * (shear - 1))**2, 0)**(1 / 2)


@calculation
def shat_gxi(calc_vars):
    '''
    Effective Magnetic Shear

    TRANSP uses a different definition for shat than we use here, so we expect
    our values to be similar, but not a direct match.
    '''

    a = calc_vars.rmin.values[-1, :]
    gxi = calc_vars.gxi.values
    shear = calc_vars.shear.values

    return np.maximum(2 * shear - 1 + ((a * gxi) * (shear - 1))**2, 0)**(1 / 2)


@calculation
def shear(calc_vars):
    '''Magnetic Shear'''
    gq = calc_vars.gq.values
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values

    return gq * rmin / rmaj


@calculation
def tau(calc_vars):
    '''Temperature Ratio te / ti'''
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    return te / ti


@calculation
def vthe(calc_vars):
    '''Thermal Velocity of Electrons'''
    zckb = constants.ZCKB
    zcme = constants.ZCME
    te = calc_vars.te.values

    return (2 * zckb * te / zcme)**(1 / 2)


@calculation
def vthe_qrmaj(calc_vars):
    '''Thermal Velocity of Electrons / (qR)'''
    vthe = calc_vars.vthe.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values

    return vthe / (q * rmaj)


@calculation
def vthi(calc_vars):
    '''Thermal Velocity of Ions'''
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    aimass = calc_vars.aimass.values
    ti = calc_vars.ti.values

    return (2 * zckb * ti / (zcmp * aimass))**(1 / 2)


@calculation
def vpar(calc_vars):
    '''Parallel Velocity'''
    bpol = calc_vars.bpol.values
    btor = calc_vars.btor.values
    vpol = calc_vars.vpol.values
    vtor = calc_vars.vtor.values

    return vtor + vpol * bpol / btor


@calculation
def xetgm_const(calc_vars):
    '''ETGM Diffusivity Factor'''
    lareunit = calc_vars.lareunit.values
    vthe = calc_vars.vthe.values
    gmaxunit = calc_vars.gmaxunit.values
    gte = np.maximum(np.minimum(calc_vars.gte.values, gmaxunit), -gmaxunit)
    rmaj = calc_vars.rmaj.values

    return lareunit ** 2 * vthe * gte / rmaj


@calculation
def zeff(calc_vars):
    '''Effective Charge'''
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    zimp = calc_vars.zimp.values

    return (nh + nf + zimp**2 * nz) / ne


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

    nh0(calc_vars)
    nh(calc_vars)
    ni(calc_vars)
    ahyd(calc_vars)
    aimass(calc_vars)
    zeff(calc_vars)
    btor(calc_vars)
    bpol(calc_vars)
    vpar(calc_vars)


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
    # Positive drmin
    gradient('gq', 'q', drmin, calc_vars)
    # Negative drmin
    gradient('gne', 'ne', -drmin, calc_vars)
    gradient('gnh', 'nh', -drmin, calc_vars)
    gradient('gni', 'ni', -drmin, calc_vars)
    gradient('gnz', 'nz', -drmin, calc_vars)
    gradient('gte', 'te', -drmin, calc_vars)
    gradient('gti', 'ti', -drmin, calc_vars)
    gradient('gvpar', 'vpar', -drmin, calc_vars)
    gradient('gvpol', 'vpol', -drmin, calc_vars)
    gradient('gvtor', 'vtor', -drmin, calc_vars)

    if calc_vars.options.use_gnezero:
        calc_vars.gne.values[:, :] = 1e-8

    if calc_vars.options.use_gtezero:
        calc_vars.gte.values[:, :] = 1e-8

    if calc_vars.options.use_gneabs:
        calc_vars.gne.values = np.absolute(calc_vars.gne.values)

    if calc_vars.options.use_gnethreshold:
        calc_vars.gne.values[:, :] = 1
        calc_vars.gte.values[:, :] = 1e-12
        calc_vars.gti.values[:, :] = 1e-12
        calc_vars.gni.values[:, :] = 1e-12

    if calc_vars.options.use_gtethreshold:
        calc_vars.gne.values = np.absolute(calc_vars.gne.values)
        calc_vars.gne.values[:, :] = 1e-12
        calc_vars.gte.values[:, :] = 1
        calc_vars.gti.values[:, :] = 1e-12
        calc_vars.gni.values[:, :] = 1e-12


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

    rhochi(calc_vars)
    bunit(calc_vars)
    bunit_btor(calc_vars)
    tau(calc_vars)
    eps(calc_vars)
    p(calc_vars)
    beta(calc_vars)
    betae(calc_vars)
    betaeunit(calc_vars)
    csound(calc_vars)
    csound_a(calc_vars)
    loge(calc_vars)
    nuei(calc_vars)
    nuei2(calc_vars)
    vthe(calc_vars)
    vthe_qrmaj(calc_vars)
    vthi(calc_vars)
    nuste(calc_vars)
    nusti(calc_vars)
    gyrfe(calc_vars)
    gyrfeunit(calc_vars)
    gyrfi(calc_vars)
    gyrfiunit(calc_vars)
    lare(calc_vars)
    lareunit(calc_vars)
    gmax(calc_vars)
    gmaxunit(calc_vars)
    gmaxunit_gte(calc_vars)
    gmaxunit_gne(calc_vars)
    gmaxunit_gnh(calc_vars)
    gmaxunit(calc_vars)
    shear(calc_vars)
    shat(calc_vars)
    shat_gxi(calc_vars)
    alphamhd(calc_vars)
    alphamhdunit(calc_vars)
    gave(calc_vars)
    gave_noss(calc_vars)
    fle_etgm(calc_vars)
    xetgm_const(calc_vars)
    etae(calc_vars)
    etai(calc_vars)


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

    calc_vars = datahelper.deepcopy_data(cdf_vars)
    calculate_base_variables(calc_vars)
    calculate_gradient_variables(calc_vars)
    calculate_additional_variables(calc_vars)

    return calc_vars


def get_calculated_vars():
    '''
    Gets list of all calculated variables

    Note: Gradients only show up here after their calculations were made

    Returns:
    * (list[str]): Names of all calculated variables
    '''

    calculations = ([
        o[0] for o in inspect.getmembers(sys.modules[__name__])
        if inspect.isfunction(o[1]) and hasattr(o[1], calculation.__name__)
    ])
    return calculations + list(_gradients)
