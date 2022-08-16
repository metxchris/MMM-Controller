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
    set_interp = interp1d(x, dxvar, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)
    dxvar = set_interp(xb)

    # take gradient
    gradient_values = rmaj * dxvar / var.values
    gvar.set(values=gradient_values, units='')

    if calc_vars.options.apply_smoothing:
        gvar.apply_smoothing()

    gvar.set_origin_to_zero()
    gvar.clamp_values(constants.MAX_GRADIENT)
    gvar.set_minvalue(ignore_exceptions=calc_vars.options.ignore_exceptions)
    gvar.check_for_nan(ignore_exceptions=calc_vars.options.ignore_exceptions)


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

        var.set_minvalue(ignore_exceptions=calc_vars.options.ignore_exceptions)
        var.check_for_nan(ignore_exceptions=calc_vars.options.ignore_exceptions)

        return func

    wrapper.calculation = True
    return wrapper


def calculation_output(func):
    '''
    Decorator function that wraps each non-gradient variable calculation

    Same as the calculation wrapper, but calculations are made with and stored
    in output variables.  Output variables will be 1D arrays compared with
    the 2D arrays of calc_vars.

    Note: The name of the variable functions must match the name of the
    variable in the OutputVariables class in order for the calculations to
    work.

    Parameters:
    * func (function): Function of the variable to calculate
    '''

    @functools.wraps(func)  # Preserves the name of functions decorated with @calculation_output
    def wrapper(calc_vars, output_vars):
        var = getattr(output_vars, func.__name__)  # Get the variable corresponding to func
        var.values = func(calc_vars, output_vars)  # Do the calculation

        if output_vars.options.apply_smoothing:
            var.apply_smoothing()

        var.set_minvalue(ignore_exceptions=calc_vars.options.ignore_exceptions)
        var.check_for_nan(ignore_exceptions=calc_vars.options.ignore_exceptions)

        return func

    wrapper.calculation = True
    return wrapper


@calculation
def ah(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values

    return (nh0 + 2 * nd) / (nh0 + nd)


@calculation
def ai(calc_vars):
    '''Mean Atomic Mass of Thermal Ions'''
    ah = calc_vars.ah.values
    az = calc_vars.az.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values

    return (ah * nh + az * nz) / (nh + nz)


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
def alphamhdu(calc_vars):
    '''Alpha MHD (Weiland Definition)'''
    betaeu = calc_vars.betaeu.values
    gne = calc_vars.gne.values
    gni = calc_vars.gni.values
    gte = calc_vars.gte.values
    gti = calc_vars.gti.values
    q = calc_vars.q.values
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    return q**2 * betaeu * (gne + gte + ti / te * (gni + gti))


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
def betapu(calc_vars):
    '''Beta Prime (Unit)'''
    alphamhdu = calc_vars.alphamhdu.values
    q = calc_vars.q.values

    return 2.54 * alphamhdu / (2 * q**2)


@calculation
def betaeu(calc_vars):
    '''Electron Beta'''
    zckb = constants.ZCKB
    zcmu0 = constants.ZCMU0
    bu = calc_vars.bu.values
    ne = calc_vars.ne.values
    te = calc_vars.te.values

    return 2 * zcmu0 * ne * te * zckb / bu**2


@calculation
def bpol(calc_vars):
    '''Poloidal Magnetic Field'''
    btor = calc_vars.btor.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values[-1, :]
    rmin = calc_vars.rmin.values

    return rmin / rmaj * btor / q


@calculation
def btor(calc_vars):
    '''Toroidal Magnetic Field'''
    bzxr = calc_vars.bzxr.values
    rmaj = calc_vars.rmaj.values

    return bzxr / rmaj


@calculation
def bu(calc_vars):
    '''Magnetic Field (unit)'''
    btor0 = calc_vars.btor.values[0, :]
    rhochi = calc_vars.rhochi.values
    rmin = calc_vars.rmin.values
    x = calc_vars.x.values[:, 0]  # same for all time values
    xb = calc_vars.xb.values[:, 0]  # same for all time values

    drho_drmin = np.diff(rhochi, axis=0) / np.diff(rmin, axis=0)

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, drho_drmin, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)
    dxrho = set_interp(xb)

    bu = np.empty_like(dxrho)
    bu[1:, :] = btor0 * rhochi[1:, :] / rmin[1:, :] * dxrho[1:, :]
    bu[0, :] = bu[1, :]

    return bu


@calculation
def bu_btor(calc_vars):
    '''Toroidal Magnetic Field'''
    bu = calc_vars.bu.values
    btor = calc_vars.btor.values

    return bu / btor


@calculation
def csound(calc_vars):
    '''Sound Speed'''
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    ai = calc_vars.ai.values
    te = calc_vars.te.values

    return (zckb * te / (zcmp * ai))**(0.5)


@calculation
def csound_a(calc_vars):
    '''Sound Frequency'''
    csound = calc_vars.csound.values
    amin = calc_vars.rmin.values[-1, :]

    return csound / amin


@calculation
def curlh(calc_vars):
    '''LH Current'''
    curdlh = calc_vars.curdlh.values
    area = calc_vars.surf.values
    rmin = calc_vars.rmin.values

    # return curdlh * area
    return curdlh * constants.PI * rmin**2


@calculation
def curoh(calc_vars):
    '''OH Current'''
    curdoh = calc_vars.curdoh.values
    area = calc_vars.darea.values

    return curdoh * area


@calculation
def dbp(calc_vars):
    '''dBpol / drho'''
    bpol = calc_vars.bpol.values
    rho = calc_vars.rho.values
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]  # includes origin

    # partial derivative along radial dimension
    dbpol_drho = np.diff(bpol, axis=0) / np.diff(rho, axis=0)

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, dbpol_drho, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)

    return set_interp(xb)


@calculation
def d2bp(calc_vars):
    '''d^2Bpol / drho^2'''
    dbp = calc_vars.dbp.values
    rho = calc_vars.rho.values
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]  # includes origin

    # partial derivative along radial dimension
    d2bpol_drho = np.diff(dbp, axis=0) / np.diff(rho, axis=0)

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, d2bpol_drho, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)

    return set_interp(xb)


@calculation
def e_r_grp(calc_vars):
    '''Radial Electric Field (Pressure Term)'''
    zce = constants.ZCE
    zckb = constants.ZCKB
    ni = calc_vars.ni.values
    rmin = calc_vars.rmin.values
    ti = calc_vars.ti.values
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]
    zeff = calc_vars.zeff.values

    p_i = ti * ni * zckb
    drmin = np.diff(rmin, axis=0)
    dpdr_x = np.diff(p_i, axis=0) / drmin

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, dpdr_x, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)
    dpdr = set_interp(xb)

    # From pt_vflows_mod.f90:
    # zE_r_grp(lcentr:lep1) =  1.0 / ( xzeffp(lcentr:lep1,1) * ze * rhoth(lcentr:lep1,2) * 1.0E6 ) * zgrp(lcentr:lep1)

    return dpdr / (zeff * zce * ni)


@calculation
def e_r_phi(calc_vars):
    '''Radial Electric Field (vtor Term)'''
    bpol = calc_vars.bpol.values
    omega = calc_vars.omega.values
    rmaj = calc_vars.rmaj.values

    # From pt_vflows_mod.f90:
    # zE_r_phi(lcentr:lep1) = omega(lcentr:lep1,1) * rmjrmp(lcentr:lep1,1) * zcm_to_m * bpol(lcentr:lep1)

    return omega * rmaj * bpol


@calculation
def e_r_tht(calc_vars):
    '''Radial Electric Field (vpol Term)'''
    vpol = calc_vars.vpol.values
    bzxr = calc_vars.bzxr.values
    rmaj = calc_vars.rmaj.values

    # From pt_vflows_mod.f90:
    # zE_r_tht(lcentr:lep1) = -1.0 * zvpol(lcentr:lep1) * zcm_to_m * bzxr / rmjrmp(lcentr:lep1,1)

    return -vpol * bzxr / rmaj


@calculation
def eps(calc_vars):
    '''Inverse Aspect Ratio'''
    arat = calc_vars.arat.values

    return 1 / arat


@calculation
def epsne(calc_vars):
    '''Pinch Term'''
    gbu = calc_vars.gbu.values
    gne = calc_vars.gne.values

    return 2 * gbu / gne


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
def gelong(calc_vars):
    '''Elongation Gradient (delong / deps)'''
    _gradients.add('gelong')
    elong = calc_vars.elong.values
    rmin = calc_vars.rmin.values
    rmaj0 = calc_vars.rmaj.values[0, :]
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]  # includes origin

    # Inverse aspect ratio
    eps = rmin / rmaj0

    # partial derivative along radial dimension
    delong_deps = np.diff(elong, axis=0) / np.diff(eps, axis=0)

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, delong_deps, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)

    return set_interp(xb)


@calculation
def gyrfi(calc_vars):
    '''Ion Gyrofrequency'''
    zce = constants.ZCE
    zcmp = constants.ZCMP
    ai = calc_vars.ai.values
    btor = calc_vars.btor.values

    return zce * btor / (zcmp * ai)


@calculation
def gyrfiu(calc_vars):
    '''Ion Gyrofrequency'''
    zce = constants.ZCE
    zcmp = constants.ZCMP
    ai = calc_vars.ai.values
    bu = calc_vars.bu.values

    return zce * bu / (zcmp * ai)


@calculation
def gyrfe(calc_vars):
    '''Electron Gyrofrequency'''
    zce = constants.ZCE
    zcme = constants.ZCME
    btor = calc_vars.btor.values

    return zce * btor / zcme


@calculation
def gyrfeu(calc_vars):
    '''Electron Gyrofrequency'''
    zce = constants.ZCE
    zcme = constants.ZCME
    bu = calc_vars.bu.values

    return zce * bu / zcme


@calculation
def gxi(calc_vars):

    rhochi = calc_vars.rhochi.values
    elong = calc_vars.elong.values
    rmin = calc_vars.rmin.values
    x = calc_vars.x.values[:, 0]
    xb = calc_vars.xb.values[:, 0]

    drmin = np.diff(calc_vars.rmin.values, axis=0)

    dxvar = np.diff(rhochi, axis=0) / drmin

    # interpolate from x grid to xb grid
    set_interp = interp1d(x, dxvar, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)
    dxvar2 = set_interp(xb)

    return dxvar2 * rmin[-1, :] * elong[-1, :]**0.5
    # return (1 + elong**2 / (2 * elong**2))**0.5



@calculation
def lare(calc_vars):
    '''Electron Gyroradius'''
    vthe = calc_vars.vthe.values
    gyrfe = calc_vars.gyrfe.values

    return vthe / gyrfe


@calculation
def lareu(calc_vars):
    '''Electron Gyroradius'''
    vthe = calc_vars.vthe.values
    gyrfeu = calc_vars.gyrfeu.values

    return vthe / gyrfeu


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
    return 39.23 - np.log(zeff * ne**(0.5) / te)


@calculation
def nh0(calc_vars):
    '''Hydrogen Ion Density'''
    nd = calc_vars.nd.values
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nz = calc_vars.nz.values
    zz = calc_vars.zz.values

    return ne - zz * nz - nf - nd


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
def ni2(calc_vars):
    '''Thermal Ion Density v2'''
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    nf = calc_vars.nf.values
    zz = calc_vars.zz.values

    # TRANSP likely uses this for taking ion density gradients
    return nh + zz**2 * nz + nf


@calculation
def nuei(calc_vars):
    '''Collision Frequency'''
    zcf = constants.ZCF
    ne = calc_vars.ne.values
    te = calc_vars.te.values
    zeff = calc_vars.zeff.values
    loge = calc_vars.loge.values

    return zcf * 2**(0.5) * ne * loge * zeff / te**(1.5)


@calculation
def nuei2(calc_vars):
    '''OLD NOTE: Not sure what to call this, but it leads to the approx the correct NUSTI'''
    zcf = constants.ZCF
    ni = calc_vars.ni.values
    ti = calc_vars.ti.values
    zeff = calc_vars.zeff.values
    loge = calc_vars.loge.values

    return zcf * 2**(0.5) * ni * loge * zeff / ti**(1.5)


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

    return nuei * eps**(-1.5) * q * rmaj / vthe


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

    return nuei2 * eps**(-3 / 2) * q * rmaj / (2 * vthi) * (zcme / zcmp)**(0.5)


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
def rhochi(calc_vars):
    '''
    Rho, derived from magnetic flux
    '''

    pi = constants.PI
    btor0 = calc_vars.btor.values[0, :]
    bftor = calc_vars.bftor.values

    return (bftor / (pi * btor0))**(0.5)


@calculation
def rhos(calc_vars):
    '''
    Rhos
    '''

    csound = calc_vars.csound.values
    gyrfi = calc_vars.gyrfi.values

    return csound / gyrfi


@calculation
def rhosu(calc_vars):
    '''
    Rhos (Unit)
    '''

    csound = calc_vars.csound.values
    gyrfiu = calc_vars.gyrfiu.values

    return csound / gyrfiu


@calculation
def shat(calc_vars):
    '''Effective Magnetic Shear'''

    elong = calc_vars.elong.values
    shear = calc_vars.shear.values

    return np.maximum(2 * shear - 1 + (elong * (shear - 1))**2, 0)**(0.5)


@calculation
def shat_gxi(calc_vars):
    '''Effective Magnetic Shear'''

    a = calc_vars.rmin.values[-1, :]
    gxi = calc_vars.gxi.values
    shear = calc_vars.shear.values

    signs = np.ones_like(shear)
    signs[shear < 0] = -1

    return np.maximum(2 * shear - 1 + ((a * gxi) * (shear - 1))**2, 0)**(0.5) * signs


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
def tauh(calc_vars):
    '''Temperature Ratio ti / te'''
    tau = calc_vars.tau.values

    return 1 / tau


@calculation
def vei_nc(calc_vars):
    '''Collision Frequency for DBM'''
    zcc = constants.ZCC
    zce = constants.ZCE
    zcme = constants.ZCME
    zceps0 = constants.ZCEPS0
    zcmu0 = constants.ZCMU0
    etanc = calc_vars.etanc.values
    ne = calc_vars.ne.values

    wpe = ((ne * zce**2) / (zcme * zceps0))**(0.5)

    return 1.96 * wpe**2 / zcc**2 * etanc / zcmu0


@calculation
def vthe(calc_vars):
    '''Thermal Velocity of Electrons'''
    zckb = constants.ZCKB
    zcme = constants.ZCME
    te = calc_vars.te.values

    return (2 * zckb * te / zcme)**(0.5)


@calculation
def vthi(calc_vars):
    '''Thermal Velocity of Ions'''
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    ai = calc_vars.ai.values
    ti = calc_vars.ti.values

    return (2 * zckb * ti / (zcmp * ai))**(0.5)


@calculation
def vpar(calc_vars):
    '''Parallel Velocity'''
    bpol = calc_vars.bpol.values
    btor = calc_vars.btor.values
    vpol = calc_vars.vpol.values
    vtor = calc_vars.vtor.values

    return vtor + vpol * bpol / btor


@calculation
def vpolx(calc_vars):
    '''Parallel Velocity'''
    rmaj = calc_vars.rmaj.values
    theta = calc_vars.omega.values

    return rmaj * theta


@calculation
def vtorx(calc_vars):
    '''Parallel Velocity'''
    rmaj = calc_vars.rmaj.values
    omega = calc_vars.omega.values

    return rmaj * omega


@calculation
def wbounce(calc_vars):
    '''Bounce Frequency'''
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values
    wtransit = calc_vars.wtransit.values

    return (rmin / (2 * rmaj))**(0.5) * wtransit


@calculation
def wtransit(calc_vars):
    '''Transit Frequency'''
    vthe = calc_vars.vthe.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values

    return vthe / (q * rmaj)


@calculation
def wexb(calc_vars):
    '''ExB Shear Rate (adapted from pt_vflows_mod.f90)'''

    def dfdr(E_r):
        # zE_r(lcp1:lep1)/ ( bpol(lcp1:lep1) * rmjrmp(lcp1:lep1,1) * zcm_to_m )

        f = E_r / (bpol * rmaj)

        drmin = np.diff(rmin, axis=0)
        dfdr_x = np.diff(f, axis=0) / drmin

        # interpolate from x grid to xb grid
        set_interp = interp1d(x, dfdr_x, kind=constants.INTERP_TYPE, fill_value="extrapolate", axis=0)
        return set_interp(xb)

    x = calc_vars.x.values[:, 0]  # same for all time values
    xb = calc_vars.xb.values[:, 0]  # same for all time values
    bpol = calc_vars.bpol.values
    bzxr = calc_vars.bzxr.values
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values

    # Electric Field Components
    E_r_phi = calc_vars.e_r_phi.values
    E_r_tht = calc_vars.e_r_tht.values
    E_r_grp = calc_vars.e_r_grp.values

    bratio = rmaj * bpol / (bzxr / rmaj[0, :])

    wexb_phi = np.minimum(bratio * dfdr(E_r_phi), 1e6)
    wexb_tht = np.minimum(bratio * dfdr(E_r_tht), 1e6) * 1e-2
    wexb_grp = np.minimum(bratio * dfdr(E_r_grp), 1e6)

    return np.absolute(wexb_phi + wexb_tht + wexb_grp)


@calculation
def xetgm_const(calc_vars):
    '''ETGM Diffusivity Factor'''
    lareu = calc_vars.lareu.values
    vthe = calc_vars.vthe.values
    gte = calc_vars.gte.values
    rmaj = calc_vars.rmaj.values

    return lareu ** 2 * vthe * gte / rmaj


@calculation
def xke(calc_vars):
    '''Total Electron Thermal Diffusivity from CDF'''
    condepr = calc_vars.condepr.values
    condewnc = calc_vars.condewnc.values
    xkepaleo = calc_vars.xkepaleo.values

    xke = condepr + condewnc + xkepaleo

    return xke


@calculation
def xki(calc_vars):
    '''Total Ion Thermal Diffusivity from CDF'''
    condipr = calc_vars.condipr.values
    condiwnc = calc_vars.condiwnc.values

    xki = condipr + condiwnc

    return xki


@calculation
def zeff(calc_vars):
    '''Effective Charge'''
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    zz = calc_vars.zz.values

    return (nh + nf + zz**2 * nz) / np.maximum(ne, 1e-16)


@calculation
def pf(calc_vars):
    '''Effective Charge'''
    nf = calc_vars.nf.values
    tf = calc_vars.ti.values

    return nf * tf


@calculation
def gpf2(calc_vars):
    '''Effective Charge'''
    gnf = calc_vars.gnf.values
    gtf = calc_vars.gti.values

    return gnf + gtf


@calculation_output
def gmanMTM(calc_vars, output_vars):
    '''MTM Growth Rate Normalized'''
    t = calc_vars.options.time_idx
    gmaMTM = output_vars.gmaMTM.values
    csound_a = calc_vars.csound_a.values

    return gmaMTM / np.maximum(csound_a[:, t], 1e-16)


@calculation_output
def gmadiffETGM(calc_vars, output_vars):
    '''ETGM Growth Rate resonance'''
    gmaETGM = output_vars.gmaETGM.values
    wdeETGM = output_vars.wdeETGM.values

    return gmaETGM - wdeETGM


@calculation_output
def nEPM(calc_vars, output_vars):
    '''EPM toroidal wave number'''
    t = calc_vars.options.time_idx
    kyrhosEPM = output_vars.kyrhosEPM.values
    rhos = calc_vars.rhos.values[:, t]
    rmin = calc_vars.rmin.values[:, t]
    q = calc_vars.q.values[:, t]

    return kyrhosEPM / rhos * rmin / q


@calculation_output
def omgdiffETGM(calc_vars, output_vars):
    '''ETGM Frequency resonance'''
    omgETGM = output_vars.omgETGM.values
    wdeETGM = output_vars.wdeETGM.values

    return omgETGM - wdeETGM


@calculation_output
def wde_gaveETGM(calc_vars, output_vars):
    '''ETGM Frequency resonance'''
    wdeETGM = output_vars.wdeETGM.values
    gaveETGM = output_vars.gaveETGM.values
    sign = output_vars.gaveETGM.get_sign()

    return wdeETGM / (sign * np.maximum(np.absolute(gaveETGM), 1e-16))


@calculation_output
def wseETGM(calc_vars, output_vars):
    '''omega_*e'''
    t = calc_vars.options.time_idx
    gne = calc_vars.gne.values[:, t]
    wde_gaveETGM = output_vars.wde_gaveETGM.values

    return 0.5 * gne * wde_gaveETGM


@calculation_output
def wsetaETGM(calc_vars, output_vars):
    '''omega_*e * (1 + eta_e)'''
    t = calc_vars.options.time_idx
    etae = calc_vars.etae.values[:, t]
    wseETGM = output_vars.wseETGM.values

    return (1 + etae) * wseETGM


@calculation_output
def wteETGM(calc_vars, output_vars):
    '''omega_Te'''
    t = calc_vars.options.time_idx
    gte = calc_vars.gte.values[:, t]
    wde_gaveETGM = output_vars.wde_gaveETGM.values

    return 0.5 * gte * wde_gaveETGM


@calculation_output
def omgnMTM(calc_vars, output_vars):
    '''MTM Growth Rate Normalized'''
    t = calc_vars.options.time_idx
    omgMTM = output_vars.omgMTM.values
    csound_a = calc_vars.csound_a.values

    return omgMTM / np.maximum(csound_a[:, t], 1e-16)


@calculation_output
def waETGM(calc_vars, output_vars):
    '''Alfven Frequency (Unit)'''
    t = calc_vars.options.time_idx
    zcmu0 = constants.ZCMU0
    zcmp = constants.ZCMP
    bu = calc_vars.bu.values[:, t]
    ai = calc_vars.ai.values[:, t]
    ni = calc_vars.ni.values[:, t]
    kpara2 = output_vars.kpara2ETGM.values

    return kpara2**(0.5) * bu / (zcmu0 * zcmp * ai * ni)**(0.5)


def calculate_output_variables(calc_vars, output_vars, controls):
    '''
    Calculations using output variables

    Output calculations can depend on any base variable, gradient, additional
    variable, or any output variable obtained from MMM.  Note that calc_vars
    will contain 2D arrays, whereas output_vars will only contain 1D arrays.

    Parameters:
    * calc_vars (InputVariables): Object containing variable data
    * output_vars (OutputVariables): Object containing variable data
    * controls (InputControls): Object containing controls data
    '''

    # Each use of the following calculate_variable functions are passing in
    # the local function of the variable to be calculated, which shares the
    # same name as the variable it calculates. Note that calculation order
    # matters here.

    if controls.cmodel_etgm.values:
        wde_gaveETGM(calc_vars, output_vars)
        wseETGM(calc_vars, output_vars)
        wsetaETGM(calc_vars, output_vars)
        wteETGM(calc_vars, output_vars)
        omgdiffETGM(calc_vars, output_vars)
        gmadiffETGM(calc_vars, output_vars)
        waETGM(calc_vars, output_vars)

    if controls.cmodel_mtm.values:
        gmanMTM(calc_vars, output_vars)
        omgnMTM(calc_vars, output_vars)

    if controls.cmodel_epm.values:
        nEPM(calc_vars, output_vars)


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

    # nh0(calc_vars)
    nh(calc_vars)
    ni(calc_vars)
    ni2(calc_vars)
    ah(calc_vars)
    ai(calc_vars)
    zeff(calc_vars)
    btor(calc_vars)
    rhochi(calc_vars)
    bu(calc_vars)
    bpol(calc_vars)
    dbp(calc_vars)
    d2bp(calc_vars)
    vpar(calc_vars)
    vtorx(calc_vars)
    vpolx(calc_vars)
    pf(calc_vars)

    # Calculations for testing
    # 
    # e_r_grp(calc_vars)
    # e_r_phi(calc_vars)
    # e_r_tht(calc_vars)
    # wexb(calc_vars)
    # gxi(calc_vars)

    if hasattr(calc_vars.options, 'use_etgm_btor') and calc_vars.options.use_etgm_btor:
        calc_vars.bu.values = calc_vars.btor.values


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
    gradient('gbu', 'bu', drmin, calc_vars)
    gradient('gbtor', 'btor', drmin, calc_vars)
    # Negative drmin
    gradient('gne', 'ne', -drmin, calc_vars)
    gradient('gnh', 'nh', -drmin, calc_vars)
    gradient('gnf', 'nf', -drmin, calc_vars)
    gradient('gpf', 'pf', -drmin, calc_vars)
    gradient('gni', 'ni', -drmin, calc_vars)
    gradient('gnz', 'nz', -drmin, calc_vars)
    gradient('gte', 'te', -drmin, calc_vars)
    gradient('gti', 'ti', -drmin, calc_vars)
    gradient('gvpar', 'vpar', -drmin, calc_vars)
    gradient('gvpol', 'vpol', -drmin, calc_vars)
    gradient('gvtor', 'vtor', -drmin, calc_vars)
    # gelong calculated differently than other normalized gradients
    gelong(calc_vars)

    gpf2(calc_vars)

    if hasattr(calc_vars.options, 'use_gnezero') and calc_vars.options.use_gnezero:
        calc_vars.gne.values[:, :] = 1e-12

    if hasattr(calc_vars.options, 'use_gtezero') and calc_vars.options.use_gtezero:
        calc_vars.gte.values[:, :] = 1e-12

    if hasattr(calc_vars.options, 'use_gneabs') and calc_vars.options.use_gneabs:
        calc_vars.gne.values = np.absolute(calc_vars.gne.values)

    if hasattr(calc_vars.options, 'use_gnethreshold') and calc_vars.options.use_gnethreshold:
        calc_vars.gne.values[:, :] = 1
        calc_vars.gte.values[:, :] = 1e-12
        calc_vars.gti.values[:, :] = 1e-12
        calc_vars.gni.values[:, :] = 1e-12

    if hasattr(calc_vars.options, 'use_gtethreshold') and calc_vars.options.use_gtethreshold:
        calc_vars.gne.values = np.absolute(calc_vars.gne.values)
        calc_vars.gne.values[:, :] = 1e-12
        calc_vars.gte.values[:, :] = 1
        calc_vars.gti.values[:, :] = 1e-12
        calc_vars.gni.values[:, :] = 1e-12

    if hasattr(calc_vars.options, 'use_etgm_btor') and calc_vars.options.use_etgm_btor:
        calc_vars.gbu.values = calc_vars.gbtor.values


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

    bu_btor(calc_vars)
    tau(calc_vars)
    tauh(calc_vars)
    eps(calc_vars)
    p(calc_vars)
    beta(calc_vars)
    betae(calc_vars)
    betaeu(calc_vars)
    csound(calc_vars)
    csound_a(calc_vars)
    loge(calc_vars)
    nuei(calc_vars)
    nuei2(calc_vars)
    vthe(calc_vars)
    vthi(calc_vars)
    wtransit(calc_vars)
    wbounce(calc_vars)
    nuste(calc_vars)
    nusti(calc_vars)
    gyrfe(calc_vars)
    gyrfeu(calc_vars)
    gyrfi(calc_vars)
    gyrfiu(calc_vars)
    lare(calc_vars)
    lareu(calc_vars)
    rhos(calc_vars)
    rhosu(calc_vars)
    shear(calc_vars)
    shat(calc_vars)
    shat_gxi(calc_vars)
    alphamhd(calc_vars)
    alphamhdu(calc_vars)
    betapu(calc_vars)
    xetgm_const(calc_vars)
    etae(calc_vars)
    etai(calc_vars)
    epsne(calc_vars)
    vei_nc(calc_vars)
    # curlh(calc_vars)
    # curoh(calc_vars)
    xke(calc_vars)
    xki(calc_vars)


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
