"""Calculates variable values

Variable Classifications:
* Base variable: A non-gradient variable that is either used as input to MMM,
  or is needed to calculate an MMM input variable.  Base variables do not
  depend on values of gradient or additional variables.
* gradient variable: A variable that is a normalized gradient, and may be used
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
from scipy.interpolate import interp1d, Akima1DInterpolator

# Local Packages
import settings
import modules.constants as constants
import modules.datahelper as datahelper


_gradients = set()  # Stores the names of calculated gradient variables


def differential(x, method=''):
    '''
    Take the radial differential of variable x

    Parameters:
    * x (np.ndarray): values to take differential of

    Raises:
    * ValueError: If settings.GRADIENT_METHOD does not match a valid gradient method
    '''
    if not method:
        method = settings.GRADIENT_METHOD

    if method == 'interpolate':
        return np.diff(x, axis=0)

    elif method == 'traditional':
        dx = np.zeros_like(x)
        tmp = np.diff(x, axis=0)
        dx[:-1, :] = tmp
        dx[1:, :] += tmp
        return dx

    elif method == 'ptsolver':
        dx = np.zeros_like(x)
        tmp = np.diff(x, axis=0)
        dx[1:, :] += tmp
        dx[0, :] = dx[1, :]
        return dx

    elif method == 'akima':
        dx = np.zeros_like(x)
        tmp = np.diff(x, axis=0)
        dx[:-1, :] = tmp
        dx[1:, :] += tmp
        return dx  # I don't think akima uses this, but we need drmin to be nonzero

    else:
        raise ValueError(
            f'\'{settings.GRADIENT_METHOD}\' is not a valid argument for the gradient method\n'
            f'Possible choices: \'interpolate\', \'traditional\', \'akima\', \'ptsolver\''
        )


def differentiate(yname, xname, calc_vars, method=''):
    '''
    Take the derivative of variable y with respect to variable x

    Parameters:
    * yname (str): name of variable to use for y
    * xname (str): name of variable to use for x
    * calc_vars (InputVariables): Object containing variable data
    '''
    if not method:
        method = settings.GRADIENT_METHOD

    if method == 'akima':
        return akima_differentiate(yname, xname, calc_vars)

    dy = differential(getattr(calc_vars, yname).values, method)

    if xname == 'rmin' and not method:
        dx = calc_vars.drmin.values
    else:
        dx = differential(getattr(calc_vars, xname).values, method)

    dy_dx = dy / dx

    if method == 'interpolate':

        rhox = calc_vars.x.values[:, 0]
        rhoxb = calc_vars.xb.values[:, 0]  # includes origin

        # interpolate from rhox grid to rhoxb grid
        set_interp = interp1d(rhox, dy_dx, kind=settings.INTERPOLATION_METHOD, fill_value="extrapolate", axis=0)
        dy_dx = set_interp(rhoxb)

    return dy_dx


def akima_differentiate(yname, xname, mmm_vars):

    yvals = getattr(mmm_vars, yname).values
    xvals = getattr(mmm_vars, xname).values
    dy_dx = np.zeros_like(yvals)
    time_count = dy_dx.shape[1]

    # Can't vectorize Akima derivative
    for i in range(time_count):
        set_interp = Akima1DInterpolator(xvals[:, i], yvals[:, i], axis=0)
        dy_dx[:, i] = set_interp(xvals[:, i], 1)
    return dy_dx

def gradient(gvar_name, var_name, gsign, calc_vars, method=''):
    '''
    Calculates the normalized gradient using interpolation

    After the gradient value is calculated, optional smoothing is applied, and
    then the gradient is checked for min and nan values.  The overall sign of
    the gradient equation is determined by the sign given for drmin.

    Parameters:
    * gvar_name (str): The name of the variable to store the gradient result in
    * var_name (str): The name of the variable to take the gradient of
    * gsign (int): sign on gradient, should be +1 or -1
    * calc_vars (InputVariables): Object containing variable data
    '''

    _gradients.add(gvar_name)
    rmaj = calc_vars.rmaj.values

    # get variables related to the gradient from variable names
    gvar = getattr(calc_vars, gvar_name)
    var = getattr(calc_vars, var_name)

    dy_dx = differentiate(var_name, 'rmin', calc_vars, method)
    gvar.set(values=gsign * rmaj * dy_dx / var.values, units='')

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
def dwtor_dr(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    wtor = calc_vars.omega.values
    gwtor = calc_vars.gwtor.values
    rmaj = calc_vars.rmaj.values

    return -gwtor * wtor / rmaj


@calculation
def dvtor_dr(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    vtor = calc_vars.vtor.values
    gvt = calc_vars.gvtor.values
    rmaj = calc_vars.rmaj.values

    return -gvt * vtor / rmaj


@calculation
def dvtor_dwtor(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    dvtor_dr = calc_vars.dvtor_dr.values
    dwtor_dr = calc_vars.dwtor_dr.values

    return dvtor_dr / dwtor_dr


@calculation
def ah(calc_vars):
    '''Mean Atomic Mass of Hydrogenic Ions (Hydrogen + Deuterium)'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values
    nt = calc_vars.nt.values

    return (nh0 + 2 * nd + 3 * nt) / (nh0 + nd + nt)


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
def agxi_1(calc_vars):
    '''Beta'''
    a = calc_vars.rmin.values[-1, :]
    gxi = calc_vars.gxi.values
    elong = calc_vars.elong.values
    gr = gxi + (1 / a - gxi) / elong

    return 1 / (a * gr)


@calculation
def agxi_1b(calc_vars):
    '''Beta'''
    # a = calc_vars.rmin.values[-1, :]
    # gxi = calc_vars.gxi.values
    # gxi2 = calc_vars.gxi2.values

    a = calc_vars.rmin.values[-1, :]
    gxi = calc_vars.gxi.values
    elong = calc_vars.elong.values
    gr = gxi + (1 / a - gxi) / elong**2

    return 1 / (a * gr)

    return 1 / (a * gxi)


@calculation
def agxi2_1(calc_vars):
    '''Beta'''
    # a = calc_vars.rmin.values[-1, :]
    # gxi = calc_vars.gxi.values
    # gxi2 = calc_vars.gxi2.values

    a = calc_vars.rmin.values[-1, :]
    gxi = calc_vars.gxi.values

    return 1 / (a * gxi)


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
def bftorsqrt(mmm_vars):
    """Poloidal Magnetic Field"""
    bftor = mmm_vars.bftor.values
    return bftor**(0.5)


@calculation
def bpol(calc_vars):
    '''Poloidal Magnetic Field'''
    # btor = calc_vars.btor.values
    # q = calc_vars.q.values
    # rmaj = calc_vars.rmaj.values
    # rmin = calc_vars.rmin.values
    # xb = calc_vars.xb.values
    # rhochi = calc_vars.rhochi.values
    # return rhochi / rmaj * btor / q
    # # return rhochi[-1, :] * xb / rmaj * btor / q
    # # return rmin[-1, :] * xb / rmaj * btor / q
    # return rmin / rmaj * btor / q

    bftor = calc_vars.bftor.values[-1, :]
    q = calc_vars.q.values
    rho = calc_vars.rho.values
    lpol = calc_vars.lpol.values
    dvol = calc_vars.dvol.values
    dvol_drho = calc_vars.dvol_drho.values
    drho = np.diff(rho, axis=0)[0, 0]

    return 2 * bftor * rho * lpol / np.maximum(dvol / drho * q, 1e-16)
    # return 2 * bftor * rho * lpol / (dvol_drho * q)


@calculation
def btor(calc_vars):
    '''Toroidal Magnetic Field'''
    bzxr = calc_vars.bzxr.values
    rmaj = calc_vars.rmaj.values

    return bzxr / rmaj


# @calculation
# def bu(calc_vars):
#     '''Magnetic Field (unit)'''
#     btor0 = calc_vars.btor.values[0, :]
#     rhochi = calc_vars.rhochi.values
#     rmin = calc_vars.rmin.values
#     x = calc_vars.x.values[:, 0]  # same for all time values
#     xb = calc_vars.xb.values[:, 0]  # same for all time values

#     drho_drmin = np.diff(rhochi, axis=0) / np.diff(rmin, axis=0)

#     # interpolate from x grid to xb grid
#     set_interp = interp1d(x, drho_drmin, kind=settings.INTERPOLATION_METHOD, fill_value="extrapolate", axis=0)
#     dxrho = set_interp(xb)

#     bu = np.empty_like(dxrho)
#     bu[1:, :] = btor0 * rhochi[1:, :] / rmin[1:, :] * dxrho[1:, :]
#     bu[0, :] = bu[1, :]

#     return bu


@calculation
def bu(mmm_vars):
    """Magnetic Field (unit)"""
    bftorsqrt = mmm_vars.bftorsqrt.values
    rmin = mmm_vars.rmin.values
    dbftorsqrt_dr = differentiate('bftorsqrt', 'rmin', mmm_vars)
    return bftorsqrt * dbftorsqrt_dr / np.maximum(np.pi * rmin, 1e-3)


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
    darea = calc_vars.darea.values
    rmin = calc_vars.rmin.values

    return curdlh * darea
    # return curdlh * constants.PI * rmin**2


@calculation
def curoh(calc_vars):
    '''OH Current'''
    curdoh = calc_vars.curdoh.values
    darea = calc_vars.darea.values

    icur = np.zeros((darea.shape[1], 1))

    for i in range(len(icur)):
        icur[i] = np.vdot(curdoh[:, i], darea[:, i])

    return icur.flatten()


@calculation
def curohrho(calc_vars):
    '''OH Current'''
    curdoh = calc_vars.curdoh.values
    darea = calc_vars.darea.values
    t = calc_vars.options.time_idx

    icur = np.zeros((darea.shape[0], 1))

    for i in range(len(icur)):
        icur[i] = np.vdot(curdoh[:i, t], darea[:i, t])

    return icur.flatten()


@calculation
def dbp(calc_vars):
    '''dBpol / drho'''
    return differentiate('bpol', 'rmin', calc_vars)
    # return akima_differentiate('bpol', 'rho', calc_vars)


@calculation
def d2bp(calc_vars):
    '''d^2Bpol / drho^2'''
    return differentiate('dbp', 'rmin', calc_vars)
    # return akima_differentiate('dbp', 'rho', calc_vars)


@calculation
def dvol_drho(calc_vars):
    ''''''
    return differentiate('dvol', 'rho', calc_vars)
    # return akima_differentiate('dbp', 'rho', calc_vars)


@calculation
def drmin(calc_vars):
    '''differential rmin'''
    return differential(calc_vars.rmin.values)


@calculation
def tmhdf(calc_vars):
    '''d^2Bpol / drho^2'''
    pmhdf = calc_vars.pmhdf.values
    nf = calc_vars.nf.values

    # partial derivative along radial dimension
    return pmhdf / nf / constants.ZCKB


@calculation
def tfpa(calc_vars):
    '''Fast Ion Temperature (parallel)'''
    ufastpa = calc_vars.ufastpa.values
    nf = calc_vars.nf.values

    return ufastpa / nf


@calculation
def tfpp(calc_vars):
    '''Fast Ion Temperature (perpendicular)'''
    ufastpp = calc_vars.ufastpp.values
    nf = calc_vars.nf.values

    return ufastpp / nf


@calculation
def tf(calc_vars):
    '''Fast Ion Temperature'''
    tfpa = calc_vars.tfpa.values
    tfpp = calc_vars.tfpp.values

    return tfpa + 0.5 * tfpp


@calculation
def tfast(calc_vars):
    '''Fast ion temp from EBEAM_D'''

    return calc_vars.ebeam.values / 1.5


@calculation
def rmajm(calc_vars):
    '''Fast Ion Temperature'''
    # return calc_vars.rmajm.values - 1
    return calc_vars.rmajm.values


@calculation
def ebeam2(calc_vars):
    '''Radial Electric Field (Pressure Term)'''
    ebeam = calc_vars.ebeam.values

    return ebeam / 1.5


@calculation
def ebeamr(calc_vars):
    '''Radial Electric Field (Pressure Term)'''
    ebeam = calc_vars.ebeam.values
    tf = calc_vars.tmhdf.values

    return ebeam / 1.5 / tf


@calculation
def ebeamsum(calc_vars):
    '''Radial Electric Field (Pressure Term)'''
    ebeampp = calc_vars.ebeampp.values
    ebeampl = calc_vars.ebeampl.values

    return ebeampp + ebeampl


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
    set_interp = interp1d(x, dpdr_x, kind=settings.INTERPOLATION_METHOD, fill_value="extrapolate", axis=0)
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
    '''Elongation gradient (delong / deps)'''
    _gradients.add('gelong')
    # rmaj0 = calc_vars.rmaj.values[0, :]
    # rmaj = calc_vars.rmaj.values
    # return rmaj0 * differentiate('elong', 'rmin', calc_vars)
    return differentiate('elong', 'eps', calc_vars)


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
def gxi(mmm_vars):
    """Grad of normalized xb (rho)"""
    return differentiate('xb', 'rmin', mmm_vars)
    # return akima_differentiate('xb', 'rmin', mmm_vars)


@calculation
def icur(calc_vars):
    '''delta darea'''
    darea = calc_vars.darea.values
    jcur = calc_vars.jcur.values

    icur = np.zeros((darea.shape[1], 1))

    for i in range(len(icur)):
        icur[i] = np.vdot(jcur[:, i], darea[:, i])

    return icur.flatten()


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
def mmmtime(calc_vars):
    '''MMM Runtime in TRANSP'''
    walltime = calc_vars.walltime.values
    cpmcfi = calc_vars.cpmcfi.values
    cpout = calc_vars.cpout.values
    cptrk = calc_vars.cptrk.values
    cpgeom = calc_vars.cpgeom.values
    cpmhdq = calc_vars.cpmhdq.values
    cpxgpl = calc_vars.cpxgpl.values
    cpsc0 = calc_vars.cpsc0.values
    cpbmax = calc_vars.cpbmax.values

    return walltime - cpmcfi - cpout - cptrk - cpgeom - cpmhdq - cpxgpl - cpsc0


@calculation
def nh0(calc_vars):
    '''Hydrogen Ion Density'''
    nd = calc_vars.nd.values
    ne = calc_vars.ne.values
    nf = calc_vars.nf.values
    nz = calc_vars.nz.values
    nt = calc_vars.nt.values
    zz = calc_vars.zz.values
    
    return ne - zz * nz - nf - nd - nt


@calculation
def nh(calc_vars):
    '''Total Hydrogenic Ion Density'''
    nh0 = calc_vars.nh0.values
    nd = calc_vars.nd.values
    nt = calc_vars.nt.values

    return nh0 + nd + nt


@calculation
def ni(calc_vars):
    '''Thermal Ion Density'''
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    nf = calc_vars.nf.values

    # TRANSP Definition
    return nh + nz + nf


@calculation
def ne(calc_vars):
    '''Thermal Ion Density'''
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    zz = calc_vars.zz.values
    nf = calc_vars.nf.values

    # return nh + zz*nz + nf
    return nh + zz*nz


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
def nuste(calc_vars):
    '''Electron Collisionality'''

    eps = calc_vars.eps.values
    nuei = calc_vars.nuei.values
    q = calc_vars.q.values
    rmaj = calc_vars.rmaj.values
    vthe = calc_vars.vthe.values

    return nuei * eps**(-1.5) * q * rmaj / vthe


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
    '''Rhos'''

    csound = calc_vars.csound.values
    gyrfi = calc_vars.gyrfi.values

    return csound / gyrfi


@calculation
def rhosu(calc_vars):
    '''Rhos (Unit)'''

    csound = calc_vars.csound.values
    gyrfiu = calc_vars.gyrfiu.values

    return csound / gyrfiu


@calculation
def shat(calc_vars):
    '''Effective Magnetic Shear'''

    elong = calc_vars.elong.values
    shear = calc_vars.shear.values

    signs = np.ones_like(shear)
    signs[shear < 0] = -1
    
    return np.maximum(2 * shear - 1 + (elong * (shear - 1))**2, shear**2)**(0.5) * signs


@calculation
def shat_gxi(calc_vars):
    '''Effective Magnetic Shear'''

    a = calc_vars.rmin.values[-1, :]
    gxi = calc_vars.gxi.values
    shear = calc_vars.shear.values

    signs = np.ones_like(shear)
    signs[shear < 0] = -1

    return np.maximum(2 * shear - 1 + ((a * gxi) * (shear - 1))**2, shear**2)**(0.5) * signs


@calculation
def shear(calc_vars):
    '''Magnetic Shear'''
    gq = calc_vars.gq.values
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values

    return gq * rmin / rmaj


@calculation
def te_ti(calc_vars):
    '''Temperature Ratio te / ti'''
    te = calc_vars.te.values
    ti = calc_vars.ti.values

    return te / ti


@calculation
def ti_te(calc_vars):
    '''Temperature Ratio ti / te'''
    te_ti = calc_vars.te_ti.values

    return 1 / te_ti


@calculation
def tf_te(calc_vars):
    '''Temperature Ratio tf / te'''
    te = calc_vars.te.values
    tf = calc_vars.tf.values

    return tf / te


@calculation
def tbtbe(calc_vars):
    '''Temperature Ratio ti / te'''
    btbe = calc_vars.btbe.values
    btor = calc_vars.bzxr.values
    nf = calc_vars.nf.values
    rmaj = calc_vars.rmaj.values
    zcmu0 = constants.ZCMU0
    zckb = constants.ZCKB
    pi = constants.PI

    gr2i = calc_vars.gr2i.values
    dvol = np.maximum(calc_vars.dvol.values, 0.001)
    bzxr = calc_vars.bzxr.values

    # print(dvol)

    return btbe * btor**2 / (2 * zcmu0 * zckb * nf) * 1.8
    # return btbe * btor**2 / (2 * zcmu0 * zckb * nf) * 8 / 15 * pi
    # return 4 / 15 * 1000 * 3.1415 * rmaj**2 / (zcmu0 * nf * gr2i * dvol)


@calculation
def vA(calc_vars):
    '''Alfven Velocity'''
    zcmp = constants.ZCMP
    zcmu0 = constants.ZCMU0
    btor = calc_vars.btor.values
    ai = calc_vars.ai.values
    ni = calc_vars.ne.values  # Assuming quasi-neutrality

    return btor / (zcmu0 * ni * ai * zcmp)**(0.5)


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
def vtor(calc_vars):
    '''Parallel Velocity'''
    rmaj = calc_vars.rmaj.values
    omega = calc_vars.omega.values

    return rmaj * omega


@calculation
def wbe(calc_vars):
    '''Bounce Frequency'''
    rmaj = calc_vars.rmaj.values
    rmin = calc_vars.rmin.values
    wte = calc_vars.wte.values

    return (rmin / (2 * rmaj))**(0.5) * wte


@calculation
def wte(calc_vars):
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
        set_interp = interp1d(x, dfdr_x, kind=settings.INTERPOLATION_METHOD, fill_value="extrapolate", axis=0)
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

    # xke = condepr + condewnc + xkepaleo
    xke = condepr

    return xke


@calculation
def fke(calc_vars):
    '''Total Electron Thermal Diffusivity from CDF'''
    xke = calc_vars.xke.values
    gte = calc_vars.gte.values
    te = calc_vars.te.values
    rmaj = calc_vars.rmaj.values

    return xke * gte * te / rmaj



@calculation
def xkeetgm(calc_vars):
    '''ETGM Electron Thermal Diffusivity from CDF'''
    xkemmm = calc_vars.xkemmm.values
    xkemtm = calc_vars.xkemtm.values
    xkew20 = calc_vars.xkew20.values
    xkedrbm = calc_vars.xkedrbm.values

    return xkemmm - xkemtm - xkew20 - xkedrbm


@calculation
def xki(calc_vars):
    '''Total Ion Thermal Diffusivity from CDF'''
    condipr = calc_vars.condipr.values
    # condiwnc = calc_vars.condiwnc.values

    # xki = condipr + condiwnc
    xki = condipr

    return xki


@calculation
def fki(calc_vars):
    '''Total Ion Thermal Diffusivity from CDF'''
    xki = calc_vars.xki.values
    gti = calc_vars.gti.values
    ti = calc_vars.ti.values
    rmaj = calc_vars.rmaj.values

    return xki * gti * ti / rmaj


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
def zave(calc_vars):
    '''Effective Charge'''
    nf = calc_vars.nf.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values
    zz = calc_vars.zz.values

    return (nh + nf + zz * nz) / (nh + nf + nz)
    # return (nh + zz * nz) / (nh + nz)


@calculation
def zni(calc_vars):
    '''Effective Charge'''
    zave = calc_vars.zave.values
    nf = calc_vars.nf.values
    nh = calc_vars.nh.values
    nz = calc_vars.nz.values

    return zave * (nh + nf + nz)
    # return zave * (nh + nz)


@calculation_output
def gmaW20(calc_vars, output_vars):
    '''Growth Rate of Most Unstable Mode'''
    gmaW20i = output_vars.gmaW20i.values
    gmaW20e = output_vars.gmaW20e.values
    gmaW20 = np.zeros_like(gmaW20i)
    gmaW20[gmaW20i >= gmaW20e] = gmaW20i[gmaW20i >= gmaW20e]
    gmaW20[gmaW20i < gmaW20e] = gmaW20e[gmaW20i < gmaW20e]

    return gmaW20


@calculation_output
def omgW20(calc_vars, output_vars):
    '''Frequency of Most Unstable Mode'''
    gmaW20i = output_vars.gmaW20i.values
    gmaW20e = output_vars.gmaW20e.values
    omgW20i = output_vars.omgW20i.values
    omgW20e = output_vars.omgW20e.values
    omgW20 = np.zeros_like(omgW20i)
    omgW20[gmaW20i >= gmaW20e] = omgW20i[gmaW20i >= gmaW20e]
    omgW20[gmaW20i < gmaW20e] = omgW20e[gmaW20i < gmaW20e]

    return omgW20


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
def gmanETGM(calc_vars, output_vars):
    '''ETGM Growth Rate resonance'''
    t = calc_vars.options.time_idx
    gmaETGM = output_vars.gmaETGM.values
    csound_a = calc_vars.csound_a.values[:, t]

    return gmaETGM / csound_a



@calculation_output
def omgnETGM(calc_vars, output_vars):
    '''ETGM Growth Rate resonance'''
    t = calc_vars.options.time_idx
    omgETGM = output_vars.omgETGM.values
    csound_a = calc_vars.csound_a.values[:, t]

    return -omgETGM / csound_a


@calculation_output
def omgdiffETGM(calc_vars, output_vars):
    '''ETGM Frequency resonance'''
    omgETGM = output_vars.omgETGM.values
    wdeETGM = output_vars.wdeETGM.values

    return omgETGM - wdeETGM


@calculation_output
def wdeETGM(calc_vars, output_vars):
    '''ETGM Frequency resonance'''
    t = calc_vars.options.time_idx
    rmaj = calc_vars.rmaj.values[:, t]
    csound = calc_vars.csound.values[:, t]
    kyrhos = output_vars.kyrhosETGM.values

    return 2 * kyrhos * csound / rmaj


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
    kpara = output_vars.kparaETGM.values

    return kpara * bu / (zcmu0 * zcmp * ai * ni)**(0.5)


@calculation_output
def wdeEPM(calc_vars, output_vars):
    '''wde'''
    t = calc_vars.options.time_idx
    csound = calc_vars.csound.values[:, t]
    rmaj = calc_vars.rmaj.values[:, t]
    kyrhos = output_vars.kyrhosEPM.values

    return 2 * csound * kyrhos / rmaj


@calculation_output
def wdfEPM(calc_vars, output_vars):
    '''wdf'''
    t = calc_vars.options.time_idx
    tf = calc_vars.tf.values[:, t]
    te = calc_vars.te.values[:, t]
    wdeEPM = output_vars.wdeEPM.values

    return -tf / te * wdeEPM


@calculation_output
def wseEPM(calc_vars, output_vars):
    '''wdf'''
    t = calc_vars.options.time_idx
    gne = calc_vars.gne.values[:, t]
    wdeEPM = output_vars.wdeEPM.values

    return 0.5 * gne * wdeEPM


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

    if calc_vars.options.save_model_outputs:
        if calc_vars.options.cmodel_etgm > 0:
            wdeETGM(calc_vars, output_vars)
            wde_gaveETGM(calc_vars, output_vars)
            wseETGM(calc_vars, output_vars)
            wsetaETGM(calc_vars, output_vars)
            wteETGM(calc_vars, output_vars)
            omgdiffETGM(calc_vars, output_vars)
            gmadiffETGM(calc_vars, output_vars)
            waETGM(calc_vars, output_vars)
            gmanETGM(calc_vars, output_vars)
            omgnETGM(calc_vars, output_vars)

        if calc_vars.options.cmodel_mtm > 0:
            if settings.SAVE_ADDITIONAL_VARIABLES:
                gmanMTM(calc_vars, output_vars)
                omgnMTM(calc_vars, output_vars)

        if calc_vars.options.cmodel_epm > 0:
            wdeEPM(calc_vars, output_vars)
            wdfEPM(calc_vars, output_vars)
            wseEPM(calc_vars, output_vars)

        if calc_vars.options.cmodel_w20 > 0:
            gmaW20(calc_vars, output_vars)
            omgW20(calc_vars, output_vars)


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

    dvol_drho(calc_vars)
    # nh0(calc_vars)
    drmin(calc_vars)
    nh(calc_vars)
    ni(calc_vars)
    # ne(calc_vars)
    ah(calc_vars)
    zeff(calc_vars)
    btor(calc_vars)
    bftorsqrt(calc_vars)
    bu(calc_vars)
    bpol(calc_vars)
    dbp(calc_vars)
    d2bp(calc_vars)
    vtor(calc_vars)
    vpar(calc_vars)
    tfpa(calc_vars)
    tfpp(calc_vars)
    tf(calc_vars)
    eps(calc_vars)  # for gelong

    # Calculations for testing
    # 
    # e_r_grp(calc_vars)
    # e_r_phi(calc_vars)
    # e_r_tht(calc_vars)
    # wexb(calc_vars)
    gxi(calc_vars)

    if hasattr(calc_vars.options, 'use_etgm_btor') and calc_vars.options.use_etgm_btor:
        calc_vars.bu.values = calc_vars.btor.values


def calculate_gradient_variables(calc_vars):
    '''
    Calculates gradient variables

    Since the input parameter calc_vars is a reference, no return value for
    calculations are needed.  gradient calculations do not depend on any
    additional variables.

    Parameters:
    * calc_vars (InputVariables): Object containing variable data
    '''

    # Each use of the following calculate_variable functions are passing in
    # the local function of the variable to be calculated, which shares the
    # same name as the variable it calculates.

    # Positive sign
    gradient('gq', 'q', 1, calc_vars)
    gradient('gbu', 'bu', 1, calc_vars)
    gradient('gbtor', 'btor', 1, calc_vars)

    # Negative sign
    gradient('gne', 'ne', -1, calc_vars)
    gradient('gnh', 'nh', -1, calc_vars)
    gradient('gnf', 'nf', -1, calc_vars)
    gradient('gni', 'ni', -1, calc_vars)
    gradient('gnz', 'nz', -1, calc_vars)
    gradient('gte', 'te', -1, calc_vars)
    gradient('gtfpa', 'tfpa', -1, calc_vars)
    gradient('gtfpp', 'tfpp', -1, calc_vars)
    gradient('gtf', 'tf', -1, calc_vars)
    gradient('gti', 'ti', -1, calc_vars)
    gradient('gvpol', 'vpol', -1, calc_vars)
    gradient('gvtor', 'vtor', -1, calc_vars)
    gradient('gvpar', 'vpar', -1, calc_vars)

    method = settings.SOLVER_GRADIENT_METHOD
    gradient('gti_solver', 'ti', -1, calc_vars, method)
    gradient('gte_solver', 'te', -1, calc_vars, method)
    gradient('gne_solver', 'ne', -1, calc_vars, method)
    gradient('gnz_solver', 'nz', -1, calc_vars, method)
    gradient('gvtor_solver', 'vtor', -1, calc_vars, method)
    gradient('gvpol_solver', 'vpol', -1, calc_vars, method)

    # gelong calculated differently than other normalized gradients
    gelong(calc_vars)

    if hasattr(calc_vars.options, 'use_gnezero') and calc_vars.options.use_gnezero:
        calc_vars.gne.values[:, :] = 1e-12

    if hasattr(calc_vars.options, 'use_gtezero') and calc_vars.options.use_gtezero:
        calc_vars.gte.values[:, :] = 1e-12

    if hasattr(calc_vars.options, 'use_gtizero') and calc_vars.options.use_gtizero:
        calc_vars.gti.values[:, :] = 1e-12

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

    if not settings.SAVE_ADDITIONAL_VARIABLES:
        return

    gradient('gwtor', 'omega', -1, calc_vars)
    dwtor_dr(calc_vars)
    dvtor_dr(calc_vars)
    dvtor_dwtor(calc_vars)
    agxi_1(calc_vars)
    agxi_1b(calc_vars)
    agxi2_1(calc_vars)
    rhochi(calc_vars)
    ai(calc_vars)
    te_ti(calc_vars)
    ti_te(calc_vars)
    tf_te(calc_vars)
    p(calc_vars)
    beta(calc_vars)
    betae(calc_vars)
    betaeu(calc_vars)
    csound(calc_vars)
    csound_a(calc_vars)
    loge(calc_vars)
    nuei(calc_vars)
    vthe(calc_vars)
    vthi(calc_vars)
    wte(calc_vars)
    wbe(calc_vars)
    nuste(calc_vars)
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
    vA(calc_vars)
    zave(calc_vars)
    zni(calc_vars)

    # curlh(calc_vars)
    # curoh(calc_vars)
    # curohrho(calc_vars)

    # Non essential calculations that may have been removed from memory
    if calc_vars.xkemmm is not None:
        xkeetgm(calc_vars)
    if calc_vars.xke is not None:
        xke(calc_vars)
        fke(calc_vars)
    if calc_vars.xki is not None:
        xki(calc_vars)
        fki(calc_vars)
    if calc_vars.walltime is not None:
        mmmtime(calc_vars)

    # icur(calc_vars)
    tmhdf(calc_vars)

    tbtbe(calc_vars)
    ebeamsum(calc_vars)
    rmajm(calc_vars)
    tfast(calc_vars)


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
