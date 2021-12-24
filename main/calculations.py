# Standard Packages
import sys
import copy
import inspect

# 3rd Party Packages
import numpy as np
from scipy.interpolate import interp1d  # TODO: use Akima1DInterpolator?

# Local Packages
import main.options as options
import main.constants as constants


def nh0(vars):
    '''Hydrogen Ion Density'''
    nd = vars.nd.values
    ne = vars.ne.values
    nf = vars.nf.values
    nz = vars.nz.values
    zimp = vars.zimp.values

    nh0 = ne - zimp * nz - nf - nd

    vars.nh0.set(values=nh0, units=vars.ne.units)


def nh(vars):
    '''Total Hydrogenic Ion Density'''
    nh0 = vars.nh0.values
    nd = vars.nd.values

    nh = nh0 + nd

    vars.nh.set(values=nh, units=vars.nd.units)


def ni(vars):
    '''Thermal Ion Density'''
    nh = vars.nh.values
    nz = vars.nz.values

    # TRANSP Definition
    ni = nh + nz

    # Alternate Definition
    # zimp = vars.zimp.values
    # ni = nh + zimp * nz

    vars.ni.set(values=ni, units=vars.ne.units)


def ahyd(vars):
    '''Mean atomic mass of hydrogenic ions (hydrogen + deuterium)'''
    nh0 = vars.nh0.values
    nd = vars.nd.values

    ahyd = (nh0 + 2 * nd) / (nh0 + nd)

    vars.ahyd.set(values=ahyd, units='')


def aimass(vars):
    '''Mean Atomic Mass of Thermal Ions'''
    ahyd = vars.ahyd.values
    aimp = vars.aimp.values
    nh = vars.nh.values
    nz = vars.nz.values

    aimass = (ahyd * nh + aimp * nz) / (nh + nz)

    vars.aimass.set(values=aimass, units='')


def tau(vars):
    '''Temperature Ratio te / ti'''
    te = vars.te.values
    ti = vars.ti.values

    tau = te / ti

    vars.tau.set(values=tau, units='')


def vpar(vars):
    '''Parallel Velocity'''
    bpol = vars.bpol.values
    btor = vars.btor.values
    vpol = vars.vpol.values
    vtor = vars.vtor.values

    vpar = vtor + vpol * bpol / btor

    vars.vpar.set(values=vpar, units=vars.vtor.units)


def zeff(vars):
    '''Effective Charge'''
    ne = vars.ne.values
    nf = vars.nf.values
    nh = vars.nh.values
    nz = vars.nz.values
    zimp = vars.zimp.values

    zeff = (nh + nf + zimp**2 * nz) / ne

    vars.zeff.set(values=zeff, units='')


def btor(vars):
    '''Toroidal Magnetic Field'''
    bz = vars.bz.values
    raxis = vars.rmaj.values[0, :]
    rmaj = vars.rmaj.values

    btor = raxis / rmaj * bz

    vars.btor.set(values=btor, units=vars.bz.units)


def bpol(vars):
    '''Poloidal Magnetic Field'''
    btor = vars.btor.values
    q = vars.q.values
    rmaj = vars.rmaj.values
    rmin = vars.rmin.values

    bpol = rmin / rmaj * btor / q

    vars.bpol.set(values=bpol, units=vars.btor.units)


def eps(vars):
    '''Inverse Aspect Ratio'''
    arat = vars.arat.values

    eps = 1 / arat

    vars.eps.set(values=eps, units='')


def p(vars):
    '''Plasma Pressure'''
    zckb = constants.ZCKB
    ne = vars.ne.values
    ni = vars.ni.values
    te = vars.te.values
    ti = vars.ti.values

    p = (ne * te + ni * ti) * zckb

    vars.p.set(values=p, units='Pa')


def beta(vars):
    '''Beta'''
    zcmu0 = constants.ZCMU0
    btor = vars.btor.values
    p = vars.p.values

    beta = 2 * zcmu0 * p / btor**2

    vars.beta.set(values=beta, units='')


def betae(vars):
    '''Electron Beta'''
    zckb = constants.ZCKB
    zcmu0 = constants.ZCMU0
    btor = vars.btor.values
    ne = vars.ne.values
    te = vars.te.values

    betae = 2 * zcmu0 * ne * te * zckb / btor**2

    vars.betae.set(values=betae, units='')


def loge(vars):
    '''Electron Coulomb Logarithm'''

    # TODO: Need to add equations for different TE ranges
    ne = vars.ne.values
    te = vars.te.values

    # NRL Plasma Formulary Definition
    loge = 37.8 - np.log(ne**(1 / 2) / te)

    # TRANSP definition (equivalent)
    # zeff = vars.zeff.values
    # loge = 39.23 - np.log(zeff*ne**(1 / 2) / te)

    vars.loge.set(values=loge, units='')


def nuei(vars):
    '''Collision Frequency (NU_{ei})'''
    zcf = constants.ZCF
    ne = vars.ne.values
    te = vars.te.values
    zeff = vars.zeff.values
    loge = vars.loge.values

    nuei = zcf * 2**(1 / 2) * ne * loge * zeff / te**(3 / 2)

    vars.nuei.set(values=nuei, units='s^-1')


def nuei2(vars):
    '''OLD NOTE: Not sure what to call this, but it leads to the approx the correct NUSTI'''
    zcf = constants.ZCF
    ni = vars.ni.values
    ti = vars.ti.values
    zeff = vars.zeff.values
    loge = vars.loge.values

    nuei2 = zcf * 2**(1 / 2) * ni * loge * zeff / ti**(3 / 2)

    vars.nuei2.set(values=nuei2, units='s^-1')


def vthe(vars):
    '''Thermal Velocity of Electrons'''
    zckb = constants.ZCKB
    zcme = constants.ZCME
    te = vars.te.values

    vthe = (2 * zckb * te / zcme)**(1 / 2)

    vars.vthe.set(values=vthe, units='m/s')


def vthi(vars):
    '''Thermal Velocity of Ions'''
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    aimass = vars.aimass.values
    ti = vars.ti.values

    vthi = (zckb * ti / (zcmp * aimass))**(1 / 2)

    vars.vthi.set(values=vthi, units='m/s')


def nuste(vars):
    '''
    Electron Collisionality (NU^{*}_{e}) TODO: units? OLD NOTE: This is in
    approximate agreement with NUSTE in TRANSP.  One source of the
    disagreement is likely because the modmmm7_1.f90 Coulomb logarithm
    (loge) does not match perfectly with the TRANSP version (CLOGE).
    '''

    eps = vars.eps.values
    nuei = vars.nuei.values
    q = vars.q.values
    rmaj = vars.rmaj.values
    vthe = vars.vthe.values

    nuste = nuei * eps**(-3 / 2) * q * rmaj / vthe

    vars.nuste.set(values=nuste, units='')


def nusti(vars):
    '''
    Ion Collisionality (NUSTI = NU^{*}_{i}) TODO: Units OLD NOTE: This is approx
    correct, but agreement is also somewhat time-dependent.  The issue is
    possibly due to the artificial AIMASS that we are using.  We likely also need
    to use the coulomb logarithm for ions as well.
    '''

    zcme = constants.ZCME
    zcmp = constants.ZCMP
    eps = vars.eps.values
    q = vars.q.values
    nuei2 = vars.nuei2.values
    rmaj = vars.rmaj.values
    vthi = vars.vthi.values

    nusti = nuei2 * eps**(-3 / 2) * q * rmaj / (2 * vthi) * (zcme / zcmp)**(1 / 2)

    vars.nusti.set(values=nusti, units='')


def gyrfi(vars):
    '''Ion Gyrofrequency'''
    zce = constants.ZCE
    zcmp = constants.ZCMP
    aimass = vars.aimass.values
    btor = vars.btor.values

    gyrfi = zce * btor / (zcmp * aimass)

    vars.gyrfi.set(values=gyrfi, units='s^-1')


def gmax(vars):
    '''Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm.f90)'''
    eps = vars.eps.values
    q = vars.q.values
    rmaj = vars.rmaj.values
    gyrfi = vars.gyrfi.values
    vthi = vars.vthi.values

    gmax = rmaj / (vthi / gyrfi * q / eps)

    vars.gmax.set(values=gmax, units='')


def shear(vars):
    '''Magnetic Shear'''
    gq = vars.gq.values
    rmaj = vars.rmaj.values
    rmin = vars.rmin.values

    shear = gq * rmin / rmaj

    vars.shear.set(values=shear, units='')


def shat(vars):
    '''Effective Magnetic Shear'''
    elong = vars.elong.values
    shear = vars.shear.values

    shat = (2 * shear - 1 + (elong * (shear - 1))**2)**(1 / 2)
    shat[shat < 0] = 0

    vars.shat.set(values=shat, units='')


def alphamhd(vars):
    '''Alpha MHD (Weiland Definition)'''
    betae = vars.betae.values
    gne = vars.gne.values
    gni = vars.gni.values
    gte = vars.gte.values
    gti = vars.gti.values
    q = vars.q.values
    te = vars.te.values
    ti = vars.ti.values

    alphamhd = q**2 * betae * (gne + gte + ti / te * (gni + gti))

    vars.alphamhd.set(values=alphamhd, units='')


def gave(vars):
    '''Average Magnetic Surface Curvature'''
    shear = vars.shear.values
    alphamhd = vars.alphamhd.values

    gave = 2 / 3 + 5 / 9 * shear - 5 / 12 * alphamhd

    vars.gave.set(values=gave, units='')


def etae(vars):
    gte = vars.gte.values
    gne = vars.gne.values

    etae = gte / gne

    vars.etae.set(values=etae, units='')


def etai(vars):
    gti = vars.gti.values
    gni = vars.gni.values

    etai = gti / gni

    vars.etai.set(values=etai, units='')


def test(vars):
    '''Test variables are just used for testing calculations, and are not sent to the MMM driver'''
    nh = vars.nh.values
    nd = vars.nd.values

    ni = nh + nd

    vars.test.set(values=ni)


def test2(vars):
    '''Test variables are just used for testing calculations, and are not sent to the MMM driver'''
    gti = vars.gti.values
    gni = vars.gtest.values

    test2 = gti / gni

    vars.test2.set(values=test2)


def calculate_gradient(gvar_name, var_name, drmin, vars):
    rmaj = vars.rmaj.values
    x = vars.x.values[:, 0]
    xb = vars.xb.values[:, 0]  # includes origin

    # get variables related to the gradient from variable names
    gvar = getattr(vars, gvar_name)
    var = getattr(vars, var_name)

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

    gvar.clamp_gradient(100)
    gvar.set_minvalue()

    if opts.reject_outliers:
        gvar.reject_outliers()

    gvar.check_for_nan()


def calculate_variable(var_function, vars):
    '''Calculate the variable specified by it's corresponding function'''

    var_function(vars)

    # Get the variable name specified by var_function
    var_name = var_function.__name__

    if options.instance.apply_smoothing:
        getattr(vars, var_name).apply_smoothing(options.instance.input_points)

    getattr(vars, var_name).set_minvalue()

    if options.instance.reject_outliers:
        getattr(vars, var_name).reject_outliers()

    getattr(vars, var_name).check_for_nan()


def calculate_additional_variables(vars):

    calculate_variable(tau, vars)
    calculate_variable(eps, vars)
    calculate_variable(p, vars)
    calculate_variable(beta, vars)
    calculate_variable(betae, vars)
    calculate_variable(loge, vars)
    calculate_variable(nuei, vars)
    calculate_variable(nuei2, vars)
    calculate_variable(vthe, vars)
    calculate_variable(vthi, vars)
    calculate_variable(nuste, vars)
    calculate_variable(nusti, vars)
    calculate_variable(gyrfi, vars)
    calculate_variable(gmax, vars)
    calculate_variable(shear, vars)
    calculate_variable(shat, vars)
    calculate_variable(alphamhd, vars)
    calculate_variable(gave, vars)
    calculate_variable(etae, vars)
    calculate_variable(etai, vars)


def calculate_inputs(cdf_vars):
    '''
    Calculates new variables needed for MMM and data display

    Note that each use of the calculate_variable function below is passing in
    the function of the variable to be calculated, which shares the same name
    as the variable it calculates.

    Parameters:
    * cdf_vars (InputVariables): Variables object containing data from a CDF

    Returns:
    * vars (InputVariables): Variables object containing calculation results
    '''

    vars = copy.deepcopy(cdf_vars)

    '''
    Base Calculations:
    * These are either MMM input variables or are needed for input calculations
    * Some calculations depend on values from previous calculations
    * Base calculations do not depend on gradient values
    '''
    calculate_variable(nh0, vars)
    calculate_variable(nh, vars)
    calculate_variable(ni, vars)
    calculate_variable(ahyd, vars)
    calculate_variable(aimass, vars)
    calculate_variable(zeff, vars)
    calculate_variable(btor, vars)
    calculate_variable(bpol, vars)
    calculate_variable(vpar, vars)

    '''
    Gradient Calculations:
    * The sign on drmin (differential rmin) sets the sign of the gradient equation
    * Note that everything but gq has a negative drmin
    '''
    drmin = np.diff(vars.rmin.values, axis=0)
    calculate_gradient('gne', 'ne', -drmin, vars)
    calculate_gradient('gnh', 'nh', -drmin, vars)
    calculate_gradient('gni', 'ni', -drmin, vars)
    calculate_gradient('gnz', 'nz', -drmin, vars)
    calculate_gradient('gq', 'q', drmin, vars)
    calculate_gradient('gte', 'te', -drmin, vars)
    calculate_gradient('gti', 'ti', -drmin, vars)
    calculate_gradient('gvpar', 'vpar', -drmin, vars)
    calculate_gradient('gvpol', 'vpol', -drmin, vars)
    calculate_gradient('gvtor', 'vtor', -drmin, vars)

    # Additional Calculations (some depend on gradient calculations)
    calculate_additional_variables(vars)

    # Test variables are just used for testing calculations, and are not sent to the MMM driver
    calculate_variable(test, vars)
    calculate_gradient('gtest', 'test', -drmin, vars)
    calculate_variable(test2, vars)

    return vars


def get_calculated_vars():
    '''Returns (list of str): function names of calculated variables in this module (other than gradient calculations)'''
    return [o[0] for o in inspect.getmembers(sys.modules[__name__]) if inspect.isfunction(o[1]) and 'calculate' not in o[0]]
