# Standard Packages
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np 
# Local Packages
from mmm_package import variables, constants

# Set VPOL using VPOLD or VPOLH, if possible # TODO: handle this better
def calculate_vpol(vars):
    vpol = np.zeros((vars.xb.values.shape[0], vars.time.values.shape[0]))
    if vars.vpold.values is not None:
        vpol = vars.vpold.values
    elif vars.vpolh.values is not None:
        vpol = vars.vpolh.values

    vars.vpol.set_variable(vpol, 'M/SEC')

# Hydrogenic Ion Density
def calculate_nh(vars):
    ne = vars.ne.values
    nf = vars.nf.values
    nz = vars.nz.values
    zimp = vars.zimp.values

    nh = ne - zimp * nz - nf

    vars.nh.set_variable(nh, vars.ne.get_units())

# AIMASS (also AHYD?)
def calculate_aimass(vars):
    nh = vars.nh.values
    nd = vars.nd.values

    aimass = (nh + 2 * nd) / (nh + nd)

    vars.aimass.set_variable(aimass, vars.nh.get_units())

# Minor Radius, and set origin value to 0
def calculate_rmin(vars):
    arat = vars.arat.values
    rmaj = vars.rmaj.values

    rmin = (rmaj / arat)
    rmin[0, :] = np.zeros(vars.time.values.shape[0])

    vars.rmin.set_variable(rmin, vars.rmaj.get_units())

# Temperature Ratio
def calculate_tau(vars):
    te = vars.te.values
    ti = vars.ti.values

    tau = te / ti

    vars.tau.set_variable(tau, None)

# VTOR
def calculate_vtor(vars):
    rmaj = vars.rmaj.values
    omega = vars.omega.values

    vtor = rmaj * omega

    vars.vtor.set_variable(vtor, 'M/SEC')

# VPAR (setting equal to VTOR for now)
def calculate_vpar(vars):
    vars.vpar.set_variable(vars.vtor.values, vars.vtor.units)

# Effective Charge
def calculate_zeff(vars):
    ne = vars.ne.values
    nf = vars.nf.values
    nh = vars.nh.values
    nz = vars.nz.values
    zimp = vars.zimp.values

    zeff = (nh + nf + zimp**2 * nz) / ne

    vars.zeff.set_variable(zeff, None)

# Calculate BTOR
def calculate_btor(vars):
    bz = vars.bz.values
    raxis = vars.rmaj.values[0, :]
    rmaj = vars.rmaj.values

    btor = rmaj / raxis * bz

    vars.btor.set_variable(btor, vars.bz.get_units())

# Calculate Inverse Aspect Ratio
def calculate_eps(vars):
    arat = vars.arat.values

    eps = arat**(-1)

    vars.eps.set_variable(eps, None)

# Plasma Pressure
def calculate_p(vars):
    zckb = constants.ZCKB
    ne = vars.ne.values
    ni = vars.ni.values
    te = vars.te.values
    ti = vars.ti.values

    p = (ne * te + ni * ti) * zckb

    vars.p.set_variable(p, 'PA')

# Beta (in %)
def calculate_beta(vars):
    zcmu0 = constants.ZCMU0
    btor = vars.btor.values
    p = vars.p.values

    beta = 2 * zcmu0 * p / btor**2 * 100

    vars.beta.set_variable(beta, '%')

# Electron Beta (in %)
def calculate_betae(vars):
    zckb = constants.ZCKB
    zcmu0 = constants.ZCMU0
    btor = vars.btor.values
    ne = vars.ne.values
    te = vars.te.values

    betae = 2 * zcmu0 * ne * te * zckb / btor**2 * 100

    vars.betae.set_variable(betae, '%')

# Coulomb Logarithm TODO: what is 37.8? what are units?
def calculate_zlog(vars):
    ne = vars.ne.values
    te = vars.te.values

    zlog = 37.8 - np.log(ne**(1/2) / te)

    vars.zlog.set_variable(zlog, None)

# Collision Frequency (NU_{ei}) TODO: units?
def calculate_nuei(vars):
    zcf = constants.ZCF
    ne = vars.ne.values
    te = vars.te.values
    zeff = vars.zeff.values
    zlog = vars.zlog.values

    nuei = zcf * 2**(1/2) * ne * zlog * zeff / te**(3/2)

    vars.nuei.set_variable(nuei, None)

# OLD NOTE: Not sure what to call this, but it leads to the approx the correct NUSTI
def calculate_nuei2(vars):
    zcf = constants.ZCF
    ni = vars.ni.values
    ti = vars.ti.values
    zeff = vars.zeff.values
    zlog = vars.zlog.values

    nuei2 = zcf * 2**(1/2) * ni * zlog * zeff / ti**(3/2)

    vars.nuei2.set_variable(nuei2, None)

# Thermal Velocity of Electrons TODO: units?
def calculate_zvthe(vars):
    zckb = constants.ZCKB
    zcme = constants.ZCME
    te = vars.te.values

    zvthe = (2 * zckb * te / zcme)**(1/2)

    vars.zvthe.set_variable(zvthe, None)

# Thermal Velocity of Ions TODO: units?
def calculate_zvthi(vars):
    zckb = constants.ZCKB
    zcmp = constants.ZCMP
    aimass = vars.aimass.values
    ti = vars.ti.values

    zvthi = (zckb * ti / (zcmp * aimass))**(1/2)

    vars.zvthi.set_variable(zvthi, None)

# Electron Collisionality (NU^{*}_{e}) TODO: units?
# OLD NOTE: This is in approximate
# agreement with NUSTE in transp.  One source of the disagreement is
# likely because the modmmm7_1.f90 Coulomb logarithm (zlog) does not
# match perfectly with the TRANSP version (CLOGE).
def calculate_nuste(vars):
    eps = vars.eps.values
    nuei = vars.nuei.values
    q = vars.q.values
    rmaj = vars.rmaj.values
    zvthe = vars.zvthe.values

    nuste = nuei * eps**(-3/2) * q * rmaj / zvthe

    vars.nuste.set_variable(nuste, None)

# Ion Collisionality (NUSTI = NU^{*}_{i}) TODO: Units
# OLD NOTE: This is approx correct, but
# agreement is also somewhat time-dependent.  The issue is possibly due
# to the artifical AIMASS that we are using.  We likely also need to
# use the coulomb logarithm for ions as well.
def calculate_nusti(vars):
    zcme = constants.ZCME
    zcmp = constants.ZCMP
    eps = vars.eps.values
    q = vars.q.values
    nuei2 = vars.nuei2.values
    rmaj = vars.rmaj.values
    zvthi = vars.zvthi.values

    nusti = nuei2 * eps**(-3/2) * q * rmaj / (2 * zvthi) * (zcme / zcmp)**(1/2)

    vars.nusti.set_variable(nusti, None)

# Ion Gyrofrequency TODO: units
def calculate_zgyrfi(vars):
    zce = constants.ZCE
    zcmp = constants.ZCMP
    aimass = vars.aimass.values
    btor = vars.btor.values

    zgyrfi = zce * btor / (zcmp * aimass)

    vars.zgyrfi.set_variable(zgyrfi, None)

# Upper bound for ne, nh, te, and ti gradients in DRBM model (modmmm7_1.f90) TODO: units
def calculate_zgmax(vars):
    eps = vars.eps.values
    q = vars.q.values
    rmaj = vars.rmaj.values
    zgyrfi = vars.zgyrfi.values
    zvthi = vars.zvthi.values

    zgmax = rmaj / (zvthi / zgyrfi * q / eps)

    vars.zgmax.set_variable(zgmax, None)

# Calculates new variables needed for MMM and data display from CDF variables
# Values are stored to vars within each function call
def calculate_inputs(vars):

    # Some calculations depend on values from previous calculations
    calculate_vpol(vars)
    calculate_nh(vars)
    calculate_aimass(vars)
    calculate_rmin(vars)
    calculate_tau(vars)
    calculate_vtor(vars)
    calculate_vpar(vars)
    calculate_zeff(vars)
    calculate_btor(vars)
    calculate_eps(vars)
    calculate_p(vars)
    calculate_beta(vars)
    calculate_betae(vars)
    calculate_zlog(vars)
    calculate_nuei(vars)
    calculate_nuei2(vars)
    calculate_zvthe(vars)
    calculate_zvthi(vars)
    calculate_nuste(vars)
    calculate_nusti(vars)
    calculate_zgyrfi(vars)
    calculate_zgmax(vars)

    # Calculate gradients

    # Calculate inputs dependent on gradients

    return vars

if __name__ == '__main__':
    pass