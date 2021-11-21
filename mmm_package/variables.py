class Variable(object):
    def __init__(self, name, desc=None, cdfvar=None, mmmvar=None, units=None, values=None):
        self.name = name
        self.desc = desc
        self.cdfvar = cdfvar
        self.mmmvar = mmmvar
        self.units = units
        self.values = values

    def __str__(self):
        return str(self.name)

class Variables(object):
    def __init__(self):
        # CDF Independent Variables
        self.time = Variable('Time', cdfvar='TIME')
        self.x = Variable('X', cdfvar='X')
        self.xb = Variable('XB', cdfvar='XB')

        # CDF Variables
        self.aimp = Variable('AIMP', cdfvar='AIMP')
        self.arat = Variable('Aspect Ratio', cdfvar='ARAT')
        self.bz = Variable('BZ', cdfvar='BZ')
        self.elong = Variable('Elongation', cdfvar='ELONG')
        self.omega = Variable('OMEGA', cdfvar='OMEGA')
        self.ne = Variable('Electron Density', cdfvar='NE')
        self.nf = Variable('NF', cdfvar='BDENS')
        self.nh2 = Variable('NH', cdfvar='NH') # TODO: why do i have this?
        self.nd2 = Variable('ND', cdfvar='ND') # TODO: why do i have this?
        self.ni = Variable('Thermal Ion Density', cdfvar='NI')
        self.nz = Variable('NZ', cdfvar='NIMP')
        self.pcur = Variable('PCUR', cdfvar='PCUR')
        self.q = Variable('Safety Factor', cdfvar='Q')
        self.rmaj = Variable('Major Radius', cdfvar='RMJMP')
        self.te = Variable('Electron Temperature', cdfvar='TE')
        self.ti = Variable('Thermal Ion Temperature', cdfvar='TI')
        self.triang = Variable('TRIANG', cdfvar='TRIANG')
        self.vpold = Variable('VPOL', cdfvar='VPOLD_NC')
        self.vpolh = Variable('VPOL', cdfvar='VPOLH_NC')
        self.wexbs = Variable('ExB Shear Rate', cdfvar='SREXBA')
        self.zimp = Variable('ZIMP', cdfvar='XZIMP')

        # Calculated Variables (some are also in the CDF)
        # TODO: Check that calculated values match CDF values
        self.aimass = Variable('AIMASS')
        self.alphamhd = Variable('Alpha_MHD')
        self.betae = Variable('Electron Beta') # cdfvar='BETAE'
        self.btor = Variable('Toroidal Magnetic Field')
        self.eps = Variable('Inverse Aspect Ratio')
        self.nd = Variable('ND') # cdfvar='ND'
        self.nh = Variable('Hydrogenic Ion Density') # cdfvar='NH'
        self.nuei = Variable('Collision Frequency')
        self.nuei2 = Variable('NUEI2')
        self.nuste = Variable('Electron Collisionality') # cdfvar='NUSTE'
        self.nusti = Variable('Ion Collisionality') # cdfvar='NUSTI'
        self.p = Variable('Plasma Pressure') # cdfvar='PPLAS'
        self.raxis = Variable('RAXIS')
        self.rmin = Variable('Minor Radius')
        self.shat = Variable('Effective Magnetic Shear') # cdfvar='SHAT'
        self.shear = Variable('Magnetic Shear')
        self.vpar = Variable('VPAR')
        self.vtor = Variable('VTOR')
        self.tau = Variable('Temperature Ratio')
        self.zcf = Variable('Collision Frequency Factor')
        self.zeff = Variable('Effective Charge') # cdfvar='ZEFF'
        self.zgmax = Variable('ZGMAX')
        self.zgyrfi = Variable('Ion Gyrofrequency')
        self.zlari = Variable('ZLARI')
        self.zlarpo = Variable('ZLARPO')
        self.zlog = Variable('Coulomb Logarithm')
        self.zvthe = Variable('Electron Thermal Velocity')
        self.zvthi = Variable('Ion Thermal Velocity')

        # Calculated Gradients
        self.gne = Variable('Electron Density Gradient')
        self.gnh = Variable('Hydrogenic Ion Density Gradient')
        self.gni = Variable('Thermal Ion Density Gradient')
        self.gnz = Variable('NZ Gradient')
        self.gq = Variable('Safety Factor Gradient')
        self.gte = Variable('Electron Temperature Gradient')
        self.gti = Variable('Thermal Ion Temperature Gradient')
        self.gvpar = Variable('VPAR Gradient')
        self.gvpol = Variable('VPOL Gradient')
        self.gvtor = Variable('VTOR Gradient')

    def get_variables(self):
        return [var for var in dir(self) if not callable(getattr(self, var)) and not var.startswith("__")]
        
    def get_nonzero_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).values is not None]

    def get_cdf_variables(self):
        vars = self.get_variables()
        return [var for var in vars if getattr(self, var).cdfvar is not None]

    def print_nonzero_variables(self):
        vars = self.get_nonzero_variables()
        for var in vars:
            print(var + ", "
                  + str(getattr(self, var).name) + ", "
                  + str(getattr(self, var).desc) + ", " 
                  + str(getattr(self, var).units) + ", "
                  + str(getattr(self, var).values.shape))

    def __str__(self):
        return str(self.get_nonzero_variables())

    def get_nzones(self):
        return self.x.values.shape[0] if self.xb.values is not None and self.xb.values.ndim > 0 else 0

    def get_ntimes(self):
        return self.x.values.shape[1] if self.xb.values is not None and self.xb.values.ndim > 1 else 0
