# Standard Packages
import sys
sys.path.insert(0, '../')

# Local Packages
from main.enums import ShotType


class InputControls:
    def __init__(self, options):
        self.npoints = Control('npoints', 'Number of radial points', options.input_points)
        # Switches for component models
        self.cmodel_weiland = Control('cmodel_weiland', 'Weiland', 1)
        self.cmodel_dribm = Control('cmodel_dribm', 'DRIBM', 1 if options.shot_type != ShotType.NSTX else 0)
        self.cmodel_etg = Control('cmodel_etg', 'ETG', 1)
        self.cmodel_etgm = Control('cmodel_etgm', 'ETGM', 1)
        self.cmodel_mtm = Control('cmodel_mtm', 'MTM', 1)
        # Weiland real options
        self.exbs_weiland = Control('exbs_weiland', 'ExB shear coefficient', 1)
        self.mpsf_weiland = Control('mpsf_weiland', 'Momentum pinch scaling factor', 1)
        self.lbetd_weiland = Control('lbetd_weiland', 'Lower bound of electron thermal diffusivity', 0)
        self.ubetd_weiland = Control('ubetd_weiland', 'Upper bound of electron thermal diffusivity', 100)
        self.lbitd_weiland = Control('lbitd_weiland', 'Lower bound of ion thermal diffusivity', 0)
        self.ubitd_weiland = Control('ubitd_weiland', 'Upper bound of ion thermal diffusivity', 100)
        # DRIBM real options
        self.exbs_dribm = Control('exbs_dribm', 'ExB shear coefficient', 1)
        self.kyrhos_dribm = Control('kyrhos_dribm', 'kyrhos', 0.1)
        # MTM real options
        self.ky_kx_mtm = Control('ky_kx_mtm', 'ky/kx for MTM', 0.2)
        self.cf_mtm = Control('cf_mtm', 'calibration factor', 1)
        # ETG integer options
        self.jenko_threshold_etg = Control('jenko_threshold_etg', 'Jenko threshold', 2)
        # ETG real options
        self.cees_scale_etg = Control('cees_scale_etg', 'CEES scale', 0.06)
        self.ceem_scale_etg = Control('ceem_scale_etg', 'CEEM scale', 0.06)
        # ETGM integer options
        self.cl_etgm = Control('cl_etgm', 'Collisionless limit', 1)
        # ETGM real options
        self.exbs_etgm = Control('exbs_etgm', 'ExB shear coefficient', 0)
        self.kyrhoe_etgm = Control('kyrhoe_etgm', 'kyrhoe', 0.25)
        self.kyrhos_etgm = Control('kyrhos_etgm', 'kyrhos', 0.33)
        # Verbose level (int)
        self.lprint = Control('lprint', 'Verbose Level', 0)

    def get_mmm_header(self):
        return MMM_HEADER.format(
            npoints=self.npoints.get_values_int(),
            cmodel_weiland=self.cmodel_weiland.get_values_real(),
            cmodel_dribm=self.cmodel_dribm.get_values_real(),
            cmodel_etg=self.cmodel_etg.get_values_real(),
            cmodel_etgm=self.cmodel_etgm.get_values_real(),
            cmodel_mtm=self.cmodel_mtm.get_values_real(),
            exbs_weiland=self.exbs_weiland.get_values_real(),
            mpsf_weiland=self.mpsf_weiland.get_values_real(),
            lbetd_weiland=self.lbetd_weiland.get_values_real(),
            ubetd_weiland=self.ubetd_weiland.get_values_real(),
            lbitd_weiland=self.lbitd_weiland.get_values_real(),
            ubitd_weiland=self.ubitd_weiland.get_values_real(),
            exbs_dribm=self.exbs_dribm.get_values_real(),
            kyrhos_dribm=self.kyrhos_dribm.get_values_real(),
            ky_kx_mtm=self.ky_kx_mtm.get_values_real(),
            cf_mtm=self.cf_mtm.get_values_real(),
            jenko_threshold_etg=self.jenko_threshold_etg.get_values_int(),
            cees_scale_etg=self.cees_scale_etg.get_values_real(),
            ceem_scale_etg=self.ceem_scale_etg.get_values_real(),
            cl_etgm=self.cl_etgm.get_values_int(),
            exbs_etgm=self.exbs_etgm.get_values_real(),
            kyrhos_etgm=self.kyrhos_etgm.get_values_real(),
            kyrhoe_etgm=self.kyrhoe_etgm.get_values_real(),
            lprint=self.lprint.get_values_int(),
        )


class Control:
    def __init__(self, name, desc, values):
        self.name = name
        self.desc = desc
        self.values = values

    def get_values_real(self):
        return f'{self.values}D0'

    def get_values_int(self):
        return str(int(self.values))


# Header for MMM input file
MMM_HEADER = (
'''&testmmm_input_control
 npoints = {npoints}    ! Number of radial points
 input_kind = 1
/
&testmmm_input_1stkind
! This is a sample input file of the first kind
! of an NSTX discharge

!.. Switches for component models
!   1D0 - ON, 0D0 - OFF
cmodel  =
   {cmodel_weiland}    ! Weiland
   {cmodel_dribm}    ! DRIBM
   {cmodel_etg}    ! ETG
   {cmodel_etgm}    ! ETGM
   {cmodel_mtm}    ! MTM  
    
!.. Weiland real options
cW20 =
   {exbs_weiland}     ! ExB shear coefficient
   {mpsf_weiland}     ! Momentum pinch scaling factor
   {lbetd_weiland}     ! Lower bound of electron thermal diffusivity
   {ubetd_weiland}     ! Upper bound of electron thermal diffusivity
   {lbitd_weiland}     ! Lower bound of ion thermal diffusivity
   {ubitd_weiland}     ! Upper bound of ion thermal diffusivity

!.. DRIBM real options
cDBM =
   {exbs_dribm}     ! ExB shear coefficient
   {kyrhos_dribm}   ! kyrhos
   
!.. MTM real options
cMTM =
   {ky_kx_mtm}   ! ky/kx for MTM
   {cf_mtm}   ! calibration factor

!.. ETG integer options
lETG =
   {jenko_threshold_etg}       ! Jenko threshold
           ! applied to both electrostatic and electromagnetic regimes

!.. ETG real options
cETG =
   {cees_scale_etg}    ! CEES scale
   {ceem_scale_etg}    ! CEEM scale
   
!.. ETGM integer options
lETGM =
   {cl_etgm}      ! Collisionless limit

!.. ETGM real options
cETGM =
   {exbs_etgm}     ! ExB shear coefficient
   {kyrhos_etgm}   ! kyrhos
   {kyrhoe_etgm}   ! kyrhoe
   
lprint   = {lprint}      ! Verbose level\n\n''')


# For testing purposes: Run to print mmm input file header
if __name__ == '__main__':
    from main.options import Options
    Options.instance.set_options(shot_type = ShotType.NSTX, input_points = 51)
    print(InputControls(Options.instance).get_mmm_header())
