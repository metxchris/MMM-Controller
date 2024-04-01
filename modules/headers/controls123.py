# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')

# Local Packages
import modules.controls as controls
import settings


VALID_VERSIONS = ['#123']


def get_mmm_header(self):
    '''
    Gets the header for the MMM input file

    Raises:
    * ValueError: If settings.MMM_HEADER_VERSION not in VALID_VERSIONS
    '''

    # Temporary EPM switch

    if settings.MMM_HEADER_VERSION not in VALID_VERSIONS:
        raise ValueError(
            'Incorrect value for settings.MMM_HEADER_VERSION'
            f'\n\tExpected: {", ".join(VALID_VERSIONS)}'
            f'\n\tReceived: {settings.MMM_HEADER_VERSION}'
        )

    epm1 = ''
    epm2 = ''
    if settings.USE_EPM:
        epm1 = f'   {self.cmodel_epm.get_input_line()}'
        epm2 = (
            '!.. EPM integer options\n'
            'iEPM =\n'
            f'   {self.epm_direction.get_input_line()}'
            f'   {self.epm_n_start.get_input_line()}'
            f'   {self.epm_n_end.get_input_line()}'
            f'   {self.epm_n_step.get_input_line()}'
            f'   {self.epm_sum_modes.get_input_line()}'
            '\n'
            '!.. EPM real options\n'
            'rEPM =\n'
            f'   {self.epm_exbs.get_input_line()}'
            f'   {self.epm_xti_min.get_input_line()}'
            f'   {self.epm_xti_max.get_input_line()}'
            f'   {self.epm_xde_min.get_input_line()}'
            f'   {self.epm_xde_max.get_input_line()}'
            f'   {self.epm_xte_min.get_input_line()}'
            f'   {self.epm_xte_max.get_input_line()}'
            f'   {self.epm_chi_cal.get_input_line()}'
            '\n'
        )

    return (
        '&testmmm_input_control\n'
        f'   npoints = {self.input_points.get_input_line()}'
        f'   input_kind = 1\n'
        '/\n'
        '&testmmm_input_1stkind\n'
        '\n'
        '!.. Switches for component models (1D0 - ON, 0D0 - OFF)\n'
        'cmodel  =\n'
        f'   {self.cmodel_weiland.get_input_line()}'
        f'   {self.cmodel_dribm.get_input_line()}'
        f'   {self.cmodel_etgm.get_input_line()}'
        f'   {self.cmodel_mtm.get_input_line()}'
        f'{epm1}'
        '\n'
        '!.. MMM integer options (All Models)\n'
        'iMMM =\n'
        f'   {self.mmm_separate_conv_vel.get_input_line()}'
        f'   {self.mmm_convert_negative_chi.get_input_line()}'
        f'   {self.mmm_allow_negative_chi.get_input_line()}'
        f'   {self.mmm_use_solver_grads.get_input_line()}'
        f'   {self.mmm_limit_small_grads.get_input_line()}'
        '\n'
        '!.. MMM real options (All Models)\n'
        'rMMM =\n'
        f'   {self.mmm_xti_max.get_input_line()}'
        f'   {self.mmm_xte_max.get_input_line()}'
        f'   {self.mmm_xde_max.get_input_line()}'
        f'   {self.mmm_xdz_max.get_input_line()}'
        f'   {self.mmm_xvt_max.get_input_line()}'
        f'   {self.mmm_xvp_max.get_input_line()}'
        f'   {self.mmm_vti_max.get_input_line()}'
        f'   {self.mmm_vte_max.get_input_line()}'
        f'   {self.mmm_vde_max.get_input_line()}'
        f'   {self.mmm_vdz_max.get_input_line()}'
        f'   {self.mmm_vvt_max.get_input_line()}'
        f'   {self.mmm_vvp_max.get_input_line()}'
        '\n'
        '!.. Weiland integer options\n'
        'iW20 =\n'
        f'   {self.w20_shear_def.get_input_line()}'
        '\n'
        '!.. Weiland real options\n'
        'rW20 =\n'
        f'   {self.w20_exbs.get_input_line()}'
        f'   {self.w20_vconv_mult.get_input_line()}'
        f'   {self.w20_diff_mult.get_input_line()}'
        f'   {self.w20_kyrhos.get_input_line()}'
        '\n'
        + epm2 +
        '!.. DRIBM integer options\n'
        'iDBM =\n'
        f'   {self.dribm_diffusivity_type.get_input_line()}'
        f'   {self.dribm_direction.get_input_line()}'
        f'   {self.dribm_sum_modes.get_input_line()}'
        f'   {self.dribm_kyrhos_scan.get_input_line()}'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   {self.dribm_kyrhos_type.get_input_line()}'
        f'   {self.dribm_kyrhos_layers.get_input_line()}'
        '\n'
        '!.. DRIBM real options\n'
        'rDBM =\n'
        f'   {self.dribm_exbs.get_input_line()}'
        f'   {self.dribm_diff_mult.get_input_line()}'
        f'   {self.dribm_kyrhos_min.get_input_line()}'
        f'   {self.dribm_kyrhos_max.get_input_line()}'
        '\n'
        '!.. MTM integer options\n'
        'iMTM =\n'
        f'   {self.mtm_kyrhos_loops.get_input_line()}'
        f'   {self.mtm_kyrhos_layer_loops.get_input_line()}'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   {self.mtm_kyrhos_type.get_input_line()}'
        f'   {self.mtm_kyrhos_layers.get_input_line()}'
        '\n'
        '!.. MTM real options\n'
        'rMTM =\n'
        f'   {self.mtm_diff_mult.get_input_line()}'
        f'   {self.mtm_kyrhos_min.get_input_line()}'
        f'   {self.mtm_kyrhos_max.get_input_line()}'
        '\n'
        '!.. ETGM integer options\n'
        'iETGM =\n'
        f'   {self.etgm_diffusivity_type.get_input_line()}'
        f'   {self.etgm_jenko_threshold.get_input_line()}'
        f'   {self.etgm_direction.get_input_line()}'
        f'   {self.etgm_sum_modes.get_input_line()}'
        f'   {self.etgm_kyrhos_scan.get_input_line()}'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   0\n'
        f'   {self.etgm_kyrhos_type.get_input_line()}'
        f'   {self.etgm_kyrhos_layers.get_input_line()}'
        '\n'
        '!.. ETGM real options\n'
        'rETGM =\n'
        f'   {self.etgm_exbs.get_input_line()}'
        f'   {self.etgm_diff_mult.get_input_line()}'
        f'   {self.etgm_kyrhos_min.get_input_line()}'
        f'   {self.etgm_kyrhos_max.get_input_line()}'
        '\n'
    )


if __name__ == '__main__':
    settings.MMM_HEADER_VERSION = '#123'

    import modules.options
    options = modules.options.Options(runid='TEST', scan_num=373, input_points=51)

    '''Print sample MMM Header from user-specified Options'''
    controls = controls.InputControls(options)
    print(f'Example MMM Header:\n{"-"*50}')
    print(controls.get_mmm_header())
