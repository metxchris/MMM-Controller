#!/usr/bin/python3

"""Handles the entire flow of data both to and from MMM

The controller is meant to be ran directly to execute one of three different
MMM runs.  The different run types the controller supports are:
* Basic Run: MMM is ran once using data obtained from TRANSP.  Various profile
  plot PDFs are created, and all used data is saved to CSVs.  This run is
  fundamental to storing all non-scanned variable data, so a basic run is
  also executed at the start of both variable and control scans.
* Variable Scan: Values of a specified variable are adjusted using a range of
  scan factors, which are multiplicative factors.  MMM is ran once for each
  factor in the scan range.  Valid variables to scan must be members of the
  InputVariables class.
* Control Scan: Values of a specified control are adjusted using a range of
  scan factors, in similar fashion to a variable scan.  Valid controls to
  scan must be members of the InputControls class.

When conducting either a variable or control scan, values of all input,
additional, and output variables are saved to CSV for each factor in the scan
range.  Factor files contain variable data as functions of rmin (or rho),
where each file corresponds to a different scan factor value (noted by the
filename). After the scan is completed, this saved data is then reshaped and
saved again into new rho files, which contain values of each variable type at
a specified point of rho.  Rho files contain variable data as function of the
scanned variable at each point of rho, where each file corresponds to a
different rho value (noted by the filename).

A new scan number is created each time a new run is executed.  This means that
new data generated by MMM will never overwrite previously generated data. The
trade-off of this is that many scan folders can be created over time, and it
is up to the user to maintain the folders as necessary.

In terms of plot generation, the controller only produces plots of base profiles
(unaltered by scan factors).  Plots from data stored in scan factor files or
rho files should be generated by directly running the various modules
provided in the plotting directory.

Example Usage:
* See commands listed at the bottom of this file
"""

# 3rd Party Packages
import numpy as np

# Local Packages
import settings
import modules.controls
import modules.constants
import modules.options
import modules.calculations as calculations
import modules.adjustments as adjustments
import modules.datahelper as datahelper
import modules.mmm as mmm
import modules.reshaper as reshaper
import modules.utils as utils
import plotting.modules.profiles as profiles
from modules.enums import ShotType, ScanType, ProfileType


def _execute_variable_scan(mmm_vars, controls):
    '''
    Executes an input variable scan, where the values of an input variable are
    varied over a specified range and are then sent to the MMM driver for
    each value of the range

    Create a copy of mmm_vars as modified_vars to keep variables that are
    modified over the course of the scan separate from base MMM input
    variables.  For each factor of the scan_range, we modify the value of the
    specified var_to_scan, and then adjust any dependent variables. The MMM
    driver is ran each time var_to_scan is adjusted, and all input and output
    variable data is saved to a subfolder named after var_to_scan.
    Afterwards, the saved CSV data is reshaped into data dependent on the
    scanned parameter, and is saved to another set of CSV within a new
    subfolder labeled rho.

    Parameter scan PDFs are not produced here, and the output data is intended
    to be plotted by a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write the MMM input file
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    runid = mmm_vars.options.runid
    scan_num = mmm_vars.options.scan_num
    scan_range = mmm_vars.options.scan_range
    var_to_scan = mmm_vars.options.var_to_scan

    for i, scan_factor in enumerate(scan_range):
        print(f'{runid}.{scan_num} {var_to_scan} scan: {i + 1} / {len(scan_range)}')
        adjusted_vars = adjustments.adjust_scanned_variable(mmm_vars, scan_factor)
        adjusted_vars.save(scan_factor)
        output_vars = mmm.run_wrapper(adjusted_vars, controls)
        calculations.calculate_output_variables(mmm_vars, output_vars, controls)
        output_vars.save(scan_factor)


def _execute_control_scan(mmm_vars, controls):
    '''
    Executes an input control scan, where the values of an input control are
    varied over a specified range and are then sent to the MMM driver for
    each value of the range

    Parameter scan PDFs are not produced here, and the output data is intended
    to be plotted by a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    # Create a reference to control being scanned in InputControls. Modifying
    # scanned_control values will modify its corresponding values in
    # adjusted_controls

    runid = mmm_vars.options.runid
    scan_num = mmm_vars.options.scan_num
    scan_range = mmm_vars.options.scan_range
    var_to_scan = mmm_vars.options.var_to_scan

    adjusted_controls = datahelper.deepcopy_data(controls)
    scanned_control = adjusted_controls.get_scanned_control()
    base_control = controls.get_scanned_control()

    for i, scan_factor in enumerate(scan_range):
        print(f'{runid}.{scan_num} {var_to_scan} scan: {i + 1} / {len(scan_range)}')
        scanned_control.values = scan_factor * base_control.values
        mmm_vars.save(scan_factor)
        adjusted_controls.save(scan_factor)
        output_vars = mmm.run_wrapper(mmm_vars, adjusted_controls)
        calculations.calculate_output_variables(mmm_vars, output_vars, controls)
        output_vars.save(scan_factor)


def _execute_time_scan(mmm_vars, controls):
    '''
    Executes an input time scan

    Parameter scan PDFs are not produced here, and the output data is intended
    to be plotted by a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    # Set the time scan range
    mmm_vars.options.set_time_ranges(mmm_vars.time.values)

    # Save options again to save the computed time ranges
    mmm_vars.options.save()

    scan_range_idxs = mmm_vars.options.scan_range_idxs
    var_to_scan = mmm_vars.options.var_to_scan

    for i, time_idx in enumerate(scan_range_idxs):
        print(f'{options.runid}.{options.scan_num} {var_to_scan} scan: {i + 1} / {len(scan_range_idxs)}')
        options.time_idx = time_idx
        options.time_str = options.scan_range[i]
        time_scan_str = f'{float(options.time_str):{modules.constants.SCAN_FACTOR_FMT}}'
        mmm_vars.save(time_scan_str)
        output_vars = mmm.run_wrapper(mmm_vars, controls)
        calculations.calculate_output_variables(mmm_vars, output_vars, controls)
        output_vars.save(time_scan_str)


def main(scanned_vars, controls):
    '''
    Runs the controller for MMM

    Needed output folders are created and a unique scan number is chosen for
    storing output data. All input variable objects are initialized and
    corresponding plot PDFs are created.  The MMM driver is then ran once,
    and then an optional variable scan can be ran afterwards.

    Parameters:
    * scanned_vars (Dict): Dictionary of variables being scanned
        - keys (str | None): The variable being scanned
        - values (np.ndarray | None): The range of factors to scan over
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    utils.init_logging()
    options = controls.options  # Creates a reference

    # Remove default option when a specific scan option is specified
    if None in scanned_vars and len(scanned_vars) > 1:
        del scanned_vars[None]

    # TODO: Add validation for all items in scanned_vars
    for adjustment_name, scan_range in scanned_vars.items():
        options.scan_num = utils.get_scan_num(options.runid)
        options.set(adjustment_name=adjustment_name, scan_range=scan_range)

        print(f'\nRunning MMM Controller for {options.runid}, scan {options.scan_num}...')

        utils.init_output_dirs(options)

        mmm_vars, cdf_vars, __ = datahelper.initialize_variables(options)
        output_vars = mmm.run_wrapper(mmm_vars, controls)
        calculations.calculate_output_variables(mmm_vars, output_vars, controls)

        options.save()
        controls.save()
        mmm_vars.save()
        output_vars.save()

        if settings.MAKE_PROFILE_PDFS:
            profiles.plot_profiles(ProfileType.INPUT, mmm_vars)
            profiles.plot_profiles(ProfileType.ADDITIONAL, mmm_vars)
            profiles.plot_profiles(ProfileType.COMPARED, mmm_vars, cdf_vars)
            profiles.plot_profiles(ProfileType.OUTPUT, output_vars)

        # Variable and control scans
        if options.scan_type.value:
            if options.scan_type is ScanType.VARIABLE:
                _execute_variable_scan(mmm_vars, controls)
            elif options.scan_type is ScanType.CONTROL:
                _execute_control_scan(mmm_vars, controls)
            elif options.scan_type is ScanType.TIME:
                _execute_time_scan(mmm_vars, controls)

            reshaper.create_rho_files(options)
            print(f'\nScan complete: {options.runid}, scan {options.scan_num}, {options.var_to_scan}\n')


# Run this file directly to plot variable profiles and run the MMM driver
if __name__ == '__main__':
    import sys
    np.set_printoptions(threshold=sys.maxsize)  # Print all array values (diagnostic purposes)

    scanned_vars = {None: None}

    '''
    TRANSP Data:
    * Uncomment the line you wish to use
    # '''
    runid, shot_type, input_time = '121123K55', ShotType.NSTU, 11.8

    runid, shot_type, input_time = '120968A02', ShotType.NSTX, 0.56  # High
    runid, shot_type, input_time = '120968A02', ShotType.NSTX, 0.752  # High
    # runid, shot_type, input_time = '120982A09', ShotType.NSTX, 0.62  # Low
    # runid, shot_type, input_time = '129016A04', ShotType.NSTX, 0.46
    # runid, shot_type, input_time = '129017A04', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '129018A02', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '129019A02', ShotType.NSTX, 0.62
    # runid, shot_type, input_time = '129020A02', ShotType.NSTX, 0.56
    # runid, shot_type, input_time = '129041A10', ShotType.NSTX, 0.49  # Medium
    runid, shot_type, input_time = '138536A01', ShotType.NSTX, 0.629
    # runid, shot_type, input_time = '141007A10', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '141031A01', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '141032A01', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '141040A01', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '141716A80', ShotType.NSTX, 0.5

    # runid, shot_type, input_time = '98777V06', ShotType.D3D, 5.6
    # runid, shot_type, input_time = '101381J05', ShotType.D3D, 5.6
    # runid, shot_type, input_time = '101381T31', ShotType.D3D, 5.6
    # runid, shot_type, input_time = '101391J08', ShotType.D3D, 5.6
    # runid, shot_type, input_time = '118341T54', ShotType.D3D, 5.6
    # runid, shot_type, input_time = '132017T01', ShotType.D3D, 2
    # runid, shot_type, input_time = '132411T02', ShotType.D3D, 0.56
    # runid, shot_type, input_time = '132498J05', ShotType.D3D, 0.56
    # runid, shot_type, input_time = '141552A01', ShotType.D3D, 2
    # runid, shot_type, input_time = '150840T02', ShotType.D3D, 2
    # runid, shot_type, input_time = '153283T50', ShotType.D3D, 2

    # runid, shot_type, input_time = '85126T02', ShotType.EAST, 2
    # runid, shot_type, input_time = '85610T01', ShotType.EAST, 2
    # runid, shot_type, input_time = '85122T04', ShotType.EAST, 2
    # runid, shot_type, input_time = '80208T04', ShotType.EAST, 2
    # runid, shot_type, input_time = '90328T01', ShotType.EAST, 2

    # runid, shot_type, input_time = '15334T03', ShotType.KSTR, 2
    # runid, shot_type, input_time = '16295T10', ShotType.KSTR, 2
    # runid, shot_type, input_time = '16297T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '16299T01', ShotType.KSTR, 2.25
    # runid, shot_type, input_time = '16325T10', ShotType.KSTR, 4.023
    # runid, shot_type, input_time = '16901T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18399T05', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18400T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18402T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18404T05', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18476T02', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18477T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18492T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18495T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18499T01', ShotType.KSTR, 2
    # runid, shot_type, input_time = '18602T01', ShotType.KSTR, 2

    # runid, shot_type, input_time = '08505Z06', ShotType.MAST, 2
    # runid, shot_type, input_time = '18696B01', ShotType.MAST, 2
    # runid, shot_type, input_time = '22341P37', ShotType.MAST, 0.25
    # runid, shot_type, input_time = '24899R05', ShotType.MAST, 2
    # runid, shot_type, input_time = '27527M34', ShotType.MAST, 2
    # runid, shot_type, input_time = '29271A01', ShotType.MAST, 2
    # runid, shot_type, input_time = '29976U69', ShotType.MAST, 2
    # runid, shot_type, input_time = '45424H01', ShotType.MAST, 2

    # runid, shot_type, input_time = '18696R06', ShotType.ITER, 2
    # runid, shot_type, input_time = '38265R80', ShotType.ITER, 2
    # runid, shot_type, input_time = '50000A10', ShotType.ITER, 2
    # runid, shot_type, input_time = '80200A13', ShotType.ITER, 1.0
    # runid, shot_type, input_time = '80300A02', ShotType.ITER, 2

    # runid, shot_type, input_time = 'TEST', ShotType.NSTX, 0.5

    # # !runid, shot_type, input_time = '144449T54', ShotType.D3D, 2
    # # !runid, shot_type, input_time = '85124T02', ShotType.EAST, 2
    # # !runid, shot_type, input_time = '16296T10', ShotType.KSTR, 2
    # # !runid, shot_type, input_time = '16949T02', ShotType.KSTR, 2
    # # runid, shot_type, input_time = '20100J26', ShotType.ITER, 2
    # runid, shot_type, input_time = '20102A12', ShotType.ITER, 2
    # runid, shot_type, input_time = '20105J18', ShotType.ITER, 2
    # runid, shot_type, input_time = '20160H15', ShotType.ITER, 2
    #     runid, shot_type, input_time = '28014T01', ShotType.ITER, 2
    # runid, shot_type, input_time = '33610E10', ShotType.ITER, 2
    # runid, shot_type, input_time = '33616N07', ShotType.ITER, 2
    # runid, shot_type, input_time = '33701K02', ShotType.ITER, 2
    # runid, shot_type, input_time = '35601B11', ShotType.ITER, 2
    # runid, shot_type, input_time = '59100A05', ShotType.ITER, 2

    '''
    Scanned Variables:
    * Uncomment the lines you wish to include in scanned_vars
    * MMM will only run once if no specific scan is specified below
    '''

    # scanned_vars['etgm_kyrhos_min'] = np.arange(start=1e-6, stop=1.01 + 1e-6, step=0.005)
    # scanned_vars['dribm_kyrhos_min'] = np.arange(start=0.225, stop=0.235 + 1e-6, step=0.001)
    # scanned_vars['mtm_kyrhos_min'] = np.arange(start=0.1, stop=3 + 1e-6, step=0.1)
    # scanned_vars['mtm_kyrhos_max'] = np.arange(start=5, stop=10.0 + 1e-6, step=0.5)
    # scanned_vars['weiland_kyrhos'] = np.arange(start=0.2, stop=0.4 + 1e-6, step=0.005)

    # scanned_vars['etgm_kyrhos_min'] = np.arange(start=1, stop=40 + 1e-6, step=0.2)
    # scanned_vars['etgm_alpha_mult'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['etgm_betae_mult'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['etgm_nuei_mult'] = np.arange(start=0.2, stop=4 + 1e-6, step=0.025)
    # scanned_vars['etgm_vthe_mult'] = np.arange(start=0.2, stop=4 + 1e-6, step=0.025)
    # scanned_vars['etgm_betaep_mult'] = np.arange(start=0, stop=3 + 1e-6, step=0.025)

    # scanned_vars['betaeu'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['bu'] = np.arange(start=0.2, stop=3 + 1e-6, step=0.02)
    # scanned_vars['gne'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)
    # scanned_vars['gte'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)
    # scanned_vars['ne'] = np.arange(start=0.2, stop=4 + 1e-6, step=0.025)
    # scanned_vars['q'] = np.arange(start=0.5, stop=2 + 1e-6, step=0.01)
    # scanned_vars['shear'] = np.arange(start=-3, stop=3 + 1e-6, step=0.05)
    # scanned_vars['te'] = np.arange(start=0.2, stop=4 + 1e-6, step=0.025)

    # scanned_vars['betaeu_alphaconst'] = np.arange(start=0.05, stop=3 + 1e-6, step=0.01)
    # scanned_vars['betaeu'] = np.arange(start=0.05, stop=3 + 1e-6, step=0.05)
    # scanned_vars['nuei_alphaconst'] = np.arange(start=0.2, stop=4 + 1e-6, step=0.02)
    # scanned_vars['nuei_lareunitconst'] = np.arange(start=0.2, stop=4 + 1e-6, step=0.02)
    # scanned_vars['ti'] = np.arange(start=0.2, stop=5 + 1e-6, step=0.02)

    # scanned_vars['zeff'] = np.arange(start=0.1, stop=4 + 1e-6, step=0.02)**2
    # scanned_vars['ai'] = np.arange(start=0.5, stop=1.5 + 1e-6, step=0.1)

    # gte = 0
    # scanned_vars['gne'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)

    # gne = 0
    # scanned_vars['gte'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)
    # scanned_vars['q'] = np.arange(start=0.5, stop=2 + 1e-6, step=0.01)

    # gneabs
    # scanned_vars['gne'] = np.arange(start=-4.0, stop=8 + 1e-6, step=0.05)
    # scanned_vars['gne_alphaconst'] = np.arange(start=-4.0, stop=8 + 1e-6, step=0.1)

    # gne threshold
    # scanned_vars['gne'] = np.arange(start=0.00, stop=201 + 1e-6, step=0.5)

    # gte threshold
    # scanned_vars['gte'] = np.arange(start=0.00, stop=201 + 1e-6, step=0.5)

    # scanned_vars['mtm_kyrhos'] = np.arange(start=0.02, stop=32 + 1e-6, step=0.02)

    # scanned_vars['etgm_kxoky_mult'] = np.arange(start=0, stop=1 + 1e-6, step=0.02)
    # scanned_vars['etgm_extra_mult'] = np.arange(start=0.1, stop=6 + 1e-6, step=0.05)

    # time scan (options.normalize_time_range = 0)
    # scanned_vars['time'] = np.arange(start=0.3, stop=0.6 + 1e-6, step=0.001)

    # normalized time scan (options.normalize_time_range = 1)
    # scanned_vars['time'] = np.linspace(start=0, stop=1, num=40)
    scanned_vars['time'] = np.linspace(start=0.0, stop=1.0, num=100)
    # scanned_vars['time'] = np.linspace(start=0.0, stop=1.0, num=300)


    # DRIBM Scans
    # scanned_vars['time'] = np.linspace(start=0.0, stop=1.0, num=100)
    # scanned_vars['ah'] = np.arange(start=0.5, stop=1.5 + 1e-6, step=0.025)
    # scanned_vars['betaeu'] = np.arange(start=0.05, stop=3 + 1e-6, step=0.05)
    # scanned_vars['gte'] = np.arange(start=0.1, stop=6 + 1e-6, step=0.1)
    # scanned_vars['gne'] = np.arange(start=0.1, stop=6 + 1e-6, step=0.1)
    # scanned_vars['q'] = np.arange(start=0.5, stop=2 + 1e-6, step=0.025)
    # scanned_vars['shear'] = np.arange(start=-3, stop=3 + 1e-6, step=0.1)
    # scanned_vars['dribm_kpsh_mult'] = np.arange(start=0.5, stop=2.0 + 1e-6, step=0.025)
    # scanned_vars['dribm_ti_te_mult'] = np.arange(start=0.25, stop=4 + 1e-6, step=0.05)
    # scanned_vars['dribm_vei_mult'] = np.arange(start=0.25, stop=4 + 1e-6, step=0.05)
    # scanned_vars['etae'] = np.arange(start=0.25, stop=4 + 1e-6, step=0.05)

    '''
    Options:
    * input_points is the number of points to use when making the MMM input file
    * Set input_points = None to match the number of points used in the CDF
    * apply_smoothing enables smoothing of all variables that have a smooth value set in the Variables class
    '''
    options = modules.options.Options(
        runid=runid,
        shot_type=shot_type,
        input_time=input_time,
        input_points=101,
        apply_smoothing=1,
        use_gtezero=0,
        use_gnezero=0,
        use_gneabs=0,
        use_gnethreshold=0,
        use_gtethreshold=0,
        use_etgm_btor=0,
        use_experimental_profiles=0,
    )

    '''
    Input Controls:
    * cmodel enables (disables) the corresponding model if set to 1 (0)
    '''
    wexb = 0
    controls = modules.controls.InputControls(
        options,
        # CMODEL
        cmodel_weiland=0,
        cmodel_dribm=0,
        cmodel_epm=1,
        cmodel_etgm=0,
        cmodel_mtm=0,
        cmodel_etg=0,
        # W20
        weiland_shear_def=0,
        weiland_extra_int=0,
        weiland_exbs=wexb,
        weiland_kyrhos=0.316,
        # EPM
        epm_n_start=1,
        epm_n_end=10,
        epm_n_step=1,
        # DRBM
        dribm_exbs=wexb,
        dribm_direction=1,
        dribm_sum_modes=0,
        dribm_kyrhos_scan=0,
        dribm_kyrhos_min=1,
        dribm_kyrhos_max=0.3,
        dribm_kxoky=1,
        dribm_sat_expo=2,
        dribm_kyrhos_type=1,
        dribm_chi_max_cal=0.15,
        dribm_chi_sum_cal=0.15 * 0.05,
        dribm_chi2_max_cal=0.45,
        dribm_chi2_sum_cal=0.45 * 0.05,
        dribm_kpsh_mult=1,
        dribm_ti_te_mult=1,
        dribm_vei_mult=1,
        dribm_diffusivity_type=1,
        # ETGM
        etgm_exbs=wexb,
        etgm_direction=0,
        etgm_sum_modes=0,
        etgm_kyrhos_scan=200,
        etgm_kyrhos_min=1,
        etgm_kyrhos_max=50,
        etgm_kxoky=0.5,
        etgm_sat_expo=2,
        etgm_kyrhos_type=1,
        etgm_xte_max_cal=0.05,
        etgm_xte_sum_cal=0.05 * 0.05,
        etgm_xte2_max_cal=1,
        etgm_xte2_sum_cal=1 * 0.05,
        etgm_alpha_mult=1,  # Extra
        etgm_betae_mult=1,  # Extra
        etgm_nuei_mult=1,  # Extra
        etgm_vthe_mult=1,  # Extra
        etgm_betaep_mult=1,  # Extra
        etgm_disable_geometry=0,  # Extra
        etgm_electrostatic=0,  # Extra
        # MTM
        mtm_kyrhos_loops=100,
        mtm_kyrhos_min=0.005,
        mtm_kyrhos_max=10,
        mtm_kyrhos_type=1,
        mtm_cf=0.001,
        # Verbose level
        lprint=0,
    )

    '''
    Output Profile Comparisons:
    '''
    # Base
    # controls.etgm_kyrhoe_scan.values = 0
    # controls.etgm_kyrhos_scan.values = 0
    # controls.etgm_exbs.values = 1
    # options.use_gnezero = 0
    # options.use_gtezero = 0
    # shafranov shift (recompile MMM)
    # shear vs shat_gxi (recompile MMM)
    # controls.etgm_use_gne_in.values = 0

    settings.AUTO_OPEN_PDFS = 0
    settings.MAKE_PROFILE_PDFS = 0
    settings.PRINT_MMM_RESPONSE = 0
    settings.SAVE_ADDITIONAL_VARIABLES = 0
    settings.STARTING_SCAN_NUMBER = 1

    main(scanned_vars, controls)
